#!/usr/bin/bash
# pegar arquivos de alinhamento gerado pelo guidance e rodar swamp - filtros para rodar paml
# Fernanda T. 08-09-2020
# Usage: bash swampPipeline.sh geneList.txt swamp_path
# Arquivos swamp_gene_stats.txt e genes_threshold_swamp.txt estao sempre sendo sobrescritos; logo, apagar se precisar recomecar a corrida

gene=($(cat $1)) # geneList.txt - primeiro parametro, com lista de genes (um por linha)
swamp_path=$2 # swamp_path - segundo parametro, dir que sera criado e onde os outputs cada gene vao estar, sem o ultimo "/"
guidance_path='/media/labgenoma/DATAPART5/vera/Selecao' # path completo onde estao as pastas de cada gene gerado pelo guidance, sem o ultimo "/"
work_dir='/media/labgenoma/DATAPART5/vera/swamp_teste' # path completo de onde tu esta rodando o pipeline, onde editCTL_v2.py e otters_v2.tre e branch_v2.txt devem estar, sem o ultimo "/"
swamp_instal_path='/home/labgenoma/programas/SWAMP' # path onde esta o script SWAMP.py, sem o ultimo "/"
threshold='50' # porcentagem maxima de codons mascarados
num_sp='13' # numero de sequencias/especies/individuos

mkdir $swamp_path

for (( i=0; i<"${#gene[@]}"; i++ )); do
	echo -e "\n####################"
        echo -e "#### Working on gene ${gene[i]}!!!"
	echo -e "####################\n"

	# Create a folder for each gene
	echo -e "### Creating ./$swamp_path/${gene[i]} directory"
	mkdir $swamp_path/${gene[i]}

	# Copy the guidance alignment
	cp $guidance_path/${gene[i]}/Seqs.Orig_DNA.fas.FIXED $swamp_path/${gene[i]}/${gene[i]}.fasta
	cd $swamp_path/${gene[i]}

	# Edit the fasta header to change Aonyx capensis name
	echo -e "### Editing ./$swamp_path/${gene[i]}/${gene[i]}.fasta file"
	sed 's/Aonyx_capensis_3349/Aocape3349/g' ${gene[i]}.fasta | sed 's/Aonyx_capensis_4194/Aocape4194/g' > ${gene[i]}.renamed.fasta

	# Sort the fasta file
	perl -pe 's/[\r\n]+/;/g; s/>/\n>/g' ${gene[i]}.renamed.fasta | sort -t"[" -k2,2V | sed 's/;/\n/g' | sed '/^$/d' > ${gene[i]}.renamed.sorted.fasta

	# Convert fasta to phylip using PRANK
	prank -d=${gene[i]}.renamed.sorted.fasta -o=${gene[i]} -f=phylips -convert

	# Edit the CTL file
	echo -e "### Running editCTL_v2.py for gene ${gene[i]} (...)"
	touch ${gene[i]}.ctl
	python $work_dir/editCTL_v2.py $PWD/ ${gene[i]}

	# Run codeml for each gene
	echo -e "### Running codeml for gene ${gene[i]} (...)"
	cp $work_dir/otters_v2.tre .
	yes "" | codeml ${gene[i]}.ctl &> codeml.log # esse yes eh para ele passar por pedido de "press enter"

	# Run SWAMP, print alignment and save to file
	echo -e "### Running swamp for gene ${gene[i]} (...)"
	python $swamp_instal_path/SWAMP.py -i $PWD -b $work_dir/branch_v2.txt -t 2 -w 20 &> swamp.log
	echo -e "swamp results for gene ${gene[i]}" >> $work_dir/$swamp_path/swamp_gene_stats.txt
	python $swamp_instal_path/SWAMP.py --print-alignment ${gene[i]}_masked.phy >> $work_dir/$swamp_path/swamp_gene_stats.txt

	# Choose genes to run PAML based on the masked threshold
	codons_total=$(awk -v i="$work_dir/$swamp_path/${gene[i]}/${gene[i]}.phy" '$2 == i {print $5}' swamp.log)
	echo -e "Number of ${gene[i]} codons: $codons_total * $num_sp"
	codons_masked=$(awk '$1 == "masked" {print $2}' swamp.log)
	echo -e "Number of ${gene[i]} masked codons: $codons_masked"

	check=$(awk '$1 == "1" {print $1}' swamp.log) # aqui eh para checar se o gene ja nao caiu no filtro de < 33 codons
	if [ -z "$check" ]; then

		if [[ "$codons_masked" -gt 0 ]]; then # aqui eh pra checar se teve algum codon mascarado
			prop_masked=$(( ($codons_masked*100)/($codons_total*$num_sp) ))
			echo -e "Percentage of ${gene[i]} masked codons: $prop_masked%"

			if (( "$prop_masked" < "$threshold" )); then
				echo -e "${gene[i]}" >> $work_dir/$swamp_path/genes_$threshold\threshold_pass_swamp.txt
			fi
		else
			echo -e "${gene[i]}" >> $work_dir/$swamp_path/genes_$threshold\threshold_pass_swamp.txt
		fi
	else
		echo -e "${gene[i]} has sequence(s) containing fewer than 33 informative codons after masking!"
	fi

	cd $work_dir
done
