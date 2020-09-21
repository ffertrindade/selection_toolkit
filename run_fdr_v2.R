# Rodar FDR para resultado de selecao com lista de valores de p
# Fernanda T 21-09-2020

## Check for the arguments
library(foreach)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Rscript run_fdr_v2.R geneList number_species sp1_pval.txt,sp2_pval.txt,(...)", call.=FALSE)
}
if(require(doMC)){
  registerDoMC(10)
} else {
  registerDoSEQ()
}

# Parameters
#setwd("C:/Users/ffert/OneDrive - PUCRS - BR/Papers/vera")
geneList <- read.table(args[1], header = FALSE) # genes_50threshold_pass_swamp.txt
numSp <- as.numeric(args[2]) # 11
fileList <- unlist(strsplit(args[3],split=",")) #Acapensis_pvalue_v2.txt,Acinereus_pvalue_v2.txt,(...)
#geneList <- read.table("genes_50threshold_pass_swamp - Copy.txt", header = FALSE) # genes_50threshold_pass_swamp.txt
#numSp <- as.numeric('2')
#fileList <- unlist(strsplit("Acapensis_pvalue_v2.txt,Acinereus_pvalue_v2.txt",split=",")) #Acapensis_pvalue_v2.txt,Acinereus_pvalue_v2.txt

# Read pvalues for geneList in all species
for (file in fileList) {
  f <- read.table(file, header = FALSE)
  for (i in 1:length(geneList$V1)) {
    gen <- paste("^", geneList$V1[i],"$", sep="", collapse="")
    pos <- print(grep(gen, f$V1))
    if (length(pos) == 0) {
      geneList$V2[i] <- NA
    } else {
      geneList$V2[i] <- f$V2[pos]
    }
  }
  names(geneList)[names(geneList) == 'V2'] <- file
}

# Run FDR for each gene correcting by numSp
df <- data.frame(t(geneList))
colnames(df) <- geneList$V1
df <- df[-1,]

qval <- foreach(i=1:length(df)) %dopar% {
  p.adjust(as.double(as.character(df[,i])), method = "fdr", n = numSp)
}
names(qval) <- geneList$V1

# Output
fout <- paste("fdr_",numSp,"sp_",length(geneList$V1),"genes.txt", sep="", collapse="")
write.table(qval, file = fout)

temp <- read.table(fout, header = TRUE, sep = " ")
row.names(temp) <- fileList
temp <- t(temp)
write.table(temp, file = fout)

