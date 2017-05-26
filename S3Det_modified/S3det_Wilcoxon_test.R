len<-length(commandArgs())
totalparametros<-len-2;
arrayParam <- commandArgs()[ (len-(totalparametros-1)) : len]
print(arrayParam[1])
setwd(arrayParam[1])

# library(exactRankTests)
matrix<-read.table("tmp_Wilcoxon_test_input.txt",sep="\t",skip=1,strip.white=TRUE)


if(length(matrix$V1)>50){mode=0}else{mode=1}
number_of_test<-dim(matrix)[2]-1

WILTEST<-numeric(number_of_test)

for(i in 1:number_of_test)
{
	WIL_TEST <- wilcox.test(matrix[,i],matrix[,i+1], alternative = "less",exact=mode)
	WILTEST[i] <-WIL_TEST$p.value 
}

write.table(WILTEST , file="tmp_Wilcoxon_test_output.txt",append = FALSE, quote=FALSE, sep="\t")
