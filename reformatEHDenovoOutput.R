input=read.table('EHDnResults.combinedCounts.filtered.bed',sep='\t',header=F,stringsAsFactors=F)
samples=read.table('manifest.txt',header=F,stringsAsFactors = F)
output=data.frame(chr=character(nrow(input)),start=numeric(nrow(input)),end=numeric(nrow(input)),motif=character(nrow(input)),stringsAsFactors=F)
sampleMatrix=data.frame(matrix(0,nrow=nrow(input),ncol=nrow(samples)+1))
colnames(sampleMatrix)=c(samples$V1)
output=cbind(output,sampleMatrix)

for (i in 1:nrow(input)) {
  output[i,1:4]=input[i,1:4]
  tmpArray=strsplit(as.character(input[i,7]),',')
  for (j in 1:length(tmpArray[[1]])) {
    val=strsplit(as.character(tmpArray[[1]][j]),':')
    output[i,val[[1]][1]]=val[[1]][2]
  }
}

output2=output[order(output$chr, output$start),]
write.table(output2,file='EHDnResults.combinedCounts.filtered.reformatted.bed',sep='\t',quote=F,row.names=F)

