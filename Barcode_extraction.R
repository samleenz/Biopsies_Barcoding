library(edgeR)

fastq1 = "Undetermined_S0_R1.fastq"  # Name of read 1 fastq file
fastq2 = "Undetermined_S0_R2.fastq"  # Name of read 2 fastq file
sampleanno = "TPH09_TrueIndexes.txt"  # Name of tab delimited file containing sample annotation information
outputfile = "TPHseq09"  # Give your experiment a name - all results files will begin with this prefix

# Count sequences in fastq file
counts = processAmplicons(readfile=fastq1, readfile2=fastq2, barcodefile=sampleanno, hairpinfile="LongBarcodes.txt", barcodeStart=1, 
                          barcodeEnd=8, barcodeStartRev=1, barcodeEndRev=8, hairpinStart=25, hairpinEnd=39, verbose=TRUE)

# Write counts out to file
raw = data.frame("ID"=rownames(counts$counts), counts$counts)
colnames(raw)[-1] =  as.character(counts$samples[,5])
write.table(raw, file=paste(outputfile, "Counts.txt", sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)

