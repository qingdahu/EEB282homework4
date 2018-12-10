# Homework 4

November 19, 2018
Qingda Hu

## Summarize partitions of a genome assembly

Get the Drosophila melanogaster genome and unzip it:

`wget ftp://ftp.flybase.net:21/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz

gzip -d dmel-all-chromosome-r6.24.fasta.gz`

### Calculate the following for all sequences ≤ 100kb and all sequences > 100kb: Total number of nucleotides, Total number of Ns, Total number of sequences

`module load jje/kent

module load jje/jjeutils/0.1a

gunzip dmel-all-chromosome-r6.24.fasta.gz

faSize <(bioawk -c fastx '{ if(length($seq) < 100001) { print ">"$name; print $seq }}'  dmel-all-chromosome-r6.24.fasta )

faSize <(bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}'  dmel-all-chromosome-r6.24.fasta)`


Less than 100kb: 6178042 bases (662593 N's 5515449 real) in 1863 sequences

Larger than 100kb: 137547960 bases (490385 N's 137057575 real) in 7 sequences 



### Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:   Sequence length distribution, Sequence GC% distribution, Cumulative genome size sorted from largest to smallest sequences

For whole genome sequence length distributions, I outputted all the length of the sequences and GC content with bioawk. I generated CDF with plotCDF tool. Then I imported the lengths and GC content into R for plotting.

`bioawk -c fastx '{ print length($seq) }' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length.txt

bioawk -c fastx '{ print gc($seq) }' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc.txt

module load R

plotCDF length.txt lengthCDF.png

R`

![lengthCDF.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF.png?raw=true)

After making sure the working directory is set to where length.txt is, load in the file generated with bioawk and plot the data. 

`getwd()

fullgenomeseqlengths <- read.table('length.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqlengths$V1) 

jpeg('rplot-full-length.jpg')

plot(d, main="whole genome sequence length distributions")

dev.off()`

![rplot-full-length.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length.jpg?raw=true)

`fullgenomeseqgc <- read.table('gc.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqgc $V1) 

jpeg('rplot-full-gc.jpg')

plot(d, main="whole genome sequence gc distributions")

dev.off()`

![rplot-full-gc.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-gc.jpg?raw=true)

Similarly for sequences ≤ 100kb:

`bioawk -c fastx '{ if(length($seq) < 100001) {print length($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length-lessthan100kb.txt
bioawk -c fastx '{ if(length($seq) < 100001) {print gc($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc-lessthan100kb.txt
plotCDF length-lessthan100kb.txt lengthCDF-lessthan100kb.png `

![lengthCDF-lessthan100kb.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF-lessthan100kb.png?raw=true)

`fullgenomeseqlengths <- read.table('length-lessthan100kb.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqlengths$V1) 


jpeg('rplot-full-length-lessthan100kb.jpg')

plot(d, main="length distributions for less than 100kb")

dev.off()`

![rplot-full-length-lessthan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length-lessthan100kb.jpg?raw=true)

`fullgenomeseqgc <- read.table('gc-lessthan100kb.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqgc $V1) 

jpeg('rplot-full-gc-lessthan100kb.jpg')

plot(d, main="gc distributions for less than 100kb")

dev.off()`

![rplot-full-gc-lessthan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-gc-lessthan100kb.jpg?raw=true)





And for all sequences > 100kb

`bioawk -c fastx '{ if(length($seq) > 100000) {print length($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length-morethan100kb.txt

bioawk -c fastx '{ if(length($seq) < 100001) {print gc($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc-morethan100kb.txt

plotCDF length-morethan100kb.txt lengthCDF-morethan100kb.png `

![lengthCDF-morethan100kb.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF-morethan100kb.png?raw=true)


`fullgenomeseqlengths <- read.table('length-morethan100kb.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqlengths$V1) 

jpeg('rplot-full-length-morethan100kb.jpg')

plot(d, main="length distributions for more than 100kb")

dev.off()`

![rplot-full-length-morethan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length-morethan100kb.jpg?raw=true)

`fullgenomeseqgc <- read.table('gc-morethan100kb.txt', sep = '\t',header = F) 

d <- density(fullgenomeseqgc $V1) 

jpeg('rplot-full-gc-morethan100kb.jpg')

plot(d, main="gc distributions for more than 100kb")

dev.off()`

![rplot-full-gc-morethan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-gc-morethan100kb.jpg?raw=true)













# Genome assembly

## Assemble a genome from MinION reads

```
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz


qrsh -q free128 -pe openmp 64

module load jje/jjeutils

minimap=$(which minimap)

miniasm=$(which miniasm)

basedir=/pub/jje/ee282/$USER

projname=nanopore_assembly

basedir=$basedir/$projname

raw=$basedir/$projname/data/raw

processed=$basedir/$projname/data/processed

figures=$basedir/$projname/output/figures

reports=$basedir/$projname/output/reports

createProject $projname $basedir

ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 64 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| fold -w 60 \
> $processed/unitigs.fa
```



### N50

obtain N50 using the code that we built during class:

`bioawk -c fastx '{ l+=length($2); print length($2)} END {print l}; ' $processed/unitigs.fa | sort -nr | awk 'NR ==1 {l=$1; } NR>1 {cumm+=$1/l; if(cumm>=0.5) {print $1; exit}}'| less -S`

Output is 4494246. Compare to Scaffold N50 25,286,936 from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4



### MUMmer plot





