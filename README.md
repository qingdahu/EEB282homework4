# Homework 4

November 19, 2018
Qingda Hu

## Summarize partitions of a genome assembly

Get the Drosophila melanogaster genome and unzip it:

```
wget ftp://ftp.flybase.net:21/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz
gzip -d dmel-all-chromosome-r6.24.fasta.gz
```

### Calculate the following for all sequences ≤ 100kb and all sequences > 100kb: Total number of nucleotides, Total number of Ns, Total number of sequences

```
module load jje/kent
module load jje/jjeutils/0.1a
gunzip dmel-all-chromosome-r6.24.fasta.gz
faSize <(bioawk -c fastx '{ if(length($seq) < 100001) { print ">"$name; print $seq }}'  dmel-all-chromosome-r6.24.fasta )
faSize <(bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}'  dmel-all-chromosome-r6.24.fasta)
```


Less than 100kb: 6178042 bases (662593 N's 5515449 real) in 1863 sequences

Larger than 100kb: 137547960 bases (490385 N's 137057575 real) in 7 sequences 



### Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:   Sequence length distribution, Sequence GC% distribution, Cumulative genome size sorted from largest to smallest sequences

For whole genome sequence length distributions, I outputted all the length of the sequences and GC content with bioawk. I generated CDF with plotCDF tool. Then I imported the lengths and GC content into R for plotting.

```
bioawk -c fastx '{ print length($seq) }' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length.txt
bioawk -c fastx '{ print gc($seq) }' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc.txt
module load R
plotCDF length.txt lengthCDF.png
R
```

![lengthCDF.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF.png?raw=true)

After making sure the working directory is set to where length.txt is, load in the file generated with bioawk and plot the data. 

```
getwd()
fullgenomeseqlengths <- read.table('length.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqlengths$V1) 
jpeg('rplot-full-length.jpg')
plot(d, main="whole genome sequence length distributions")
dev.off()
```

![rplot-full-length.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length.jpg?raw=true)

```
fullgenomeseqgc <- read.table('gc.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqgc $V1) 
jpeg('rplot-full-gc.jpg')
plot(d, main="whole genome sequence gc distributions")
dev.off()
```

![rplot-full-gc.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-gc.jpg?raw=true)

Similarly for sequences ≤ 100kb:

```
bioawk -c fastx '{ if(length($seq) < 100001) {print length($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length-lessthan100kb.txt
bioawk -c fastx '{ if(length($seq) < 100001) {print gc($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc-lessthan100kb.txt
plotCDF length-lessthan100kb.txt lengthCDF-lessthan100kb.png 
```

![lengthCDF-lessthan100kb.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF-lessthan100kb.png?raw=true)

```
fullgenomeseqlengths <- read.table('length-lessthan100kb.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqlengths$V1) 
jpeg('rplot-full-length-lessthan100kb.jpg')
plot(d, main="length distributions for less than 100kb")
dev.off()
```

![rplot-full-length-lessthan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length-lessthan100kb.jpg?raw=true)

```
fullgenomeseqgc <- read.table('gc-lessthan100kb.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqgc $V1) 
jpeg('rplot-full-gc-lessthan100kb.jpg')
plot(d, main="gc distributions for less than 100kb")
dev.off()
```

![rplot-full-gc-lessthan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-gc-lessthan100kb.jpg?raw=true)





And for all sequences > 100kb

```
bioawk -c fastx '{ if(length($seq) > 100000) {print length($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > length-morethan100kb.txt
bioawk -c fastx '{ if(length($seq) < 100001) {print gc($seq) }}' dmel-all-chromosome-r6.24.fasta |sort -rnk 1,1 > gc-morethan100kb.txt
plotCDF length-morethan100kb.txt lengthCDF-morethan100kb.png 
```

![lengthCDF-morethan100kb.png](https://github.com/qingdahu/EEB282homework4/blob/master/lengthCDF-morethan100kb.png?raw=true)


```
fullgenomeseqlengths <- read.table('length-morethan100kb.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqlengths$V1) 
jpeg('rplot-full-length-morethan100kb.jpg')
plot(d, main="length distributions for more than 100kb")
dev.off()
```

![rplot-full-length-morethan100kb.jpg](https://github.com/qingdahu/EEB282homework4/blob/master/rplot-full-length-morethan100kb.jpg?raw=true)

```
fullgenomeseqgc <- read.table('gc-morethan100kb.txt', sep = '\t',header = F) 
d <- density(fullgenomeseqgc $V1) 
jpeg('rplot-full-gc-morethan100kb.jpg')
plot(d, main="gc distributions for more than 100kb")
dev.off()
```

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

```
bioawk -c fastx '{ l+=length($2); print length($2)} END {print l}; ' $processed/unitigs.fa | sort -nr | awk 'NR ==1 {l=$1; } NR>1 {cumm+=$1/l; if(cumm>=0.5) {print $1; exit}}'| less -S
```

Output is 4494246. Compare to Scaffold N50 25,286,936 from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4




### MUMmer plot


I started by using faSplitByN to break the scaffold into contigs. 

```
module load jje/jjeutils
faSplitByN dmel-all-chromosome-r6.24.fasta dmel-contig.fa 10
module unload jje/jjeutils
```

For easy of coding I moved the unitig.fa and dmel-contig.fa to the same folder. 


```
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot/4.6.0


REF="dmel-contig.fa"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)


nucmer -l 100 -c 150 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
  -R ${REF} -Q ${QRY} --png

```

Here is the result.

![MUMmer plot](https://github.com/qingdahu/EEB282homework4/blob/master/flybase_unitigs.png?raw=true)


I also accidentally ran it for the scaffold. See below. 

![MUMmer with scaffold](https://github.com/qingdahu/EEB282homework4/blob/master/flybase_unitigs.scaffold.png?raw=true)


### Contiguity plot


------------------------
fifo example 
#!/usr/bin/env bash
module load perl
module load jje/jjeutils/0.1a
module load rstudio/0.99.9.9

r6url="https://drive.google.com/uc?export=download&id=1sSStvljlakBAjGhON2FJe_tCj6JLE1rZ"
trusequrl="https://drive.google.com/uc?export=download&id=0B89qjFzcZ81rdXdadmhRRkxwamc"

scratchdir=/pub/jje/ee282/$USER
mkdir -p $scratchdir; cd $scratchdir

mkfifo {r6scaff,r6ctg,truseq}_fifo

wget -O - -q $r6url \
| tee >( \
  bioawk -c fastx ' { print length($seq) } ' \
  | sort -rn \
  | awk ' BEGIN { print "Assembly\tLength\nFB_Scaff\t0" } { print "FB_Scaff\t" $1 } ' \
  > r6scaff_fifo & ) \
| faSplitByN /dev/stdin /dev/stdout 10 \
| bioawk -c fastx ' { print length($seq) } ' \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nFB_Ctg\t0" } { print "FB_Ctg\t" $1 } ' \
> r6ctg_fifo &

wget -O - -q $trusequrl \
| bioawk -c fastx ' { print length($seq) } ' \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nTruSeq_Ctg\t0" } { print "TruSeq_Ctg\t" $1 } ' \
> truseq_fifo &

plotCDF2 {r6scaff,r6ctg,truseq}_fifo /dev/stdout \
| tee r6_v_truseq.png \
| display 

rm {r6scaff,r6ctg,truseq}_fifo





### BUSCO score


Here is the script for busco run with unitigs. 

```
#!/bin/bash
#
#$ -N busco_qh
#$ -q free128,abio128,bio,abio
#$ -pe openmp 32
#$ -R Y
### -m beas
### -M qingdah@uci.edu

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="unitigs.fa"
MYEXT=".fa"

#my busco run
BUSCO.py -c ${NSLOTS} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}
```

The results were:

```
# BUSCO version is: 2.0 
# The lineage dataset is: diptera_odb9 (Creation date: 2016-10-21, number of species: 25, number of BUSCOs: 2799)
# To reproduce this run: python /pub/jje/ee282/bin/busco/BUSCO.py -i unitigs.fa -o unitigs_diptera_odb9 -l /pub/jje/ee282/bin/busco/lineages/diptera_odb9/ -m genome -c 32 -sp fly
#
# Summarized benchmarking in BUSCO notation for file unitigs.fa
# BUSCO was run in mode: genome

        C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799

        13      Complete BUSCOs (C)
        13      Complete and single-copy BUSCOs (S)
        0       Complete and duplicated BUSCOs (D)
        32      Fragmented BUSCOs (F)
        2754    Missing BUSCOs (M)
        2799    Total BUSCO groups searched
```

I then ran busco with the flybase contig. 

```
#!/bin/bash
#
#$ -N busco_qh
#$ -q free128,abio128,bio,abio
#$ -pe openmp 32
#$ -R Y
### -m beas
### -M qingdah@uci.edu

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="dmel-contig.fa"
MYEXT=".fa"

#my busco run
BUSCO.py -c ${NSLOTS} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}
```

The results for 
















