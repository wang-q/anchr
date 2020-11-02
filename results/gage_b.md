# Assemble genomes from GAGE-B data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from GAGE-B data sets](#assemble-genomes-from-gage-b-data-sets)
- [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
  - [bcer_mi: reference](#bcer_mi-reference)
  - [bcer_mi: download](#bcer_mi-download)
  - [bcer_mi: template](#bcer_mi-template)
  - [bcer_mi: run](#bcer_mi-run)


# *Bacillus cereus* ATCC 10987

* Taxonomy ID: [222523](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=222523)
* Assembly: [GCF_000008005.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000008005.1)
* Proportion of paralogs (> 1000 bp): 0.0344

B_cereus_MiSeq


## bcer_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/008/005/GCF_000008005.1_ASM800v1/ \
    Bcer/

```

```shell script
mkdir -p ~/data/anchr/bcer_mi/1_genome
cd ~/data/anchr/bcer_mi/1_genome

find ~/data/anchr/ref/Bcer/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    xargs gzip -dcf |
    faops filter -N -s stdin genome.fa

cat ~/data/anchr/paralogs/model/Results/Bcer/Bcer.multi.fas |
    faops filter -N -d stdin stdout \
    > paralogs.fa

```

## bcer_mi: download

* Illumina

```shell script
cd ~/data/anchr/bcer_mi

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/datasets/B_cereus_MiSeq.tar.gz

# NOT gzipped tar
tar xvf B_cereus_MiSeq.tar.gz raw/frag_1__cov100x.fastq
tar xvf B_cereus_MiSeq.tar.gz raw/frag_2__cov100x.fastq

cat raw/frag_1__cov100x.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/frag_2__cov100x.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* GAGE-B assemblies

```shell script
cd ${HOME}/data/anchr/bcer_mi

mkdir -p 8_competitor
cd 8_competitor

aria2c -x 9 -s 3 -c http://ccb.jhu.edu/gage_b/genomeAssemblies/B_cereus_MiSeq.tar.gz

tar xvfz B_cereus_MiSeq.tar.gz abyss_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz cabog_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz mira_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz msrca_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz sga_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz soap_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz spades_ctg.fasta
tar xvfz B_cereus_MiSeq.tar.gz velvet_ctg.fasta

```

## bcer_mi: template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=bcer_mi

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 5432652  \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --tile" \
    --qual "20 25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    --cov "40 50 60 all" \
    --statp 2 \
    --redo \
    \
    --extend

```

## bcer_mi: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=bcer_mi

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

#bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

```
