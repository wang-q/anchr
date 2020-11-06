# Assemble genomes from GAGE-B data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from GAGE-B data sets](#assemble-genomes-from-gage-b-data-sets)
- [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
  - [bcer_mi: reference](#bcer_mi-reference)
  - [bcer_mi: download](#bcer_mi-download)
  - [bcer_mi: template](#bcer_mi-template)
  - [bcer_mi: run](#bcer_mi-run)
- [*Mycobacterium abscessus* 6G-0125-R](#mycobacterium-abscessus-6g-0125-r)
  - [Mabs_mi: reference](#mabs_mi-reference)
- [*Rhodobacter sphaeroides* 2.4.1](#rhodobacter-sphaeroides-241)
  - [Rsph_mi: reference](#rsph_mi-reference)
- [*Vibrio cholerae* CP1032(5)](#vibrio-cholerae-cp10325)
  - [Vcho_mi: reference](#vcho_mi-reference)


# *Bacillus cereus* ATCC 10987

## bcer_mi: reference

* Reference genome

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

rm -fr 9_quast_competitor
quast --no-check --threads 16 \
    -R 1_genome/genome.fa \
    8_competitor/abyss_ctg.fasta \
    8_competitor/cabog_ctg.fasta \
    8_competitor/mira_ctg.fasta \
    8_competitor/msrca_ctg.fasta \
    8_competitor/sga_ctg.fasta \
    8_competitor/soap_ctg.fasta \
    8_competitor/spades_ctg.fasta \
    8_competitor/velvet_ctg.fasta \
    7_merge_anchors/anchor.merge.fasta \
    7_merge_anchors/others.non-contained.fasta \
    1_genome/paralogs.fa \
    --label "abyss,cabog,mira,msrca,sga,soap,spades,velvet,merge,others,paralogs" \
    -o 9_quast_competitor

```


# *Mycobacterium abscessus* 6G-0125-R

## Mabs_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Mabs_mi/1_genome
cd ~/data/anchr/Mabs_mi/1_genome

cp ~/data/anchr/ref/Mabs/genome.fa .
cp ~/data/anchr/ref/Mabs/paralogs.fa .

```

# *Rhodobacter sphaeroides* 2.4.1

## Rsph_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Rsph_mi/1_genome
cd ~/data/anchr/Rsph_mi/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .

```


# *Vibrio cholerae* CP1032(5)

* *Vibrio cholerae* O1 biovar El Tor str. N16961
  * Taxid: [243277](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=243277)
  * Assembly: [GCF_000006745.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000006745.1)
  * Proportion of paralogs (> 1000 bp): 0.0216
* *Vibrio cholerae* CP1032(5)
  * Assembly: [GCF_000279305.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000279305.1)

## Vcho_mi: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Vcho_mi/1_genome
cd ~/data/anchr/Vcho_mi/1_genome

cp ~/data/anchr/ref/Vcho/genome.fa .
cp ~/data/anchr/ref/Vcho/paralogs.fa .

```

