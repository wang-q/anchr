# Assemble genomes from FDA-ARGOS data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from FDA-ARGOS data sets](#assemble-genomes-from-fda-argos-data-sets)
- [Francisella tularensis FDAARGOS_247](#francisella-tularensis-fdaargos_247)
  - [Ftul: reference](#ftul-reference)
  - [Ftul: download](#ftul-download)
  - [Ftul: template](#ftul-template)
  - [Ftul: run](#ftul-run)
- [Haemophilus influenzae FDAARGOS_199](#haemophilus-influenzae-fdaargos_199)
  - [Hinf: reference](#hinf-reference)
  - [Hinf: download](#hinf-download)
  - [Hinf: template](#hinf-template)
  - [Hinf: run](#hinf-run)
- [Campylobacter jejuni subsp. jejuni ATCC 700819](#campylobacter-jejuni-subsp-jejuni-atcc-700819)
  - [Cjej: reference](#cjej-reference)
  - [Cjej: download](#cjej-download)
  - [Cjej: template](#cjej-template)
  - [Cjej: run](#cjej-run)


* Rsync to hpcc

```shell script
for D in Ftul Hinf Cjej; do
    rsync -avP \
        ~/data/anchr/${D}/ \
        wangq@202.119.37.251:data/anchr/${D}
done

```


# Francisella tularensis FDAARGOS_247

## Ftul: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Ftul/1_genome
cd ~/data/anchr/Ftul/1_genome

cp ~/data/anchr/ref/Ftul/genome.fa .
cp ~/data/anchr/ref/Ftul/paralogs.fa .

```

## Ftul: download

```shell script
cd ~/data/anchr/Ftul

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX2105481,Ftul,HiSeq 2500 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name | srx        | platform | layout | ilength | srr        | spot     | base |
|:-----|:-----------|:---------|:-------|:--------|:-----------|:---------|:-----|
| Ftul | SRX2105481 | ILLUMINA | PAIRED | 560     | SRR4124773 | 10615135 | 2G   |


* Illumina

```shell script
cd ~/data/anchr/Ftul

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR4124773_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR4124773_2.fastq.gz R2.fq.gz

```

## Ftul: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Ftul

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1892775 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --sample 300 \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --redo \
    \
    --extend

```

## Ftul: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Ftul

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Ftul bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh
#rm -fr 4_down_sampling 6_down_sampling

```


# Haemophilus influenzae FDAARGOS_199

## Hinf: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Hinf/1_genome
cd ~/data/anchr/Hinf/1_genome

cp ~/data/anchr/ref/Hinf/genome.fa .
cp ~/data/anchr/ref/Hinf/paralogs.fa .

```

## Hinf: download

```shell script
cd ~/data/anchr/Hinf

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX2104758,Hinf,HiSeq 2500 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name | srx        | platform | layout | ilength | srr        | spot    | base  |
|:-----|:-----------|:---------|:-------|:--------|:-----------|:--------|:------|
| Hinf | SRX2104758 | ILLUMINA | PAIRED | 516     | SRR4123928 | 6115624 | 1.15G |


* Illumina

```shell script
cd ~/data/anchr/Hinf

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR4123928_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR4123928_2.fastq.gz R2.fq.gz

```

## Hinf: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Hinf

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1830138 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --redo \
    \
    --extend

```

## Hinf: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Hinf

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Hinf bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh
#rm -fr 4_down_sampling 6_down_sampling

```

# Campylobacter jejuni subsp. jejuni ATCC 700819

## Cjej: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Cjej/1_genome
cd ~/data/anchr/Cjej/1_genome

cp ~/data/anchr/ref/Cjej/genome.fa .
cp ~/data/anchr/ref/Cjej/paralogs.fa .

```

## Cjej: download

```shell script
cd ~/data/anchr/Cjej

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX2107012,Cjej,HiSeq 2500 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name | srx        | platform | layout | ilength | srr        | spot    | base  |
|:-----|:-----------|:---------|:-------|:--------|:-----------|:--------|:------|
| Cjej | SRX2107012 | ILLUMINA | PAIRED | 562     | SRR4125016 | 7696800 | 1.45G |


* Illumina

```shell script
cd ~/data/anchr/Cjej

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR4125016_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR4125016_2.fastq.gz R2.fq.gz

```

## Cjej: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Cjej

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1830138 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --redo \
    \
    --extend

```

## Cjej: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Cjej

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Cjej bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh
#rm -fr 4_down_sampling 6_down_sampling

```

