# Assemble genomes of model organisms by `anchr`

<!-- TOC -->
* [Assemble genomes of model organisms by `anchr`](#assemble-genomes-of-model-organisms-by-anchr)
  * [*Mycoplasma genitalium* G37](#mycoplasma-genitalium-g37)
    * [g37: reference](#g37-reference)
    * [g37: download](#g37-download)
    * [g37: template](#g37-template)
    * [g37: run](#g37-run)
  * [*E. coli* str. K-12 substr. MG1655](#e-coli-str-k-12-substr-mg1655)
    * [mg1655: reference](#mg1655-reference)
    * [mg1655: download](#mg1655-download)
    * [mg1655: template](#mg1655-template)
    * [mg1655: run](#mg1655-run)
  * [*E. coli* str. K-12 substr. DH5alpha](#e-coli-str-k-12-substr-dh5alpha)
    * [dh5alpha: reference](#dh5alpha-reference)
    * [dh5alpha: download](#dh5alpha-download)
    * [dh5alpha: template](#dh5alpha-template)
    * [dh5alpha: run](#dh5alpha-run)
  * [*Bacillus cereus* ATCC 10987](#bacillus-cereus-atcc-10987)
    * [bcer: reference](#bcer-reference)
    * [bcer: download](#bcer-download)
    * [bcer: template](#bcer-template)
    * [bcer: run](#bcer-run)
  * [*Rhodobacter sphaeroides* 2.4.1](#rhodobacter-sphaeroides-241)
    * [rsph: reference](#rsph-reference)
    * [rsph: download](#rsph-download)
    * [rsph: template](#rsph-template)
    * [rsph: run](#rsph-run)
<!-- TOC -->

> **Note**: `asm multik` 门禁（multik 单次调用 / 性能优化 / end-multiplicity
> 门控全链门禁，以及 bcer/rsph 现代流程门禁待跑）已迁移至
> `scripts/asm_gate.md` §gate history。本文件只保留数据集定义、download、
> template、run 与质量基线表格。引用门禁基线请查 `scripts/asm_gate.md`。

## *Mycoplasma genitalium* G37

### g37: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/g37/1_genome
cd ~/data/anchr/g37/1_genome

cp ~/data/anchr/ref/g37/genome.fa .
cp ~/data/anchr/ref/g37/paralogs.fa .

```

### g37: download

```shell
mkdir -p ~/data/anchr/g37/ena
cd ~/data/anchr/g37/ena

cat << EOF > source.csv
ERX452667,G37,MiSeq
EOF

anchr ena meta source.csv > ena_info.json
anchr ena manifest ena_info.json

tva to md ena_info.tsv --fmt

aria2c -j 4 -x 4 -s 2 -c --file-allocation=none -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name | srx       | platform | layout | ilength | srr       |   spots | bases  |
|------|-----------|----------|--------|--------:|-----------|--------:|--------|
| G37  | ERX452667 | ILLUMINA | PAIRED |     447 | ERR486835 | 680,644 | 97.37M |

* Illumina

```shell
cd ~/data/anchr/g37

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/ERR486835_1.fastq.gz R1.fq.gz
ln -s ../ena/ERR486835_2.fastq.gz R2.fq.gz

```

### g37: template

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=g37

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 580076 \
    --parallel 16 \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### g37: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=g37

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group    |  Mean | Median | STDev | Pairs% | Orientation |
| -------- | ----: | -----: | ----: | -----: | ----------- |
| R.genome | 466.5 |    451 | 127.8 | 97.74% | FR          |
| R.contig | 462.1 |    447 | 126.1 | 92.63% | FR          |

Table: statInsertSize

| K    | property              |        min |        max |
| ---- | --------------------- | ---------: | ---------: |
| R.21 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 577,874 bp |
|      | Genome Repeat Length  |   3,639 bp |   3,645 bp |
|      | Genome Unique Length  | 573,825 bp | 574,639 bp |
|      | Model Fit             |   92.5827% |   93.3366% |
|      | Read Error Rate       |            |  0.137189% |
|      | Kmer Cov              |            |      148.8 |
| R.51 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 578,025 bp |
|      | Genome Repeat Length  |            |       0 bp |
|      | Genome Unique Length  |            | 578,025 bp |
|      | Model Fit             |   95.6385% |   95.7255% |
|      | Read Error Rate       |            | 0.0942288% |
|      | Kmer Cov              |            |      112.2 |
| R.81 | Homozygous (a)        |            |       100% |
|      | Genome Haploid Length |            | 578,301 bp |
|      | Genome Repeat Length  |            |       0 bp |
|      | Genome Unique Length  |            | 578,301 bp |
|      | Model Fit             |   97.0931% |   97.1094% |
|      | Read Error Rate       |            | 0.0778153% |
|      | Kmer Cov              |            |       77.3 |

Table: statFastK

| chr       | chrLength |  size | coverage |
| --------- | --------: | ----: | -------: |
| NC_000908 |    580076 |  6905 |   0.0119 |
| all       |    580076 |  6905 |   0.0119 |

Table: statRepetitive

| Name       |    N50 |     Sum |      # |
| ---------- | -----: | ------: | -----: |
| genome     | 580076 | 580.08K |      1 |
| paralogs   |   1567 |  11.53K |      8 |
| repetitive |    184 |   6.91K |     41 |
| Illumina.R |    150 |  102.1M | 680644 |
| trim.R     |    150 | 101.44M | 676644 |
| Q0L0       |    150 | 101.44M | 676644 |
| Q25L60     |    150 |  98.23M | 657747 |
| Q30L60     |    150 |  95.03M | 638160 |

Table: statReads

| Name     |  N50 |     Sum |      # |
| -------- | ---: | ------: | -----: |
| clumpify |  150 | 102.01M | 680076 |
| highpass |  150 | 101.54M | 676954 |
| trim     |  150 | 101.44M | 676644 |
| filter   |  150 | 101.44M | 676644 |
| R1       |  150 |  50.74M | 338322 |
| R2       |  150 |   50.7M | 338322 |
| Rs       |    0 |       0 |      0 |

Table: statTrimReads

```text
#R.trim
#Matched	589	0.08701%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

| Name          |  N50 |     Sum |      # |
| ------------- | ---: | ------: | -----: |
| clumped       |  150 | 101.44M | 676644 |
| ecco          |  150 | 101.44M | 676644 |
| ecct          |  150 |  97.96M | 653378 |
| extended      |  190 | 123.74M | 653378 |
| merged.raw    |  442 |  77.19M | 181762 |
| unmerged.raw  |  190 |  54.74M | 289854 |
| unmerged.trim |  190 |   54.7M | 289808 |
| M1            |  442 |  77.19M | 181762 |
| U1            |  190 |  27.46M | 144904 |
| U2            |  190 |  27.24M | 144904 |
| Us            |    0 |       0 |      0 |
| M.cor         |  362 | 132.08M | 653332 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
| ------------------ | ----: | -----: | ----: | -----: |
| M.ihist.merge1.txt | 258.3 |    264 |  24.7 |  4.23% |
| M.ihist.merge.txt  | 424.6 |    430 |  64.8 | 55.68% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% |   RealG | RunTime |
| -------- | ----: | -----: | -------: | ------: | ------: |
| Q0L0.R   | 174.9 |  160.8 |    8.07% | 580.08K | 0:00:10 |
| Q25L60.R | 169.4 |  160.5 |    5.26% | 580.08K | 0:00:09 |
| Q30L60.R | 163.9 |  157.3 |    4.01% | 580.08K | 0:00:08 |

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.972 |     37567 | 565.13K |   21 |       874 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.972 |     34592 | 578.44K |   29 |      1404 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P002   |   40.0 |  0.972 |     38759 | 576.46K |   22 |       895 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.966 |     14761 | 567.41K |   59 |      1193 |     79 |    9 |  17.3 | 212.0 |
| Q0L0X80P001   |   80.0 |  0.966 |     21883 | 566.02K |   54 |      1237 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X40P000 |   40.0 |  0.972 |     34622 | 569.35K |   25 |       973 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P001 |   40.0 |  0.972 |     38774 | 565.36K |   24 |      1865 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P002 |   40.0 |  0.971 |     37594 | 564.99K |   28 |      1122 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X80P000 |   80.0 |  0.964 |     19064 |  570.5K |   58 |      1662 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X80P001 |   80.0 |  0.966 |     17836 | 566.39K |   53 |      1000 |     79 |    9 |  17.3 | 212.0 |
| Q30L60X40P000 |   40.0 |  0.972 |     48771 | 565.12K |   23 |      1673 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P001 |   40.0 |  0.971 |     31848 | 570.28K |   29 |      1810 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P002 |   40.0 |  0.973 |     38030 | 565.11K |   20 |       995 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X80P000 |   80.0 |  0.968 |     17586 | 567.54K |   52 |      1403 |     79 |    9 |  17.3 | 212.0 |

Table: statUnitigsMultik.md

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.972 |     37567 | 565.13K |   21 |       874 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.972 |     34592 | 578.44K |   29 |      1404 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P002   |   40.0 |  0.972 |     38759 | 576.46K |   22 |       895 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.966 |     14761 | 567.41K |   59 |      1193 |     79 |    9 |  17.3 | 212.0 |
| Q0L0X80P001   |   80.0 |  0.966 |     21883 | 566.02K |   54 |      1237 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X40P000 |   40.0 |  0.972 |     34622 | 569.35K |   25 |       973 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P001 |   40.0 |  0.972 |     38774 | 565.36K |   24 |      1865 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X40P002 |   40.0 |  0.971 |     37594 | 564.99K |   28 |      1122 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X80P000 |   80.0 |  0.964 |     19064 |  570.5K |   58 |      1662 |     79 |    9 |  17.3 | 212.0 |
| Q25L60X80P001 |   80.0 |  0.966 |     17836 | 566.39K |   53 |      1000 |     79 |    9 |  17.3 | 212.0 |
| Q30L60X40P000 |   40.0 |  0.972 |     48771 | 565.12K |   23 |      1673 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P001 |   40.0 |  0.971 |     31848 | 570.28K |   29 |      1810 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X40P002 |   40.0 |  0.973 |     38030 | 565.11K |   20 |       995 |     40 |    5 |   8.3 | 110.0 |
| Q30L60X80P000 |   80.0 |  0.968 |     17586 | 567.54K |   52 |      1403 |     79 |    9 |  17.3 | 212.0 |

Table: statUnitigsUnitig.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.709 |     81602 | 582.23K |   25 |      6561 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.709 |     77931 | 581.81K |   22 |      7925 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.710 |     81650 | 580.96K |   18 |      5780 |     39 |    5 |   8.0 | 108.0 |
| MRX80P000 |   80.0 |  0.708 |     40795 | 580.65K |   30 |      7209 |     78 |    9 |  17.0 | 210.0 |
| MRX80P001 |   80.0 |  0.709 |     47580 | 580.05K |   25 |      6116 |     79 |    9 |  17.3 | 212.0 |

Table: statMRUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.709 |     81602 | 603.76K |   24 |      6175 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.709 |     81221 | 582.83K |   22 |      7708 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.710 |    119297 |    581K |   16 |      5291 |     39 |    5 |   8.0 | 108.0 |
| MRX80P000 |   80.0 |  0.708 |     40795 | 581.69K |   31 |      7684 |     78 |    9 |  17.0 | 210.0 |
| MRX80P001 |   80.0 |  0.709 |     47580 | 581.25K |   25 |      5784 |     79 |    9 |  17.3 | 212.0 |

Table: statMRUnitigsUnitig.md

| Name                      | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.709 |    121443 | 580.83K |    9 |       322 |    225 |   22 |  53.0 | 582.0 |
| 7_merge_mr_unitigs_multik |  0.709 |    121369 | 586.24K |   17 |      1052 |    224 |   23 |  51.7 | 586.0 |
| 7_merge_mr_unitigs_unitig |  0.708 |     93018 | 580.51K |   11 |       547 |    225 |   22 |  53.0 | 582.0 |
| 7_merge_unitigs_multik    |  0.692 |     48860 | 567.51K |   17 |      1444 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_unitig    |  0.692 |     40559 | 591.06K |   21 |      1482 |    227 |   23 |  52.7 | 592.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.988 |     83431 | 573.62K |    9 |      7382 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_spades  |  0.710 |    580031 | 580.03K |    1 |       601 |    225 |   22 |  53.0 | 582.0 |
| 8_megahit    |  0.978 |     54983 | 564.83K |   19 |     13431 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_megahit |  0.710 |    319000 | 579.65K |    8 |     15840 |    224 |   23 |  51.7 | 586.0 |

Table: statOtherAnchors.md

| Name                    |    N50 |     Sum |    # |
| ----------------------- | -----: | ------: | ---: |
| Genome                  | 580076 | 580.08K |    1 |
| Paralogs                |   1567 |  11.53K |    8 |
| repetitive              |    184 |   6.91K |   41 |
| 7_merge_anchors.anchors | 121443 | 580.83K |    9 |
| spades.contig           | 163847 |    581K |   38 |
| spades.scaffold         | 163847 |    581K |   38 |
| mr_spades.contig        | 580506 | 580.63K |    2 |
| mr_spades.scaffold      | 580506 | 580.63K |    2 |
| megahit.contig          |  55053 | 578.27K |   48 |
| mr_megahit.contig       | 319186 | 595.49K |   33 |

Table: statFinal

* Assembly quality by QUAST

| Assembly         | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| ---------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik     |        17 |  130886 | 567510 |  48860 |     0 | 97.338 | 1.002 |   0.00 |  112.24 |      31.11 |
| merge_mr_multik  |        17 |  187458 | 586244 | 121369 |     0 | 98.666 | 1.005 |   0.00 |  241.36 |      70.08 |
| merge_unitig     |        21 |  106959 | 591059 |  40559 |     0 | 97.350 | 1.044 |   0.00 |  111.48 |      29.52 |
| merge_mr_unitig  |        11 |  235838 | 580508 |  93018 |     0 | 99.384 | 1.000 |   0.00 |  267.71 |      80.11 |
| merge_anchors    |         9 |  236389 | 580827 | 121443 |     0 | 99.464 | 1.000 |   0.00 |  272.94 |      81.45 |
| spades           |        35 |  236302 | 580752 | 163847 |     1 | 99.056 | 1.002 |   0.00 |  222.74 |      66.49 |
| mr_spades        |         2 |  580506 | 580632 | 580506 |     0 | 100.000 | 1.001 |   0.00 |  300.77 |      91.64 |
| megahit          |        48 |  179893 | 578265 |  55053 |     2 | 97.731 | 1.002 |   0.00 |  104.60 |      23.95 |
| mr_megahit       |        33 |  319186 | 595492 | 319186 |     4 | 99.714 | 1.022 |   0.00 |  304.04 |      88.44 |

Table: statQuast

## *E. coli* str. K-12 substr. MG1655

### mg1655: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/mg1655/1_genome
cd ~/data/anchr/mg1655/1_genome

cp ~/data/anchr/ref/mg1655/genome.fa .
cp ~/data/anchr/ref/mg1655/paralogs.fa .

```

### mg1655: download

```shell
mkdir -p ~/data/anchr/mg1655/ena
cd ~/data/anchr/mg1655/ena

aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz

```

* Illumina

```shell
cd ~/data/anchr/mg1655

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz R1.fq.gz
ln -s ../ena/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz R2.fq.gz

```

### mg1655: template

* template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4641652 \
    --parallel 16 \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --merge \
    \
    --bwa "Q25L60" \
    \
    --cov "40 80" \
    --unitigger "multik unitig" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### mg1655: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=mg1655

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group    |  Mean | Median | STDev | Pairs% | Orientation |
| -------- | ----: | -----: | ----: | -----: | ----------- |
| R.genome | 298.1 |    298 |  17.7 | 99.40% | FR          |
| R.contig | 297.8 |    297 |  17.7 | 97.75% | FR          |

Table: statInsertSize

| K    | property              |          min |          max |
| ---- | --------------------- | -----------: | -----------: |
| R.21 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,477,162 bp |
|      | Genome Repeat Length  |   135,854 bp |   135,942 bp |
|      | Genome Unique Length  | 4,339,855 bp | 4,342,674 bp |
|      | Model Fit             |     95.3903% |     97.3934% |
|      | Read Error Rate       |              |    0.531924% |
|      | Kmer Cov              |              |        299.7 |
| R.51 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,384,614 bp |
|      | Genome Repeat Length  |    90,815 bp |    90,901 bp |
|      | Genome Unique Length  | 4,291,722 bp | 4,295,793 bp |
|      | Model Fit             |     95.9309% |     97.6265% |
|      | Read Error Rate       |              |    0.326151% |
|      | Kmer Cov              |              |        223.4 |
| R.81 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,303,893 bp |
|      | Genome Repeat Length  |    59,727 bp |    59,810 bp |
|      | Genome Unique Length  | 4,241,157 bp | 4,247,095 bp |
|      | Model Fit             |     95.9005% |     97.4791% |
|      | Read Error Rate       |              |    0.272518% |
|      | Kmer Cov              |              |        151.5 |

Table: statFastK

| chr       | chrLength |   size | coverage |
| --------- | --------: | -----: | -------: |
| NC_000913 |   4641652 | 124422 |   0.0268 |
| all       |   4641652 | 124422 |   0.0268 |

Table: statRepetitive

| Name       |     N50 |     Sum |        # |
| ---------- | ------: | ------: | -------: |
| genome     | 4641652 |   4.64M |        1 |
| paralogs   |    1737 |    193K |      112 |
| repetitive |    1265 | 124.42K |      155 |
| Illumina.R |     151 |   1.73G | 11458940 |
| trim.R     |     149 |   1.46G | 10636412 |
| Q0L0       |     149 |   1.46G | 10636412 |
| Q25L60     |     148 |   1.35G | 10185038 |
| Q30L60     |     128 |   1.12G |  9506942 |

Table: statReads

| Name     |  N50 |     Sum |        # |
| -------- | ---: | ------: | -------: |
| clumpify |  151 |   1.73G | 11439000 |
| highpass |  151 |    1.7G | 11272376 |
| trim     |  149 |   1.46G | 10636412 |
| filter   |  149 |   1.46G | 10636412 |
| R1       |  150 | 753.89M |  5318206 |
| R2       |  144 | 706.08M |  5318206 |
| Rs       |    0 |       0 |        0 |

Table: statTrimReads

```text
#R.trim
#Matched	17107	0.15176%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

| Name          |  N50 |    Sum |        # |
| ------------- | ---: | -----: | -------: |
| clumped       |  149 |  1.46G | 10635230 |
| ecco          |  149 |  1.46G | 10635230 |
| ecct          |  149 |  1.45G | 10566976 |
| extended      |  189 |  1.87G | 10566976 |
| merged.raw    |  339 |  1.76G |  5210725 |
| unmerged.raw  |  172 | 21.12M |   145526 |
| unmerged.trim |  170 | 20.22M |   140094 |
| M1            |  339 |  1.76G |  5210571 |
| U1            |  175 | 10.62M |    70047 |
| U2            |  165 |   9.6M |    70047 |
| Us            |    0 |      0 |        0 |
| M.cor         |  338 |  1.78G | 10561236 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
| ------------------ | ----: | -----: | ----: | -----: |
| M.ihist.merge1.txt | 271.3 |    277 |  24.3 | 10.45% |
| M.ihist.merge.txt  | 337.6 |    338 |  19.3 | 98.77% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real | RunTime |
| -------- | ----: | -----: | -------: | ---: | ----: | ----: | -------: | ------: |
| Q0L0.R   | 314.5 |  291.7 |    7.27% | "24" | 4.64M | 4.47M |     0.96 | 0:02:46 |
| Q25L60.R | 289.9 |  278.4 |    3.95% | "24" | 4.64M | 4.42M |     0.95 | 0:02:31 |
| Q30L60.R | 242.2 |  237.1 |    2.12% | "24" | 4.64M | 4.37M |     0.94 | 0:02:05 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.976 |     15337 | 4.81M |  520 |     32162 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P001   |   40.0 |  0.976 |     14493 | 4.81M |  509 |     28910 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P002   |   40.0 |  0.976 |     14481 | 4.77M |  510 |     30287 |     40 |    7 |   6.3 | 122.0 |
| Q0L0X80P000   |   80.0 |  0.976 |     16422 | 4.78M |  470 |     29446 |     79 |   11 |  15.3 | 224.0 |
| Q0L0X80P001   |   80.0 |  0.976 |     15544 | 4.79M |  482 |     27038 |     79 |   11 |  15.3 | 224.0 |
| Q0L0X80P002   |   80.0 |  0.977 |     16106 | 4.73M |  465 |     28339 |     79 |   11 |  15.3 | 224.0 |
| Q25L60X40P000 |   40.0 |  0.978 |     10267 | 4.81M |  658 |     33214 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X40P001 |   40.0 |  0.977 |     11022 | 4.79M |  633 |     32558 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X40P002 |   40.0 |  0.978 |     11976 | 4.79M |  622 |     31944 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X80P000 |   80.0 |  0.979 |     14114 | 4.74M |  517 |     28806 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P001 |   80.0 |  0.979 |     14209 | 4.75M |  509 |     25773 |     80 |   12 |  14.7 | 232.0 |
| Q25L60X80P002 |   80.0 |  0.979 |     14093 | 4.74M |  500 |     26831 |     80 |   12 |  14.7 | 232.0 |
| Q30L60X40P000 |   40.0 |  0.970 |      7025 | 4.77M |  907 |     37334 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X40P001 |   40.0 |  0.972 |      7036 | 4.72M |  898 |     36172 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X40P002 |   40.0 |  0.969 |      7280 | 4.73M |  889 |     35073 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X80P000 |   80.0 |  0.979 |      9458 | 4.79M |  700 |     26388 |     79 |   14 |  12.3 | 242.0 |
| Q30L60X80P001 |   80.0 |  0.979 |      9694 | 4.82M |  703 |     25628 |     79 |   14 |  12.3 | 242.0 |

Table: statUnitigsMultik.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.976 |     14647 | 4.81M |  529 |     32383 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P001   |   40.0 |  0.976 |     14389 | 4.84M |  525 |     29268 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P002   |   40.0 |  0.976 |     14029 | 4.81M |  524 |     30647 |     40 |    7 |   6.3 | 122.0 |
| Q0L0X80P000   |   80.0 |  0.976 |     16168 | 4.78M |  475 |     29502 |     79 |   11 |  15.3 | 224.0 |
| Q0L0X80P001   |   80.0 |  0.976 |     15533 | 4.79M |  485 |     27099 |     79 |   11 |  15.3 | 224.0 |
| Q0L0X80P002   |   80.0 |  0.977 |     15975 | 4.74M |  472 |     28454 |     79 |   11 |  15.3 | 224.0 |
| Q25L60X40P000 |   40.0 |  0.978 |     10077 | 4.84M |  677 |     34189 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X40P001 |   40.0 |  0.977 |     10757 | 4.82M |  647 |     32940 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X40P002 |   40.0 |  0.978 |     11761 | 4.81M |  635 |     32296 |     40 |    7 |   6.3 | 122.0 |
| Q25L60X80P000 |   80.0 |  0.979 |     13956 | 4.78M |  531 |     28928 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P001 |   80.0 |  0.979 |     14103 | 4.75M |  515 |     25831 |     80 |   12 |  14.7 | 232.0 |
| Q25L60X80P002 |   80.0 |  0.979 |     13823 | 4.74M |  509 |     26717 |     79 |   12 |  14.3 | 230.0 |
| Q30L60X40P000 |   40.0 |  0.970 |      6947 | 4.81M |  922 |     37681 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X40P001 |   40.0 |  0.972 |      6869 | 4.76M |  916 |     36580 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X40P002 |   40.0 |  0.969 |      7092 | 4.76M |  903 |     35290 |     40 |    8 |   5.3 | 128.0 |
| Q30L60X80P000 |   80.0 |  0.979 |      9162 | 4.84M |  722 |     27156 |     79 |   14 |  12.3 | 242.0 |
| Q30L60X80P001 |   80.0 |  0.979 |      9219 | 4.84M |  716 |     26352 |     79 |   14 |  12.3 | 242.0 |

Table: statUnitigsUnitig.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.498 |     80369 | 4.61M |  115 |     52747 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.498 |    112151 | 4.57M |  108 |     51317 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.498 |    102005 | 4.59M |  112 |     54891 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.498 |     95584 | 4.56M |  110 |     50530 |     79 |   10 |  16.3 | 218.0 |
| MRX80P001 |   80.0 |  0.498 |     88489 | 4.56M |  107 |     52335 |     79 |   10 |  16.3 | 218.0 |
| MRX80P002 |   80.0 |  0.498 |     77856 | 4.61M |  118 |     54033 |     79 |   10 |  16.3 | 218.0 |

Table: statMRUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.497 |     61058 | 4.62M |  147 |     54805 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.497 |     80391 |  4.6M |  138 |     53084 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.497 |     68769 | 4.61M |  142 |     56884 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.498 |     85548 | 4.56M |  115 |     50789 |     79 |   10 |  16.3 | 218.0 |
| MRX80P001 |   80.0 |  0.497 |     78749 | 4.61M |  120 |     52825 |     79 |   10 |  16.3 | 218.0 |
| MRX80P002 |   80.0 |  0.497 |     63333 | 4.61M |  128 |     54501 |     79 |   10 |  16.3 | 218.0 |

Table: statMRUnitigsUnitig.md

| Name                      | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.496 |    125881 | 4.69M |   85 |      3159 |    371 |   36 |  87.7 | 958.0 |
| 7_merge_mr_unitigs_multik |  0.494 |    117787 | 4.58M |  102 |      2051 |    371 |   36 |  87.7 | 958.0 |
| 7_merge_mr_unitigs_unitig |  0.495 |    117841 | 4.62M |   81 |      2260 |    371 |   36 |  87.7 | 958.0 |
| 7_merge_unitigs_multik    |  0.491 |     41393 | 5.19M |  246 |      7515 |    370 |   37 |  86.3 | 962.0 |
| 7_merge_unitigs_unitig    |  0.491 |     58785 | 5.09M |  200 |      5246 |    370 |   36 |  87.3 | 956.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.986 |     95433 | 4.54M |   96 |     17.3K |    285 |   35 |    60 |   780 |
| 8_mr_spades  |  0.983 |    148463 | 4.55M |   71 |    22.68K |    371 |   36 |  87.7 |   958 |
| 8_megahit    |  0.982 |     67324 | 4.53M |  125 |    18.69K |    285 |   35 |    60 |   780 |
| 8_mr_megahit |  0.988 |    132782 | 4.56M |   68 |    17.91K |    372 |   36 |    88 |   960 |

Table: statOtherAnchors.md

| Name                    |     N50 |     Sum |    # |
| ----------------------- | ------: | ------: | ---: |
| Genome                  | 4641652 |   4.64M |    1 |
| Paralogs                |    2003 | 260.35K |  131 |
| repetitive              |    1235 |  91.99K |  169 |
| 7_merge_anchors.anchors |  125881 |   4.69M |   85 |
| spades.contig           |  125603 |   4.57M |  129 |
| spades.scaffold         |  132608 |   4.57M |  125 |
| mr_spades.contig        |  148607 |   4.59M |  148 |
| mr_spades.scaffold      |  148607 |   4.59M |  146 |
| megahit.contig          |   82825 |   4.57M |  149 |
| mr_megahit.contig       |  132896 |   4.61M |  124 |

Table: statFinal

* Assembly quality by QUAST

| Assembly        | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| --------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik    |       246 |  104347 | 5194161 |  41393 |     6 | 97.963 | 1.142 |   0.00 |    0.27 |       0.12 |
| merge_mr_multik |       102 |  268282 | 4584693 | 117787 |     0 | 98.360 | 1.004 |   0.00 |    0.00 |       0.04 |
| merge_unitig    |       200 |  132778 | 5089742 |  58785 |     1 | 97.933 | 1.120 |   0.00 |    0.39 |       0.04 |
| merge_mr_unitig |        81 |  315530 | 4620846 | 117841 |     0 | 98.428 | 1.011 |   0.00 |    0.26 |       0.02 |
| merge_anchors   |        85 |  315530 | 4686121 | 125881 |     6 | 98.531 | 1.025 |   0.00 |    0.49 |       0.02 |
| spades          |       124 |  224028 | 4573092 | 125603 |     0 | 98.462 | 1.000 |   0.00 |    0.92 |       0.20 |
| mr_spades       |       148 |  284843 | 4587655 | 148607 |     0 | 98.637 | 1.002 |   0.00 |    0.74 |       0.17 |
| megahit         |       149 |  236116 | 4569032 |  82825 |     2 | 98.368 | 1.001 |   0.00 |    2.91 |       0.57 |
| mr_megahit      |       124 |  313220 | 4608030 | 132896 |     3 | 98.912 | 1.004 |   0.00 |    2.24 |       0.43 |

Table: statQuast

## *E. coli* str. K-12 substr. DH5alpha

### dh5alpha: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/dh5alpha/1_genome
cd ~/data/anchr/dh5alpha/1_genome

cp ~/data/anchr/ref/dh5alpha/genome.fa .
cp ~/data/anchr/ref/dh5alpha/paralogs.fa .

```

### dh5alpha: download

```shell
cd ~/data/anchr/dh5alpha

mkdir -p ena
cd ena

cat << EOF > source.csv
SRP251726,dh5alpha,HiSeq 2500 PE125
EOF

anchr ena meta source.csv > ena_info.json
anchr ena manifest ena_info.json

tva to md ena_info.tsv --fmt

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name     | srx        | platform        | layout | ilength | srr         | spot    | base  |
|:---------|:-----------|:----------------|:-------|:--------|:------------|:--------|:------|
| dh5alpha | SRX7856678 | ILLUMINA        | PAIRED |         | SRR11245239 | 5881654 | 1.37G |
| dh5alpha | SRX7856679 | OXFORD_NANOPORE | SINGLE |         | SRR11245238 | 346489  | 3.35G |

* Illumina

```shell
cd ~/data/anchr/dh5alpha

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR11245239_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR11245239_2.fastq.gz R2.fq.gz

```

### dh5alpha: template

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4583637 \
    --parallel 16 \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

### dh5alpha: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=dh5alpha

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

| Group             |  Mean | Median |  STDev | Pairs%/Orientation |
|-------------------|------:|-------:|-------:|-------------------:|
| R.genome.bbtools  | 470.7 |    346 | 2460.3 |             99.98% |
| R.tadpole.bbtools | 389.0 |    340 |  206.1 |             94.97% |
| R.genome.picard   | 394.8 |    346 |  208.3 |                 FR |
| R.tadpole.picard  | 389.0 |    340 |  205.8 |                 FR |

Table: statInsertSize

| K    | property              |          min |          max |
|------|-----------------------|-------------:|-------------:|
| R.21 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,422,294 bp |
|      | Genome Repeat Length  |   120,434 bp |   120,514 bp |
|      | Genome Unique Length  | 4,300,398 bp | 4,303,243 bp |
|      | Model Fit             |     97.6629% |     98.0363% |
|      | Read Error Rate       |              |    0.209622% |
|      | Kmer Cov              |              |        267.1 |
| R.51 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,441,231 bp |
|      | Genome Repeat Length  |    94,823 bp |    94,884 bp |
|      | Genome Unique Length  | 4,344,991 bp | 4,347,764 bp |
|      | Model Fit             |     97.9127% |      98.422% |
|      | Read Error Rate       |              |    0.131634% |
|      | Kmer Cov              |              |        185.6 |
| R.81 | Homozygous (a)        |              |         100% |
|      | Genome Haploid Length |              | 4,466,733 bp |
|      | Genome Repeat Length  |    87,251 bp |    87,306 bp |
|      | Genome Unique Length  | 4,378,088 bp | 4,380,823 bp |
|      | Model Fit             |     97.4061% |      98.803% |
|      | Read Error Rate       |              |    0.109848% |
|      | Kmer Cov              |              |        108.3 |

Table: statFastK

| chr       | chrLength |   size | coverage |
| --------- | --------: | -----: | -------: |
| NZ_CP017100 | 4583637 | 110586 |   0.0241 |
| all       |   4583637 | 110586 |   0.0241 |

Table: statRepetitive

| Name       |     N50 |     Sum |        # |
|------------|--------:|--------:|---------:|
| genome     | 4583637 |   4.58M |        1 |
| paralogs   |    1737 | 188.16K |      111 |
| repetitive |    1175 | 110.59K |      190 |
| Illumina.R |     125 |   1.47G | 11763308 |
| trim.R     |     125 |   1.37G | 10962178 |
| Q0L0       |     125 |   1.37G | 10962178 |
| Q25L60     |     125 |   1.25G | 10280852 |
| Q30L60     |     125 |   1.13G |  9405463 |

Table: statReads

| Name     | N50 |     Sum |        # |
|----------|----:|--------:|---------:|
| clumpify | 125 |   1.37G | 10970448 |
| highpass | 125 |   1.37G | 10966054 |
| trim     | 125 |   1.37G | 10962178 |
| filter   | 125 |   1.37G | 10962178 |
| R1       | 125 | 682.99M |  5481089 |
| R2       | 125 | 683.63M |  5481089 |
| Rs       |   0 |       0 |        0 |

Table: statTrimReads

```text
#R.trim
#Matched	5620	0.05125%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

| Name          | N50 |     Sum |        # |
|---------------|----:|--------:|---------:|
| clumped       | 125 |   1.37G | 10959360 |
| ecco          | 125 |   1.37G | 10959360 |
| eccc          | 125 |   1.37G | 10959360 |
| ecct          | 125 |   1.37G | 10952518 |
| extended      | 165 |    1.8G | 10952518 |
| merged.raw    | 343 |    1.1G |  3510900 |
| unmerged.raw  | 165 | 646.07M |  3930718 |
| unmerged.trim | 165 | 646.07M |  3930718 |
| M1            | 343 |   1.06G |  3403638 |
| U1            | 165 | 322.94M |  1965359 |
| U2            | 165 | 323.12M |  1965359 |
| Us            |   0 |       0 |        0 |
| M.cor         | 250 |   1.71G | 10737994 |

Table: statMergeReads

| Group              |  Mean | Median | STDev | Pairs% |
|--------------------|------:|-------:|------:|-------:|
| M.ihist.merge1.txt | 172.9 |    173 |  27.8 | 21.04% |
| M.ihist.merge.txt  | 312.5 |    310 |  87.7 | 64.11% |

Table: statMergeInsert

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real | RunTime |
|----------|------:|-------:|---------:|-----:|------:|------:|---------:|--------:|
| Q0L0.R   | 298.2 |  262.0 |   12.13% | "87" | 4.58M |  4.6M |     1.00 | 0:02:07 |
| Q25L60.R | 273.7 |  253.4 |    7.43% | "87" | 4.58M | 4.53M |     0.99 | 0:01:57 |
| Q30L60.R | 246.5 |  232.2 |    5.80% | "87" | 4.58M | 4.52M |     0.99 | 0:01:44 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.976 |     54931 | 4.54M |  179 |     17199 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.977 |     54863 | 4.57M |  177 |     16924 |     40 |    6 |   7.3 | 116.0 |
| Q0L0X40P002   |   40.0 |  0.978 |     57804 | 4.58M |  169 |     18787 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.978 |     63676 |  4.5M |  142 |     21256 |     80 |    9 |  17.7 | 214.0 |
| Q0L0X80P001   |   80.0 |  0.978 |     61232 | 4.54M |  134 |     25687 |     80 |    9 |  17.7 | 214.0 |
| Q0L0X80P002   |   80.0 |  0.978 |     63676 | 4.49M |  127 |     22936 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X40P000 |   40.0 |  0.978 |     58783 | 4.58M |  176 |     15789 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X40P001 |   40.0 |  0.978 |     55933 |  4.6M |  189 |     19875 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X40P002 |   40.0 |  0.978 |     54846 | 4.54M |  176 |     19225 |     40 |    5 |   8.3 | 110.0 |
| Q25L60X80P000 |   80.0 |  0.979 |     64186 |  4.5M |  123 |     21229 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X80P001 |   80.0 |  0.979 |     78657 | 4.48M |  125 |     23786 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X80P002 |   80.0 |  0.979 |     67297 | 4.49M |  123 |     22328 |     80 |    9 |  17.7 | 214.0 |
| Q30L60X40P000 |   40.0 |  0.978 |     57814 | 4.57M |  192 |     18780 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X40P001 |   40.0 |  0.978 |     57025 |  4.6M |  181 |     15701 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X40P002 |   40.0 |  0.978 |     58430 | 4.53M |  175 |     18761 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X80P000 |   80.0 |  0.980 |     66134 | 4.53M |  129 |     22429 |     80 |    9 |  17.7 | 214.0 |
| Q30L60X80P001 |   80.0 |  0.980 |     78665 | 4.48M |  121 |     22051 |     80 |    9 |  17.7 | 214.0 |

Table: statUnitigsMultik.md

| Name          | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.976 |     57110 | 4.54M |  176 |     17644 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X40P001   |   40.0 |  0.977 |     54083 |  4.6M |  184 |     15663 |     40 |    6 |   7.3 | 116.0 |
| Q0L0X40P002   |   40.0 |  0.978 |     55676 | 4.55M |  172 |     19514 |     40 |    5 |   8.3 | 110.0 |
| Q0L0X80P000   |   80.0 |  0.978 |     54886 | 4.48M |  171 |     20953 |     80 |    9 |  17.7 | 214.0 |
| Q0L0X80P001   |   80.0 |  0.977 |     42172 | 4.48M |  182 |     23536 |     80 |    9 |  17.7 | 214.0 |
| Q0L0X80P002   |   80.0 |  0.978 |     45621 | 4.48M |  158 |     22443 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X40P000 |   40.0 |  0.978 |     59716 | 4.59M |  169 |     14758 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X40P001 |   40.0 |  0.978 |     58342 | 4.58M |  174 |     17884 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X40P002 |   40.0 |  0.978 |     63662 | 4.54M |  161 |     14683 |     40 |    6 |   7.3 | 116.0 |
| Q25L60X80P000 |   80.0 |  0.979 |     59740 | 4.55M |  143 |     20322 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X80P001 |   80.0 |  0.979 |     58517 | 4.48M |  139 |     22769 |     80 |    9 |  17.7 | 214.0 |
| Q25L60X80P002 |   80.0 |  0.979 |     59738 | 4.54M |  140 |     21125 |     80 |    9 |  17.7 | 214.0 |
| Q30L60X40P000 |   40.0 |  0.978 |     60389 | 4.57M |  173 |     15659 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X40P001 |   40.0 |  0.978 |     63102 | 4.58M |  174 |     16476 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X40P002 |   40.0 |  0.978 |     58430 | 4.53M |  160 |     18163 |     40 |    6 |   7.3 | 116.0 |
| Q30L60X80P000 |   80.0 |  0.979 |     63680 | 4.48M |  140 |     21304 |     80 |    9 |  17.7 | 214.0 |
| Q30L60X80P001 |   80.0 |  0.979 |     58854 | 4.48M |  138 |     21643 |     80 |    9 |  17.7 | 214.0 |

Table: statUnitigsUnitig.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.674 |     79994 | 4.49M |  134 |     63615 |     40 |    5 |   8.3 | 110.0 |
| MRX40P001 |   40.0 |  0.674 |     86851 |  4.5M |  135 |     59661 |     40 |    5 |   8.3 | 110.0 |
| MRX40P002 |   40.0 |  0.674 |     84687 |  4.5M |  135 |     60077 |     40 |    5 |   8.3 | 110.0 |
| MRX80P000 |   80.0 |  0.674 |     95175 | 4.48M |  123 |     62681 |     80 |    8 |  18.7 | 208.0 |
| MRX80P001 |   80.0 |  0.674 |     95148 | 4.48M |  116 |     60606 |     80 |    8 |  18.7 | 208.0 |
| MRX80P002 |   80.0 |  0.674 |     95105 | 4.48M |  127 |     65572 |     80 |    8 |  18.7 | 208.0 |

Table: statMRUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.670 |     85670 |  4.7M |  153 |     19796 |     40 |    5 |   8.3 | 110.0 |
| MRX40P001 |   40.0 |  0.670 |     85700 | 4.57M |  144 |     18496 |     40 |    5 |   8.3 | 110.0 |
| MRX40P002 |   40.0 |  0.670 |     83066 | 4.59M |  149 |     18890 |     40 |    5 |   8.3 | 110.0 |
| MRX80P000 |   80.0 |  0.670 |     82967 | 4.53M |  106 |     19455 |     80 |    8 |  18.7 | 208.0 |
| MRX80P001 |   80.0 |  0.670 |     86404 | 4.56M |  106 |     19654 |     80 |    8 |  18.7 | 208.0 |
| MRX80P002 |   80.0 |  0.670 |     85779 | 4.51M |  119 |     19555 |     80 |    8 |  18.7 | 208.0 |

Table: statMRUnitigsUnitig.md

| Name                      | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.669 |    112871 | 4.71M |   96 |      6980 |    376 |   29 |  96.3 | 926.0 |
| 7_merge_mr_unitigs_multik |  0.665 |    102446 | 4.56M |  116 |      7650 |    376 |   29 |  96.3 | 926.0 |
| 7_merge_mr_unitigs_unitig |  0.669 |    112871 | 4.69M |  111 |      6447 |    376 |   29 |  96.3 | 926.0 |
| 7_merge_unitigs_multik    |  0.664 |     80614 | 4.59M |  117 |      7562 |    375 |   29 |  96.0 | 924.0 |
| 7_merge_unitigs_unitig    |  0.663 |     80614 |  4.5M |  114 |      6697 |    375 |   29 |  96.0 | 924.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |  Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ---: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.984 |    112448 | 4.47M |   76 |   18.41K |    263 |   22 |  65.7 |   658 |
| 8_mr_spades  |  0.679 |    132652 |  4.5M |   73 |     25274 |    376 |   29 |  96.3 | 926.0 |
| 8_megahit    |  0.985 |     67322 | 4.46M |  117 |   25.03K |    263 |   22 |  65.7 |   658 |
| 8_mr_megahit |  0.680 |    132628 |  4.5M |   79 |     51938 |    376 |   29 |  96.3 | 926.0 |

Table: statOtherAnchors.md

| Name                    |     N50 |     Sum |    # |
| ----------------------- | ------: | ------: | ---: |
| Genome                  | 4583637 |   4.58M |    1 |
| Paralogs                |    1737 | 188.16K |  111 |
| repetitive              |    1175 | 110.59K |  190 |
| 7_merge_anchors.anchors |  112871 |   4.71M |   96 |
| spades.contig           |  114710 |   4.52M |  171 |
| spades.scaffold         |  143522 |   4.52M |  163 |
| mr_spades.contig        |  178373 |   4.52M |   87 |
| mr_spades.scaffold      |  203812 |   4.52M |   84 |
| megahit.contig          |   85613 |   4.51M |  175 |
| mr_megahit.contig       |  133730 |   4.56M |  154 |

Table: statFinal

* Assembly quality by QUAST

| Assembly        | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| --------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik    |       117 |  258515 | 4587460 |  80614 |     2 | 97.787 | 1.023 |   0.00 |    0.54 |       0.04 |
| merge_mr_multik |       116 |  258601 | 4559551 | 102446 |     0 | 97.962 | 1.015 |   0.00 |    0.07 |       0.09 |
| merge_unitig    |       114 |  258516 | 4503384 |  80614 |     0 | 97.752 | 1.005 |   0.00 |    0.67 |       0.02 |
| merge_mr_unitig |       111 |  310161 | 4687516 | 112871 |     1 | 98.296 | 1.040 |   0.00 |    0.51 |       0.15 |
| merge_anchors   |        96 |  310394 | 4714287 | 112871 |     3 | 98.365 | 1.046 |   0.00 |    0.55 |       0.15 |
| spades          |       154 |  268214 | 4513903 | 114710 |     1 | 98.381 | 1.001 |   0.00 |    1.15 |       0.38 |
| mr_spades       |        87 |  402458 | 4523351 | 178373 |     2 | 98.612 | 1.001 |   0.00 |    1.41 |       0.20 |
| megahit         |       175 |  258343 | 4510262 |  85613 |    11 | 98.262 | 1.001 |   0.00 |    2.91 |       0.38 |
| mr_megahit      |       154 |  359335 | 4556625 | 133730 |     2 | 98.906 | 1.005 |   0.00 |    2.41 |       0.44 |

Table: statQuast

## *Bacillus cereus* ATCC 10987

> GAGE-B MiSeq 数据集（100× 子集），reads 与参考**同株**（见 reference
> 节实测）。老流程基线（2025-06）来自 `results/gage_b.md`；现代流程门禁
> 待跑（见 `scripts/asm_gate.md` §gate history）。

### bcer: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/Bcer_100x/1_genome
cd ~/data/anchr/Bcer_100x/1_genome

cp ~/data/anchr/ref/Bcer/genome.fa .
cp ~/data/anchr/ref/Bcer/paralogs.fa .
cp ~/data/anchr/ref/Bcer/repetitive.fa .

```

* 参考构成：NC_003909（染色体 5,224,283 bp）+ NC_005707（质粒 pBc10987，
  208,369 bp）= 5,432,652 bp / 2 复制子；GC ~38%（低 GC）。
* **reads 与参考同株**：GAGE-B 数据集即 ATCC 10987。本地实测（50k 原始
  R1，2026-08-18）：reads 的 solid 31-mer 99.66% 出现在参考中、85.7%
  整条 perfect 回贴——菌株差异≈0，QUAST 判据可直接使用（对比 Mabs/Vcho
  的跨株情况，见下）。

### bcer: download

* Illumina（GAGE-B 100× 子集）

```shell
cd ~/data/anchr/Bcer_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c https://ccb.jhu.edu/gage_b/datasets/B_cereus_MiSeq.tar.gz

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

* 本地已就位：`~/data/anchr/Bcer_100x/2_illumina/` 有 R1.fq.gz + R2.fq.gz
  （2024-12 下载）；`1_genome/` 尚空，按 reference 节 cp 即可。

### bcer: template

* 现代流程模板（与 g37/dh5alpha 同口径；待跑）

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 5432652 \
    --parallel 24 \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig bcalm" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

* 老流程模板（2025-06 GAGE-B 运行，基线数据见下节）：`--genome 5432652
  --fastqc --insertsize --kat --trim "--dedupe --tile --cutoff 5 --cutk 31"
  --qual "20 25 30" --len "60" --filter "adapter artifact" --quorum --merge
  --ecphase "1 2 3" --cov "40 50 60 all" --unitigger "superreads bcalm
  tadpole" --statp 2 --readl 250 --uscale 2 --lscale 3 --redo --extend`
  （完整命令见 `results/gage_b.md`）。

### bcer: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Bcer_100x

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

现代流程门禁**待跑**（本地数据已齐，可直接执行；见
`scripts/asm_gate.md` §gate history）；老流程 2025-06 运行见
`results/gage_b.md`。

### bcer: 老流程基线（2025-06 GAGE-B）

reads：2,080,000 对 × 250 bp（~100×，481.02 Mb）；接头污染极低
（#R.trim 匹配 0.06267%）。老流程 quorum（Q0L0）丢弃 12.98%。

Table: statReads（老流程）

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 5224283 | 5432652 |       2 |
| Paralogs    |    2295 |  220468 |     101 |
| repetitive |    2461 |  113050 |     173 |
| Illumina.R  |     251 | 481.02M | 2080000 |
| trim.R      |     250 | 404.36M | 1807384 |
| Q20L60      |     250 | 396.79M | 1758525 |
| Q25L60      |     250 | 379.57M | 1706208 |
| Q30L60      |     250 | 344.25M | 1611233 |

Table: statTrimReads（老流程）

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 251 | 480.99M | 2079856 |
| filteredbytile | 251 | 463.51M | 2005704 |
| highpass       | 251 | 459.84M | 1989878 |
| trim           | 250 | 404.41M | 1807576 |
| filter         | 250 | 404.36M | 1807384 |
| R1             | 250 | 209.26M |  903692 |
| R2             | 247 |  195.1M |  903692 |
| Rs             |   0 |       0 |       0 |

```text
#R.trim
#Matched	1247	0.06267%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	99	0.00548%
#Name	Reads	ReadsPct
```

Table: statMergeReads（老流程）

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 250 | 404.36M | 1807382 |
| ecco          | 250 | 404.36M | 1807382 |
| eccc          | 250 | 404.36M | 1807382 |
| ecct          | 250 | 400.39M | 1786476 |
| extended      | 290 | 470.96M | 1786476 |
| merged.raw    | 586 | 316.57M |  584078 |
| unmerged.raw  | 285 | 149.75M |  618320 |
| unmerged.trim | 285 | 149.74M |  618292 |
| M1            | 586 | 316.54M |  584024 |
| U1            | 290 |  79.96M |  309146 |
| U2            | 270 |  69.78M |  309146 |
| Us            |   0 |       0 |       0 |
| M.cor         | 518 | 466.86M | 1786340 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 362.0 |    388 |  97.6 |         19.50% |
| M.ihist.merge.txt  | 542.0 |    564 | 120.0 |         65.39% |

Table: statQuorum（老流程）

| Name     | CovIn | CovOut | Discard% |  Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|------:|------:|------:|---------:|----------:|
| Q0L0.R   |  74.4 |   64.8 |   12.98% | "127" | 5.43M | 5.35M |     0.98 | 0:01'00'' |
| Q20L60.R |  73.0 |   64.7 |   11.48% | "127" | 5.43M | 5.35M |     0.98 | 0:00'52'' |
| Q25L60.R |  69.9 |   63.8 |    8.71% | "127" | 5.43M | 5.34M |     0.98 | 0:00'50'' |
| Q30L60.R |  63.4 |   59.7 |    5.75% | "127" | 5.43M | 5.34M |     0.98 | 0:00'49'' |

Table: statMRUnitigsSuperreads（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.62% |     37925 | 5.31M | 240 |        81 | 22.62K | 519 |   39.0 |  7.0 |   6.0 | 120.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:34 |
| MRX40P001  |   40.0 |  97.54% |     40827 | 5.31M | 237 |        64 | 18.73K | 454 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:32 |
| MRX50P000  |   50.0 |  97.58% |     37545 | 5.31M | 249 |        90 | 22.64K | 547 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:36 |
| MRX60P000  |   60.0 |  97.53% |     37952 | 5.31M | 250 |        92 | 20.41K | 527 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:37 |
| MRXallP000 |   85.9 |  97.46% |     36340 | 5.33M | 253 |        90 |  18.9K | 537 |   84.0 | 15.0 |  13.0 | 258.0 | "31,41,51,61,71,81" |   0:02:00 |   0:00:35 |

Table: statMRUnitigsBcalm（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.81% |     39857 | 5.34M | 232 |        87 | 18.83K | 442 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:30 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     42781 | 5.31M | 229 |        86 | 18.97K | 412 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:31 |
| MRX50P000  |   50.0 |  97.44% |     39857 | 5.34M | 236 |        80 | 16.34K | 466 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:01:39 |   0:00:33 |
| MRX60P000  |   60.0 |  97.42% |     41648 | 5.38M | 238 |        98 | 17.86K | 463 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:01:40 |   0:00:35 |
| MRXallP000 |   85.9 |  97.38% |     39857 | 5.31M | 246 |        97 | 16.29K | 489 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:01:54 |   0:00:35 |

Table: statMRUnitigsTadpole（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000  |   40.0 |  97.83% |     44440 | 5.32M | 217 |       107 |  18.4K | 404 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:32 |
| MRX40P001  |   40.0 |  97.77% |     44376 | 5.36M | 219 |        83 | 16.18K | 378 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:33 |
| MRX50P000  |   50.0 |  97.76% |     42799 | 5.34M | 224 |       131 | 19.98K | 444 |   49.0 |  9.0 |   7.3 | 152.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:32 |
| MRX60P000  |   60.0 |  97.77% |     42132 | 5.34M | 228 |       131 | 19.13K | 443 |   59.0 | 11.0 |   8.7 | 184.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:34 |
| MRXallP000 |   85.9 |  97.70% |     41648 | 5.32M | 239 |       120 | 18.37K | 484 |   85.0 | 15.0 |  13.3 | 260.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:35 |

Table: statMergeAnchors（老流程）

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.76% |     36737 |  5.3M | 254 |     41729 | 291.14K | 43 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:42 |
| 7_merge_mr_unitigs_bcalm      |  98.15% |     37087 | 5.21M | 252 |     41793 |  70.11K |  9 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_mr_unitigs_superreads |  97.83% |     34002 | 5.15M | 259 |      1091 |    6.8K |  7 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:43 |
| 7_merge_mr_unitigs_tadpole    |  98.14% |     42769 | 5.21M | 232 |     41634 |  53.78K | 12 |   63.0 | 13.0 |   8.0 | 204.0 |   0:00:43 |
| 7_merge_unitigs_bcalm         |  98.04% |     32860 | 5.31M | 273 |     36517 |  63.09K | 19 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:49 |
| 7_merge_unitigs_superreads    |  98.10% |     32344 | 5.31M | 273 |     41729 |  64.41K | 24 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:46 |
| 7_merge_unitigs_tadpole       |  98.02% |     32710 | 5.31M | 277 |     62295 |  81.41K | 21 |   64.0 | 12.0 |   9.3 | 200.0 |   0:00:45 |

Table: statOtherAnchors（老流程）

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.96% |     49586 | 2.14M |  80 |       899 |  9.35K | 106 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:33 |
| 8_mr_spades  |  98.81% |     78041 | 5.35M | 130 |       217 | 10.11K | 212 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |
| 8_megahit    |  98.69% |     37648 | 4.85M | 218 |       125 | 21.73K | 356 |   64.0 | 13.0 |   8.3 | 206.0 |   0:00:40 |
| 8_mr_megahit |  98.85% |     65403 | 5.36M | 153 |       451 | 13.53K | 265 |   85.0 | 15.0 |  13.3 | 260.0 |   0:00:36 |

Table: statFinal（老流程）

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 5224283 | 5432652 |   2 |
| Paralogs                 |    2295 |  220468 | 101 |
| repetitive              |    2461 |  113050 | 173 |
| 7_merge_anchors.anchors  |   36737 | 5304696 | 254 |
| 7_merge_anchors.others   |   41729 |  291138 |  43 |
| glue_anchors             |   36737 | 5304175 | 253 |
| fill_anchors             |   61075 | 5317582 | 159 |
| spades.contig            |  207470 | 5366804 | 153 |
| spades.scaffold          |  285416 | 5367163 | 139 |
| spades.non-contained     |  207470 | 5349666 |  58 |
| mr_spades.contig         |  100015 | 5367895 | 128 |
| mr_spades.scaffold       |  284294 | 5374592 |  66 |
| mr_spades.non-contained  |  100015 | 5361433 | 105 |
| megahit.contig           |   59732 | 5360219 | 204 |
| megahit.non-contained    |   59732 | 5341409 | 158 |
| mr_megahit.contig        |   75019 | 5388594 | 186 |
| mr_megahit.non-contained |   75019 | 5369027 | 141 |

## *Rhodobacter sphaeroides* 2.4.1

> GAGE-B MiSeq 数据集（100× 子集）。参考即 2.4.1（名义同株，无本地
> reads 实测）；结构复杂（2 染色体 + 5 质粒 = 7 复制子、重复 12.8%、
> GC ~68%）正是纳入理由——对"无 N 染色体"目标的多复制子压力测试。
> 老流程基线（2025-06）来自 `results/gage_b.md`；现代流程门禁待跑（见
> `scripts/asm_gate.md` §gate history）。

### rsph: reference

* Reference genome

```shell
mkdir -p ~/data/anchr/Rsph_100x/1_genome
cd ~/data/anchr/Rsph_100x/1_genome

cp ~/data/anchr/ref/Rsph/genome.fa .
cp ~/data/anchr/ref/Rsph/paralogs.fa .
cp ~/data/anchr/ref/Rsph/repetitive.fa .

```

* 参考构成：NC_007493（染色体 1，3,188,524 bp）+ NC_007494（染色体 2，
  1,314,453 bp）+ 5 个质粒（NC_009007/NC_007488/NC_007489/NC_007490/
  NC_009008）= 4,602,977 bp / 7 复制子；GC ~68%（极高 GC），重复含量
  12.8%（GAGE-B 四菌株中最高）。

### rsph: download

* Illumina（GAGE-B 100× 子集）

```shell
cd ~/data/anchr/Rsph_100x

mkdir -p 2_illumina
cd 2_illumina

aria2c -x 4 -s 2 -c http://ccb.jhu.edu/gage_b/datasets/R_sphaeroides_MiSeq.tar.gz

# NOT gzipped tar
tar xvf R_sphaeroides_MiSeq.tar.gz raw/insert_540_1__cov100x.fastq
tar xvf R_sphaeroides_MiSeq.tar.gz raw/insert_540_2__cov100x.fastq

cat raw/insert_540_1__cov100x.fastq |
    pigz -p 8 -c \
    > R1.fq.gz
cat raw/insert_540_2__cov100x.fastq |
    pigz -p 8 -c \
    > R2.fq.gz

rm -fr raw

```

* 本地暂无 reads（`~/data/anchr/Rsph_100x/` 为空），需先下载。

### rsph: template

* 现代流程模板（与 g37/dh5alpha 同口径；待跑）

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}

rm 0_script/*
anchr template \
    --genome 4602977 \
    --parallel 24 \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --merge \
    \
    --cov "40 80" \
    --unitigger "multik unitig bcalm" \
    --statp 2 \
    --uscale 2 \
    --lscale 3

```

* 老流程模板（2025-06 GAGE-B 运行，基线数据见下节）：`--genome 4602977
  --fastqc --insertsize --kat --trim "--dedupe" --qual "20 25 30" --len
  "60" --filter "adapter artifact" --quorum --merge --ecphase "1 3" --cov
  "30 all" --unitigger "superreads bcalm tadpole" --statp 2 --readl 250
  --uscale 2 --lscale 3 --redo --extend`（完整命令见
  `results/gage_b.md`）。

### rsph: run

```shell
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Rsph_100x

cd ${WORKING_DIR}/${BASE_NAME}

bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh

# bash 0_script/0_cleanup.sh

```

现代流程门禁**待跑**（需先按 download 节取数；见
`scripts/asm_gate.md` §gate history）；老流程 2025-06 运行见
`results/gage_b.md`（QUAST 用 `--no-check`，多复制子参考）。

### rsph: 老流程基线（2025-06 GAGE-B）

reads：1,800,000 对 × 250 bp（~100×，451.8 Mb）；接头匹配 6.39%
（Reverse_adapter 4.58% 为主）。老流程 quorum（Q0L0）丢弃 11.10%。

Table: statReads（老流程）

| Name        |     N50 |     Sum |       # |
|:------------|--------:|--------:|--------:|
| Genome      | 3188524 | 4602977 |       7 |
| Paralogs    |    2337 |  146789 |      66 |
| repetitive |     572 |   57281 |     165 |
| Illumina.R  |     251 |  451.8M | 1800000 |
| trim.R      |     148 |  200.1M | 1452706 |
| Q20L60      |     148 | 193.66M | 1401466 |
| Q25L60      |     139 | 169.12M | 1304628 |
| Q30L60      |     119 | 125.02M | 1123194 |

Table: statTrimReads（老流程）

| Name     | N50 |     Sum |       # |
|:---------|----:|--------:|--------:|
| clumpify | 251 | 447.53M | 1782994 |
| trim     | 148 |  200.1M | 1452706 |
| filter   | 148 |  200.1M | 1452706 |
| R1       | 164 | 100.23M |  655190 |
| R2       | 133 |  81.52M |  655190 |
| Rs       | 141 |  18.34M |  142326 |

```text
#R.trim
#Matched	113970	6.39206%
#Name	Reads	ReadsPct
Reverse_adapter	81598	4.57646%
pcr_dimer	14481	0.81217%
PCR_Primers	8081	0.45323%
TruSeq_Universal_Adapter	5665	0.31772%
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

Table: statMergeReads（老流程）

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 148 | 200.09M | 1452579 |
| ecco          | 148 | 199.84M | 1452578 |
| ecct          | 148 | 198.72M | 1444000 |
| extended      | 186 | 255.79M | 1444000 |
| merged.raw    | 455 | 197.38M |  475527 |
| unmerged.raw  | 172 |  80.09M |  492946 |
| unmerged.trim | 172 |  80.07M |  492605 |
| M1            | 455 |  197.2M |  475127 |
| U1            | 172 |  19.67M |  121605 |
| U2            | 151 |  17.53M |  121605 |
| Us            | 182 |  42.86M |  249395 |
| M.cor         | 443 | 277.99M | 1692254 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 184.3 |    179 |  66.0 |         10.54% |
| M.ihist.merge.txt  | 415.1 |    452 |  89.0 |         65.86% |

Table: statQuorum（老流程）

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   |  43.5 |   38.7 |   11.10% | "39" |  4.6M | 4.55M |     0.99 | 0:00'37'' |
| Q20L60.R |  42.1 |   37.9 |    9.98% | "39" |  4.6M | 4.55M |     0.99 | 0:00'32'' |
| Q25L60.R |  36.8 |   34.9 |    5.03% | "35" |  4.6M | 4.54M |     0.99 | 0:00'29'' |
| Q30L60.R |  27.2 |   26.6 |    2.20% | "31" |  4.6M | 4.52M |     0.98 | 0:00'26'' |

Table: statMRUnitigsSuperreads（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.39% |     20508 | 4.35M | 364 |      4332 |  140.5K | 739 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRX30P001  |   30.0 |  97.40% |     21605 | 4.34M | 357 |      3817 | 154.53K | 743 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| MRXallP000 |   60.4 |  97.36% |     22279 | 4.34M | 345 |      5800 | 167.85K | 749 |   55.0 | 9.0 |   9.3 | 164.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:29 |

Table: statMRUnitigsBcalm（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.41% |     18774 | 4.35M | 404 |      5004 | 173.16K | 899 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:30 |
| MRX30P001  |   30.0 |  97.35% |     19071 | 4.34M | 389 |      6101 | 196.85K | 857 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:28 |
| MRXallP000 |   60.4 |  97.44% |     20747 | 4.33M | 353 |      6101 |  199.4K | 792 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:30 |

Table: statMRUnitigsTadpole（老流程）

| Name       | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:-----------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX30P000  |   30.0 |  97.50% |     19948 | 4.34M | 379 |      5224 | 161.26K | 801 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:28 |
| MRX30P001  |   30.0 |  97.55% |     19862 | 4.34M | 372 |      5475 | 173.65K | 799 |   27.0 | 5.0 |   5.0 |  84.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:29 |
| MRXallP000 |   60.4 |  97.54% |     21368 | 4.33M | 344 |      6100 | 177.53K | 759 |   55.0 | 8.0 |  10.3 | 158.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:29 |

Table: statMergeAnchors（老流程）

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  91.17% |     36740 | 4.41M | 246 |      4689 | 360.25K | 133 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:33 |
| 7_merge_mr_unitigs_bcalm      |  92.62% |     21353 | 4.33M | 349 |      5224 | 203.81K |  65 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_mr_unitigs_superreads |  92.48% |     22986 | 4.39M | 328 |      5487 | 188.38K |  62 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:34 |
| 7_merge_mr_unitigs_tadpole    |  92.68% |     22194 | 4.33M | 335 |      5819 | 199.44K |  57 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:35 |
| 7_merge_unitigs_bcalm         |  90.95% |     17117 | 4.33M | 424 |      5224 | 311.57K | 102 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:30 |
| 7_merge_unitigs_superreads    |  91.58% |     30174 | 4.35M | 271 |      4689 | 308.94K | 114 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:31 |
| 7_merge_unitigs_tadpole       |  91.26% |     22998 | 4.34M | 319 |      5018 | 335.15K | 115 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:30 |

Table: statOtherAnchors（老流程）

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.20% |     45176 | 1.94M |  83 |      7794 |  86.43K | 125 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:29 |
| 8_mr_spades  |  99.13% |     30528 | 3.93M | 216 |     16123 |  151.9K | 304 |   55.0 |  8.0 |  10.3 | 158.0 |   0:00:31 |
| 8_megahit    |  98.54% |     37964 |  4.1M | 230 |      5525 | 153.99K | 387 |   35.0 |  7.0 |   5.0 | 112.0 |   0:00:28 |
| 8_mr_megahit |  99.34% |     26589 | 4.38M | 287 |     16123 | 164.89K | 523 |   55.0 |  8.0 |  10.3 | 158.0 |   0:00:31 |

Table: statFinal（老流程）

| Name                     |     N50 |     Sum |    # |
|:-------------------------|--------:|--------:|-----:|
| Genome                   | 3188524 | 4602977 |    7 |
| Paralogs                 |    2337 |  146789 |   66 |
| repetitive              |     572 |   57281 |  165 |
| 7_merge_anchors.anchors  |   36740 | 4407623 |  246 |
| 7_merge_anchors.others   |    4689 |  360250 |  133 |
| glue_anchors             |   37591 | 4404044 |  240 |
| fill_anchors             |   48535 | 4406444 |  183 |
| spades.contig            |  150729 | 4576779 |  136 |
| spades.scaffold          |  172916 | 4577123 |  131 |
| spades.non-contained     |  150729 | 4562257 |   71 |
| mr_spades.contig         |   55603 | 4566224 |  170 |
| mr_spades.scaffold       |   89512 | 4567395 |  121 |
| mr_spades.non-contained  |   55603 | 4555133 |  149 |
| megahit.contig           |   52830 | 4572904 |  245 |
| megahit.non-contained    |   52830 | 4541309 |  182 |
| mr_megahit.contig        |   31157 | 4576803 |  282 |
| mr_megahit.non-contained |   31157 | 4563775 |  255 |
