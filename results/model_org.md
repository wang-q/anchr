# Assemble genomes of model organisms by `anchr`

<!-- TOC -->
* [Assemble genomes of model organisms by `anchr`](#assemble-genomes-of-model-organisms-by-anchr)
  * [*Mycoplasma genitalium* G37](#mycoplasma-genitalium-g37)
    * [g37: reference](#g37-reference)
    * [g37: download](#g37-download)
    * [g37: template](#g37-template)
    * [g37: run](#g37-run)
<!-- TOC -->

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
    --parallel 8 \
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

| Group             |   Mean | Median |   STDev | Pairs%/Orientation |
| ----------------- | -----: | -----: | ------: | -----------------: |
| R.genome.bbtools  | 2035.3 |    452 | 21883.0 |             97.74% |
| R.tadpole.bbtools |  462.4 |    447 |   130.2 |             92.63% |
| R.genome.picard   |  466.5 |    452 |   127.8 |                 FR |
| R.tadpole.picard  |  462.1 |    447 |   126.1 |                 FR |

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
| NC_000908 |    580076 | 16875 |   0.0291 |
| all       |    580076 | 16875 |   0.0291 |

Table: statRepetitive

| Name       |    N50 |     Sum |      # |
| ---------- | -----: | ------: | -----: |
| genome     | 580076 | 580.08K |      1 |
| paralogs   |   1567 |  11.53K |      8 |
| repetitive |    499 |  16.88K |     53 |
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

| Name     | CovIn | CovOut | Discard% | Kmer |   RealG |    EstG | Est/Real | RunTime |
| -------- | ----: | -----: | -------: | ---: | ------: | ------: | -------: | ------: |
| Q0L0.R   | 174.9 |  160.8 |    8.07% | "24" | 580.08K |  578.3K |     1.00 | 0:00:11 |
| Q25L60.R | 169.4 |  160.5 |    5.26% | "24" | 580.08K | 576.79K |     0.99 | 0:00:10 |
| Q30L60.R | 163.9 |  157.3 |    4.01% | "24" | 580.08K |  578.3K |     1.00 | 0:00:10 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.958 |     21302 | 555.35K |   44 |     32953 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P001   |   40.0 |  0.959 |     22949 | 555.44K |   38 |     32005 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P002   |   40.0 |  0.959 |     26649 | 554.88K |   38 |     32809 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X80P000   |   80.0 |  0.954 |     24026 | 552.55K |   55 |     39733 |     78 |   10 |  16.0 | 216.0 |
| Q0L0X80P001   |   80.0 |  0.958 |     26589 | 555.31K |   42 |     36453 |     79 |   10 |  16.3 | 218.0 |
| Q25L60X40P000 |   40.0 |  0.958 |     27487 | 555.77K |   45 |     30906 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P001 |   40.0 |  0.959 |     30447 | 554.97K |   35 |     31552 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P002 |   40.0 |  0.958 |     22705 | 556.24K |   45 |     31053 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X80P000 |   80.0 |  0.954 |     17826 | 552.93K |   55 |     39782 |     78 |   10 |  16.0 | 216.0 |
| Q25L60X80P001 |   80.0 |  0.955 |     27545 | 553.61K |   52 |     37777 |     78 |   10 |  16.0 | 216.0 |
| Q30L60X40P000 |   40.0 |  0.959 |     18901 | 555.32K |   42 |     31349 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P001 |   40.0 |  0.959 |     22543 |  555.7K |   43 |     31661 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P002 |   40.0 |  0.957 |     26675 | 555.18K |   43 |     32509 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X80P000 |   80.0 |  0.956 |     26602 | 553.76K |   51 |     38388 |     79 |   10 |  16.3 | 218.0 |

Table: statUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.682 |     33961 | 555.56K |   23 |     28166 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.683 |     30348 | 554.67K |   30 |     29843 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.681 |     33076 | 555.03K |   29 |     29646 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.678 |     26566 | 554.15K |   45 |     33866 |     78 |   10 |  16.0 | 216.0 |
| MRX80P001 |   80.0 |  0.678 |     33930 | 553.93K |   37 |     33683 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsMultik.md

| Name                      | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.682 |     54928 | 556.64K |   15 |      2119 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_multik |  0.680 |     54770 | 562.64K |   17 |      2232 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_multik    |  0.682 |     74328 | 609.46K |   18 |      2433 |    230 |   23 |  53.7 | 598.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.988 |     83398 | 572.73K |    8 |      8276 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_spades  |  0.710 |    580031 | 580.03K |    1 |       601 |    225 |   22 |  53.0 | 582.0 |
| 8_megahit    |  0.978 |     39418 | 565.37K |   25 |     11901 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_megahit |  0.711 |    319000 | 579.38K |    8 |     20462 |    223 |   23 |  51.3 | 584.0 |

Table: statOtherAnchors.md

| Name                    |    N50 |     Sum |    # |
| ----------------------- | -----: | ------: | ---: |
| Genome                  | 580076 | 580.08K |    1 |
| Paralogs                |   1567 |  11.53K |    8 |
| repetitive              |    499 |  16.88K |   53 |
| 7_merge_anchors.anchors |  54928 | 556.64K |   15 |
| spades.contig           | 163847 |    581K |   38 |
| spades.scaffold         | 163847 |    581K |   38 |
| mr_spades.contig        | 580506 | 580.63K |    2 |
| mr_spades.scaffold      | 580506 | 580.63K |    2 |
| megahit.contig          |  39526 | 577.27K |   46 |
| mr_megahit.contig       | 319187 | 599.85K |   47 |

Table: statFinal

* Assembly quality by QUAST

| Assembly       | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik   |        18 |   80539 | 609463 |  74328 |     0 | 95.92 | 1.095 |   0.00 |   30.86 |       3.45 |
| merge_mr_multik |       17 |  179634 | 562642 |  54770 |     0 | 95.62 | 1.014 |   0.00 |   28.62 |       3.02 |
| merge_anchors  |        15 |  179675 | 556641 |  54928 |     0 | 95.94 | 1.000 |   0.00 |   30.73 |       3.41 |
| spades         |        35 |  236302 | 580752 | 163847 |     1 | 99.06 | 1.002 |   0.00 |  222.74 |      66.49 |
| mr_spades      |         2 |  580506 | 580632 | 580506 |     0 | 100.0 | 1.001 |   0.00 |  300.77 |      91.64 |
| megahit        |        46 |  136243 | 577269 |  39526 |     6 | 97.45 | 1.003 |   0.00 |   92.26 |      20.46 |
| mr_megahit     |        47 |  319187 | 599845 | 319187 |     5 | 99.91 | 1.026 |   0.00 |  318.69 |      93.84 |

Table: statQuast
