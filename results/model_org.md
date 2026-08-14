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
    --xmx 12g \
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
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --statp 2 \
    --readl 125 \
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
| Q0L0X40P000   |   40.0 |  0.957 |     24867 | 554.54K |   46 |     33745 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P001   |   40.0 |  0.960 |     33897 | 555.99K |   32 |     30186 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P002   |   40.0 |  0.959 |     30455 | 554.76K |   36 |     32563 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X80P000   |   80.0 |  0.956 |     17233 | 553.23K |   49 |     38995 |     79 |   10 |  16.3 | 218.0 |
| Q0L0X80P001   |   80.0 |  0.952 |     22275 | 551.97K |   56 |     40648 |     78 |   10 |  16.0 | 216.0 |
| Q25L60X40P000 |   40.0 |  0.953 |     22058 |    552K |   42 |     31882 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P001 |   40.0 |  0.959 |     21733 | 555.72K |   38 |     31003 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P002 |   40.0 |  0.958 |     31333 | 555.76K |   38 |     30528 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X80P000 |   80.0 |  0.957 |     26608 | 554.63K |   46 |     37262 |     79 |   10 |  16.3 | 218.0 |
| Q25L60X80P001 |   80.0 |  0.953 |     18035 | 553.53K |   53 |     38514 |     78 |   10 |  16.0 | 216.0 |
| Q30L60X40P000 |   40.0 |  0.958 |     23164 | 555.28K |   45 |     33037 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P001 |   40.0 |  0.960 |     31455 |  555.7K |   39 |     32190 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P002 |   40.0 |  0.959 |     27384 |  555.5K |   36 |     30467 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X80P000 |   80.0 |  0.953 |     18873 | 551.91K |   54 |     39008 |     79 |   10 |  16.3 | 218.0 |

Table: statUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.681 |     33954 | 554.85K |   33 |     29783 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.682 |     28070 |  555.3K |   30 |     29235 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.681 |     39079 | 555.35K |   29 |     28209 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.680 |     33904 | 553.86K |   37 |     33008 |     78 |   10 |  16.0 | 216.0 |
| MRX80P001 |   80.0 |  0.678 |     33337 | 554.08K |   43 |     34396 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsMultik.md

| Name                      | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.682 |     54929 | 556.87K |   16 |      2271 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_multik |  0.680 |     54917 |  554.9K |   16 |      2182 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_multik    |  0.682 |     79554 |  619.2K |   18 |      2532 |    230 |   23 |  53.7 | 598.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.988 |     83398 | 572.73K |    8 |      8276 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_spades  |  0.710 |    580031 | 580.03K |    1 |       601 |    225 |   22 |  53.0 | 582.0 |
| 8_megahit    |  0.978 |     42279 | 564.75K |   25 |     13478 |    159 |   17 |  36.0 | 420.0 |
| 8_mr_megahit |  0.710 |    319000 | 579.38K |    8 |     17621 |    224 |   23 |  51.7 | 586.0 |

Table: statOtherAnchors.md

| Name                    |    N50 |     Sum |    # |
| ----------------------- | -----: | ------: | ---: |
| Genome                  | 580076 | 580.08K |    1 |
| Paralogs                |   1567 |  11.53K |    8 |
| repetitive              |    499 |  16.88K |   53 |
| 7_merge_anchors.anchors |  54929 | 556.87K |   16 |
| spades.contig           | 163847 |    581K |   38 |
| spades.scaffold         | 163847 |    581K |   38 |
| mr_spades.contig        | 580506 | 580.63K |    2 |
| mr_spades.scaffold      | 580506 | 580.63K |    2 |
| megahit.contig          |  42359 | 578.23K |   49 |
| mr_megahit.contig       | 319186 |    597K |   39 |

Table: statFinal

* Assembly quality by QUAST

| Assembly       | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik   |        18 |  135827 | 619203 |  79554 |     0 | 95.93 | 1.112 |   0.00 |   32.00 |       3.23 |
| merge_mr_multik |       16 |  179608 | 554898 |  54917 |     0 | 95.64 | 1.000 |   0.00 |   28.66 |       3.06 |
| merge_anchors  |        16 |  179741 | 556866 |  54929 |     0 | 95.95 | 1.000 |   0.00 |   31.45 |       3.23 |
| spades         |        35 |  236302 | 580752 | 163847 |     1 | 99.06 | 1.002 |   0.00 |  222.74 |      66.49 |
| mr_spades      |         2 |  580506 | 580632 | 580506 |     0 | 100.0 | 1.001 |   0.00 |  300.77 |      91.64 |
| megahit        |        49 |   97809 | 578231 |  42359 |     6 | 97.49 | 1.003 |   0.00 |   99.76 |      22.03 |
| mr_megahit     |        39 |  319186 | 597002 | 319186 |     3 | 99.88 | 1.025 |   0.00 |  318.39 |      94.96 |

Table: statQuast
