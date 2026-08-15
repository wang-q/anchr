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
| Q25L60.R | 169.4 |  160.5 |    5.26% | "24" | 580.08K | 576.79K |     0.99 | 0:00:11 |
| Q30L60.R | 163.9 |  157.3 |    4.01% | "24" | 580.08K |  578.3K |     1.00 | 0:00:10 |

Table: statQuorum

| Name          | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| Q0L0X40P000   |   40.0 |  0.960 |     31796 |  556.4K |   34 |     28902 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P001   |   40.0 |  0.960 |     39117 | 556.69K |   31 |     28179 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X40P002   |   40.0 |  0.957 |     22406 |  556.1K |   46 |     31784 |     39 |    6 |   7.0 | 114.0 |
| Q0L0X80P000   |   80.0 |  0.953 |     19609 | 554.92K |   55 |     36497 |     78 |   10 |  16.0 | 216.0 |
| Q0L0X80P001   |   80.0 |  0.952 |     19364 | 554.53K |   61 |     38693 |     78 |   10 |  16.0 | 216.0 |
| Q25L60X40P000 |   40.0 |  0.958 |     34062 | 556.64K |   39 |     28935 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P001 |   40.0 |  0.959 |     26560 | 556.15K |   42 |     30516 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X40P002 |   40.0 |  0.960 |     31991 | 556.06K |   34 |     29897 |     39 |    6 |   7.0 | 114.0 |
| Q25L60X80P000 |   80.0 |  0.953 |     15777 | 555.95K |   57 |     36840 |     78 |   10 |  16.0 | 216.0 |
| Q25L60X80P001 |   80.0 |  0.956 |     22431 | 556.07K |   47 |     36356 |     79 |   10 |  16.3 | 218.0 |
| Q30L60X40P000 |   40.0 |  0.959 |     23060 | 556.92K |   40 |     30302 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P001 |   40.0 |  0.960 |     31356 | 556.07K |   38 |     30698 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X40P002 |   40.0 |  0.959 |     45204 | 556.36K |   35 |     28947 |     39 |    6 |   7.0 | 114.0 |
| Q30L60X80P000 |   80.0 |  0.953 |     20810 | 555.69K |   59 |     37352 |     78 |   10 |  16.0 | 216.0 |

Table: statUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.682 |     39116 | 555.91K |   28 |     28142 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.679 |     35869 | 555.49K |   35 |     28529 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.682 |     33952 | 556.79K |   27 |     27070 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.677 |     32672 | 555.88K |   44 |     31191 |     78 |   10 |  16.0 | 216.0 |
| MRX80P001 |   80.0 |  0.677 |     36626 | 555.32K |   41 |     31821 |     78 |    9 |  17.0 | 210.0 |

Table: statMRUnitigsMultik.md

| Name                      | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.683 |     54945 | 557.01K |   16 |       577 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_mr_unitigs_multik |  0.681 |     54916 | 555.45K |   16 |       413 |    226 |   22 |  53.3 | 584.0 |
| 7_merge_unitigs_multik    |  0.682 |     54945 | 556.67K |   18 |       526 |    226 |   22 |  53.3 | 584.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |     Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ------: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.988 |     83431 | 573.62K |    9 |      7382 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_spades  |  0.710 |    580031 | 580.03K |    1 |       601 |    225 |   22 |  53.0 | 582.0 |
| 8_megahit    |  0.976 |     47946 | 563.64K |   21 |     13611 |    159 |   16 |  37.0 | 414.0 |
| 8_mr_megahit |  0.710 |    319000 | 580.18K |   10 |     19598 |    223 |   23 |  51.3 | 584.0 |

Table: statOtherAnchors.md

| Name                    |    N50 |     Sum |    # |
| ----------------------- | -----: | ------: | ---: |
| Genome                  | 580076 | 580.08K |    1 |
| Paralogs                |   1567 |  11.53K |    8 |
| repetitive              |    499 |  16.88K |   53 |
| 7_merge_anchors.anchors |  54945 | 557.01K |   16 |
| spades.contig           | 163847 |    581K |   38 |
| spades.scaffold         | 163847 |    581K |   38 |
| mr_spades.contig        | 580506 | 580.63K |    2 |
| mr_spades.scaffold      | 580506 | 580.63K |    2 |
| megahit.contig          |  42958 | 577.25K |   48 |
| mr_megahit.contig       | 319187 | 599.78K |   46 |

Table: statFinal

* Assembly quality by QUAST

| Assembly       | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik   |        18 |  179605 | 556672 |  54945 |     0 | 95.92 | 1.000 |   0.00 |   32.17 |       2.88 |
| merge_mr_multik |       16 |  179794 | 555445 |  54916 |     0 | 95.73 | 1.000 |   0.00 |   33.68 |       3.78 |
| merge_anchors  |        16 |  179802 | 557007 |  54945 |     0 | 95.99 | 1.000 |   0.00 |   37.18 |       3.95 |
| spades         |        35 |  236302 | 580752 | 163847 |     1 | 99.06 | 1.002 |   0.00 |  222.74 |      66.49 |
| mr_spades      |         2 |  580506 | 580632 | 580506 |     0 | 100.0 | 1.001 |   0.00 |  300.77 |      91.64 |
| megahit        |        48 |  179705 | 577254 |  42958 |     5 | 97.61 | 1.002 |   0.00 |   93.92 |      21.85 |
| mr_megahit     |        46 |  319187 | 599775 | 319187 |     2 | 99.81 | 1.027 |   0.00 |  313.64 |      90.36 |

Table: statQuast

### mg1655: multik vs bcalm 对照（2026-08-15）

同一批 `6_down_sampling` reads（MRX40P000/P001/P002 + MRX80P000/P001，
5 组），唯一变量是 unitigger：`asm multik` vs `asm unitig`（自研 BCALM
语义）vs 外部 bcalm（每 k 31..81 独立 unitigs + `asm olc --unitigs` 跨 k
合并），下游 anchor/OLC 合并参数完全相同。

#### unitig / anchor 阶段（N50 / 条数 / Sum）

| 组        | multik unitigs | unitig unitigs | bcalm unitigs |
| --------- | -------------: | -------------: | ------------: |
| MRX40P000 |  21238 / 1455  |  61235 / 131   |  56477 / 142  |
| MRX40P001 |  19872 / 1455  |  61235 / 129   |  54908 / 147  |
| MRX40P002 |  21238 / 1455  |  59716 / 130   |  57961 / 143  |
| MRX80P000 |  19331 / 1525  |  53834 / 158   |  42412 / 194  |
| MRX80P001 |  20068 / 1525  |  50795 / 163   |  42234 / 195  |

multik 的 unitigs 比 bcalm/unitig 短约 2.5-3 倍、条数多约 9-11 倍；
`asm unitig`（自研）与 bcalm 等价且略优；两链 anchor Sum ≈ 4.53M，比
multik 的 4.47M 更接近基因组。

#### 最终合并（`asm olc --unitigs`，5 组合并）与 QUAST

| 指标      | unitig 链（自研） | bcalm 链（现代） | multik 链（现代，同输入） | legacy bcalm | legacy merge_anchors |
| --------- | ----------------: | ---------------: | ------------------------: | ------------: | --------------------: |
| # contigs |               102 |              108 |                       317 |           101 |                   103 |
| N50       |           105,719 |           95,478 |                    23,403 |        78,596 |                95,484 |
| Largest   |          246,019 |          202,937 |                    73,030 |       174,107 |               204,605 |
| # mis     |                 1 |                1 |                         4 |             0 |                     0 |
| GF%       |            97.77 |            97.77 |                     96.44 |         97.26 |                 97.67 |
| Dup       |           1.000 |           1.001 |                      1.003 |         1.000 |                  1.000 |
| mm/100k   |            0.27 |            0.04 |                      0.18 |          0.00 |                  0.40 |
| indel/100k|            0.29 |            0.23 |                      0.02 |          0.02 |                  0.13 |

结论：N50 差距基本全部来自 unitig 阶段——`asm unitig` 链（105.7K）与
bcalm 链（95.5K）都追平并超过 legacy（95.5K），multik 只有 23.4K；mis
从 multik 的 4 降到 unitig/bcalm 链的 1（legacy 0）；`asm unitig` 与外部
bcalm 等价（自研可替代外部依赖）。

#### Dup 修复（2026-08-15 追加）

原两链 Dup 1.07-1.08 偏高：anchor 阶段无近似重复 contig，重复由
`asm olc --unitigs` 跨组合并产生（跨组 anchors 边界不一致，exact overlap
检测连不上，残留同区域不同边界的 contig 对）。`asm olc` 现在在 consensus
后增加近似 overlap 合并（`consensus::merge_overlapping_contigs`：31-mer
定位主导 offset + 头部锚定 + 重叠区 ≥99% 一致才合并，嵌合 contig 的
多块对齐被拒绝）：unitig 链 Dup **1.079 → 1.000**、bcalm 链 **1.068 →
1.001**，GF 97.77% 不变，unitig 链 N50 105.7K → 110.2K（93 contigs），
bcalm 链 N50 95.5K → 88.1K（102 contigs）。mis 仍为 1（contig_26
relocation，嵌合修复属另一任务）。

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
    --bwa "Q25L60" \
    \
    --cov "40 80" \
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
| Q0L0X40P000   |   40.0 |  0.964 |     13614 | 4.48M |  545 |    125438 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P001   |   40.0 |  0.964 |     13663 | 4.48M |  546 |    127750 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X40P002   |   40.0 |  0.965 |     14305 | 4.48M |  511 |    119182 |     39 |    7 |   6.0 | 120.0 |
| Q0L0X80P000   |   80.0 |  0.959 |     10077 | 4.48M |  761 |    152824 |     79 |   12 |  14.3 | 230.0 |
| Q0L0X80P001   |   80.0 |  0.959 |      9538 | 4.48M |  756 |    149797 |     79 |   12 |  14.3 | 230.0 |
| Q0L0X80P002   |   80.0 |  0.959 |      9004 | 4.47M |  769 |    156952 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X40P000 |   40.0 |  0.969 |     14007 | 4.47M |  513 |    115498 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X40P001 |   40.0 |  0.970 |     14179 | 4.48M |  508 |    116620 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X40P002 |   40.0 |  0.969 |     13648 | 4.47M |  513 |    118120 |     39 |    7 |   6.0 | 120.0 |
| Q25L60X80P000 |   80.0 |  0.968 |     12815 | 4.48M |  582 |    118105 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P001 |   80.0 |  0.968 |     12039 | 4.48M |  596 |    122840 |     79 |   12 |  14.3 | 230.0 |
| Q25L60X80P002 |   80.0 |  0.969 |     12705 | 4.48M |  583 |    121945 |     79 |   12 |  14.3 | 230.0 |
| Q30L60X40P000 |   40.0 |  0.967 |     11397 | 4.44M |  603 |    114763 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X40P001 |   40.0 |  0.968 |     11784 | 4.45M |  590 |    113624 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X40P002 |   40.0 |  0.966 |     11228 | 4.44M |  598 |    114871 |     39 |    8 |   5.0 | 126.0 |
| Q30L60X80P000 |   80.0 |  0.972 |     13557 | 4.48M |  530 |    106155 |     79 |   14 |  12.3 | 242.0 |
| Q30L60X80P001 |   80.0 |  0.972 |     13097 | 4.48M |  529 |    107238 |     79 |   14 |  12.3 | 242.0 |

Table: statUnitigsMultik.md

| Name      | CovCor | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| --------- | -----: | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| MRX40P000 |   40.0 |  0.481 |     21714 | 4.47M |  351 |    123635 |     39 |    6 |   7.0 | 114.0 |
| MRX40P001 |   40.0 |  0.481 |     22510 | 4.47M |  352 |    123687 |     39 |    6 |   7.0 | 114.0 |
| MRX40P002 |   40.0 |  0.481 |     21169 | 4.47M |  356 |    122936 |     39 |    6 |   7.0 | 114.0 |
| MRX80P000 |   80.0 |  0.480 |     20006 | 4.48M |  373 |    115009 |     78 |   10 |  16.0 | 216.0 |
| MRX80P001 |   80.0 |  0.480 |     19961 | 4.48M |  375 |    112524 |     78 |   10 |  16.0 | 216.0 |
| MRX80P002 |   80.0 |  0.480 |     20779 | 4.48M |  373 |    113540 |     78 |   10 |  16.0 | 216.0 |

Table: statMRUnitigsMultik.md

| Name                      | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------------------- | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 7_merge_anchors           |  0.479 |     28019 | 4.61M |  278 |      5841 |    381 |   38 |  89.0 | 990.0 |
| 7_merge_mr_unitigs_multik |  0.477 |     23594 | 4.53M |  317 |      4097 |    380 |   38 |  88.7 | 988.0 |
| 7_merge_unitigs_multik    |  0.479 |     23956 | 4.73M |  318 |      6618 |    380 |   38 |  88.7 | 988.0 |

Table: statMergeAnchors.md

| Name         | Mapped | N50Anchor |   Sum |    # | SumOthers | median |  MAD | lower | upper |
| ------------ | -----: | --------: | ----: | ---: | --------: | -----: | ---: | ----: | ----: |
| 8_spades     |  0.989 |     82215 | 4.55M |  121 |     23296 |    292 |   35 |  62.3 | 794.0 |
| 8_mr_spades  |  0.500 |    145104 | 4.56M |   84 |     25914 |    382 |   37 |  90.3 | 986.0 |
| 8_megahit    |  0.985 |     35210 | 4.54M |  210 |     46173 |    291 |   36 |  61.0 | 798.0 |
| 8_mr_megahit |  0.502 |    126306 | 4.58M |   87 |     31206 |    382 |   37 |  90.3 | 986.0 |

Table: statOtherAnchors.md

| Name                    |     N50 |     Sum |    # |
| ----------------------- | ------: | ------: | ---: |
| Genome                  | 4641652 |   4.64M |    1 |
| Paralogs                |    1737 |    193K |  112 |
| repetitive              |    1265 | 124.42K |  155 |
| 7_merge_anchors.anchors |   28019 |   4.61M |  278 |
| spades.contig           |  125607 |   4.57M |  135 |
| spades.scaffold         |  132608 |   4.57M |  131 |
| mr_spades.contig        |  148607 |   4.59M |  152 |
| mr_spades.scaffold      |  148607 |   4.59M |  149 |
| megahit.contig          |   43891 |   4.59M |  273 |
| mr_megahit.contig       |  126312 |   4.61M |  138 |

Table: statFinal

* Assembly quality by QUAST

| Assembly       | # contigs | Largest |  Total |    N50 | # mis |   GF% |  Dup | N/100k | mm/100k | indel/100k |
| -------------- | --------: | ------: | -----: | -----: | ----: | ----: | ---: | -----: | ------: | ---------: |
| merge_multik   |       318 |   81358 | 4726538 |  23956 |    12 | 96.78 | 1.052 |   0.00 |    0.63 |       0.19 |
| merge_mr_multik |      317 |   73037 | 4527602 |  23594 |     5 | 96.49 | 1.011 |   0.00 |    0.24 |       0.02 |
| merge_anchors  |       278 |   83824 | 4606970 |  28019 |    12 | 96.84 | 1.025 |   0.00 |    0.63 |       0.20 |
| spades         |       125 |  224028 | 4572988 | 125607 |     0 | 98.47 | 1.000 |   0.00 |    0.74 |       0.17 |
| mr_spades      |       152 |  284843 | 4588125 | 148607 |     0 | 98.64 | 1.002 |   0.00 |    0.92 |       0.17 |
| megahit        |       273 |  175838 | 4588105 |  43891 |    92 | 98.26 | 1.004 |   0.00 |    3.73 |       0.50 |
| mr_megahit     |       138 |  311797 | 4611104 | 126312 |     2 | 98.89 | 1.004 |   0.00 |    2.13 |       0.26 |

Table: statQuast
