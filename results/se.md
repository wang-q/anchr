# Single End

[TOC levels=1-3]: # ""

- [Single End](#single-end)
- [SE](#se)
  - [se: reference](#se-reference)
  - [se: download](#se-download)
  - [se: template](#se-template)
  - [se: run](#se-run)
- [SE2](#se2)
  - [se2: reference](#se2-reference)
  - [se2: download](#se2-download)
  - [se2: template](#se2-template)
  - [se2: run](#se2-run)
- [Results](#results)
  - [Reads](#reads)
  - [Anchors](#anchors)
  - [Comparison](#comparison)


# SE

## se: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/se/1_genome
cd ~/data/anchr/se/1_genome

cp ../../mg1655/1_genome/genome.fa .
cp ../../mg1655/1_genome/paralogs.fa .

```

## se: download

* Illumina

```shell script
mkdir -p ~/data/anchr/se/2_illumina
cd ~/data/anchr/se/2_illumina

cp ../../mg1655/2_illumina/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz R1.fq.gz

```

## se: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=se

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --se \
    --genome 4641652 \
    --parallel 24 \
    --xmx 80g \
    \
    --trim "--dedupe --tile --cutoff 5 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --cov "40 80" \
    --statp 2 \
    --redo \
    \
    --extend

```

## se: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=se

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

#bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

```

# SE2

## se2: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/se2/1_genome
cd ~/data/anchr/se2/1_genome

cp ../../mg1655/1_genome/genome.fa .
cp ../../mg1655/1_genome/paralogs.fa .

```

## se2: download

* Illumina

```shell script
mkdir -p ~/data/anchr/se2/2_illumina
cd ~/data/anchr/se2/2_illumina

cp ../../mg1655/2_illumina/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz R1.fq.gz

```

## se2: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=se2

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --se \
    --genome 4641652 \
    --parallel 24 \
    --xmx 80g \
    \
    --trim "--dedupe --tile --cutoff 5 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --cov "40 80" \
    --statp 2 \
    --redo \
    \
    --extend

```

## se2: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=se2

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

#bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh

```

# Results

## Reads

Table: statReads

| Name       |     N50 |     Sum |       # |
|:-----------|--------:|--------:|--------:|
| Genome     | 4641652 | 4641652 |       1 |
| Paralogs   |    1937 |  187300 |     106 |
| Illumina.R |     151 | 865.15M | 5729470 |
| trim.R     |     150 | 624.87M | 4416799 |
| Q25L60     |     150 | 585.96M | 4261791 |
| Q30L60     |     139 | 506.27M | 4004763 |

| Name       |     N50 |     Sum |       # |
|:-----------|--------:|--------:|--------:|
| Genome     | 4641652 | 4641652 |       1 |
| Paralogs   |    1937 |  187300 |     106 |
| Illumina.R |     151 | 865.15M | 5729470 |
| trim.R     |     144 | 616.67M | 4616321 |
| Q25L60     |     138 | 557.07M | 4398522 |
| Q30L60     |     117 | 448.71M | 4060083 |

Table: statTrimReads

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 151 | 717.62M | 4752462 |
| filteredbytile | 151 | 696.69M | 4613843 |
| highpass       | 151 | 671.82M | 4449143 |
| trim           | 150 |  624.9M | 4417025 |
| filter         | 150 | 624.87M | 4416799 |
| R1             | 150 | 624.87M | 4416799 |

| Name           | N50 |     Sum |       # |
|:---------------|----:|--------:|--------:|
| clumpify       | 151 | 778.63M | 5156497 |
| filteredbytile | 151 | 753.91M | 4992807 |
| highpass       | 151 | 703.74M | 4660551 |
| trim           | 144 | 616.69M | 4616521 |
| filter         | 144 | 616.67M | 4616321 |
| R1             | 144 | 616.67M | 4616321 |

Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 135.6 |  125.9 |    7.16% | "31" | 4.64M | 4.58M |     0.99 | 0:01'20'' |
| Q25L60.R | 127.2 |  121.5 |    4.48% | "31" | 4.64M | 4.56M |     0.98 | 0:01'19'' |
| Q30L60.R | 109.9 |  107.0 |    2.67% | "31" | 4.64M | 4.55M |     0.98 | 0:01'11'' |

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 133.8 |  121.5 |    9.25% | "31" | 4.64M | 4.57M |     0.99 | 0:01'21'' |
| Q25L60.R | 121.0 |  115.2 |    4.74% | "31" | 4.64M | 4.56M |     0.98 | 0:01'17'' |
| Q30L60.R |  97.5 |   95.5 |    2.13% | "31" | 4.64M | 4.55M |     0.98 | 0:01'07'' |

## Anchors


Table: statUnitigsAnchors.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  96.68% |      9871 | 4.44M | 670 |        81 |  97.89K | 1496 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  96.66% |      9887 | 4.44M | 650 |        78 |  93.27K | 1427 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  96.37% |     10018 | 4.43M | 680 |        75 |  89.15K | 1469 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  94.91% |      7136 | 4.38M | 883 |        69 | 106.02K | 1849 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:35 |
| Q25L60X40P000 |   40.0 |  98.10% |     16618 | 4.46M | 457 |        80 |  73.61K | 1174 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q25L60X40P001 |   40.0 |  98.11% |     14596 | 4.48M | 462 |        90 |  83.58K | 1218 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q25L60X40P002 |   40.0 |  98.26% |     15219 | 4.49M | 474 |        82 |  78.11K | 1235 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:40 |
| Q25L60X80P000 |   80.0 |  97.61% |     14054 | 4.47M | 511 |        72 |  70.46K | 1223 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:11 |   0:00:39 |
| Q30L60X40P000 |   40.0 |  98.85% |     17944 | 4.39M | 412 |       530 | 106.19K | 1319 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:42 |
| Q30L60X40P001 |   40.0 |  98.78% |     15822 | 4.26M | 411 |       670 | 110.87K | 1266 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:42 |
| Q30L60X80P000 |   80.0 |  98.79% |     21605 | 4.43M | 335 |       174 |  77.09K | 1086 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:10 |   0:00:45 |


| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.71% |     13078 | 4.48M | 548 |        84 |  94.38K | 1422 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:38 |
| Q0L0X40P001   |   40.0 |  97.52% |     12485 | 4.46M | 569 |        97 | 101.76K | 1429 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q0L0X40P002   |   40.0 |  97.53% |     11801 | 4.47M | 574 |        80 |  92.92K | 1450 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:40 |
| Q0L0X80P000   |   80.0 |  96.46% |      9178 | 4.45M | 694 |        67 |  87.64K | 1573 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  98.42% |     15079 |  4.4M | 442 |       328 | 108.36K | 1339 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:41 |
| Q25L60X40P001 |   40.0 |  98.45% |     15550 | 4.42M | 427 |       516 | 119.88K | 1416 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:42 |
| Q25L60X80P000 |   80.0 |  98.10% |     16911 | 4.48M | 417 |       269 |  87.62K | 1156 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:09 |   0:00:42 |
| Q30L60X40P000 |   40.0 |  98.78% |     13759 | 3.98M | 470 |       747 | 184.74K | 1585 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:42 |
| Q30L60X40P001 |   40.0 |  98.73% |     12818 |  3.9M | 474 |       879 | 203.27K | 1524 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:43 |
| Q30L60X80P000 |   80.0 |  98.83% |     19657 |  4.2M | 325 |       108 |  76.02K | 1255 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:01:04 |   0:00:47 |

Table: statMergeAnchors.md

| Name                    | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:------------------------|--------:|----------:|------:|----:|----------:|--------:|---:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors         |  97.19% |     40612 | 4.51M | 187 |      1135 | 112.47K | 89 |  127.0 | 4.0 |  38.3 | 208.5 |   0:00:37 |
| 7_merge_unitigs_anchors |  99.04% |     40794 | 4.52M | 184 |      1135 | 112.47K | 89 |  126.0 | 4.0 |  38.0 | 207.0 |   0:01:10 |

| Name                    | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors         |  97.01% |     47142 |  4.5M | 190 |      1116 | 135.73K | 122 |  122.0 | 4.0 |  36.7 | 201.0 |   0:00:37 |
| 7_merge_unitigs_anchors |  98.96% |     47183 | 4.51M | 191 |      1134 | 134.71K | 121 |  122.0 | 4.0 |  36.7 | 201.0 |   0:01:09 |


Table: statOtherAnchors.md

| Name       | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-----------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|----------:|
| 8_megahit  |  98.52% |     56469 | 4.31M | 150 |      1133 | 34.82K | 271 |  127.0 | 4.0 |  38.3 | 208.5 |   0:00:37 |
| 8_platanus |  98.08% |     26990 | 4.53M | 258 |        59 | 28.41K | 516 |  127.0 | 4.0 |  38.3 | 208.5 |   0:00:39 |

| Name       | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-----------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|----------:|
| 8_megahit  |  98.62% |     48888 | 4.25M | 172 |      1054 | 35.85K | 290 |  122.0 | 4.0 |  36.7 | 201.0 |   0:00:36 |
| 8_platanus |  97.78% |     21800 | 4.51M | 353 |        57 |  38.5K | 708 |  122.0 | 5.0 |  35.7 | 205.5 |   0:00:37 |

## Comparison

| Name                    |     N50 |     Sum |   # |
|:------------------------|--------:|--------:|----:|
| Genome                  | 4641652 | 4641652 |   1 |
| Paralogs                |    1937 |  187300 | 106 |
| 7_merge_anchors.anchors |   40612 | 4505933 | 187 |
| 7_merge_anchors.others  |    1135 |  112469 |  89 |
| glue_anchors            |   40612 | 4505933 | 187 |
| fill_anchors            |   40612 | 4505933 | 187 |
| spades.non-contained    |       0 |       0 |   0 |
| megahit.contig          |   67336 | 4562142 | 149 |
| megahit.non-contained   |   67336 | 4549720 | 122 |
| platanus.contig         |   27038 | 4642481 | 653 |
| platanus.non-contained  |   27762 | 4555356 | 258 |

| Name                    |     N50 |     Sum |   # |
|:------------------------|--------:|--------:|----:|
| Genome                  | 4641652 | 4641652 |   1 |
| Paralogs                |    1937 |  187300 | 106 |
| 7_merge_anchors.anchors |   47142 | 4497778 | 190 |
| 7_merge_anchors.others  |    1116 |  135725 | 122 |
| glue_anchors            |   47142 | 4497778 | 190 |
| fill_anchors            |   47142 | 4497778 | 190 |
| spades.non-contained    |       0 |       0 |   0 |
| megahit.contig          |   63432 | 4562785 | 157 |
| megahit.non-contained   |   63432 | 4550225 | 125 |
| platanus.contig         |   21710 | 4658040 | 879 |
| platanus.non-contained  |   21925 | 4553224 | 355 |

