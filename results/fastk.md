# FastK, GENESCOPE.FK, and MERQURY.FK

[FastK](https://github.com/thegenemyers/FASTK)

[GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

[MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

## Installation

```shell
brew install --HEAD wang-q/tap/fastk
brew install --HEAD wang-q/tap/merquryfk

# GeneScope.FK is bundled

```

## GeneScope

```shell
cd ~/data/anchr/mg1655/2_illumina

time FastK -v -t1 -k21 R1.fq.gz R2.fq.gz -NTable-21
#real	0m50.853s
#user	1m27.788s
#sys	0m4.377s

time FastK -v -t1 -k51 R1.fq.gz R2.fq.gz -NTable-51
#real	1m3.061s
#user	1m36.951s
#sys	0m6.808s

Histex -G Table-21 | Rscript ~/Scripts/rust/anchr/templates/genescopefk.R -k 21 -p 1 -o GeneScope-21

Histex -G Table-51 | Rscript ~/Scripts/rust/anchr/templates/genescopefk.R -k 51 -p 1 -o GeneScope-51

Fastrm Table-21 Table-51

for K in 21 51; do
    echo ${K}
    cat GeneScope-${K}/model.txt |
        grep '^kmercov' |
        tr -s ' ' '\t' |
        cut -f 2 |
        perl -nl -e 'printf qq{%.1f\n}, $_'
done |
    paste - - |
    ( echo -e 'K\tCov' && cat ) |
    mlr --itsv --omd cat

for K in 21 51; do
    cat GeneScope-${K}/summary.txt |
        sed '1,6 d' |
        sed '1 s/^/K\t/' |
        sed "2 s/^/${K}\t/" |
        sed "3,7 s/^/\t/" |
        perl -nlp -e 's/\s{2,}/\t/g; s/\s+$//g;'
done |
    tsv-uniq |
    mlr --itsv --omd cat

ll |
    grep Table |
    tr -s ' ' '\t' |
    cut -f 5,9 |
    sed 's/\*$//' |
    ( echo -e 'Size\tName' && cat ) |
    mlr --itsv --omd cat

```

| K   | Cov   |
|-----|-------|
| 21  | 299.7 |
| 51  | 223.4 |

| K   | property              | min          | max          |
|-----|-----------------------|--------------|--------------|
| 21  | Homozygous (a)        | 100%         | 100%         |
|     | Genome Haploid Length | NA bp        | 4,478,083 bp |
|     | Genome Repeat Length  | 136,755 bp   | 136,883 bp   |
|     | Genome Unique Length  | 4,339,242 bp | 4,343,288 bp |
|     | Model Fit             | 97.2749%     | 97.3934%     |
|     | Read Error Rate       | 0.531821%    | 0.531821%    |
| 51  | Homozygous (a)        | 100%         | 100%         |
|     | Genome Haploid Length | NA bp        | 4,385,263 bp |
|     | Genome Repeat Length  | 91,444 bp    | 91,569 bp    |
|     | Genome Unique Length  | 4,290,813 bp | 4,296,704 bp |
|     | Model Fit             | 97.3222%     | 97.6265%     |
|     | Read Error Rate       | 0.326106%    | 0.326106%    |

| Size      | Name             |
|-----------|------------------|
| 145354802 | .Table-21.ktab.1 |
| 128966832 | .Table-21.ktab.2 |
| 138372242 | .Table-21.ktab.3 |
| 133722077 | .Table-21.ktab.4 |
| 403680936 | .Table-51.ktab.1 |
| 383925912 | .Table-51.ktab.2 |
| 424657164 | .Table-51.ktab.3 |
| 394298388 | .Table-51.ktab.4 |
| 262164    | Table-21.hist    |
| 134217744 | Table-21.ktab    |
| 262164    | Table-51.hist    |
| 134217744 | Table-51.ktab    |


