# FastK, GENESCOPE.FK, and MERQURY.FK

[FastK](https://github.com/thegenemyers/FASTK)

[GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

[MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

## Installation

```shell
brew install --HEAD wang-q/tap/fastk
brew install --HEAD wang-q/tap/merquryfk

# GeneScope.FK is bundled

parallel -j 1 -k --line-buffer '
    Rscript -e '\'' if (!requireNamespace("{}", quietly = FALSE)) { install.packages("{}", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN") } '\''
    ' ::: \
        argparse minpack.lm \
        ggplot2 scales viridis

```

## GeneScope

```shell
cd ~/data/anchr/g37/ena/

time FastK -v -t1 -k21 *.fastq.gz -NTable-21
#real	0m50.853s
#user	1m27.788s
#sys	0m4.377s

time FastK -v -t1 -k51 R1.fq.gz R2.fq.gz -NTable-51
#real	1m3.061s
#user	1m36.951s
#sys	0m6.808s

Histex -G Table-21 | Rscript ~/Scripts/rust/anchr/templates/genescopefk.R -k 21 -p 1 -o GeneScope-21

Histex -G Table-51 | Rscript ~/Scripts/rust/anchr/templates/genescopefk.R -k 51 -p 1 -o GeneScope-51

# disk usages
ll |
    grep Table |
    tr -s ' ' '\t' |
    cut -f 5,9 |
    sed 's/\*$//' |
    ( echo -e 'Size\tName' && cat ) |
    mlr --itsv --omd cat

Fastrm Table-21 Table-51

# Reports
for K in 21 51; do
    COV=$(
        cat GeneScope-${K}/model.txt |
            grep '^kmercov' |
            tr -s ' ' '\t' |
            cut -f 2 |
            perl -nl -e 'printf qq{%.1f\n}, $_'
    )

    cat GeneScope-${K}/summary.txt |
        sed '1,6 d' |
        sed '1 s/^/K\t/' |
        sed "2 s/^/${K}\t/" |
        sed "3,7 s/^/\t/" |
        perl -nlp -e 's/\s{2,}/\t/g; s/\s+$//g;' |
        perl -nla -F'\t' -e '
            @fields = map {/\bNA\b/ ? q{} : $_ } @F;        # Remove NA fields
            $fields[2] = q{} if $fields[2] eq $fields[3];   # Remove identical fields
            print join qq{\t}, @fields
        '

    printf "\tKmer Cov\t\t${COV}\n"
done |
    keep-header -- grep -v '^K' |
    mlr --itsv --omd cat

```

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

| K   | property              | min          | max          |
|-----|-----------------------|--------------|--------------|
| 21  | Homozygous (a)        |              | 100%         |
|     | Genome Haploid Length |              | 4,478,083 bp |
|     | Genome Repeat Length  | 136,755 bp   | 136,883 bp   |
|     | Genome Unique Length  | 4,339,242 bp | 4,343,288 bp |
|     | Model Fit             | 97.2749%     | 97.3934%     |
|     | Read Error Rate       |              | 0.531821%    |
|     | Kmer Cov              |              | 299.7        |
| 51  | Homozygous (a)        |              | 100%         |
|     | Genome Haploid Length |              | 4,385,263 bp |
|     | Genome Repeat Length  | 91,444 bp    | 91,569 bp    |
|     | Genome Unique Length  | 4,290,813 bp | 4,296,704 bp |
|     | Model Fit             | 97.3222%     | 97.6265%     |
|     | Read Error Rate       |              | 0.326106%    |
|     | Kmer Cov              |              | 223.4        |

## Merqury

> these tables must be produced with the option -t1 set

```shell
cd ~/data/anchr/mg1655/2_illumina

mkdir "Merqury"

# KatGC
FastK -v -t1 -k21 R*.fq.gz -NTable-21

FastK -v -t1 -k51 R*.fq.gz -NTable-51

KatGC -T4 -x1.9 -s Table-21 Merqury/KatGC-21
KatGC -T4 -x1.9 -s Table-51 Merqury/KatGC-51

Fastrm Table-21 Table-51

# KatComp
FastK -v -t1 -k21 R1.fq.gz -NR1-21
FastK -v -t1 -k21 R2.fq.gz -NR2-21

KatComp -T4 -s R1-21 R2-21 Merqury/KatComp-21

Fastrm R1-21 R2-21

```

## Find repetitive regions

```shell
cd ~/data/anchr/mg1655/1_genome/

faops size genome.fa > chr.sizes

FastK -v -p -k21 genome

cat chr.sizes |
    number-lines |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        Profex genome {1} |
            sed "1,2 d" |
            perl -nl -e '\''/(\d+).+(\d+)/ and print qq{$1\t$2}'\'' |
            tsv-filter --ge 2:2 |
            cut -f 1 |
            sed "s/^/{2}:/" |
            spanr cover stdin |
            spanr span --op fill -n 2 stdin |
            spanr span --op excise -n 100 stdin |
            spanr span --op fill -n 10 stdin
    ' |
    spanr combine stdin \
    > repetitive.yml

spanr stat chr.sizes repetitive.yml
#chr,chrLength,size,coverage
#NC_000913,4641652,91989,0.0198
#all,4641652,91989,0.0198

faops size repetitives.fa | tsv-summarize --sum 2
#110138

faops size paralogs.fa | tsv-summarize --sum 2
#187300

```
