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

Fastrm R1 Table

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

ll |
    grep Table |
    tr -s ' ' '\t' |
    cut -f 5,9 |
    sed 's/\*$//' |
    ( echo -e 'Size\tName' && cat ) |
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


