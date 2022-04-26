# basecov hist

`basecov.txt` from `~/data/anchr/g37/4_unitigs_bcalm/Q30L60X80P000/anchor`

* The header is `#RefName	Pos	Coverage`
* Pos is 0-based

## Test data

```shell
gzip -dcf tests/G37/basecov.txt.gz | head -n 1
# #RefName        Pos     Coverage

gzip -dcf tests/G37/basecov.txt.gz |
    tsv-select -f 1 |
    tsv-uniq -H \
    > tests/G37/refname.txt

gzip -dcf tests/G37/basecov.txt.gz |
    grep -Fw -f <(head -n 2 tests/G37/refname.txt) \
    > tests/G37/basecov.1.txt

gzip -dcf tests/G37/basecov.txt.gz |
    grep -Fw -f <(head -n 11 tests/G37/refname.txt) \
    > tests/G37/basecov.10.txt

```

## bins

```shell
plotr hist --col 3 tests/G37/basecov.1.txt --device png -o tests/G37/hist.1.png

plotr hist --col 3 tests/G37/basecov.10.txt --device png -o tests/G37/hist.10.png

gzip -dcf tests/G37/basecov.txt.gz > tests/G37/basecov.txt.tmp
plotr hist --col 3 tests/G37/basecov.txt.tmp --device png -o tests/G37/hist.png

```

## hist

```shell
cat tests/G37/basecov.1.txt |
    tsv-summarize -H --group-by 3 --count |
    keep-header -- tsv-sort -k1,1n

cat tests/G37/basecov.10.txt |
    tsv-summarize -H --group-by 3 --count |
    keep-header -- tsv-sort -k1,1n

gzip -dcf tests/G37/basecov.txt.gz |
    tsv-summarize -H --group-by 3 --count |
    keep-header -- tsv-sort -k1,1n

```
