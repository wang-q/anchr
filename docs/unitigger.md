# Comparisons of unitiggers

```shell
cd ~/data/anchr/mg1655

Bifrost build -t 4 -k 31 --clip-tips --del-isolated --fasta -s 2_illumina/Q25L60/pe.cor.fa.gz -o Bf-31

Bifrost build -t 4 -k 81 --clip-tips --del-isolated --fasta -s 2_illumina/Q25L60/pe.cor.fa.gz -o Bf-81

bcalm -in 2_illumina/Q25L60/pe.cor.fa.gz \
    -kmer-size 31 -abundance-min 3 -verbose 0 \
    -nb-cores 4 -out bcalm-31

bcalm -in 2_illumina/Q25L60/pe.cor.fa.gz \
    -kmer-size 81 -abundance-min 3 -verbose 0 \
    -nb-cores 4 -out bcalm-81

tadpole.sh \
    in=2_illumina/Q25L60/pe.cor.fa.gz \
    out=tadpole-31.fasta \
    threads=4 \
    k=31 \
    overwrite

tadpole.sh \
    in=2_illumina/Q25L60/pe.cor.fa.gz \
    out=tadpole-81.fasta \
    threads=4 \
    k=81 \
    overwrite

# Add masurca to $PATH
export PATH="$(readlink -f "$(which masurca)" | xargs dirname):$PATH"
create_k_unitigs_large_k \
    -c 30 -t 4 \
    -m 31 -n 4641652 -l 31 -f 0.000001 \
    <(gzip -dcf 2_illumina/Q25L60/pe.cor.fa.gz) \
    > superreads-31.fasta

create_k_unitigs_large_k \
    -c 80 -t 4 \
    -m 81 -n 4641652 -l 81 -f 0.000001 \
    <(gzip -dcf 2_illumina/Q25L60/pe.cor.fa.gz) \
    > superreads-81.fasta


hnsm n50 -S -C -g 4641652 Bf-31.fasta
hnsm n50 -S -C -g 4641652 Bf-81.fasta

hnsm n50 -S -C -g 4641652 bcalm-31.unitigs.fa
hnsm n50 -S -C -g 4641652 bcalm-81.unitigs.fa

hnsm n50 -S -C -g 4641652 tadpole-31.fasta
hnsm n50 -S -C -g 4641652 tadpole-81.fasta

hnsm n50 -S -C -g 4641652 superreads-31.fasta
hnsm n50 -S -C -g 4641652 superreads-81.fasta

```
