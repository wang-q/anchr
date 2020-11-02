# QA resources: genome, annotations, paralogs, and repeats

End users of [Anchr](https://github.com/wang-q/anchr) don't need to run the following codes. We use
these data just for quality assessments.

These steps require another project: [App::Egaz](https://github.com/wang-q/App-Egaz).

Paralogs detected here **may** overlap with transposons/retrotransposons.

## Download

```shell script
mkdir -p ~/data/anchr/assembly
cd ~/data/anchr/assembly

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/ \
    lambda/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/ \
    mg1655/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/001/723/505/GCF_001723505.1_ASM172350v1/ \
    dh5alpha/

```

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

for STRAIN in \
    lambda mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    ; do
    echo >&2 ${STRAIN};
    mkdir -p ${STRAIN}

    if [ ! -d ../assembly/${STRAIN} ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

    find ../assembly/${STRAIN}/ -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs gzip -dcf |
        faops filter -N -s stdin ${STRAIN}/genome.fa

done

```

## Prepare

```shell script
mkdir -p ~/data/anchr/paralogs/genomes
cd ~/data/anchr/paralogs/genomes

for STRAIN in \
    lambda mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    ; do
    if [ -d ${STRAIN} ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

    if [ ! -e ../ref/${STRAIN}/genome.fa ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

    egaz prepseq \
        ../ref/${STRAIN}/genome.fa -o ${STRAIN} \
        --repeatmasker '--parallel 16' -v
done

for STRAIN in Bper Cdif Cdip Cjej Ftul Hinf Lmon Lpne Ngon Nmen Sfle Vpar; do
    if [ -d ${STRAIN} ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

done

```

## Self-alignments

```shell script
cd ~/data/anchr/paralogs

egaz template \
    genomes/lambda \
    genomes/mg1655 genomes/dh5alpha \
    --self -o e_coli/ \
    --circos \
    --length 1000 --parallel 16 -v

bash e_coli/1_self.sh
bash e_coli/3_proc.sh
bash e_coli/4_circos.sh

```

```shell script
cd ~/data/anchr/paralogs

egaz template \
    genomes/Bcer genomes/Mabs genomes/Rsph genomes/Vcho \
    --self -o gage/ \
    --circos \
    --length 1000 --parallel 16 -v

bash gage/1_self.sh
bash gage/3_proc.sh
bash gage/4_circos.sh

```

