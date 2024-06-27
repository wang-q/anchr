# QA resources: genome, annotations, paralogs, and repeats

End users of [Anchr](https://github.com/wang-q/anchr) don't need to run the following codes. We use
these data just for quality assessments.

These steps require another project: [App::Egaz](https://github.com/wang-q/App-Egaz).

Paralogs detected here **may** overlap with transposons/retrotransposons.

## Download

### Model organisms

* *Mycoplasma genitalium* G37

    * Taxonomy ID: [243273](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=243273)
    * Ref. Assembly: [GCF_000027325.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000027325.1)

* *E. coli* str. K-12 substr. MG1655

    * Taxonomy ID: [511145](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145)
    * Assembly: [GCF_000005845.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)
    * Proportion of paralogs (> 1000 bp): 0.0325

* *E. coli* str. K-12 substr. DH5alpha

    * Taxonomy ID: [83333](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=83333)
    * Assembly: [GCF_001723515.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001723515.1)
    * Proportion of paralogs (> 1000 bp): 0.0342

```shell script
mkdir -p ~/data/anchr/assembly
cd ~/data/anchr/assembly

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/ \
    g37/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/ \
    mg1655/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/001/723/505/GCF_001723505.1_ASM172350v1/ \
    dh5alpha/

```

### GAGE-B

* *Bacillus cereus* ATCC 10987; 蜡样芽胞杆菌

    * Taxonomy ID: [222523](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=222523)
    * Assembly: [GCF_000008005.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000008005.1)
    * Proportion of paralogs (> 1000 bp): 0.0344

* *Mycobacterium abscessus* 6G-0125-R; 脓肿分枝杆菌

    * *Mycobacterium abscessus* ATCC 19977
        * Taxonomy ID: [561007](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=561007)
        * Assembly: [GCF_000069185.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000069185.1)
        * Proportion of paralogs (> 1000 bp): 0.0344

    * *Mycobacterium abscessus* 6G-0125-R
        * Assembly: [GCF_000270985.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000270985.1)
        * Illumina and 454, 40x

* *Rhodobacter sphaeroides* 2.4.1; 类球红细菌

    * Taxonomy ID: [272943](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=272943)
    * Assembly: [GCF_000012905.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000012905.2)
    * Proportion of paralogs (> 1000 bp): 0.0293

* *Vibrio cholerae* CP1032(5); 霍乱弧菌

    * *Vibrio cholerae* O1 biovar El Tor str. N16961
        * Taxonomy ID: [243277](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=243277)
        * Assembly: [GCF_000006745.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000006745.1)
        * Proportion of paralogs (> 1000 bp): 0.0216
    * *Vibrio cholerae* CP1032(5)
        * Assembly: [GCF_000279305.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000279305.1)

```shell script
cd ~/data/anchr/assembly

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/008/005/GCF_000008005.1_ASM800v1/ \
    Bcer/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/069/185/GCF_000069185.1_ASM6918v1 \
    Mabs/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/012/905/GCF_000012905.2_ASM1290v2 \
    Rsph/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/006/745/GCF_000006745.1_ASM674v1 \
    Vcho/

```

### Yeast

* *Saccharomyces cerevisiae* S288c
    * Taxonomy ID: [559292](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=559292)
    * Assembly: [GCF_000146045.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2)
    * Proportion of paralogs (> 1000 bp): 0.058

```shell script
cd ~/data/anchr/assembly

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/146/045/GCF_000146045.2_R64 \
    s288c/

```

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

for STRAIN in \
    g37 mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    s288c \
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
    g37 mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    s288c \
    ; do
    if [ -d ${STRAIN} ]; then
        echo >&2 "==> ${STRAIN} already be processed";
        continue;
    fi

    if [ ! -e ../../ref/${STRAIN}/genome.fa ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

    egaz prepseq \
        ../../ref/${STRAIN}/genome.fa -o ${STRAIN} \
        --repeatmasker '--parallel 8' -v
done

```

## Self-alignments

```shell script
cd ~/data/anchr/paralogs

egaz template \
    genomes/g37 \
    genomes/mg1655 genomes/dh5alpha \
    --self -o e_coli/ \
    --length 1000 --parallel 4 -v

bash e_coli/1_self.sh
bash e_coli/3_proc.sh

for STRAIN in \
    g37 mg1655 dh5alpha \
    ; do
    cat e_coli/Results/${STRAIN}/${STRAIN}.multi.fas |
        faops filter -N -d stdin stdout \
        > ../ref/${STRAIN}/paralogs.fa
done

```

```shell script
cd ~/data/anchr/paralogs

egaz template \
    genomes/Bcer genomes/Mabs genomes/Rsph genomes/Vcho \
    --self -o gage_b/ \
    --length 1000 --parallel 4 -v

bash gage_b/1_self.sh
bash gage_b/3_proc.sh

for STRAIN in \
    Bcer Mabs Rsph Vcho \
    ; do
    cat gage_b/Results/${STRAIN}/${STRAIN}.multi.fas |
        faops filter -N -d stdin stdout \
        > ../ref/${STRAIN}/paralogs.fa
done

```

```shell script
cd ~/data/anchr/paralogs

egaz template \
    genomes/s288c \
    --self -o yeast/ \
    --length 1000 --parallel 4 -v

bash yeast/1_self.sh
bash yeast/3_proc.sh

for STRAIN in \
    s288c \
    ; do
    cat yeast/Results/${STRAIN}/${STRAIN}.multi.fas |
        faops filter -N -d stdin stdout \
        > ../ref/${STRAIN}/paralogs.fa
done

```

## repetitive

```shell script
mkdir -p ~/data/anchr/repetitive
cd ~/data/anchr/repetitive

for STRAIN in \
    g37 mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    s288c \
    ; do
    echo >&2 ${STRAIN};

    kat sect -t 4 -o ${STRAIN} -F ../ref/${STRAIN}/genome.fa ../ref/${STRAIN}/genome.fa

    cat ${STRAIN}-repetitive.fa |
        faops filter -N -d -a 100 stdin stdout \
        > ../ref/${STRAIN}/repetitive.fa

done

```
