# QA resources: genome, annotations, paralogs, and repeats

End users of [Anchr](https://github.com/wang-q/anchr) don't need to run the following codes. We use
these data just for quality assessments.

These steps require another project: [App::Egaz](https://github.com/wang-q/App-Egaz).

Paralogs detected here **may** overlap with transposons/retrotransposons.

## Download

### E. coli

* *Escherichia* virus Lambda

  * Taxonomy ID: [10710](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10710)
  * Ref. Assembly: [GCF_000840245.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000840245.1/)

* *Escherichia coli* str. K-12 substr. MG1655

  * Taxonomy ID: [511145](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145)
  * Assembly: [GCF_000005845.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)
  * Proportion of paralogs (> 1000 bp): 0.0325

* *Escherichia coli* str. K-12 substr. DH5alpha

  * Taxonomy ID: [83333](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=83333)
  * Assembly: [GCF_001723515.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001723515.1)
  * Proportion of paralogs (> 1000 bp): 0.0342

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

### FDA-ARGOS

* *Francisella tularensis* subsp. tularensis SCHU S4; 土拉热弗朗西斯氏菌
    * Taxonomy ID: [177416](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=177416)
    * Assembly: [GCF_000008985.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000008985.1)
    * Proportion of paralogs (> 1000 bp): 0.0438

* *Haemophilus influenzae* Rd KW20; 流感嗜血杆菌
    * Taxonomy ID: [71421](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=71421)
    * Assembly: [GCF_000027305.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000027305.1)
    * Proportion of paralogs (> 1000 bp): 0.0324

* *Campylobacter jejuni* subsp. jejuni NCTC 11168; ATCC 700819, 空肠弯曲杆菌
    * Taxonomy ID: [192222](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=192222)
    * Assembly: [GCF_000009085.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000009085.1)
    * Proportion of paralogs (> 1000 bp): 0.0196

* *Legionella pneumophila* subsp. pneumophila str. Philadelphia 1; 嗜肺军团菌
    * Taxonomy ID: [272624](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=272624)
    * Assembly: [GCF_000008485.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000008485.1)
    * Proportion of paralogs (> 1000 bp): 0.0264

* *Corynebacterium diphtheriae* NCTC 13129; 白喉杆菌
    * Taxonomy ID: [257309](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=257309)
    * Assembly: [GCF_000195815.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000195815.1)
    * Proportion of paralogs (> 1000 bp): 0.0180

* *Clostridioides difficile* 630; 艰难梭菌
    * Taxonomy ID: [272563](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=272563)
    * Assembly: [GCF_000009205.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000009205.1)
    * Proportion of paralogs (> 1000 bp): 0.0661

```shell script
cd ~/data/anchr/assembly

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/008/985/GCF_000008985.1_ASM898v1 \
    Ftul/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/027/305/GCF_000027305.1_ASM2730v1 \
    Hinf/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1 \
    Cjej/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/008/485/GCF_000008485.1_ASM848v1 \
    Lpne/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/195/815/GCF_000195815.1_ASM19581v1 \
    Cdip/

rsync -avP \
    ftp.ncbi.nlm.nih.gov::genomes/all/GCF/000/009/205/GCF_000009205.1_ASM920v1 \
    Cdif/

```

```shell script
mkdir -p ~/data/anchr/ref
cd ~/data/anchr/ref

for STRAIN in \
    lambda mg1655 dh5alpha \
    Bcer Mabs Rsph Vcho \
    Ftul Hinf Cjej Lpne Cdip Cdif \
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
    Ftul Hinf Cjej Lpne Cdip Cdif \
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
        --repeatmasker '--parallel 4' -v
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
    --length 1000 --parallel 4 -v

bash e_coli/1_self.sh
bash e_coli/3_proc.sh
bash e_coli/4_circos.sh

for STRAIN in \
    lambda mg1655 dh5alpha \
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
    --circos \
    --length 1000 --parallel 4 -v

bash gage_b/1_self.sh
bash gage_b/3_proc.sh
bash gage_b/4_circos.sh

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
    genomes/Ftul genomes/Hinf genomes/Cjej genomes/Lpne genomes/Cdip genomes/Cdif \
    --self -o fda_argos/ \
    --circos \
    --length 1000 --parallel 4 -v

bash fda_argos/1_self.sh
bash fda_argos/3_proc.sh
bash fda_argos/4_circos.sh

for STRAIN in \
    Ftul Hinf Cjej Lpne Cdip Cdif \
    ; do
    cat fda_argos/Results/${STRAIN}/${STRAIN}.multi.fas |
        faops filter -N -d stdin stdout \
        > ../ref/${STRAIN}/paralogs.fa
done

```
