# Assemble genomes from FDA-ARGOS data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from FDA-ARGOS data sets](#assemble-genomes-from-fda-argos-data-sets)
- [Download](#download)
  - [Good assemblies before 2010](#good-assemblies-before-2010)
  - [Reference Assemblies](#reference-assemblies)
  - [Illumina Reads](#illumina-reads)
  - [Reference Genomes](#reference-genomes)
  - [Paralogs](#paralogs)
  - [Repetitives](#repetitives)
  - [Copy/link files](#copylink-files)
- [Ca_jej_jejuni_NCTC_11168_ATCC_700819](#ca_jej_jejuni_nctc_11168_atcc_700819)
- [Clostridio_dif_630](#clostridio_dif_630)
- [Co_dip_NCTC_13129](#co_dip_nctc_13129)
- [Fr_tul_tularensis_SCHU_S4](#fr_tul_tularensis_schu_s4)
- [All strains](#all-strains)


* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/anchr/fda_argos/ \
    wangq@202.119.37.251:data/anchr/fda_argos

# rsync -avP wangq@202.119.37.251:data/anchr/fda_argos/ ~/data/anchr/fda_argos

```

# Download

## Good assemblies before 2010

* Get BioSample from [PRJNA231221](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA231221)
* Save `Full (text)` to `fda_argos.samples.txt` with a web browser

```shell script
mkdir -p ~/data/anchr/fda_argos/ar
cd ~/data/anchr/fda_argos/ar

cat fda_argos.samples.txt |
    grep "^Organism: " |
    sed 's/Organism: //' |
    grep -v "\[" |
    grep -v "virus" |
    sort -u \
    > organism.lst

mysql -uroot ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level, released_date
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
    " |
    perl -nla -F"\t" -e '
        $. == 1 and print and next;
        ($year) = split q(-), $F[5];
        print if $year <= 2010;
    ' |
    tsv-filter -H --istr-in-fld assembly_level:complete --or --istr-in-fld assembly_level:genome \
    > good_assembly.tsv

mysql -uroot ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level, released_date
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
    " |
    perl -nla -F"\t" -e '
        $. == 1 and print and next;
        ($year) = split q(-), $F[5];
        print if $year <= 2010;
    ' |
    tsv-filter -H --istr-in-fld assembly_level:complete --or --istr-in-fld assembly_level:genome \
    >> good_assembly.tsv

mysql -uroot ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level, released_date
    FROM ar
    WHERE 1=1
        AND taxonomy_id in (
            262316, # Mycobacterium avium subsp. paratuberculosis K-10
            83332,  # Mycobacterium tuberculosis H37Rv
            224326, # Borreliella burgdorferi B31
            272563, # Clostridioides difficile 630
            267671  # Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130
        )
    " \
    >> good_assembly.tsv

cat good_assembly.tsv |
    keep-header -- grep -Fw -f organism.lst |
    keep-header -- tsv-sort -k1,1 \
    > cross.tsv

cat cross.tsv |
    grep -v '^#' |
    tsv-filter --not-regex "1:\[.+\]" |
    tsv-filter --not-regex "1:," |
    tsv-select -f 1-5 |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[5]}++;
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > cross.assembly.tsv

cat fda_argos.samples.txt |
    sed '/identification method/d' |
    sed '/latitude and longitude/d' |
    sed '/\/host/d' |
    sed '/\/geographic/d' |
    sed '/\/collected/d' |
    sed '/\/collection/d' |
    sed '/\/culture/d' |
    sed '/\/isolation/d' |
    sed '/\/isolate/d' |
    sed '/\/Laboratory/d' |
    sed '/\/reference/d' |
    sed '/\/passage/d' |
    sed '/\/Passage/d' |
    sed '/\/type-material/d' |
    sed '/\/Extraction/d' |
    sed '/^Accession/d' \
    > samples.txt

```

## Reference Assemblies

```shell script
cd ~/data/anchr/fda_argos/

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ar/cross.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/cross.assembly.rsync.sh

bash ASSEMBLY/cross.assembly.collect.sh

```

## Illumina Reads

* Cross match strain names between `ar/samples.txt` and `ar/cross.assembly.tsv` manually

```shell script
mkdir -p ~/data/anchr/fda_argos/ena
cd ~/data/anchr/fda_argos/ena

cat << EOF > source.csv
SAMN16357368,Bact_the_VPI_5482,Bacteroides thetaiotaomicron VPI-5482
SAMN03996253,Bar_bac_KC583,Bartonella bacilliformis KC583
SAMN03996256,Bar_hen_Houston_1,Bartonella henselae str. Houston-1
SAMN03996258,Bord_bronchis_RB50,Bordetella bronchiseptica RB50
SAMN04875532,Bord_pert_Tohama_I,Bordetella pertussis Tohama I
SAMN04875533,Borr_bur_B31,Borreliella burgdorferi B31
SAMN10228570,Bu_mall_ATCC_23344,Burkholderia mallei ATCC 23344
SAMN05004753,Bu_tha_E264,Burkholderia thailandensis E264
SAMN04875589,Ca_jej_jejuni_NCTC_11168_ATCC_700819,Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819
SAMN16357198,Ci_kos_ATCC_BAA_895,Citrobacter koseri ATCC BAA-895
SAMN04875594,Clostridio_dif_630,Clostridioides difficile 630
SAMN04875534,Co_dip_NCTC_13129,Corynebacterium diphtheriae NCTC 13129
SAMN16357334,Co_kro_DSM_44385,Corynebacterium kroppenstedtii DSM 44385
SAMN16357163,Co_ure_DSM_7109,Corynebacterium urealyticum DSM 7109
SAMN11056390,Cup_met_CH34,Cupriavidus metallidurans CH34
SAMN10228557,Cut_acn_SK137,Cutibacterium acnes SK137
SAMN16357201,Es_fer_ATCC_35469,Escherichia fergusonii ATCC 35469
SAMN04875573,Fr_tul_tularensis_SCHU_S4,Francisella tularensis subsp. tularensis SCHU S4
SAMN04875536,Ha_inf_Rd_KW20,Haemophilus influenzae Rd KW20
SAMN06173312,He_pyl_26695,Helicobacter pylori 26695
SAMN06173313,He_pyl_J99,Helicobacter pylori J99
SAMN16357582,J_deni_DSM_20603,Jonesia denitrificans DSM 20603
SAMN06173315,Ko_rhi_DC2201,Kocuria rhizophila DC2201
SAMN16357230,Ky_sed_DSM_20547,Kytococcus sedentarius DSM 20547
SAMN04875539,Leg_pneumop_pneumophila_Philadelphia_1,Legionella pneumophila subsp. pneumophila str. Philadelphia 1
SAMN04875540,Lep_int_Copenhageni_Fiocruz_L1_130,Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130
SAMN16357202,Leu_mes_mesenteroides_ATCC_8293,Leuconostoc mesenteroides subsp. mesenteroides ATCC 8293
SAMN06173318,Mycobacteri_avi_paratuberculosis_K_10,Mycobacterium avium subsp. paratuberculosis K-10
SAMN07312468,Mycobacteri_tub_H37Ra,Mycobacterium tuberculosis H37Ra
SAMN11056472,Mycobacteri_tub_H37Rv,Mycobacterium tuberculosis H37Rv
SAMN11056394,Mycol_sme_MC2_155,Mycolicibacterium smegmatis MC2 155
SAMN04875544,N_gon_FA_1090,Neisseria gonorrhoeae FA 1090
SAMN04875545,N_men_FAM18,Neisseria meningitidis FAM18
SAMN04875546,N_men_MC58,Neisseria meningitidis MC58
SAMN16357208,O_anthro_ATCC_49188,Ochrobactrum anthropi ATCC 49188
#SAMN16357376,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
#SAMN16357375,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
SAMN06173357,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
SAMN06173319,Pre_mel_ATCC_25845,Prevotella melaninogenica ATCC 25845
SAMN06173320,Pse_pro_Pf_5,Pseudomonas protegens Pf-5
SAMN06173321,Psy_cry_K5,Psychrobacter cryohalolentis K5
SAMN06173322,Ros_deni_OCh_114,Roseobacter denitrificans OCh 114
SAMN11056396,She_putr_CN_32,Shewanella putrefaciens CN-32
SAMN03255464,Sta_aure_aureus_Mu50,Staphylococcus aureus subsp. aureus Mu50
SAMN03255448,Sta_aure_aureus_N315,Staphylococcus aureus subsp. aureus N315
SAMN03255470,Sta_aure_aureus_NCTC_8325,Staphylococcus aureus subsp. aureus NCTC 8325
#SAMN03255450,Sta_aure_aureus_NCTC_8325,Staphylococcus aureus subsp. aureus NCTC 8325
SAMN13450443,Sta_epi_ATCC_12228,Staphylococcus epidermidis ATCC 12228
SAMN06173368,Sta_sap_saprophyticus_ATCC_15305_NCTC_7292,Staphylococcus saprophyticus subsp. saprophyticus ATCC 15305 = NCTC 7292
SAMN06173323,Streptob_moni_DSM_12112,Streptobacillus moniliformis DSM 12112
SAMN10228564,Y_pseudot_YPIII,Yersinia pseudotuberculosis YPIII
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - -p illumina ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name                                       | srx        | platform | layout | ilength | srr         | spot     | base      |
|:-------------------------------------------|:-----------|:---------|:-------|:--------|:------------|:---------|:----------|
| Bar_bac_KC583                              | SRX1385802 | ILLUMINA | PAIRED | 443     | SRR2823707  | 7131080  | 1.34G     |
| Bar_hen_Houston_1                          | SRX1385804 | ILLUMINA | PAIRED | 391     | SRR2823712  | 6669509  | 1.25G     |
| Bord_bronchis_RB50                         | SRX1385806 | ILLUMINA | PAIRED | 397     | SRR2823715  | 2960592  | 570.33M   |
| Bord_bronchis_RB50                         | SRX1385807 | ILLUMINA | PAIRED | 397     | SRR2823716  | 3059385  | 589.37M   |
| Bord_pert_Tohama_I                         | SRX2179101 | ILLUMINA | PAIRED | 528     | SRR4271511  | 4333095  | 834.74M   |
| Bord_pert_Tohama_I                         | SRX2179104 | ILLUMINA | PAIRED | 528     | SRR4271510  | 3949224  | 760.79M   |
| Borr_bur_B31                               | SRX2179106 | ILLUMINA | PAIRED | 577     | SRR4271513  | 7284107  | 1.37G     |
| Bu_mall_ATCC_23344                         | SRX4900909 | ILLUMINA | PAIRED | 529     | SRR8072938  | 4949942  | 1.39G     |
| Bu_tha_E264                                | SRX2110113 | ILLUMINA | PAIRED | 501     | SRR4125425  | 7184151  | 1.35G     |
| Ca_jej_jejuni_NCTC_11168_ATCC_700819       | SRX2107012 | ILLUMINA | PAIRED | 562     | SRR4125016  | 7696800  | 1.45G     |
| Ci_kos_ATCC_BAA_895                        | SRX9293511 | ILLUMINA | PAIRED | 534     | SRR12825872 | 15181966 | 4.27G     |
| Clostridio_dif_630                         | SRX2107163 | ILLUMINA | PAIRED | 523     | SRR4125185  | 6595393  | 1.24G     |
| Co_dip_NCTC_13129                          | SRX2179108 | ILLUMINA | PAIRED | 593     | SRR4271515  | 5564406  | 1.05G     |
| Co_ure_DSM_7109                            | SRX9295022 | ILLUMINA | PAIRED | 592     | SRR12827384 | 7515280  | 2.11G     |
| Cup_met_CH34                               | SRX5934689 | ILLUMINA | PAIRED | 449     | SRR9161655  | 8291696  | 2.33G     |
| Cut_acn_SK137                              | SRX4900879 | ILLUMINA | PAIRED | 622     | SRR8072906  | 6793960  | 1.91G     |
| Es_fer_ATCC_35469                          | SRX9293517 | ILLUMINA | PAIRED | 537     | SRR12825878 | 7254756  | 2.04G     |
| Fr_tul_tularensis_SCHU_S4                  | SRX2105481 | ILLUMINA | PAIRED | 560     | SRR4124773  | 10615135 | 2G        |
| Ha_inf_Rd_KW20                             | SRX2104758 | ILLUMINA | PAIRED | 516     | SRR4123928  | 6115624  | 1.15G     |
| He_pyl_26695                               | SRX2742168 | ILLUMINA | PAIRED | 511     | SRR5454002  | 1440239  | 277.45M   |
| He_pyl_26695                               | SRX2742169 | ILLUMINA | PAIRED | 511     | SRR5454001  | 4696741  | 904.79M   |
| He_pyl_J99                                 | SRX2705206 | ILLUMINA | PAIRED | 530     | SRR5413257  | 5637532  | 1.06G     |
| Ko_rhi_DC2201                              | SRX2737867 | ILLUMINA | PAIRED | 524     | SRR5449077  | 1313108  | 252.96M   |
| Ko_rhi_DC2201                              | SRX2737868 | ILLUMINA | PAIRED | 524     | SRR5449076  | 4139863  | 797.51M   |
| Ky_sed_DSM_20547                           | SRX9293604 | ILLUMINA | PAIRED | 414     | SRR12825966 | 3483084  | 1,003.16M |
| Leg_pneumop_pneumophila_Philadelphia_1     | SRX2179279 | ILLUMINA | PAIRED | 570     | SRR4272054  | 5249241  | 1,011.23M |
| Lep_int_Copenhageni_Fiocruz_L1_130         | SRX2179272 | ILLUMINA | PAIRED | 466     | SRR4272049  | 6556612  | 1.23G     |
| Leu_mes_mesenteroides_ATCC_8293            | SRX9293521 | ILLUMINA | PAIRED | 531     | SRR12825883 | 13118133 | 3.69G     |
| Mycobacteri_avi_paratuberculosis_K_10      | SRX2705222 | ILLUMINA | PAIRED | 505     | SRR5413272  | 5359484  | 1.01G     |
| Mycobacteri_tub_H37Ra                      | SRX3046555 | ILLUMINA | PAIRED | 522     | SRR5879398  | 5473713  | 1.03G     |
| Mycobacteri_tub_H37Rv                      | SRX6389030 | ILLUMINA | PAIRED | 604     | SRR9626912  | 7874573  | 2.21G     |
| Mycol_sme_MC2_155                          | SRX5935584 | ILLUMINA | PAIRED | 522     | SRR9162626  | 5316609  | 1.5G      |
| N_gon_FA_1090                              | SRX2179294 | ILLUMINA | PAIRED | 533     | SRR4272072  | 7384079  | 1.39G     |
| N_men_FAM18                                | SRX2179296 | ILLUMINA | PAIRED | 519     | SRR4272074  | 12948421 | 2.44G     |
| N_men_MC58                                 | SRX2179304 | ILLUMINA | PAIRED | 536     | SRR4272082  | 6907195  | 1.3G      |
| O_anthro_ATCC_49188                        | SRX9292947 | ILLUMINA | PAIRED | 468     | SRR12825308 | 8125004  | 2.29G     |
| Par_dis_ATCC_8503                          | SRX2701505 | ILLUMINA | PAIRED | 429     | SRR5409214  | 12052578 | 3.39G     |
| Pre_mel_ATCC_25845                         | SRX2705223 | ILLUMINA | PAIRED | 563     | SRR5413273  | 5958667  | 1.12G     |
| Pse_pro_Pf_5                               | SRX2705227 | ILLUMINA | PAIRED | 527     | SRR5413277  | 4966775  | 956.81M   |
| Pse_pro_Pf_5                               | SRX2705228 | ILLUMINA | PAIRED | 527     | SRR5413278  | 1288497  | 248.22M   |
| Psy_cry_K5                                 | SRX2705242 | ILLUMINA | PAIRED | 537     | SRR5413293  | 5721091  | 1.08G     |
| Ros_deni_OCh_114                           | SRX2737869 | ILLUMINA | PAIRED | 512     | SRR5449080  | 1739986  | 335.19M   |
| Ros_deni_OCh_114                           | SRX2737871 | ILLUMINA | PAIRED | 512     | SRR5449079  | 4276643  | 823.86M   |
| She_putr_CN_32                             | SRX5936044 | ILLUMINA | PAIRED | 536     | SRR9163108  | 9195551  | 2.59G     |
| Sta_aure_aureus_Mu50                       | SRX981321  | ILLUMINA | PAIRED | 563     | SRR1955819  | 5759114  | 1.08G     |
| Sta_aure_aureus_Mu50                       | SRX981322  | ILLUMINA | PAIRED | 563     | SRR1955820  | 5484223  | 1.03G     |
| Sta_aure_aureus_N315                       | SRX981099  | ILLUMINA | PAIRED | 547     | SRR1955595  | 6533728  | 1.23G     |
| Sta_aure_aureus_NCTC_8325                  | SRX981348  | ILLUMINA | PAIRED | 536     | SRR1955845  | 6994744  | 1.32G     |
| Sta_epi_ATCC_12228                         | SRX8576561 | ILLUMINA | PAIRED | 573     | SRR12047979 | 7608242  | 2.14G     |
| Sta_sap_saprophyticus_ATCC_15305_NCTC_7292 | SRX2730648 | ILLUMINA | PAIRED | 424     | SRR5440743  | 17328248 | 4.87G     |
| Streptob_moni_DSM_12112                    | SRX2705245 | ILLUMINA | PAIRED | 524     | SRR5413295  | 5272442  | 1,015.69M |
| Y_pseudot_YPIII                            | SRX4900894 | ILLUMINA | PAIRED | 594     | SRR8072924  | 4239591  | 1.19G     |
| Y_pseudot_YPIII                            | SRX4900896 | ILLUMINA | PAIRED | 594     | SRR8072925  | 4106622  | 1.16G     |


## Reference Genomes

```shell script
mkdir -p ~/data/anchr/fda_argos/ref
cd ~/data/anchr/fda_argos/ref

for STRAIN in \
    $(
        cat ../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u

    ) \
    ; do
    echo >&2 ${STRAIN};
    mkdir -p ${STRAIN}

    if [ ! -d ../ASSEMBLY/${STRAIN} ]; then
        echo >&2 Skip ${STRAIN};
        continue;
    fi

    find ../ASSEMBLY/${STRAIN}/ -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs gzip -dcf |
        faops filter -N -s stdin ${STRAIN}/genome.fa

done

```

## Paralogs

* RepeatMasker

```shell script
mkdir -p ~/data/anchr/fda_argos/paralogs/genomes
cd ~/data/anchr/fda_argos/paralogs/genomes

for STRAIN in \
    $(
        cat ../../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u

    ) \
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

* Self-alignments

```shell script
cd ~/data/anchr/fda_argos/paralogs

egaz template \
    $(
        cat ../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u |
            xargs -I{} echo "genomes/{}"
    ) \
    --self -o ./ \
    --circos \
    --length 1000 --parallel 4 -v

bash ./1_self.sh
bash ./3_proc.sh
bash ./4_circos.sh

# paralogs.fa
for STRAIN in \
    $(
        cat ../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u

    ) \
    ; do
    cat ./Results/${STRAIN}/${STRAIN}.multi.fas |
        faops filter -N -d stdin stdout \
        > ../ref/${STRAIN}/paralogs.fa
done

# stats
for STRAIN in \
    $(
        cat ../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u

    ) \
    ; do
    cat ./Results/${STRAIN}/${STRAIN}.cover.csv |
        grep "^all" |
        sed "s/^all/${STRAIN}/"
done |
    (echo -e '#name,chrLength,size,coverage' && cat ) \
    > cover.csv

mlr --icsv --omd cat cover.csv

```

| #name                                      | chrLength | size   | coverage |
|:-------------------------------------------|:----------|:-------|:---------|
| Bact_the_VPI_5482                          | 6260361   | 286736 | 0.0458   |
| Bar_bac_KC583                              | 1445021   | 98593  | 0.0682   |
| Bar_hen_Houston_1                          | 1931047   | 146599 | 0.0759   |
| Bord_bronchis_RB50                         | 5339179   | 166533 | 0.0312   |
| Bord_pert_Tohama_I                         | 4086189   | 81815  | 0.0200   |
| Borr_bur_B31                               | 1467551   | 305695 | 0.2083   |
| Bu_mall_ATCC_23344                         | 5835527   | 223898 | 0.0384   |
| Bu_tha_E264                                | 6723972   | 274890 | 0.0409   |
| Ca_jej_jejuni_NCTC_11168_ATCC_700819       | 1641481   | 32217  | 0.0196   |
| Ci_kos_ATCC_BAA_895                        | 4720462   | 147120 | 0.0312   |
| Clostridio_dif_630                         | 4274782   | 278659 | 0.0652   |
| Co_dip_NCTC_13129                          | 2488635   | 44783  | 0.0180   |
| Co_kro_DSM_44385                           | 2446804   | 17352  | 0.0071   |
| Co_ure_DSM_7109                            | 2369219   | 83469  | 0.0352   |
| Cup_met_CH34                               | 6913352   | 362651 | 0.0525   |
| Cut_acn_SK137                              | 2495334   | 27313  | 0.0109   |
| Es_fer_ATCC_35469                          | 4643861   | 129186 | 0.0278   |
| Fr_tul_tularensis_SCHU_S4                  | 1892775   | 82812  | 0.0438   |
| Ha_inf_Rd_KW20                             | 1830138   | 59315  | 0.0324   |
| He_pyl_26695                               | 1667867   | 72577  | 0.0435   |
| He_pyl_J99                                 | 1643831   | 45987  | 0.0280   |
| J_deni_DSM_20603                           | 2749646   | 52925  | 0.0192   |
| Ko_rhi_DC2201                              | 2697540   | 29304  | 0.0109   |
| Ky_sed_DSM_20547                           | 2785024   | 95421  | 0.0343   |
| Leg_pneumop_pneumophila_Philadelphia_1     | 3397754   | 89782  | 0.0264   |
| Lep_int_Copenhageni_Fiocruz_L1_130         | 4627366   | 105815 | 0.0229   |
| Leu_mes_mesenteroides_ATCC_8293            | 2075763   | 56816  | 0.0274   |
| Mycobacteri_avi_paratuberculosis_K_10      | 4829781   | 239309 | 0.0495   |
| Mycobacteri_tub_H37Ra                      | 4419977   | 187555 | 0.0424   |
| Mycobacteri_tub_H37Rv                      | 4411532   | 176345 | 0.0400   |
| Mycol_sme_MC2_155                          | 6988209   | 350227 | 0.0501   |
| N_gon_FA_1090                              | 2153922   | 110972 | 0.0515   |
| N_men_FAM18                                | 2194961   | 123508 | 0.0563   |
| N_men_MC58                                 | 2272360   | 183408 | 0.0807   |
| O_anthro_ATCC_49188                        | 5205777   | 142165 | 0.0273   |
| Par_dis_ATCC_8503                          | 4811379   | 167145 | 0.0347   |
| Pre_mel_ATCC_25845                         | 3168282   | 72306  | 0.0228   |
| Pse_pro_Pf_5                               | 7074893   | 175507 | 0.0248   |
| Psy_cry_K5                                 | 3101097   | 81488  | 0.0263   |
| Ros_deni_OCh_114                           | 4133097   | 39434  | 0.0095   |
| She_putr_CN_32                             | 4659220   | 148175 | 0.0318   |
| Sta_aure_aureus_Mu50                       | 2903636   | 86252  | 0.0297   |
| Sta_aure_aureus_N315                       | 2814816   | 96559  | 0.0343   |
| Sta_aure_aureus_NCTC_8325                  | 2821361   | 57079  | 0.0202   |
| Sta_epi_ATCC_12228                         | 2532762   | 85141  | 0.0336   |
| Sta_sap_saprophyticus_ATCC_15305_NCTC_7292 | 2577899   | 44095  | 0.0171   |
| Streptob_moni_DSM_12112                    | 1673280   | 241650 | 0.1444   |
| Y_pseudot_YPIII                            | 4689441   | 114669 | 0.0245   |


## Repetitives

```shell script
mkdir -p ~/data/anchr/fda_argos/repetitives
cd ~/data/anchr/fda_argos/repetitives

for STRAIN in \
    $(
        cat ../ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u
    ) \
    ; do
    echo >&2 ${STRAIN};

    kat sect -t 4 -o ${STRAIN} -F ../ref/${STRAIN}/genome.fa ../ref/${STRAIN}/genome.fa

    cat ${STRAIN}-repetitive.fa |
        faops filter -N -d -a 100 stdin stdout \
        > ../ref/${STRAIN}/repetitives.fa

done

```

## Copy/link files

```shell script
cd ~/data/anchr/fda_argos

for STRAIN in \
    $(
        cat ena/source.csv |
            grep -v "^#" |
            cut -d, -f 2 |
            sort -u

    ) \
    ; do
    mkdir -p ${STRAIN}/1_genome

    cp ref/${STRAIN}/genome.fa ${STRAIN}/1_genome
    cp ref/${STRAIN}/paralogs.fa ${STRAIN}/1_genome
    cp ref/${STRAIN}/repetitives.fa ${STRAIN}/1_genome

done

# Clean symlinks
find . -name "R1.fq.gz" -or \
       -name "R2.fq.gz" -or \
       -name "S1.fq.gz" -or \
       -name "S2.fq.gz" -or \
       -name "T1.fq.gz" -or \
       -name "T2.fq.gz" |
    xargs rm
cat ena/ena_info.csv |
    grep -v "^#" |
    grep -v "^name" |
    parallel -j 1 -k --colsep "," '
        mkdir -p {1}/2_illumina
        cd {1}/2_illumina

        if [[ -L R1.fq.gz ]]; then
            if [[ -L S1.fq.gz ]]; then
                if [[ -L T1.fq.gz ]]; then
                    echo "Symlinks exist"
                else
                    ln -s ../../ena/{6}_1.fastq.gz T1.fq.gz
                    ln -s ../../ena/{6}_2.fastq.gz T2.fq.gz
                fi
            else
                ln -s ../../ena/{6}_1.fastq.gz S1.fq.gz
                ln -s ../../ena/{6}_2.fastq.gz S2.fq.gz
            fi
        else
            ln -s ../../ena/{6}_1.fastq.gz R1.fq.gz
            ln -s ../../ena/{6}_2.fastq.gz R2.fq.gz
        fi
    '

```

# Ca_jej_jejuni_NCTC_11168_ATCC_700819

```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Ca_jej_jejuni_NCTC_11168_ATCC_700819

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1641481 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --sample 300 \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 101 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md

# rm -fr 4_down_sampling 6_down_sampling

bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"

# bash 0_master.sh
# bash 0_cleanup.sh

```


Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 434.6 |    300 | 2204.8 |                         49.43% |
| R.tadpole.bbtools | 310.3 |    300 |  104.6 |                         48.23% |
| R.genome.picard   | 310.8 |    300 |   86.6 |                             FR |
| R.tadpole.picard  | 309.9 |    299 |   86.5 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 0         | 0               | 0.0000       | 33.93   |
| R.31 | 0         | 0               | 0.0000       | 33.22   |
| R.41 | 0         | 0               | 0.0000       | 32.85   |
| R.51 | 451       | 1268552         | 0.0000       | 32.60   |
| R.61 | 360       | 1589497         | 0.0000       | 32.41   |
| R.71 | 268       | 1589092         | 0.0000       | 32.24   |
| R.81 | 183       | 1622277         | 0.0000       | 32.07   |


Table: statReads

| Name        |     N50 |     Sum |        # |
|:------------|--------:|--------:|---------:|
| Genome      | 1641481 | 1641481 |        1 |
| Paralogs    |    6079 |   33258 |       13 |
| Repetitives |    2339 |   41610 |       55 |
| Illumina.R  |     101 |   1.55G | 15393600 |
| trim.R      |     100 | 475.32M |  4775662 |
| Q25L60      |     100 | 462.82M |  4660316 |
| Q30L60      |     100 | 439.09M |  4464993 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.54G | 15283902 |
| highpass | 101 |   1.54G | 15227304 |
| sample   | 101 | 492.44M |  4875688 |
| trim     | 100 | 475.42M |  4776718 |
| filter   | 100 | 475.32M |  4775662 |
| R1       | 100 | 237.93M |  2387831 |
| R2       | 100 | 237.38M |  2387831 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3726	0.07642%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	1056	0.02211%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	1056	0.02211%
```

```text
#R.peaks
#k	31
#unique_kmers	7511056
#error_kmers	5895969
#genomic_kmers	1615087
#main_peak	195
#genome_size_in_peaks	1639571
#genome_size	1640629
#haploid_genome_size	1640629
#fold_coverage	195
#haploid_fold_coverage	195
#ploidy	1
#percent_repeat_in_peaks	1.497
#percent_repeat	1.530
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 475.31M | 4775596 |
| ecco          | 100 | 475.31M | 4775596 |
| eccc          | 100 | 475.31M | 4775596 |
| ecct          | 100 | 472.64M | 4748436 |
| extended      | 140 | 662.19M | 4748436 |
| merged.raw    | 333 | 617.95M | 1930935 |
| unmerged.raw  | 140 | 123.01M |  886566 |
| unmerged.trim | 140 | 123.01M |  886556 |
| M1            | 333 | 617.43M | 1929331 |
| U1            | 140 |  61.72M |  443278 |
| U2            | 140 |  61.29M |  443278 |
| Us            |   0 |       0 |       0 |
| M.cor         | 317 | 742.37M | 4745218 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 159.5 |    164 |  20.8 |          4.68% |
| M.ihist.merge.txt  | 320.0 |    321 |  58.6 |         81.33% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 289.6 |  272.9 |    5.74% | "71" | 1.64M | 1.62M |     0.99 | 0:01'02'' |
| Q25L60.R | 282.0 |  268.8 |    4.67% | "71" | 1.64M | 1.62M |     0.99 | 0:01'02'' |
| Q30L60.R | 267.6 |  257.4 |    3.80% | "71" | 1.64M | 1.62M |     0.99 | 0:00'59'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.87% |     60522 |  1.6M | 54 |      6071 | 11.44K | 241 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| Q0L0X40P001   |   40.0 |  98.78% |     75189 |  1.6M | 50 |      6071 | 11.86K | 260 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:22 |
| Q0L0X40P002   |   40.0 |  98.81% |     55716 |  1.6M | 56 |      1019 | 13.23K | 256 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:22 |
| Q0L0X80P000   |   80.0 |  98.36% |     46692 |  1.6M | 63 |      6071 |  9.53K | 147 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q0L0X80P001   |   80.0 |  98.40% |     47742 | 1.59M | 57 |      6071 | 10.56K | 137 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q0L0X80P002   |   80.0 |  98.31% |     47680 | 1.59M | 69 |      6071 | 11.01K | 163 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q25L60X40P000 |   40.0 |  98.83% |     79980 |  1.6M | 48 |      1034 | 12.36K | 235 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:21 |
| Q25L60X40P001 |   40.0 |  98.83% |     67014 |  1.6M | 50 |      6071 | 11.67K | 234 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:21 |
| Q25L60X40P002 |   40.0 |  98.81% |     70720 |  1.6M | 52 |      6071 | 11.45K | 239 |   39.0 | 5.0 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:20 |
| Q25L60X80P000 |   80.0 |  98.50% |     56890 |  1.6M | 55 |      6071 |  9.37K | 135 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q25L60X80P001 |   80.0 |  98.23% |     57904 | 1.59M | 56 |      5734 | 11.05K | 135 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q25L60X80P002 |   80.0 |  98.47% |     52044 | 1.59M | 61 |      6071 | 11.87K | 147 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q30L60X40P000 |   40.0 |  98.86% |     70653 |  1.6M | 46 |      6071 | 11.75K | 223 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:23 |
| Q30L60X40P001 |   40.0 |  98.86% |     70723 | 1.59M | 49 |      1041 | 12.81K | 220 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:22 |
| Q30L60X40P002 |   40.0 |  98.80% |     71584 |  1.6M | 47 |      6071 | 10.67K | 215 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:23 |
| Q30L60X80P000 |   80.0 |  98.45% |     56908 | 1.59M | 53 |      6071 | 11.56K | 132 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q30L60X80P001 |   80.0 |  98.48% |     67712 | 1.59M | 52 |      6071 | 10.55K | 145 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q30L60X80P002 |   80.0 |  98.64% |     70558 | 1.59M | 50 |      6071 | 11.77K | 145 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.06% |     75105 |  1.6M | 48 |       985 |  9.63K | 114 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:21 |
| MRX40P001 |   40.0 |  97.92% |     76549 |  1.6M | 46 |      6071 | 10.29K | 105 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:21 |   0:00:20 |
| MRX40P002 |   40.0 |  98.00% |     79947 | 1.59M | 45 |      6071 | 11.43K | 109 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:21 |   0:00:24 |
| MRX80P000 |   80.0 |  97.78% |     48762 | 1.59M | 58 |      6071 | 11.04K | 129 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:22 |
| MRX80P001 |   80.0 |  97.79% |     66433 | 1.59M | 55 |      6071 |  11.9K | 123 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:21 |
| MRX80P002 |   80.0 |  97.64% |     69747 | 1.59M | 59 |      6071 | 11.29K | 132 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:21 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.84% |     71583 | 1.59M | 47 |      6071 | 10.35K | 217 |   39.0 | 5.0 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:23 |
| Q0L0X40P001   |   40.0 |  98.83% |     70844 |  1.6M | 51 |      1063 | 14.58K | 238 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:23 |
| Q0L0X40P002   |   40.0 |  98.87% |     56896 |  1.6M | 51 |      1008 | 12.85K | 244 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:22 |
| Q0L0X80P000   |   80.0 |  98.77% |     79964 |  1.6M | 49 |      6071 | 10.58K | 181 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  98.61% |     56889 |  1.6M | 51 |      6071 | 11.47K | 171 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:22 |
| Q0L0X80P002   |   80.0 |  98.77% |     70561 | 1.59M | 50 |      2340 | 12.61K | 170 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.77% |     71568 |  1.6M | 45 |      1036 | 13.47K | 229 |   39.0 | 5.0 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:21 |
| Q25L60X40P001 |   40.0 |  98.78% |     56914 |  1.6M | 49 |      6071 | 11.55K | 249 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:23 |
| Q25L60X40P002 |   40.0 |  98.74% |     60507 |  1.6M | 52 |      1056 | 12.94K | 235 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:22 |
| Q25L60X80P000 |   80.0 |  98.65% |     76595 |  1.6M | 48 |      6071 | 10.14K | 171 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:22 |
| Q25L60X80P001 |   80.0 |  98.62% |     70632 | 1.59M | 48 |      6071 | 12.01K | 158 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:21 |
| Q25L60X80P002 |   80.0 |  98.78% |     76543 |  1.6M | 49 |      2340 |  12.7K | 182 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:22 |
| Q30L60X40P000 |   40.0 |  98.90% |     71593 |  1.6M | 44 |      1087 |  13.7K | 270 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:22 |
| Q30L60X40P001 |   40.0 |  98.66% |     70696 | 1.59M | 50 |      6071 | 12.08K | 233 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:21 |
| Q30L60X40P002 |   40.0 |  98.88% |     71574 |  1.6M | 52 |      1048 | 13.36K | 251 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:22 |
| Q30L60X80P000 |   80.0 |  98.80% |     70837 |  1.6M | 47 |      6071 | 11.47K | 193 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  98.75% |     78387 | 1.59M | 46 |       775 | 12.62K | 197 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:23 |
| Q30L60X80P002 |   80.0 |  98.89% |     71569 |  1.6M | 44 |      6071 |  11.9K | 185 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.05% |     81869 |  1.6M | 39 |       985 |  8.17K |  84 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:23 |
| MRX40P001 |   40.0 |  98.04% |     79893 |  1.6M | 40 |      6071 |  9.55K |  91 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:23 |
| MRX40P002 |   40.0 |  98.00% |     81878 | 1.59M | 39 |      6071 | 11.08K |  82 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:21 |
| MRX80P000 |   80.0 |  97.76% |     75131 | 1.59M | 48 |      6071 |   9.6K |  97 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:21 |
| MRX80P001 |   80.0 |  97.70% |     79919 | 1.59M | 50 |      6071 | 10.77K | 101 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:21 |
| MRX80P002 |   80.0 |  97.72% |     79901 | 1.59M | 48 |      6071 |  9.69K |  96 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:21 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.83% |     71589 |  1.6M | 46 |      6069 | 10.22K | 177 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q0L0X40P001   |   40.0 |  98.70% |     75176 |  1.6M | 46 |      6069 | 10.67K | 197 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q0L0X40P002   |   40.0 |  98.81% |     60499 |  1.6M | 49 |      1018 | 15.28K | 217 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| Q0L0X80P000   |   80.0 |  98.58% |     79992 |  1.6M | 44 |      6069 |     9K | 115 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| Q0L0X80P001   |   80.0 |  98.53% |     79952 | 1.59M | 43 |      6069 | 10.35K | 117 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:21 |
| Q0L0X80P002   |   80.0 |  98.75% |     79985 | 1.59M | 44 |      6069 |  11.8K | 132 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q25L60X40P000 |   40.0 |  98.85% |     79894 |  1.6M | 45 |      6069 | 11.04K | 205 |   39.0 | 5.0 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q25L60X40P001 |   40.0 |  98.85% |     79967 |  1.6M | 44 |      1027 | 13.73K | 207 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q25L60X40P002 |   40.0 |  98.88% |     71572 |  1.6M | 46 |      1029 |  13.6K | 217 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q25L60X80P000 |   80.0 |  98.63% |     80753 |  1.6M | 41 |      6069 |  9.29K | 129 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q25L60X80P001 |   80.0 |  98.66% |     70814 | 1.59M | 47 |      6069 | 11.66K | 132 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| Q25L60X80P002 |   80.0 |  98.81% |     79930 | 1.59M | 43 |      6069 | 12.08K | 148 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.91% |     71593 |  1.6M | 44 |      1086 | 12.65K | 220 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q30L60X40P001 |   40.0 |  98.73% |     70694 |  1.6M | 48 |      6069 | 11.69K | 205 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q30L60X40P002 |   40.0 |  98.83% |     71572 |  1.6M | 49 |      6069 |  11.6K | 210 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q30L60X80P000 |   80.0 |  98.88% |     75187 |  1.6M | 44 |      6069 | 11.79K | 159 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| Q30L60X80P001 |   80.0 |  98.75% |     79937 | 1.59M | 43 |      6069 | 11.47K | 149 |   80.0 | 8.0 |  18.7 | 208.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q30L60X80P002 |   80.0 |  98.78% |     71572 |  1.6M | 44 |      6069 | 11.12K | 143 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.07% |     83912 |  1.6M | 38 |      1029 |  9.21K | 83 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| MRX40P001 |   40.0 |  98.12% |     80671 |  1.6M | 38 |      6069 |  9.72K | 89 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:21 |
| MRX40P002 |   40.0 |  97.88% |     81915 | 1.59M | 38 |      6069 | 12.04K | 79 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:20 |
| MRX80P000 |   80.0 |  97.83% |     75137 |  1.6M | 45 |      6069 |  9.45K | 91 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:22 |
| MRX80P001 |   80.0 |  97.83% |     80667 | 1.59M | 43 |      1101 | 12.44K | 89 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:21 |
| MRX80P002 |   80.0 |  97.81% |     80622 |  1.6M | 44 |      6069 |  9.52K | 92 |   79.0 | 8.0 |  18.3 | 206.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |  Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|-----:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  98.62% |    104866 | 1.6M | 28 |      2433 | 48.77K | 27 |  274.0 | 20.0 |  71.3 | 668.0 |   0:00:35 |
| 7_merge_mr_unitigs_bcalm      |  99.07% |    104192 | 1.6M | 29 |     41864 | 51.19K |  5 |  272.0 | 21.0 |  69.7 | 670.0 |   0:00:44 |
| 7_merge_mr_unitigs_superreads |  99.16% |     98667 | 1.6M | 32 |      6071 |  9.44K |  3 |  273.0 | 20.0 |  71.0 | 666.0 |   0:00:47 |
| 7_merge_mr_unitigs_tadpole    |  99.08% |     98667 | 1.6M | 30 |      6069 |  9.34K |  4 |  271.0 | 21.0 |  69.3 | 668.0 |   0:00:45 |
| 7_merge_unitigs_bcalm         |  99.18% |    104187 | 1.6M | 30 |      2414 |    25K | 15 |  271.0 | 20.0 |  70.3 | 662.0 |   0:00:50 |
| 7_merge_unitigs_superreads    |  99.09% |     90266 | 1.6M | 30 |      2340 | 14.33K |  8 |  272.0 | 19.0 |  71.7 | 658.0 |   0:00:46 |
| 7_merge_unitigs_tadpole       |  99.19% |    104870 | 1.6M | 31 |      1039 | 31.85K | 20 |  269.0 | 21.0 |  68.7 | 664.0 |   0:00:50 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  99.54% |    153913 | 1.61M | 15 |      6071 | 11.41K | 32 |  272.0 | 20.0 |  70.7 |  664.0 |   0:00:33 |
| 8_mr_spades  |  99.51% |    189337 | 1.61M | 12 |      6071 | 10.58K | 27 |  449.0 | 30.0 | 119.7 | 1078.0 |   0:00:30 |
| 8_megahit    |  98.70% |    115645 |  1.6M | 22 |      5957 | 10.52K | 44 |  272.0 | 20.0 |  70.7 |  664.0 |   0:00:26 |
| 8_mr_megahit |  99.31% |    154046 | 1.61M | 24 |      5850 |  9.82K | 46 |  449.0 | 30.0 | 119.7 | 1078.0 |   0:00:28 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1641481 | 1641481 |   1 |
| Paralogs                 |    6079 |   33258 |  13 |
| Repetitives              |    2339 |   41610 |  55 |
| 7_merge_anchors.anchors  |  104866 | 1602430 |  28 |
| 7_merge_anchors.others   |    2433 |   48765 |  27 |
| glue_anchors             |  112521 | 1601571 |  25 |
| fill_anchors             |  115643 | 1601689 |  23 |
| spades.contig            |  153957 | 1623445 |  28 |
| spades.scaffold          |  189387 | 1623583 |  26 |
| spades.non-contained     |  153957 | 1620303 |  17 |
| mr_spades.contig         |  189483 | 1624653 |  20 |
| mr_spades.scaffold       |  189483 | 1624653 |  20 |
| mr_spades.non-contained  |  189483 | 1623035 |  15 |
| megahit.contig           |  115696 | 1622423 |  56 |
| megahit.non-contained    |  115696 | 1607271 |  22 |
| mr_megahit.contig        |  174583 | 1631360 |  40 |
| mr_megahit.non-contained |  174583 | 1622195 |  22 |


# Clostridio_dif_630

```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Clostridio_dif_630

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4274782 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --sample 300 \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 101 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md

# rm -fr 4_down_sampling 6_down_sampling

bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"

# bash 0_master.sh
# bash 0_cleanup.sh

```


Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 692.5 |    320 | 3766.8 |                         48.50% |
| R.tadpole.bbtools | 328.8 |    317 |   97.8 |                         45.98% |
| R.genome.picard   | 330.2 |    319 |   87.9 |                             FR |
| R.tadpole.picard  | 328.5 |    317 |   87.9 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 234       | 4124061         | 0.0000       | 34.16   |
| R.31 | 203       | 4243015         | 0.0000       | 33.33   |
| R.41 | 172       | 4257860         | 0.0000       | 32.84   |
| R.51 | 146       | 4281513         | 0.0000       | 32.49   |
| R.61 | 116       | 4264133         | 0.0000       | 32.22   |
| R.71 | 87        | 4371307         | 0.0000       | 31.97   |
| R.81 | 58        | 4372831         | 0.0000       | 31.69   |


Table: statReads

| Name        |     N50 |     Sum |        # |
|:------------|--------:|--------:|---------:|
| Genome      | 4274782 | 4274782 |        1 |
| Paralogs    |    3253 |  381668 |      135 |
| Repetitives |    2524 |  165378 |      216 |
| Illumina.R  |     101 |   1.33G | 13190786 |
| trim.R      |     100 |   1.22G | 12274774 |
| Q25L60      |     100 |   1.19G | 11965165 |
| Q30L60      |     100 |   1.13G | 11482401 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.32G | 13028932 |
| highpass | 101 |   1.31G | 12944878 |
| sample   | 101 |   1.28G | 12697372 |
| trim     | 100 |   1.22G | 12275340 |
| filter   | 100 |   1.22G | 12274774 |
| R1       | 100 | 612.43M |  6137387 |
| R2       | 100 | 607.73M |  6137387 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6602	0.05200%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	284	0.00231%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	23096683
#error_kmers	18905712
#genomic_kmers	4190971
#main_peak	186
#genome_size_in_peaks	4319560
#genome_size	4350907
#haploid_genome_size	4350907
#fold_coverage	186
#haploid_fold_coverage	186
#ploidy	1
#percent_repeat_in_peaks	2.978
#percent_repeat	3.163
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |        # |
|:--------------|----:|--------:|---------:|
| clumped       | 100 |   1.22G | 12274196 |
| ecco          | 100 |   1.22G | 12274196 |
| eccc          | 100 |   1.22G | 12274196 |
| ecct          | 100 |   1.21G | 12211720 |
| extended      | 140 |    1.7G | 12211720 |
| merged.raw    | 343 |   1.52G |  4580735 |
| unmerged.raw  | 140 | 421.03M |  3050250 |
| unmerged.trim | 140 | 421.03M |  3050202 |
| M1            | 343 |   1.51G |  4559069 |
| U1            | 140 |  212.5M |  1525101 |
| U2            | 140 | 208.53M |  1525101 |
| Us            |   0 |       0 |        0 |
| M.cor         | 321 |   1.94G | 12168340 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 165.4 |    170 |  18.1 |          1.95% |
| M.ihist.merge.txt  | 331.1 |    333 |  54.3 |         75.02% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 285.4 |  266.1 |    6.76% | "71" | 4.27M | 4.22M |     0.99 | 0:02'26'' |
| Q25L60.R | 278.0 |  263.2 |    5.35% | "71" | 4.27M |  4.2M |     0.98 | 0:02'33'' |
| Q30L60.R | 264.5 |  253.7 |    4.10% | "71" | 4.27M |  4.2M |     0.98 | 0:02'23'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.37% |     80154 | 4.17M | 133 |        32 | 21.38K | 721 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  98.47% |     76246 | 4.17M | 132 |        36 | 23.18K | 740 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:38 |
| Q0L0X40P002   |   40.0 |  98.33% |     75199 | 4.17M | 128 |        33 | 21.91K | 714 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:38 |
| Q0L0X80P000   |   80.0 |  97.92% |     58767 | 4.16M | 152 |      1048 | 17.67K | 432 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:34 |
| Q0L0X80P001   |   80.0 |  97.82% |     65837 | 4.16M | 144 |        44 | 15.96K | 413 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:33 |
| Q0L0X80P002   |   80.0 |  97.78% |     63587 | 4.16M | 150 |        42 | 15.63K | 408 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:34 |
| Q25L60X40P000 |   40.0 |  98.44% |     78939 | 4.17M | 126 |      1006 | 23.63K | 669 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:36 |
| Q25L60X40P001 |   40.0 |  98.47% |     85271 | 4.17M | 125 |        41 | 20.56K | 639 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q25L60X40P002 |   40.0 |  98.44% |     85269 | 4.17M | 127 |        43 | 22.21K | 675 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:38 |
| Q25L60X80P000 |   80.0 |  97.89% |     73824 | 4.16M | 131 |      1239 | 16.76K | 379 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  97.95% |     69165 | 4.16M | 137 |        40 | 15.09K | 396 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:37 |
| Q25L60X80P002 |   80.0 |  98.05% |     80158 | 4.16M | 137 |      1048 | 17.66K | 395 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:37 |
| Q30L60X40P000 |   40.0 |  98.58% |     69400 | 4.17M | 124 |        32 | 20.69K | 706 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:39 |
| Q30L60X40P001 |   40.0 |  98.54% |     80184 | 4.17M | 126 |       976 |  25.4K | 686 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:39 |
| Q30L60X40P002 |   40.0 |  98.56% |     85262 | 4.16M | 119 |        33 | 20.79K | 686 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:41 |
| Q30L60X80P000 |   80.0 |  98.41% |     69419 | 4.16M | 125 |       754 | 18.27K | 403 |   77.0 | 10.0 |  15.7 | 214.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:44 |
| Q30L60X80P001 |   80.0 |  98.36% |     85687 | 4.16M | 119 |      1038 | 17.71K | 404 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:45 |
| Q30L60X80P002 |   80.0 |  98.38% |     78918 | 4.16M | 126 |      1038 | 18.89K | 419 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:38 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.46% |     98203 | 4.16M | 103 |      1052 | 21.67K | 302 |   38.0 | 6.0 |   6.7 | 112.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:39 |
| MRX40P001 |   40.0 |  96.85% |     85682 | 4.16M | 108 |        55 | 11.03K | 307 |   39.0 | 6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:38 |
| MRX40P002 |   40.0 |  96.90% |     85688 | 4.16M | 112 |        59 | 13.26K | 328 |   39.0 | 6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:35 |
| MRX80P000 |   80.0 |  97.25% |     80170 | 4.15M | 120 |      1052 |  26.1K | 334 |   77.0 | 9.0 |  16.7 | 208.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:42 |
| MRX80P001 |   80.0 |  97.31% |     80114 | 4.15M | 124 |      1012 | 25.95K | 352 |   78.0 | 9.0 |  17.0 | 210.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:36 |
| MRX80P002 |   80.0 |  97.29% |     71946 | 4.15M | 122 |       403 | 24.93K | 339 |   77.0 | 9.0 |  16.7 | 208.0 | "31,41,51,61,71,81" |   0:01:36 |   0:00:39 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.24% |     80126 | 4.16M | 144 |      7921 | 67.14K | 738 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:18 |   0:00:37 |
| Q0L0X40P001   |   40.0 |  98.26% |     76236 | 4.16M | 146 |      7921 | 64.48K | 789 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:37 |
| Q0L0X40P002   |   40.0 |  98.28% |     80057 | 4.16M | 143 |      7921 | 64.56K | 731 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  98.12% |    105486 | 4.16M | 112 |      7941 | 35.35K | 481 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  98.15% |     85670 | 4.16M | 111 |        63 | 19.85K | 470 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:35 |   0:00:39 |
| Q0L0X80P002   |   80.0 |  98.14% |     88764 | 4.16M | 108 |      1053 | 22.69K | 453 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  98.22% |     80067 |  4.2M | 142 |        33 | 21.58K | 722 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:37 |
| Q25L60X40P001 |   40.0 |  98.31% |     85648 | 4.17M | 141 |      7951 | 31.81K | 731 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:42 |
| Q25L60X40P002 |   40.0 |  98.26% |     78845 | 4.16M | 142 |      7921 | 66.16K | 700 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:44 |
| Q25L60X80P000 |   80.0 |  98.35% |     96039 | 4.16M | 109 |      1038 | 19.91K | 452 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:27 |   0:00:40 |
| Q25L60X80P001 |   80.0 |  98.23% |     88773 | 4.16M | 114 |      7951 | 27.13K | 473 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:43 |
| Q25L60X80P002 |   80.0 |  98.28% |    108211 | 4.16M | 109 |      7941 | 44.19K | 496 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:22 |   0:00:42 |
| Q30L60X40P000 |   40.0 |  98.36% |     82456 | 4.16M | 151 |      7931 | 61.11K | 749 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:14 |   0:00:38 |
| Q30L60X40P001 |   40.0 |  98.34% |     80156 | 4.17M | 145 |      7931 |  41.9K | 751 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:38 |
| Q30L60X40P002 |   40.0 |  98.25% |     85211 | 4.15M | 143 |      1031 | 33.43K | 719 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:36 |
| Q30L60X80P000 |   80.0 |  98.48% |     85663 | 4.16M | 115 |      7941 | 37.04K | 522 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:26 |   0:00:43 |
| Q30L60X80P001 |   80.0 |  98.46% |    104331 | 4.16M | 113 |        56 | 20.77K | 528 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:42 |
| Q30L60X80P002 |   80.0 |  98.42% |     85261 | 4.16M | 118 |      7951 | 31.57K | 519 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:24 |   0:00:40 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.81% |    105458 | 4.16M |  97 |      7931 | 57.44K | 230 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:32 |
| MRX40P001 |   40.0 |  96.82% |    105455 | 4.16M |  99 |      7931 | 56.67K | 234 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:31 |
| MRX40P002 |   40.0 |  96.82% |    105399 | 4.16M | 101 |      7931 | 50.07K | 249 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:15 |   0:00:31 |
| MRX80P000 |   80.0 |  96.80% |    102976 | 4.16M | 101 |      7921 | 51.63K | 244 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:33 |
| MRX80P001 |   80.0 |  96.81% |     85706 | 4.16M | 108 |      7931 | 58.51K | 250 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:37 |   0:00:35 |
| MRX80P002 |   80.0 |  96.81% |    105454 | 4.16M | 106 |      7921 | 52.13K | 246 |   79.0 |  9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:32 |   0:00:34 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.21% |     88772 | 4.17M | 127 |      7932 | 61.11K | 625 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  98.39% |     85665 | 4.16M | 128 |      7932 |  62.2K | 657 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q0L0X40P002   |   40.0 |  98.25% |     85208 | 4.17M | 129 |      7922 | 64.64K | 610 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:35 |
| Q0L0X80P000   |   80.0 |  98.22% |    108216 | 4.18M | 103 |      7932 | 57.76K | 352 |   78.0 |  9.0 |  17.0 | 210.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:38 |
| Q0L0X80P001   |   80.0 |  98.30% |     96060 | 4.16M | 100 |      7932 | 47.48K | 366 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:36 |
| Q0L0X80P002   |   80.0 |  98.32% |    105477 | 4.16M | 100 |      7932 | 56.04K | 365 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  98.26% |     88667 | 4.17M | 128 |      7929 | 61.64K | 633 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:34 |
| Q25L60X40P001 |   40.0 |  98.41% |     85648 | 4.16M | 124 |      7932 | 54.56K | 611 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:36 |
| Q25L60X40P002 |   40.0 |  98.42% |     88680 | 4.16M | 124 |      7922 | 65.97K | 620 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:37 |
| Q25L60X80P000 |   80.0 |  98.19% |    105496 | 4.16M | 100 |      7932 | 56.98K | 344 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  98.43% |    105489 | 4.16M | 102 |      7932 | 58.24K | 382 |   79.0 |  9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:40 |
| Q25L60X80P002 |   80.0 |  98.38% |    108211 | 4.16M |  98 |      7932 | 56.48K | 383 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |
| Q30L60X40P000 |   40.0 |  98.44% |     85664 | 4.17M | 139 |      7932 | 61.41K | 654 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:36 |
| Q30L60X40P001 |   40.0 |  98.45% |     82496 | 4.16M | 139 |      7922 | 66.63K | 659 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:36 |
| Q30L60X40P002 |   40.0 |  98.41% |     85234 | 4.16M | 136 |      7922 | 64.58K | 641 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:35 |
| Q30L60X80P000 |   80.0 |  98.49% |    104270 | 4.16M | 105 |      7932 | 59.58K | 433 |   79.0 |  9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:46 |
| Q30L60X80P001 |   80.0 |  98.55% |    105480 | 4.16M | 106 |      7932 | 51.47K | 433 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:47 |
| Q30L60X80P002 |   80.0 |  98.51% |     88768 | 4.16M | 105 |      7932 | 63.28K | 423 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.98% |    108206 | 4.16M |  93 |      7932 | 57.72K | 235 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:31 |
| MRX40P001 |   40.0 |  96.84% |    108181 | 4.16M |  96 |      7932 | 56.73K | 234 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| MRX40P002 |   40.0 |  96.95% |    108217 | 4.16M |  93 |      7932 | 57.04K | 224 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:32 |
| MRX80P000 |   80.0 |  96.82% |    102976 | 4.16M |  98 |      7932 | 59.39K | 243 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:37 |
| MRX80P001 |   80.0 |  96.85% |    108141 | 4.16M |  99 |      7932 |  59.4K | 243 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:35 |
| MRX80P002 |   80.0 |  96.90% |    105447 | 4.16M | 100 |      7932 | 60.91K | 246 |   79.0 |  9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:34 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.91% |    108213 | 4.17M | 95 |      7942 |   3.02M | 415 |  259.0 | 26.0 |  60.3 | 674.0 |   0:01:07 |
| 7_merge_mr_unitigs_bcalm      |  98.64% |    108241 | 4.16M | 97 |      7941 | 584.31K |  80 |  258.0 | 28.0 |  58.0 | 684.0 |   0:01:27 |
| 7_merge_mr_unitigs_superreads |  98.65% |    108240 | 4.16M | 97 |      7941 |  46.87K |  21 |  258.0 | 29.0 |  57.0 | 690.0 |   0:01:30 |
| 7_merge_mr_unitigs_tadpole    |  98.60% |    108233 | 4.16M | 96 |      7942 | 606.58K |  85 |  258.0 | 28.0 |  58.0 | 684.0 |   0:01:25 |
| 7_merge_unitigs_bcalm         |  98.42% |    108230 | 4.17M | 96 |      7951 | 509.54K |  92 |  264.0 | 28.0 |  60.0 | 696.0 |   0:01:18 |
| 7_merge_unitigs_superreads    |  98.70% |    108236 | 4.17M | 94 |      7951 | 113.45K |  31 |  260.0 | 28.0 |  58.7 | 688.0 |   0:01:45 |
| 7_merge_unitigs_tadpole       |  98.73% |    108225 | 4.19M | 97 |      7942 |   1.37M | 191 |  262.0 | 27.0 |  60.3 | 686.0 |   0:01:32 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  98.93% |    227375 | 4.18M | 47 |      3005 |    19K |  94 |  261.0 | 28.0 |  59.0 |  690.0 |   0:00:45 |
| 8_mr_spades  |  98.72% |    225840 | 4.19M | 38 |      2925 | 16.03K |  79 |  445.0 | 45.0 | 103.3 | 1160.0 |   0:00:48 |
| 8_megahit    |  98.61% |    108264 | 4.17M | 89 |      3049 | 17.69K | 179 |  261.0 | 28.0 |  59.0 |  690.0 |   0:00:47 |
| 8_mr_megahit |  99.29% |    225939 |  4.2M | 51 |      4720 | 22.05K | 104 |  445.0 | 45.0 | 103.3 | 1160.0 |   0:00:50 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 4274782 | 4274782 |   1 |
| Paralogs                 |    3253 |  381668 | 135 |
| Repetitives              |    2524 |  165378 | 216 |
| 7_merge_anchors.anchors  |  108213 | 4170223 |  95 |
| 7_merge_anchors.others   |    7942 | 3023532 | 415 |
| glue_anchors             |  211988 | 4168349 |  71 |
| fill_anchors             |  213521 | 4165517 |  62 |
| spades.contig            |  227405 | 4223207 | 207 |
| spades.scaffold          |  227405 | 4223264 | 205 |
| spades.non-contained     |  227405 | 4194192 |  47 |
| mr_spades.contig         |  225584 | 4221666 |  85 |
| mr_spades.scaffold       |  225921 | 4221766 |  84 |
| mr_spades.non-contained  |  225584 | 4204871 |  41 |
| megahit.contig           |  108291 | 4208909 | 148 |
| megahit.non-contained    |  108291 | 4184464 |  90 |
| mr_megahit.contig        |  226020 | 4274703 | 177 |
| mr_megahit.non-contained |  226020 | 4226796 |  53 |


# Co_dip_NCTC_13129

```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Co_dip_NCTC_13129

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 2488635 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --sample 300 \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 101 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md

# rm -fr 4_down_sampling 6_down_sampling

bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"

# bash 0_master.sh
# bash 0_cleanup.sh

```

Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 515.0 |    371 | 2326.7 |                         49.30% |
| R.tadpole.bbtools | 372.2 |    366 |  105.4 |                         43.35% |
| R.genome.picard   | 377.4 |    371 |  105.0 |                             FR |
| R.tadpole.picard  | 372.0 |    366 |  104.5 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 334       | 2415175         | 0.0000       | 54.27   |
| R.31 | 282       | 2415575         | 0.0000       | 54.29   |
| R.41 | 236       | 2415448         | 0.0000       | 54.28   |
| R.51 | 193       | 2412947         | 0.0000       | 54.25   |
| R.61 | 153       | 2419302         | 0.0000       | 54.21   |
| R.71 | 114       | 2423650         | 0.0000       | 54.16   |
| R.81 | 77        | 2450423         | 0.0000       | 54.11   |


Table: statReads

| Name        |     N50 |     Sum |        # |
|:------------|--------:|--------:|---------:|
| Genome      | 2488635 | 2488635 |        1 |
| Paralogs    |    5627 |   56034 |       18 |
| Repetitives |    2784 |   51421 |       70 |
| Illumina.R  |     101 |   1.12G | 11128812 |
| trim.R      |     100 |  655.6M |  6697056 |
| Q25L60      |     100 | 576.46M |  5975527 |
| Q30L60      |     100 | 445.37M |  4823284 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.12G | 11094550 |
| highpass | 101 |    1.1G | 10890256 |
| sample   | 101 | 746.59M |  7391986 |
| trim     | 100 |  655.6M |  6697056 |
| filter   | 100 |  655.6M |  6697056 |
| R1       | 100 | 328.82M |  3348528 |
| R2       | 100 | 326.78M |  3348528 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6007	0.08126%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	0	0.00000%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	18230195
#error_kmers	15783162
#genomic_kmers	2447033
#main_peak	174
#genome_size_in_peaks	2484697
#genome_size	2494569
#haploid_genome_size	2494569
#fold_coverage	174
#haploid_fold_coverage	174
#ploidy	1
#percent_repeat_in_peaks	1.518
#percent_repeat	1.561
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 655.59M | 6696926 |
| ecco          | 100 | 655.58M | 6696926 |
| eccc          | 100 | 655.58M | 6696926 |
| ecct          | 100 | 652.94M | 6669166 |
| extended      | 140 | 918.39M | 6669166 |
| merged.raw    | 365 | 648.18M | 1886058 |
| unmerged.raw  | 140 | 395.18M | 2897050 |
| unmerged.trim | 140 | 395.16M | 2896906 |
| M1            | 365 | 647.14M | 1883068 |
| U1            | 140 | 198.53M | 1448453 |
| U2            | 140 | 196.62M | 1448453 |
| Us            |   0 |       0 |       0 |
| M.cor         | 310 |   1.04G | 6663042 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 153.8 |    158 |  23.5 |          2.86% |
| M.ihist.merge.txt  | 343.7 |    355 |  60.4 |         56.56% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 263.4 |  234.7 |   10.91% | "51" | 2.49M |  2.7M |     1.08 | 0:01'23'' |
| Q25L60.R | 231.8 |  213.6 |    7.87% | "51" | 2.49M | 2.47M |     0.99 | 0:01'17'' |
| Q30L60.R | 179.3 |  169.1 |    5.69% | "49" | 2.49M | 2.45M |     0.98 | 0:01'06'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.25% |     19831 | 2.49M | 211 |        23 |  21.9K | 1043 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  97.96% |     18639 | 2.48M | 206 |        22 | 19.58K |  965 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:24 |
| Q0L0X40P002   |   40.0 |  98.19% |     19009 | 2.47M | 199 |        25 | 24.38K |  995 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  97.30% |     13400 | 2.44M | 261 |        26 | 18.88K |  707 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  96.86% |     14762 | 2.44M | 256 |        27 | 18.61K |  658 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  99.15% |     46819 | 2.46M | 109 |        32 | 17.87K |  623 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  98.97% |     46852 | 2.45M | 104 |        30 | 16.22K |  597 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:25 |
| Q25L60X40P002 |   40.0 |  99.14% |     55175 | 2.45M |  95 |        40 |  19.4K |  597 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:25 |
| Q25L60X80P000 |   80.0 |  98.70% |     35735 | 2.44M | 106 |        40 | 13.72K |  404 |   79.0 | 11.0 |  15.3 | 224.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q25L60X80P001 |   80.0 |  98.98% |     35264 | 2.44M | 118 |        37 | 13.77K |  442 |   78.0 | 11.0 |  15.0 | 222.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  99.30% |     55143 | 2.44M |  79 |        41 | 13.06K |  542 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:19 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.33% |     57179 | 2.44M |  81 |        58 | 18.77K |  609 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:19 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  99.27% |     55138 | 2.44M |  84 |        40 | 13.23K |  556 |   39.0 |  8.0 |   5.0 | 126.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  99.37% |     67871 | 2.44M |  64 |        37 | 13.75K |  468 |   79.0 | 15.0 |  11.3 | 248.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:32 |
| Q30L60X80P001 |   80.0 |  99.32% |     65436 | 2.44M |  66 |      1023 | 18.16K |  428 |   78.0 | 15.0 |  11.0 | 246.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:31 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.59% |     67898 | 2.44M | 62 |        92 | 11.14K | 171 |   39.0 | 5.0 |   8.0 | 108.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| MRX40P001 |   40.0 |  98.74% |     62967 | 2.44M | 64 |       105 | 12.66K | 177 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| MRX40P002 |   40.0 |  98.71% |     71789 | 2.44M | 61 |      1016 | 11.62K | 175 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:25 |
| MRX80P000 |   80.0 |  98.36% |     39636 | 2.44M | 92 |        82 | 12.87K | 221 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| MRX80P001 |   80.0 |  98.45% |     48468 | 2.43M | 83 |        96 | 13.51K | 217 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| MRX80P002 |   80.0 |  98.53% |     47725 | 2.43M | 84 |        81 | 12.88K | 213 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.63% |     34458 | 2.45M | 114 |        31 | 21.28K | 740 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.63% |     34420 | 2.46M | 115 |        35 | 21.76K | 725 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  98.64% |     40334 | 2.45M | 116 |        42 | 23.52K | 697 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  98.50% |     37648 | 2.45M | 113 |        32 |  20.3K | 617 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:28 |
| Q0L0X80P001   |   80.0 |  98.37% |     34355 | 2.45M | 110 |        39 | 23.48K | 621 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:27 |
| Q25L60X40P000 |   40.0 |  99.23% |     51095 | 2.44M |  85 |        33 | 19.62K | 611 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  99.18% |     54214 | 2.44M |  93 |        37 | 20.83K | 625 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  99.19% |     57368 | 2.45M |  89 |        38 | 20.92K | 622 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  99.11% |     58196 | 2.45M |  85 |        39 |  18.1K | 516 |   79.0 | 11.0 |  15.3 | 224.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:28 |
| Q25L60X80P001 |   80.0 |  99.16% |     57866 | 2.45M |  79 |        34 | 16.27K | 523 |   79.0 | 11.0 |  15.3 | 224.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  99.12% |     43257 | 2.44M | 113 |      1000 | 19.81K | 664 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:01:16 |   0:00:25 |
| Q30L60X40P001 |   40.0 |  99.15% |     40031 | 2.44M | 124 |        45 | 20.97K | 706 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:25 |
| Q30L60X40P002 |   40.0 |  99.17% |     31199 | 2.44M | 132 |        28 | 16.57K | 704 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:01:01 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  99.36% |     57593 | 2.44M |  78 |        31 | 16.28K | 594 |   79.0 | 15.0 |  11.3 | 248.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  99.30% |     60991 | 2.44M |  77 |        36 | 15.21K | 550 |   80.0 | 15.0 |  11.7 | 250.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:30 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |   Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.55% |    101766 | 2.44M | 49 |      2470 | 9.24K | 125 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:24 |
| MRX40P001 |   40.0 |  98.47% |     75127 | 2.44M | 53 |      2470 | 9.33K | 117 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:23 |
| MRX40P002 |   40.0 |  98.43% |    107957 | 2.44M | 48 |      2470 | 8.47K | 104 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:23 |
| MRX80P000 |   80.0 |  98.31% |     62965 | 2.44M | 62 |      2470 | 9.77K | 137 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:24 |
| MRX80P001 |   80.0 |  98.32% |     63996 | 2.44M | 65 |      2470 | 9.28K | 139 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:24 |
| MRX80P002 |   80.0 |  98.39% |     97931 | 2.44M | 55 |      2470 | 9.87K | 124 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:24 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.91% |     63763 | 2.44M |  71 |        34 | 15.09K | 477 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q0L0X40P001   |   40.0 |  98.80% |     63747 | 2.44M |  72 |        38 | 14.12K | 432 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  98.82% |     57320 | 2.44M |  75 |        44 | 16.11K | 470 |   39.0 |  6.0 |   7.0 | 114.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  98.65% |     95729 | 2.44M |  58 |        53 | 13.21K | 308 |   78.0 | 10.0 |  16.0 | 216.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  98.81% |     69302 | 2.44M |  59 |        46 | 15.12K | 357 |   79.0 | 10.0 |  16.3 | 218.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:29 |
| Q25L60X40P000 |   40.0 |  99.35% |     97903 | 2.44M |  63 |        57 | 20.43K | 508 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q25L60X40P001 |   40.0 |  99.26% |     61645 | 2.44M |  77 |        46 | 18.89K | 504 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  99.14% |     61719 | 2.44M |  73 |        37 | 14.77K | 447 |   40.0 |  6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  99.27% |     65328 | 2.44M |  60 |        44 | 14.45K | 363 |   80.0 | 11.0 |  15.7 | 226.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:29 |
| Q25L60X80P001 |   80.0 |  99.26% |     88667 | 2.44M |  61 |        41 | 13.67K | 378 |   79.0 | 11.0 |  15.3 | 224.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  99.21% |     48292 | 2.44M |  93 |        34 |  15.8K | 580 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.31% |     54019 | 2.44M | 100 |        61 | 20.86K | 656 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  99.26% |     46841 | 2.45M |  99 |        29 | 14.82K | 611 |   40.0 |  8.0 |   5.3 | 128.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:27 |
| Q30L60X80P000 |   80.0 |  99.40% |     65321 | 2.44M |  66 |        33 | 14.51K | 514 |   79.0 | 15.0 |  11.3 | 248.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  99.31% |     65329 | 2.44M |  68 |        42 | 13.98K | 458 |   79.0 | 15.0 |  11.3 | 248.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:31 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |   Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.50% |     97904 | 2.44M | 48 |      2468 | 8.94K | 114 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:24 |
| MRX40P001 |   40.0 |  98.47% |     75167 | 2.44M | 52 |      2468 | 9.28K | 113 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| MRX40P002 |   40.0 |  98.44% |    115898 | 2.44M | 45 |      2468 | 8.27K |  98 |   40.0 | 6.0 |   7.3 | 116.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:23 |
| MRX80P000 |   80.0 |  98.38% |     75181 | 2.44M | 52 |      2468 | 8.61K | 111 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:24 |
| MRX80P001 |   80.0 |  98.37% |     97182 | 2.44M | 54 |      2468 | 8.79K | 117 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:23 |
| MRX80P002 |   80.0 |  98.45% |    104194 | 2.44M | 50 |      2468 | 9.29K | 112 |   79.0 | 9.0 |  17.3 | 212.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:25 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  98.16% |    115935 | 2.44M | 46 |     15531 | 89.57K | 23 |  232.0 | 26.0 |  51.3 | 620.0 |   0:00:36 |
| 7_merge_mr_unitigs_bcalm      |  98.98% |    115936 | 2.44M | 46 |      2530 |     5K |  2 |  229.0 | 25.0 |  51.3 | 608.0 |   0:00:55 |
| 7_merge_mr_unitigs_superreads |  98.91% |    115933 | 2.44M | 46 |      8066 | 38.38K |  7 |  232.0 | 26.0 |  51.3 | 620.0 |   0:00:53 |
| 7_merge_mr_unitigs_tadpole    |  98.87% |    115932 | 2.44M | 46 |      2528 |     5K |  2 |  229.0 | 26.0 |  50.3 | 614.0 |   0:00:50 |
| 7_merge_unitigs_bcalm         |  98.96% |    115924 | 2.44M | 49 |      8923 |    50K | 15 |  230.0 | 25.0 |  51.7 | 610.0 |   0:00:59 |
| 7_merge_unitigs_superreads    |  98.84% |    115926 | 2.44M | 46 |      2530 | 28.54K | 15 |  231.0 | 25.0 |  52.0 | 612.0 |   0:00:49 |
| 7_merge_unitigs_tadpole       |  98.94% |    115927 | 2.44M | 48 |     40450 | 55.66K | 13 |  233.0 | 26.0 |  51.7 | 622.0 |   0:00:54 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  98.87% |    177182 | 2.44M | 23 |      5490 |  9.93K | 44 |  232.0 | 25.0 |  52.3 |  614.0 |   0:00:30 |
| 8_mr_spades  |  99.43% |    168252 | 2.45M | 20 |      5258 |  9.59K | 39 |  418.0 | 39.0 | 100.3 | 1070.0 |   0:00:33 |
| 8_megahit    |  98.80% |    115908 | 2.44M | 47 |      5490 | 10.98K | 91 |  232.0 | 25.0 |  52.3 |  614.0 |   0:00:29 |
| 8_mr_megahit |  99.61% |    177840 | 2.45M | 22 |      5615 |  9.97K | 43 |  418.0 | 39.0 | 100.3 | 1070.0 |   0:00:33 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 2488635 | 2488635 |   1 |
| Paralogs                 |    5627 |   56034 |  18 |
| Repetitives              |    2784 |   51421 |  70 |
| 7_merge_anchors.anchors  |  115935 | 2441270 |  46 |
| 7_merge_anchors.others   |   15531 |   89571 |  23 |
| glue_anchors             |  130277 | 2440700 |  35 |
| fill_anchors             |  165289 | 2441378 |  24 |
| spades.contig            |  310298 | 2456423 |  28 |
| spades.scaffold          |  310298 | 2456623 |  26 |
| spades.non-contained     |  310298 | 2453467 |  21 |
| mr_spades.contig         |  214636 | 2457095 |  22 |
| mr_spades.scaffold       |  214636 | 2457095 |  22 |
| mr_spades.non-contained  |  214636 | 2455082 |  19 |
| megahit.contig           |  115950 | 2455914 |  54 |
| megahit.non-contained    |  115950 | 2451627 |  44 |
| mr_megahit.contig        |  309887 | 2468408 |  42 |
| mr_megahit.non-contained |  309887 | 2459257 |  21 |


# Fr_tul_tularensis_SCHU_S4

```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Fr_tul_tularensis_SCHU_S4

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1892775 \
    --parallel 24 \
    --xmx 80g \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --sample 300 \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "superreads bcalm tadpole" \
    --statp 2 \
    --readl 101 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend

# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md

# rm -fr 4_down_sampling 6_down_sampling

bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"

# bash 0_master.sh
# bash 0_cleanup.sh

```


Table: statInsertSize

| Group             |  Mean | Median |  STDev | PercentOfPairs/PairOrientation |
|:------------------|------:|-------:|-------:|-------------------------------:|
| R.genome.bbtools  | 550.6 |    365 | 2622.0 |                         48.72% |
| R.tadpole.bbtools | 373.9 |    363 |  120.1 |                         47.64% |
| R.genome.picard   | 375.6 |    365 |  108.9 |                             FR |
| R.tadpole.picard  | 373.9 |    363 |  109.1 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 0         | 0               | 0.0000       | 35.05   |
| R.31 | 0         | 0               | 0.0000       | 34.36   |
| R.41 | 0         | 0               | 0.0000       | 34.04   |
| R.51 | 0         | 0               | 0.0000       | 33.83   |
| R.61 | 442       | 1598469         | 0.0000       | 33.67   |
| R.71 | 334       | 1751732         | 0.0000       | 33.53   |
| R.81 | 229       | 1840179         | 0.0000       | 33.39   |


Table: statReads

| Name        |     N50 |     Sum |        # |
|:------------|--------:|--------:|---------:|
| Genome      | 1892775 | 1892775 |        1 |
| Paralogs    |   33912 |   93528 |       10 |
| Repetitives |    5357 |  139352 |       72 |
| Illumina.R  |     101 |   2.14G | 21230270 |
| trim.R      |     100 | 549.71M |  5517640 |
| Q25L60      |     100 |  538.7M |  5416004 |
| Q30L60      |     100 | 516.22M |  5233886 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   2.12G | 21017824 |
| highpass | 101 |   2.12G | 20974018 |
| sample   | 101 | 567.83M |  5622104 |
| trim     | 100 | 549.79M |  5518442 |
| filter   | 100 | 549.71M |  5517640 |
| R1       | 100 | 275.12M |  2758820 |
| R2       | 100 | 274.59M |  2758820 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3641	0.06476%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	802	0.01453%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	802	0.01453%
```

```text
#R.peaks
#k	31
#unique_kmers	7538588
#error_kmers	5745386
#genomic_kmers	1793202
#main_peak	201
#genome_size_in_peaks	1840501
#genome_size	1884136
#haploid_genome_size	1884136
#fold_coverage	201
#haploid_fold_coverage	201
#ploidy	1
#percent_repeat_in_peaks	2.633
#percent_repeat	4.442
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 549.71M | 5517600 |
| ecco          | 100 | 549.71M | 5517600 |
| eccc          | 100 | 549.71M | 5517600 |
| ecct          | 100 | 547.35M | 5493616 |
| extended      | 140 | 766.63M | 5493616 |
| merged.raw    | 361 |  552.3M | 1616623 |
| unmerged.raw  | 140 |  314.8M | 2260370 |
| unmerged.trim | 140 |  314.8M | 2260366 |
| M1            | 361 | 551.61M | 1614623 |
| U1            | 140 | 157.68M | 1130183 |
| U2            | 140 | 157.12M | 1130183 |
| Us            |   0 |       0 |       0 |
| M.cor         | 311 | 868.02M | 5489612 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 154.2 |    159 |  23.6 |          2.73% |
| M.ihist.merge.txt  | 341.6 |    350 |  59.3 |         58.85% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 290.4 |  275.8 |    5.04% | "71" | 1.89M | 1.81M |     0.95 | 0:01'11'' |
| Q25L60.R | 284.6 |  272.8 |    4.15% | "71" | 1.89M |  1.8M |     0.95 | 0:01'11'' |
| Q30L60.R | 272.8 |  263.7 |    3.35% | "71" | 1.89M |  1.8M |     0.95 | 0:01'07'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.85% |     32367 | 1.79M | 82 |        31 | 13.58K | 401 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  97.85% |     32713 | 1.79M | 75 |        43 | 14.31K | 367 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  97.98% |     32716 | 1.79M | 76 |        29 | 12.32K | 395 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  97.11% |     32699 | 1.79M | 78 |      1712 |  9.57K | 177 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  96.99% |     32712 | 1.79M | 77 |      1070 |  7.85K | 174 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| Q0L0X80P002   |   80.0 |  97.16% |     32353 | 1.79M | 78 |      1712 |  9.85K | 181 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.16% |     32700 | 1.79M | 74 |        30 | 12.59K | 401 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  97.89% |     32345 |  1.8M | 78 |        31 | 12.57K | 367 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:23 |
| Q25L60X40P002 |   40.0 |  97.97% |     32701 |  1.8M | 76 |        26 | 10.69K | 385 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:23 |
| Q25L60X80P000 |   80.0 |  97.15% |     32350 | 1.79M | 77 |      1712 |  9.59K | 177 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  97.08% |     32698 | 1.79M | 75 |        39 |  7.75K | 174 |   81.0 | 7.0 |  20.0 | 204.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  97.08% |     32704 | 1.79M | 78 |      1362 |  9.04K | 177 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.05% |     32763 | 1.79M | 72 |        25 |  10.3K | 401 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q30L60X40P001 |   40.0 |  97.98% |     32693 |  1.8M | 76 |        26 | 10.24K | 378 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q30L60X40P002 |   40.0 |  97.98% |     32714 | 1.79M | 77 |        26 | 10.33K | 393 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q30L60X80P000 |   80.0 |  97.18% |     32700 | 1.79M | 73 |        31 |  6.84K | 170 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  97.09% |     32707 | 1.79M | 75 |       802 |  8.12K | 174 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q30L60X80P002 |   80.0 |  97.25% |     32708 | 1.79M | 73 |        34 |   6.9K | 170 |   81.0 | 7.0 |  20.0 | 204.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.66% |     32363 | 1.79M | 73 |        74 | 10.53K | 166 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| MRX40P001 |   40.0 |  96.61% |     32646 | 1.79M | 72 |        62 |  8.39K | 165 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:23 |
| MRX40P002 |   40.0 |  96.69% |     32677 | 1.79M | 72 |        69 |  9.06K | 164 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:23 |
| MRX80P000 |   80.0 |  96.61% |     32335 | 1.79M | 74 |        68 | 11.25K | 168 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:23 |
| MRX80P001 |   80.0 |  96.59% |     32683 | 1.79M | 75 |        65 |  9.29K | 167 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:23 |
| MRX80P002 |   80.0 |  96.60% |     32327 | 1.79M | 74 |        76 |    12K | 169 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:24 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.05% |     32675 | 1.79M | 77 |        31 |  15.5K | 451 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:23 |
| Q0L0X40P001   |   40.0 |  97.97% |     32697 | 1.79M | 77 |        32 | 15.38K | 437 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:23 |
| Q0L0X40P002   |   40.0 |  97.79% |     32693 | 1.79M | 75 |        31 | 13.98K | 398 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:23 |
| Q0L0X80P000   |   80.0 |  97.69% |     32699 | 1.79M | 72 |      1030 | 11.92K | 232 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  97.65% |     32712 | 1.79M | 73 |        34 |  9.25K | 230 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  97.58% |     32713 | 1.79M | 73 |       521 | 11.35K | 239 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:25 |
| Q25L60X40P000 |   40.0 |  98.26% |     32680 | 1.79M | 75 |        30 | 14.69K | 452 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  98.00% |     32694 | 1.79M | 74 |        30 | 13.04K | 393 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  97.82% |     32683 | 1.79M | 76 |        26 | 10.81K | 391 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:23 |
| Q25L60X80P000 |   80.0 |  97.72% |     32716 | 1.79M | 73 |      1033 | 12.08K | 234 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  97.80% |     32691 | 1.79M | 71 |        33 | 11.26K | 271 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  97.66% |     32704 | 1.79M | 71 |       349 | 10.89K | 235 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q30L60X40P000 |   40.0 |  98.04% |     32748 |  1.8M | 76 |        26 | 12.24K | 438 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:24 |
| Q30L60X40P001 |   40.0 |  97.91% |     32657 | 1.79M | 73 |        26 | 11.58K | 412 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.03% |     32701 |  1.8M | 78 |        27 |  12.5K | 440 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:23 |
| Q30L60X80P000 |   80.0 |  97.97% |     32695 | 1.79M | 73 |        30 |  9.98K | 282 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:27 |
| Q30L60X80P001 |   80.0 |  97.88% |     32707 | 1.79M | 72 |        33 | 11.11K | 280 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:25 |
| Q30L60X80P002 |   80.0 |  97.94% |     32708 | 1.79M | 74 |        34 | 12.02K | 294 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:24 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.27% |     32357 | 1.79M | 73 |        78 |  9.36K | 146 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:21 |
| MRX40P001 |   40.0 |  96.27% |     32646 | 1.79M | 71 |        66 |  7.24K | 143 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:21 |
| MRX40P002 |   40.0 |  96.29% |     32662 | 1.79M | 72 |        72 |  8.31K | 144 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:20 |
| MRX80P000 |   80.0 |  96.23% |     32335 | 1.79M | 73 |        69 |  9.84K | 146 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:22 |
| MRX80P001 |   80.0 |  96.34% |     32683 | 1.79M | 73 |        66 |  8.29K | 145 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:23 |
| MRX80P002 |   80.0 |  96.33% |     32321 | 1.79M | 72 |        83 | 10.97K | 145 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:22 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.83% |     32676 |  1.8M | 82 |        36 | 13.61K | 361 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q0L0X40P001   |   40.0 |  97.80% |     32713 |  1.8M | 78 |       831 | 15.98K | 357 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:23 |
| Q0L0X40P002   |   40.0 |  97.72% |     32699 | 1.79M | 75 |      1024 | 15.47K | 325 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q0L0X80P000   |   80.0 |  97.43% |     32699 | 1.79M | 71 |      1712 |  9.39K | 169 |   81.0 | 7.0 |  20.0 | 204.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| Q0L0X80P001   |   80.0 |  97.22% |     32712 | 1.79M | 72 |      1069 |  8.12K | 168 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| Q0L0X80P002   |   80.0 |  97.44% |     32713 | 1.79M | 71 |      1712 | 10.04K | 187 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  97.92% |     32707 | 1.79M | 74 |        46 | 14.83K | 364 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:23 |
| Q25L60X40P001 |   40.0 |  97.80% |     32716 | 1.83M | 81 |       353 | 15.08K | 340 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| Q25L60X40P002 |   40.0 |  97.69% |     32700 |  1.8M | 76 |        32 |  11.6K | 335 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| Q25L60X80P000 |   80.0 |  97.17% |     32716 | 1.79M | 72 |      1712 |   9.3K | 165 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| Q25L60X80P001 |   80.0 |  97.46% |     32698 | 1.79M | 71 |      1712 |  9.33K | 177 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  97.16% |     32704 | 1.79M | 71 |      1712 |  8.97K | 155 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:23 |
| Q30L60X40P000 |   40.0 |  97.95% |     32741 |  1.8M | 78 |        45 | 14.75K | 366 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:24 |
| Q30L60X40P001 |   40.0 |  97.58% |     32693 |  1.8M | 77 |        42 | 12.68K | 327 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q30L60X40P002 |   40.0 |  97.79% |     32693 |  1.8M | 78 |        31 | 11.32K | 339 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q30L60X80P000 |   80.0 |  97.41% |     32759 | 1.83M | 73 |        39 |  7.41K | 178 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  97.51% |     32707 | 1.79M | 73 |       801 |  9.09K | 194 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  97.58% |     32708 | 1.79M | 73 |        33 |   7.7K | 196 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.28% |     32627 | 1.79M | 72 |        99 |    11K | 145 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| MRX40P001 |   40.0 |  96.27% |     32646 | 1.79M | 71 |        66 |  7.14K | 143 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:21 |
| MRX40P002 |   40.0 |  96.29% |     32667 | 1.79M | 72 |        72 |  8.34K | 144 |   40.0 | 5.0 |   8.3 | 110.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:21 |
| MRX80P000 |   80.0 |  96.28% |     32623 | 1.79M | 72 |        69 |  9.82K | 144 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:22 |
| MRX80P001 |   80.0 |  96.34% |     32683 | 1.79M | 73 |        66 |  8.31K | 145 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| MRX80P002 |   80.0 |  96.34% |     32649 | 1.79M | 71 |        84 | 10.97K | 143 |   80.0 | 7.0 |  19.7 | 202.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:22 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.83% |     35228 | 1.79M | 70 |     14504 | 140.28K | 59 |  278.0 | 16.0 |  76.7 | 652.0 |   0:00:37 |
| 7_merge_mr_unitigs_bcalm      |  98.89% |     35224 | 1.79M | 70 |      3224 |  11.44K |  4 |  277.0 | 17.0 |  75.3 | 656.0 |   0:00:52 |
| 7_merge_mr_unitigs_superreads |  98.66% |     35224 | 1.79M | 70 |     17567 |     29K |  5 |  277.0 | 17.0 |  75.3 | 656.0 |   0:00:45 |
| 7_merge_mr_unitigs_tadpole    |  98.89% |     35227 | 1.79M | 70 |     41567 | 147.29K |  9 |  277.0 | 17.0 |  75.3 | 656.0 |   0:00:52 |
| 7_merge_unitigs_bcalm         |  98.74% |     35221 | 1.79M | 71 |      1712 |  54.47K | 31 |  278.0 | 17.0 |  75.7 | 658.0 |   0:00:46 |
| 7_merge_unitigs_superreads    |  98.82% |     35224 | 1.79M | 70 |      1045 |  27.67K | 18 |  278.0 | 17.0 |  75.7 | 658.0 |   0:00:49 |
| 7_merge_unitigs_tadpole       |  98.92% |     35222 | 1.79M | 70 |      1036 |  70.37K | 48 |  277.0 | 17.0 |  75.3 | 656.0 |   0:00:55 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median |  MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|-----:|------:|-------:|----------:|
| 8_spades     |  97.52% |     38151 | 1.79M | 65 |      2158 | 10.13K | 131 |  278.0 | 16.0 |  76.7 |  652.0 |   0:00:27 |
| 8_mr_spades  |  97.16% |     38196 | 1.79M | 63 |      2176 | 13.85K | 127 |  462.0 | 21.0 | 133.0 | 1050.0 |   0:00:32 |
| 8_megahit    |  97.20% |     35198 | 1.79M | 69 |      1714 | 10.17K | 139 |  278.0 | 16.0 |  76.7 |  652.0 |   0:00:30 |
| 8_mr_megahit |  97.44% |     35527 | 1.79M | 66 |       159 | 24.75K | 135 |  462.0 | 22.0 | 132.0 | 1056.0 |   0:00:29 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1892775 | 1892775 |   1 |
| Paralogs                 |   33912 |   93528 |  10 |
| Repetitives              |    5357 |  139352 |  72 |
| 7_merge_anchors.anchors  |   35228 | 1791016 |  70 |
| 7_merge_anchors.others   |   14504 |  140280 |  59 |
| glue_anchors             |   35228 | 1790942 |  68 |
| fill_anchors             |   38197 | 1791749 |  65 |
| spades.contig            |   37811 | 1803171 |  80 |
| spades.scaffold          |   37811 | 1803225 |  78 |
| spades.non-contained     |   37811 | 1799490 |  66 |
| mr_spades.contig         |   37904 | 1812095 |  98 |
| mr_spades.scaffold       |   37904 | 1812095 |  98 |
| mr_spades.non-contained  |   37904 | 1804927 |  64 |
| megahit.contig           |   35250 | 1802645 |  77 |
| megahit.non-contained    |   35250 | 1798804 |  70 |
| mr_megahit.contig        |   35781 | 1823752 |  94 |
| mr_megahit.non-contained |   35781 | 1815182 |  69 |


# All strains

```shell script
cd ~/data/anchr/fda_argos

cat paralogs/cover.csv |
    grep -v "^#" |
    grep -v "^name" | head -n 50 | tail -n 20 |
    parallel -j 1 -k --colsep "," '
        WORKING_DIR=${HOME}/data/anchr/fda_argos
        BASE_NAME={1}

        cd ${WORKING_DIR}/${BASE_NAME}

        rm *.sh
        anchr template \
            --genome {2} \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe --cutoff 30 --cutk 31" \
            --sample 300 \
            --qual "25 30" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --quorum \
            --merge \
            --ecphase "1 2 3" \
            \
            --cov "40 80" \
            --unitigger "superreads bcalm tadpole" \
            --statp 2 \
            --readl 101 \
            --uscale 2 \
            --lscale 3 \
            --redo \
            \
            --extend

        bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
    '

```

* Cleanup

```shell script
cd ~/data/anchr/fda_argos

find . -name "4_down_sampling" -type d | xargs rm -fr
find . -name "6_down_sampling" -type d | xargs rm -fr

```

