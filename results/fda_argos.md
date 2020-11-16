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

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

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
| R.genome.bbtools  | 431.2 |    300 | 2174.5 |                         98.84% |
| R.tadpole.bbtools | 310.1 |    299 |  102.5 |                         96.44% |
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

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 1641481 | 1641481 |        1 |
| Paralogs   |    6079 |   33258 |       13 |
| Illumina.R |     101 |   1.55G | 15393600 |
| trim.R     |     100 | 475.33M |  4775680 |
| Q25L60     |     100 | 462.88M |  4660898 |
| Q30L60     |     100 | 439.09M |  4465020 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.54G | 15283902 |
| highpass | 101 |   1.54G | 15227304 |
| sample   | 101 | 492.44M |  4875686 |
| trim     | 100 | 475.43M |  4776690 |
| filter   | 100 | 475.33M |  4775680 |
| R1       | 100 | 237.94M |  2387840 |
| R2       | 100 | 237.39M |  2387840 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3657	0.07500%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	1010	0.02114%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	1010	0.02114%
```

```text
#R.peaks
#k	31
#unique_kmers	7526149
#error_kmers	5909838
#genomic_kmers	1616311
#main_peak	199
#genome_size_in_peaks	1640491
#genome_size	1641369
#haploid_genome_size	1641369
#fold_coverage	199
#haploid_fold_coverage	199
#ploidy	1
#percent_repeat_in_peaks	1.474
#percent_repeat	1.511
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 475.32M | 4775604 |
| ecco          | 100 | 475.32M | 4775604 |
| eccc          | 100 | 475.32M | 4775604 |
| ecct          | 100 | 472.69M | 4748814 |
| extended      | 140 | 662.26M | 4748814 |
| merged.raw    | 333 | 618.32M | 1932215 |
| unmerged.raw  | 140 | 122.72M |  884384 |
| unmerged.trim | 140 | 122.72M |  884376 |
| M1            | 333 | 617.79M | 1930601 |
| U1            | 140 |  61.57M |  442188 |
| U2            | 140 |  61.15M |  442188 |
| Us            |   0 |       0 |       0 |
| M.cor         | 317 | 742.44M | 4745578 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 159.6 |    164 |  20.7 |          4.66% |
| M.ihist.merge.txt  | 320.0 |    321 |  58.6 |         81.38% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 289.6 |  272.9 |    5.75% | "71" | 1.64M | 1.62M |     0.99 | 0:01'02'' |
| Q25L60.R | 282.0 |  268.8 |    4.68% | "71" | 1.64M | 1.62M |     0.99 | 0:01'01'' |
| Q30L60.R | 267.6 |  257.4 |    3.81% | "71" | 1.64M | 1.62M |     0.99 | 0:00'59'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.84% |     64778 | 1.59M | 50 |      2340 |  14.9K | 242 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.88% |     70825 |  1.6M | 55 |      2340 | 15.42K | 254 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  98.79% |     60743 | 1.59M | 53 |      2340 | 16.52K | 230 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  98.31% |     53067 | 1.59M | 60 |      2340 | 11.31K | 138 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q0L0X80P001   |   80.0 |  98.34% |     57877 | 1.59M | 56 |      2340 | 12.17K | 136 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q0L0X80P002   |   80.0 |  98.31% |     45895 | 1.59M | 66 |      2340 | 11.94K | 161 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.77% |     66691 | 1.59M | 46 |      1577 | 18.29K | 225 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q25L60X40P001 |   40.0 |  98.81% |     70549 |  1.6M | 50 |      2340 | 15.72K | 235 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  98.93% |     71564 |  1.6M | 51 |      2340 | 16.18K | 237 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q25L60X80P000 |   80.0 |  98.57% |     70561 | 1.59M | 54 |      2340 | 12.63K | 146 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q25L60X80P001 |   80.0 |  98.55% |     50166 | 1.59M | 54 |      6071 | 11.93K | 130 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  98.40% |     56892 | 1.59M | 56 |      2340 | 12.33K | 138 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.90% |     70831 | 1.59M | 49 |      2340 |  15.6K | 257 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q30L60X40P001 |   40.0 |  99.03% |     79315 |  1.6M | 45 |      1010 |  19.2K | 273 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  98.87% |     71569 | 1.59M | 45 |      2340 | 16.42K | 222 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:23 |
| Q30L60X80P000 |   80.0 |  98.52% |     71556 | 1.59M | 49 |      2340 | 12.82K | 136 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  98.58% |     75181 | 1.59M | 47 |      6071 | 11.93K | 134 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  98.45% |     70814 | 1.59M | 44 |      6071 | 11.52K | 126 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.93% |     70799 | 1.59M | 48 |      2340 | 14.76K | 113 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:22 |   0:00:22 |
| MRX40P001 |   40.0 |  97.85% |     71554 | 1.59M | 50 |      2340 | 14.78K | 113 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:21 |   0:00:21 |
| MRX40P002 |   40.0 |  97.83% |     71496 | 1.59M | 49 |      2340 | 14.93K | 114 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:21 |   0:00:23 |
| MRX80P000 |   80.0 |  97.62% |     48742 | 1.59M | 62 |      2340 |    16K | 136 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:23 |
| MRX80P001 |   80.0 |  97.72% |     50189 | 1.59M | 58 |      2340 | 15.63K | 132 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:24 |
| MRX80P002 |   80.0 |  97.84% |     56888 | 1.59M | 58 |      2340 |  15.2K | 130 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:23 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.95% |     71562 | 1.59M | 43 |      2340 | 16.75K | 238 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:26 |
| Q0L0X40P001   |   40.0 |  98.91% |     79899 | 1.59M | 43 |      1035 | 16.88K | 240 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  98.87% |     70828 | 1.59M | 50 |      1015 | 21.19K | 265 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  98.75% |     71575 |  1.6M | 46 |      2340 |  15.9K | 200 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  98.75% |     79886 | 1.59M | 42 |      2340 | 14.88K | 190 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  98.90% |     70811 |  1.6M | 40 |      2340 | 13.79K | 200 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  99.11% |     71567 | 1.59M | 46 |      1042 | 19.91K | 260 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  98.86% |     71568 | 1.59M | 43 |      2340 | 15.25K | 230 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| Q25L60X40P002 |   40.0 |  98.87% |     71541 | 1.59M | 48 |      2340 | 15.54K | 234 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  98.80% |     79945 |  1.6M | 41 |      2340 | 14.95K | 209 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:25 |
| Q25L60X80P001 |   80.0 |  98.76% |     71551 | 1.59M | 42 |      2340 | 14.22K | 212 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  98.89% |     79928 |  1.6M | 42 |      2340 | 15.98K | 211 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:26 |
| Q30L60X40P000 |   40.0 |  98.90% |     71566 | 1.59M | 48 |      2340 | 15.68K | 245 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:26 |
| Q30L60X40P001 |   40.0 |  98.76% |     70862 | 1.59M | 47 |      1036 | 19.95K | 235 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.84% |     70821 | 1.59M | 49 |      1051 | 17.01K | 262 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q30L60X80P000 |   80.0 |  99.10% |     71573 |  1.6M | 47 |      1012 | 18.42K | 255 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:27 |
| Q30L60X80P001 |   80.0 |  98.89% |     70822 | 1.59M | 44 |      1052 |  16.9K | 222 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:25 |
| Q30L60X80P002 |   80.0 |  98.89% |     71575 | 1.59M | 43 |      2340 | 15.73K | 220 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:26 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.14% |     81837 | 1.59M | 36 |      2340 | 13.66K | 94 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:22 |
| MRX40P001 |   40.0 |  98.36% |     80585 | 1.59M | 38 |      2340 | 15.08K | 97 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:22 |
| MRX40P002 |   40.0 |  98.14% |    104097 | 1.59M | 30 |      2340 | 13.91K | 78 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:22 |
| MRX80P000 |   80.0 |  97.90% |     80598 | 1.59M | 40 |      2340 | 13.23K | 86 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:23 |
| MRX80P001 |   80.0 |  97.82% |     79878 | 1.59M | 40 |      2340 | 12.92K | 88 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:22 |
| MRX80P002 |   80.0 |  97.89% |     80648 | 1.59M | 39 |      2340 | 12.52K | 80 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:22 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.94% |     70942 |  1.6M | 46 |      2340 | 14.89K | 199 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q0L0X40P001   |   40.0 |  99.00% |     79930 |  1.6M | 45 |      2340 | 15.03K | 213 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  98.75% |     79888 | 1.59M | 42 |      2340 | 13.54K | 188 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q0L0X80P000   |   80.0 |  98.83% |     75179 |  1.6M | 45 |      2340 | 12.29K | 140 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q0L0X80P001   |   80.0 |  98.75% |     80743 | 1.59M | 41 |      6071 | 12.12K | 125 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  98.78% |     75175 | 1.59M | 43 |      2340 | 12.46K | 136 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  98.90% |     60468 | 1.59M | 46 |      1018 | 22.29K | 201 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  98.87% |     79908 |  1.6M | 43 |      1047 | 17.15K | 214 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  98.87% |     56450 | 1.59M | 49 |      2340 | 15.96K | 207 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  98.67% |     79942 | 1.59M | 43 |      6071 | 12.04K | 124 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  98.80% |     81893 | 1.59M | 36 |      6071 | 12.06K | 128 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  98.64% |     80001 | 1.59M | 39 |      6071 |    12K | 126 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.81% |     79887 | 1.59M | 46 |      2340 | 13.66K | 191 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q30L60X40P001 |   40.0 |  98.88% |     71579 |  1.6M | 47 |      1067 | 17.85K | 216 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.85% |     51956 | 1.59M | 59 |      1048 | 20.79K | 218 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| Q30L60X80P000 |   80.0 |  98.83% |     80006 | 1.59M | 41 |      2340 | 13.18K | 146 |   80.0 | 4.0 |  22.7 | 138.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:25 |
| Q30L60X80P001 |   80.0 |  98.77% |     71559 | 1.59M | 44 |      2340 | 13.59K | 162 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:25 |
| Q30L60X80P002 |   80.0 |  98.70% |     79952 | 1.59M | 41 |      6071 | 11.93K | 134 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:25 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.88% |     79886 | 1.59M | 39 |      2340 | 13.12K | 86 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| MRX40P001 |   40.0 |  97.98% |     80579 | 1.59M | 39 |      2340 | 13.65K | 80 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| MRX40P002 |   40.0 |  97.85% |     79875 | 1.59M | 38 |      2340 | 12.91K | 84 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:21 |
| MRX80P000 |   80.0 |  97.84% |     79819 | 1.59M | 43 |      2340 | 13.24K | 89 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:33 |   0:00:23 |
| MRX80P001 |   80.0 |  97.79% |     75106 | 1.59M | 43 |      2340 | 12.99K | 88 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:24 |
| MRX80P002 |   80.0 |  97.87% |     79877 | 1.59M | 41 |      2340 | 12.67K | 86 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:22 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  98.08% |    112492 | 1.59M | 24 |      1325 | 92.49K | 54 |  274.0 | 17.5 |  73.8 | 489.8 |   0:00:33 |
| 7_merge_mr_unitigs_bcalm      |  98.59% |    104598 | 1.59M | 25 |     12393 | 21.93K |  4 |  270.0 | 14.0 |  76.0 | 468.0 |   0:00:40 |
| 7_merge_mr_unitigs_superreads |  98.79% |     90150 | 1.59M | 35 |      6071 |   9.1K |  3 |  273.0 | 12.0 |  79.0 | 463.5 |   0:00:44 |
| 7_merge_mr_unitigs_tadpole    |  98.80% |    104596 | 1.59M | 26 |      6071 | 10.19K |  4 |  272.0 | 13.0 |  77.7 | 466.5 |   0:00:43 |
| 7_merge_unitigs_bcalm         |  98.77% |    104620 |  1.6M | 28 |      1042 | 43.57K | 33 |  269.0 | 13.0 |  76.7 | 462.0 |   0:00:43 |
| 7_merge_unitigs_superreads    |  99.04% |    104143 |  1.6M | 29 |      6071 | 37.84K | 18 |  275.0 | 17.0 |  74.7 | 489.0 |   0:00:49 |
| 7_merge_unitigs_tadpole       |  99.04% |    112443 | 1.59M | 26 |      1067 | 45.75K | 31 |  272.0 | 14.0 |  76.7 | 471.0 |   0:00:49 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.35% |    104590 | 1.23M | 20 |      2420 | 16.14K | 36 |  272.0 | 18.0 |  72.7 | 489.0 |   0:00:27 |
| 8_mr_spades  |  99.48% |    153958 | 1.61M | 13 |      2340 | 13.99K | 29 |  449.0 | 36.0 | 113.7 | 835.5 |   0:00:31 |
| 8_megahit    |  98.68% |    112619 |  1.6M | 24 |      6014 | 10.79K | 47 |  272.0 | 15.0 |  75.7 | 475.5 |   0:00:29 |
| 8_mr_megahit |  99.31% |    153897 | 1.61M | 22 |      2340 | 12.76K | 44 |  449.0 | 32.0 | 117.7 | 817.5 |   0:00:30 |
| 8_platanus   |  98.90% |    153842 |  1.6M | 24 |      5981 | 11.05K | 45 |  272.0 | 15.0 |  75.7 | 475.5 |   0:00:29 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1641481 | 1641481 |   1 |
| Paralogs                 |    6079 |   33258 |  13 |
| 7_merge_anchors.anchors  |  112492 | 1594189 |  24 |
| 7_merge_anchors.others   |    1325 |   92485 |  54 |
| glue_anchors             |  112492 | 1594189 |  24 |
| fill_anchors             |  153854 | 1597199 |  19 |
| spades.contig            |  153957 | 1622961 |  34 |
| spades.scaffold          |  189386 | 1623099 |  32 |
| spades.non-contained     |  153957 | 1616693 |  17 |
| mr_spades.contig         |  189483 | 1624619 |  21 |
| mr_spades.scaffold       |  189483 | 1624619 |  21 |
| mr_spades.non-contained  |  189483 | 1623055 |  16 |
| megahit.contig           |  112661 | 1622748 |  59 |
| megahit.non-contained    |  112661 | 1607256 |  23 |
| mr_megahit.contig        |  174584 | 1632240 |  43 |
| mr_megahit.non-contained |  174584 | 1622196 |  22 |
| platanus.contig          |  112552 | 1629056 | 110 |
| platanus.scaffold        |  153893 | 1622434 |  64 |
| platanus.non-contained   |  153893 | 1612616 |  21 |


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
| R.genome.bbtools  | 688.7 |    320 | 3748.9 |                         97.02% |
| R.tadpole.bbtools | 328.6 |    317 |   99.9 |                         91.93% |
| R.genome.picard   | 330.2 |    319 |   87.9 |                             FR |
| R.tadpole.picard  | 328.6 |    317 |   88.0 |                             FR |


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

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 4274782 | 4274782 |        1 |
| Paralogs   |    3253 |  381668 |      135 |
| Illumina.R |     101 |   1.33G | 13190786 |
| trim.R     |     100 |   1.22G | 12274928 |
| Q25L60     |     100 |   1.19G | 11965201 |
| Q30L60     |     100 |   1.13G | 11482311 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.32G | 13028936 |
| highpass | 101 |   1.31G | 12944880 |
| sample   | 101 |   1.28G | 12697372 |
| trim     | 100 |   1.22G | 12275498 |
| filter   | 100 |   1.22G | 12274928 |
| R1       | 100 | 612.43M |  6137464 |
| R2       | 100 | 607.74M |  6137464 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6605	0.05202%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	286	0.00233%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	23093893
#error_kmers	18902958
#genomic_kmers	4190935
#main_peak	193
#genome_size_in_peaks	4316257
#genome_size	4331329
#haploid_genome_size	4331329
#fold_coverage	193
#haploid_fold_coverage	193
#ploidy	1
#percent_repeat_in_peaks	2.904
#percent_repeat	2.998
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |        # |
|:--------------|----:|--------:|---------:|
| clumped       | 100 |   1.22G | 12274348 |
| ecco          | 100 |   1.22G | 12274348 |
| eccc          | 100 |   1.22G | 12274348 |
| ecct          | 100 |   1.21G | 12211996 |
| extended      | 140 |    1.7G | 12211996 |
| merged.raw    | 343 |   1.52G |  4580337 |
| unmerged.raw  | 140 | 421.18M |  3051322 |
| unmerged.trim | 140 | 421.17M |  3051280 |
| M1            | 343 |   1.51G |  4558638 |
| U1            | 140 | 212.57M |  1525640 |
| U2            | 140 |  208.6M |  1525640 |
| Us            |   0 |       0 |        0 |
| M.cor         | 321 |   1.94G | 12168556 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 165.7 |    170 |  17.8 |          1.94% |
| M.ihist.merge.txt  | 331.1 |    333 |  54.3 |         75.01% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 285.4 |  266.1 |    6.76% | "71" | 4.27M | 4.22M |     0.99 | 0:02'29'' |
| Q25L60.R | 278.0 |  263.2 |    5.35% | "71" | 4.27M |  4.2M |     0.98 | 0:02'22'' |
| Q30L60.R | 264.5 |  253.7 |    4.09% | "71" | 4.27M |  4.2M |     0.98 | 0:02'18'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.32% |     64752 | 4.15M | 131 |      2384 | 45.92K | 695 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:41 |
| Q0L0X40P001   |   40.0 |  98.31% |     68249 | 4.15M | 131 |      7659 | 46.36K | 691 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:43 |   0:00:41 |
| Q0L0X40P002   |   40.0 |  98.32% |     61647 | 4.14M | 136 |      1367 | 51.82K | 685 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:40 |
| Q0L0X80P000   |   80.0 |  97.85% |     53308 | 4.14M | 145 |      5174 | 39.59K | 404 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  97.86% |     64975 | 4.14M | 132 |      5461 | 39.59K | 372 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:37 |
| Q0L0X80P002   |   80.0 |  98.10% |     85261 | 4.14M | 126 |      2867 | 39.71K | 384 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:40 |
| Q25L60X40P000 |   40.0 |  98.47% |     80141 | 4.15M | 127 |      7659 | 44.84K | 673 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:37 |   0:00:41 |
| Q25L60X40P001 |   40.0 |  98.39% |     63451 | 4.14M | 138 |      1486 | 51.99K | 631 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:40 |
| Q25L60X40P002 |   40.0 |  98.35% |     66523 | 3.91M | 133 |      1215 | 54.51K | 649 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:39 |
| Q25L60X80P000 |   80.0 |  98.07% |     73815 | 4.14M | 128 |      7659 | 39.02K | 365 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:40 |
| Q25L60X80P001 |   80.0 |  98.06% |     65830 | 4.14M | 132 |      7659 |  45.4K | 371 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:40 |
| Q25L60X80P002 |   80.0 |  97.98% |     65650 | 4.14M | 132 |      2793 | 40.74K | 381 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:38 |
| Q30L60X40P000 |   40.0 |  98.57% |     80158 | 4.13M | 131 |      1367 | 51.99K | 694 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:37 |   0:00:43 |
| Q30L60X40P001 |   40.0 |  98.48% |     75860 | 4.14M | 144 |      7659 | 49.43K | 647 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:42 |
| Q30L60X40P002 |   40.0 |  98.47% |     61648 | 4.09M | 135 |      7659 | 50.68K | 669 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:39 |
| Q30L60X80P000 |   80.0 |  98.39% |     80150 | 4.14M | 128 |      7659 | 41.73K | 401 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:42 |
| Q30L60X80P001 |   80.0 |  98.39% |     75726 | 4.14M | 118 |      5174 | 40.08K | 389 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:42 |
| Q30L60X80P002 |   80.0 |  98.36% |     80155 | 4.14M | 121 |      6801 | 38.46K | 356 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:42 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.74% |     88693 | 4.13M | 102 |      1215 | 39.36K | 278 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:36 |
| MRX40P001 |   40.0 |  97.45% |     80142 | 4.13M | 109 |      7659 | 50.76K | 294 |   39.0 | 3.5 |   9.5 |  74.2 | "31,41,51,61,71,81" |   0:00:51 |   0:00:35 |
| MRX40P002 |   40.0 |  97.47% |     85152 | 4.13M | 108 |      7659 | 50.49K | 293 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:37 |
| MRX80P000 |   80.0 |  97.23% |     69168 | 4.13M | 122 |      7659 | 49.27K | 320 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:23 |   0:00:37 |
| MRX80P001 |   80.0 |  97.06% |     70686 | 4.13M | 119 |      7659 | 50.15K | 313 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:01:23 |   0:00:34 |
| MRX80P002 |   80.0 |  96.43% |     80103 | 4.13M | 120 |      1094 | 43.51K | 320 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:23 |   0:00:36 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.20% |     60053 | 4.14M | 171 |      7921 |  95.44K | 903 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:40 |
| Q0L0X40P001   |   40.0 |  98.02% |     67549 | 4.13M | 159 |      7921 |  96.96K | 841 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:16 |   0:00:39 |
| Q0L0X40P002   |   40.0 |  98.12% |     59894 | 3.96M | 165 |      7911 |  107.1K | 906 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:41 |
| Q0L0X80P000   |   80.0 |  98.26% |    104043 | 4.14M | 106 |      7931 |  82.68K | 530 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:18 |   0:00:44 |
| Q0L0X80P001   |   80.0 |  98.23% |     85658 | 4.14M | 109 |      7931 |  80.51K | 510 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:01:20 |   0:00:43 |
| Q0L0X80P002   |   80.0 |  98.34% |     85241 | 4.14M | 115 |      7659 |  55.78K | 549 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:01:19 |   0:00:43 |
| Q25L60X40P000 |   40.0 |  98.08% |     66518 | 4.13M | 161 |      7921 |  98.48K | 889 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:08 |   0:00:40 |
| Q25L60X40P001 |   40.0 |  98.10% |     58658 | 3.99M | 172 |      7911 |  96.14K | 871 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:39 |
| Q25L60X40P002 |   40.0 |  98.22% |     58743 | 4.02M | 176 |      7911 | 102.72K | 896 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:41 |
| Q25L60X80P000 |   80.0 |  98.23% |    104044 | 4.14M | 113 |      7931 |  83.07K | 509 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:19 |   0:00:41 |
| Q25L60X80P001 |   80.0 |  98.37% |     85660 | 4.14M | 113 |      7921 |  84.36K | 548 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:45 |
| Q25L60X80P002 |   80.0 |  98.42% |     85649 | 4.15M | 113 |      7659 |  55.57K | 536 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:01:19 |   0:00:46 |
| Q30L60X40P000 |   40.0 |  98.04% |     59913 | 4.11M | 168 |      7921 |  99.13K | 900 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:09 |   0:00:40 |
| Q30L60X40P001 |   40.0 |  98.18% |     76088 | 4.11M | 165 |      7921 |  99.18K | 942 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:01:08 |   0:00:43 |
| Q30L60X40P002 |   40.0 |  98.13% |     61621 | 4.11M | 169 |      7911 | 106.27K | 911 |   40.0 | 3.5 |   9.8 |  75.8 | "31,41,51,61,71,81" |   0:01:09 |   0:00:43 |
| Q30L60X80P000 |   80.0 |  98.52% |     85645 | 4.11M | 115 |      7921 |  91.39K | 597 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:01:19 |   0:00:45 |
| Q30L60X80P001 |   80.0 |  98.51% |     85221 | 4.11M | 116 |      7659 |  59.11K | 587 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:01:19 |   0:00:44 |
| Q30L60X80P002 |   80.0 |  98.55% |     85227 | 4.15M | 122 |      7921 |  87.94K | 607 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:18 |   0:00:46 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.04% |    108137 | 4.13M |  91 |      7931 | 82.25K | 208 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:31 |
| MRX40P001 |   40.0 |  96.94% |    105449 | 4.13M |  89 |      7931 | 73.93K | 194 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:01:12 |   0:00:32 |
| MRX40P002 |   40.0 |  97.07% |     88644 | 4.13M | 102 |      7921 | 86.39K | 223 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:01:13 |   0:00:33 |
| MRX80P000 |   80.0 |  96.78% |    105439 | 4.13M |  92 |      7931 |  81.7K | 200 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:32 |
| MRX80P001 |   80.0 |  96.76% |    105432 | 4.13M |  90 |      7931 | 82.64K | 195 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:01:29 |   0:00:33 |
| MRX80P002 |   80.0 |  96.76% |    105400 | 4.13M |  90 |      7931 | 81.51K | 196 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:01:28 |   0:00:33 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.32% |     89747 | 4.14M | 119 |      7883 | 82.47K | 580 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:37 |   0:00:41 |
| Q0L0X40P001   |   40.0 |  98.27% |     88671 | 4.15M | 125 |      7883 | 82.22K | 575 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:40 |
| Q0L0X40P002   |   40.0 |  98.38% |     66541 | 3.98M | 137 |      7883 | 89.63K | 606 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:40 |
| Q0L0X80P000   |   80.0 |  98.25% |    105486 | 4.14M |  92 |      7883 | 76.62K | 321 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:42 |
| Q0L0X80P001   |   80.0 |  98.28% |    108210 | 4.14M |  94 |      7883 | 78.25K | 335 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:40 |
| Q0L0X80P002   |   80.0 |  98.34% |     96060 | 4.14M |  99 |      7883 | 77.04K | 337 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:43 |
| Q25L60X40P000 |   40.0 |  98.27% |     85244 | 4.14M | 119 |      7883 | 86.79K | 566 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:37 |   0:00:40 |
| Q25L60X40P001 |   40.0 |  98.47% |     80148 | 3.93M | 133 |      7883 | 91.65K | 605 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:42 |
| Q25L60X40P002 |   40.0 |  98.30% |     66523 | 3.87M | 130 |      7883 | 95.51K | 604 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:39 |
| Q25L60X80P000 |   80.0 |  98.15% |    105474 | 4.14M |  94 |      7883 | 77.37K | 320 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:43 |   0:00:37 |
| Q25L60X80P001 |   80.0 |  98.28% |    108209 | 4.14M |  97 |      7883 | 77.65K | 334 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:41 |
| Q25L60X80P002 |   80.0 |  98.33% |    108196 | 4.14M |  93 |      7883 | 79.05K | 325 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:41 |
| Q30L60X40P000 |   40.0 |  98.29% |     80135 | 4.11M | 131 |      7883 | 87.38K | 595 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:34 |   0:00:39 |
| Q30L60X40P001 |   40.0 |  98.51% |     80144 | 4.13M | 142 |      7883 | 94.57K | 687 |   39.0 | 3.5 |   9.5 |  74.2 | "31,41,51,61,71,81" |   0:00:35 |   0:00:43 |
| Q30L60X40P002 |   40.0 |  98.41% |     77294 | 4.12M | 134 |      7883 | 94.98K | 631 |   39.0 | 3.5 |   9.5 |  74.2 | "31,41,51,61,71,81" |   0:00:37 |   0:00:41 |
| Q30L60X80P000 |   80.0 |  98.53% |    105477 | 4.11M |  98 |      7883 | 85.64K | 395 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:44 |
| Q30L60X80P001 |   80.0 |  98.47% |     88778 | 4.11M |  97 |      7883 | 82.32K | 385 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:42 |
| Q30L60X80P002 |   80.0 |  98.52% |     88674 | 4.14M | 100 |      7883 | 81.19K | 393 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:44 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.77% |    108137 | 4.13M |  90 |      7883 | 81.05K | 194 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:30 |
| MRX40P001 |   40.0 |  96.80% |    105449 | 4.13M |  90 |      7883 |  81.6K | 195 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| MRX40P002 |   40.0 |  96.76% |     85597 | 4.13M | 101 |      7883 | 84.38K | 207 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:40 |   0:00:30 |
| MRX80P000 |   80.0 |  96.77% |    102931 | 4.13M |  93 |      7883 | 82.42K | 201 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:32 |
| MRX80P001 |   80.0 |  96.75% |    101351 | 4.13M |  98 |      7883 | 83.43K | 211 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:32 |
| MRX80P002 |   80.0 |  96.75% |    105400 | 4.13M |  92 |      7883 | 81.33K | 200 |   78.0 | 7.0 |  19.0 | 148.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:32 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  96.52% |    108120 |  4.2M | 87 |      7911 |   4.36M | 590 |  263.0 | 16.0 |  71.7 | 466.5 |   0:00:53 |
| 7_merge_mr_unitigs_bcalm      |  98.26% |    108157 | 4.13M | 84 |      7951 | 800.67K |  94 |  259.0 | 23.0 |  63.3 | 492.0 |   0:01:34 |
| 7_merge_mr_unitigs_superreads |  98.32% |    108154 | 4.13M | 86 |      7941 | 120.51K |  23 |  258.0 | 20.0 |  66.0 | 477.0 |   0:01:34 |
| 7_merge_mr_unitigs_tadpole    |  98.25% |    108153 | 4.13M | 84 |      7883 | 897.64K | 106 |  259.0 | 20.0 |  66.3 | 478.5 |   0:01:33 |
| 7_merge_unitigs_bcalm         |  98.04% |    105130 | 4.14M | 92 |      7941 |   1.28M | 207 |  257.0 | 21.0 |  64.7 | 480.0 |   0:01:24 |
| 7_merge_unitigs_superreads    |  98.40% |    108181 | 4.14M | 87 |      7961 | 220.47K |  51 |  257.0 | 21.0 |  64.7 | 480.0 |   0:01:38 |
| 7_merge_unitigs_tadpole       |  98.34% |    108171 | 4.14M | 91 |      7883 |   1.45M | 230 |  264.0 | 16.0 |  72.0 | 468.0 |   0:01:31 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.94% |    227280 | 4.15M | 46 |      7659 | 41.78K |  88 |  261.0 | 26.0 |  61.0 | 508.5 |   0:00:46 |
| 8_mr_spades  |  98.73% |    225840 | 4.17M | 35 |      7659 | 38.82K |  75 |  445.0 | 45.0 | 103.3 | 870.0 |   0:00:52 |
| 8_megahit    |  98.60% |    108259 | 4.14M | 81 |      7659 | 43.43K | 171 |  261.0 | 22.0 |  65.0 | 490.5 |   0:00:45 |
| 8_mr_megahit |  99.19% |    225227 | 4.18M | 45 |      7382 |  48.1K |  99 |  445.0 | 49.5 |  98.8 | 890.0 |   0:00:52 |
| 8_platanus   |  97.85% |    224880 | 3.91M | 47 |      7970 | 34.51K |  94 |  261.0 | 24.5 |  62.5 | 501.8 |   0:00:47 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 4274782 | 4274782 |   1 |
| Paralogs                 |    3253 |  381668 | 135 |
| 7_merge_anchors.anchors  |  108120 | 4201443 |  87 |
| 7_merge_anchors.others   |    7911 | 4359974 | 590 |
| glue_anchors             |  137639 | 4119823 |  78 |
| fill_anchors             |  211859 | 4121629 |  57 |
| spades.contig            |  227405 | 4222976 | 204 |
| spades.scaffold          |  227405 | 4223033 | 202 |
| spades.non-contained     |  227405 | 4194324 |  47 |
| mr_spades.contig         |  225542 | 4221145 |  82 |
| mr_spades.scaffold       |  225921 | 4221245 |  81 |
| mr_spades.non-contained  |  225542 | 4204714 |  40 |
| megahit.contig           |  108291 | 4209286 | 150 |
| megahit.non-contained    |  108291 | 4184244 |  90 |
| mr_megahit.contig        |  225563 | 4272423 | 177 |
| mr_megahit.non-contained |  225563 | 4224024 |  54 |
| platanus.contig          |  122109 | 4271880 | 622 |
| platanus.scaffold        |  225740 | 4239538 | 412 |
| platanus.non-contained   |  225740 | 4184556 |  48 |


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
| R.genome.bbtools  | 515.8 |    371 | 2333.8 |                         98.59% |
| R.tadpole.bbtools | 371.9 |    366 |  105.0 |                         86.75% |
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

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 2488635 | 2488635 |        1 |
| Paralogs   |    5627 |   56034 |       18 |
| Illumina.R |     101 |   1.12G | 11128812 |
| trim.R     |     100 | 655.73M |  6698550 |
| Q25L60     |     100 | 576.57M |  5976695 |
| Q30L60     |     100 | 445.55M |  4824739 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.12G | 11094542 |
| highpass | 101 |    1.1G | 10890248 |
| sample   | 101 | 746.59M |  7391984 |
| trim     | 100 | 655.73M |  6698550 |
| filter   | 100 | 655.73M |  6698550 |
| R1       | 100 | 328.89M |  3349275 |
| R2       | 100 | 326.84M |  3349275 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6021	0.08145%
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
#unique_kmers	18226849
#error_kmers	15779812
#genomic_kmers	2447037
#main_peak	177
#genome_size_in_peaks	2482735
#genome_size	2490711
#haploid_genome_size	2490711
#fold_coverage	177
#haploid_fold_coverage	177
#ploidy	1
#percent_repeat_in_peaks	1.438
#percent_repeat	1.521
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 655.72M | 6698442 |
| ecco          | 100 | 655.71M | 6698442 |
| eccc          | 100 | 655.71M | 6698442 |
| ecct          | 100 | 653.04M | 6670428 |
| extended      | 140 | 918.52M | 6670428 |
| merged.raw    | 365 | 647.93M | 1885564 |
| unmerged.raw  | 140 | 395.45M | 2899300 |
| unmerged.trim | 140 | 395.43M | 2899186 |
| M1            | 365 | 646.88M | 1882547 |
| U1            | 140 | 198.68M | 1449593 |
| U2            | 140 | 196.76M | 1449593 |
| Us            |   0 |       0 |       0 |
| M.cor         | 310 |   1.04G | 6664280 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 154.3 |    158 |  23.4 |          2.81% |
| M.ihist.merge.txt  | 343.6 |    355 |  60.4 |         56.53% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 263.5 |  234.8 |   10.91% | "51" | 2.49M |  2.7M |     1.08 | 0:01'22'' |
| Q25L60.R | 231.9 |  213.6 |    7.87% | "51" | 2.49M | 2.47M |     0.99 | 0:01'16'' |
| Q30L60.R | 179.3 |  169.1 |    5.70% | "49" | 2.49M | 2.45M |     0.98 | 0:01'03'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.31% |     18236 | 2.47M | 237 |        40 | 45.32K | 1113 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:28 |
| Q0L0X40P001   |   40.0 |  98.41% |     19377 | 2.46M | 209 |        46 |  43.1K |  993 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  98.25% |     17998 | 2.47M | 216 |        38 | 38.23K |  989 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  97.09% |     13315 | 2.44M | 276 |        29 | 20.14K |  668 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  96.97% |     14739 | 2.44M | 252 |        34 | 21.93K |  612 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  99.17% |     31389 | 2.44M | 141 |       399 | 36.99K |  681 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  99.02% |     30365 | 2.44M | 138 |       152 | 33.15K |  637 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:29 |
| Q25L60X40P002 |   40.0 |  99.07% |     31706 | 2.39M | 123 |       271 | 32.13K |  606 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:28 |
| Q25L60X80P000 |   80.0 |  98.40% |     32338 | 2.44M | 118 |        36 | 14.61K |  406 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q25L60X80P001 |   80.0 |  98.75% |     35696 | 2.44M | 101 |        54 | 17.08K |  382 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  99.34% |     31697 | 2.43M | 133 |       898 | 43.71K |  638 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:32 |
| Q30L60X40P001 |   40.0 |  99.28% |     27338 | 2.43M | 142 |       908 | 39.45K |  604 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:31 |
| Q30L60X40P002 |   40.0 |  99.35% |     28763 |  2.4M | 158 |       780 | 59.06K |  658 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:32 |
| Q30L60X80P000 |   80.0 |  99.35% |     46317 | 2.43M |  95 |       365 | 27.95K |  464 |   79.0 | 9.0 |  17.3 | 158.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:36 |
| Q30L60X80P001 |   80.0 |  99.40% |     44459 | 2.43M | 112 |       334 | 30.44K |  466 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:35 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.68% |     62856 | 2.43M |  73 |        97 | 17.14K | 184 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:26 |
| MRX40P001 |   40.0 |  98.62% |     57499 | 2.43M |  86 |        93 | 16.33K | 170 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:26 |
| MRX40P002 |   40.0 |  98.52% |     61541 | 2.43M |  71 |        95 | 15.26K | 179 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:27 |
| MRX80P000 |   80.0 |  98.26% |     44385 | 2.43M | 102 |        89 | 18.79K | 235 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:26 |
| MRX80P001 |   80.0 |  98.26% |     39904 | 2.43M | 100 |        83 | 18.33K | 235 |   80.0 | 5.0 |  21.7 | 142.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| MRX80P002 |   80.0 |  98.20% |     37440 | 2.43M | 104 |        91 | 19.96K | 241 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:26 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.87% |     41909 | 2.44M | 106 |        52 | 30.78K | 670 |   39.0 |  3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:32 |
| Q0L0X40P001   |   40.0 |  98.67% |     46277 | 2.44M | 101 |        36 | 22.59K | 602 |   39.0 |  3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  98.76% |     50689 | 2.43M |  93 |       298 | 32.32K | 588 |   40.0 |  3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:53 |   0:00:32 |
| Q0L0X80P000   |   80.0 |  98.79% |     61703 | 2.45M |  71 |        42 | 20.39K | 465 |   79.0 |  6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:58 |   0:00:32 |
| Q0L0X80P001   |   80.0 |  98.75% |     69788 | 2.44M |  67 |        46 | 18.66K | 396 |   79.0 |  6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:58 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  99.09% |     43224 | 2.43M | 110 |        65 | 31.23K | 676 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  99.05% |     48075 | 2.44M | 108 |        34 | 22.42K | 635 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:31 |
| Q25L60X40P002 |   40.0 |  99.11% |     44575 | 2.43M | 101 |      1026 | 34.95K | 668 |   40.0 |  4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:30 |
| Q25L60X80P000 |   80.0 |  99.35% |     61764 | 2.44M |  82 |        38 | 21.07K | 521 |   79.0 |  6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:01:00 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  99.32% |     52452 | 2.44M |  78 |        57 | 26.12K | 508 |   79.0 |  6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:59 |   0:00:34 |
| Q30L60X40P000 |   40.0 |  98.72% |     19460 | 2.42M | 226 |       773 | 42.37K | 845 |   41.0 |  5.0 |   8.7 |  82.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:28 |
| Q30L60X40P001 |   40.0 |  98.50% |     18601 | 2.41M | 234 |       824 | 46.55K | 818 |   40.0 |  5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  98.66% |     16128 | 2.42M | 252 |       720 | 44.17K | 880 |   41.0 |  6.0 |   7.7 |  82.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  99.27% |     45133 | 2.44M | 106 |       935 | 33.03K | 638 |   80.0 | 10.0 |  16.7 | 160.0 | "31,41,51,61,71,81" |   0:01:00 |   0:00:34 |
| Q30L60X80P001 |   80.0 |  99.28% |     46294 | 2.35M | 102 |       781 | 38.44K | 628 |   80.0 |  9.0 |  17.7 | 160.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:33 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.59% |     97841 | 2.44M | 53 |        90 | 11.27K | 122 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:54 |   0:00:26 |
| MRX40P001 |   40.0 |  98.55% |     47865 | 1.96M | 69 |       149 | 12.39K | 113 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:00:54 |   0:00:26 |
| MRX40P002 |   40.0 |  98.86% |     61784 | 2.17M | 67 |      1035 | 14.53K | 114 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:00:55 |   0:00:28 |
| MRX80P000 |   80.0 |  98.50% |    102435 | 2.44M | 54 |       107 | 11.81K | 117 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:01:07 |   0:00:26 |
| MRX80P001 |   80.0 |  98.45% |    113833 | 2.44M | 50 |       206 | 11.37K | 106 |   79.0 | 4.5 |  21.8 | 138.8 | "31,41,51,61,71,81" |   0:01:03 |   0:00:26 |
| MRX80P002 |   80.0 |  98.37% |    107896 | 2.43M | 48 |      1339 | 11.74K |  96 |   79.0 | 4.5 |  21.8 | 138.8 | "31,41,51,61,71,81" |   0:01:04 |   0:00:26 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.82% |     47362 | 2.44M |  91 |        51 | 23.92K | 522 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:29 |
| Q0L0X40P001   |   40.0 |  98.82% |     49127 | 2.44M |  96 |        53 | 24.89K | 520 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:29 |
| Q0L0X40P002   |   40.0 |  98.82% |     55305 | 2.44M |  83 |        71 | 23.92K | 479 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:30 |
| Q0L0X80P000   |   80.0 |  98.66% |     61653 | 2.44M |  68 |        47 | 15.53K | 328 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:30 |
| Q0L0X80P001   |   80.0 |  98.80% |     65492 | 2.44M |  65 |        59 | 18.26K | 325 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  99.31% |     45824 | 2.36M | 103 |       402 | 34.22K | 541 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:30 |
| Q25L60X40P001 |   40.0 |  99.27% |     40083 | 2.31M | 112 |        59 | 27.45K | 567 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:30 |
| Q25L60X40P002 |   40.0 |  99.22% |     57292 | 2.43M |  75 |      1000 | 25.46K | 463 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:32 |
| Q25L60X80P000 |   80.0 |  99.31% |     61712 | 2.44M |  69 |        48 | 17.16K | 381 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q25L60X80P001 |   80.0 |  99.25% |     57383 | 2.44M |  69 |        70 | 19.95K | 370 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:31 |
| Q30L60X40P000 |   40.0 |  99.26% |     30987 | 2.43M | 136 |       903 |  44.2K | 632 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:29 |
| Q30L60X40P001 |   40.0 |  99.28% |     30238 | 2.38M | 144 |       824 | 38.28K | 639 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:30 |
| Q30L60X40P002 |   40.0 |  99.29% |     40039 | 2.39M | 129 |       784 | 45.71K | 648 |   40.0 | 5.0 |   8.3 |  80.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:33 |
| Q30L60X80P000 |   80.0 |  99.36% |     46318 | 2.44M |  99 |       573 |  29.2K | 522 |   78.0 | 9.0 |  17.0 | 156.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q30L60X80P001 |   80.0 |  99.32% |     44807 | 2.32M | 105 |       611 | 34.06K | 505 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:33 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.53% |     57539 |  1.9M | 63 |       171 | 13.12K | 120 |   39.0 | 2.0 |  11.0 |  67.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:24 |
| MRX40P001 |   40.0 |  98.51% |     80632 | 2.44M | 56 |       216 | 11.05K | 111 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:24 |
| MRX40P002 |   40.0 |  98.42% |     75987 | 2.44M | 50 |       197 | 10.18K |  98 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:26 |
| MRX80P000 |   80.0 |  98.42% |     64236 | 2.44M | 60 |       114 | 11.79K | 119 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:25 |
| MRX80P001 |   80.0 |  98.38% |     61771 | 2.43M | 60 |       102 | 12.01K | 119 |   79.0 | 4.0 |  22.3 | 136.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:26 |
| MRX80P002 |   80.0 |  98.36% |     70048 | 2.43M | 55 |      1339 | 12.62K | 114 |   79.0 | 5.0 |  21.3 | 141.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:26 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.62% |    107902 | 2.44M | 47 |     59625 | 369.84K | 81 |  233.0 | 13.0 |  64.7 | 408.0 |   0:00:35 |
| 7_merge_mr_unitigs_bcalm      |  98.90% |    115863 | 2.44M | 45 |     97854 | 117.32K | 13 |  232.0 | 14.0 |  63.3 | 411.0 |   0:00:55 |
| 7_merge_mr_unitigs_superreads |  98.96% |    107932 | 2.44M | 47 |     59625 |  99.77K | 14 |  233.0 | 13.5 |  64.2 | 410.2 |   0:00:59 |
| 7_merge_mr_unitigs_tadpole    |  98.77% |    115861 | 2.44M | 46 |      1667 |  12.07K |  9 |  233.0 | 15.0 |  62.7 | 417.0 |   0:00:51 |
| 7_merge_unitigs_bcalm         |  98.63% |    107951 | 2.44M | 58 |      1038 |  53.51K | 46 |  233.0 | 15.0 |  62.7 | 417.0 |   0:00:55 |
| 7_merge_unitigs_superreads    |  99.01% |    115897 | 2.44M | 47 |      1109 |  53.44K | 40 |  232.0 | 14.0 |  63.3 | 411.0 |   0:01:00 |
| 7_merge_unitigs_tadpole       |  99.00% |    115897 | 2.44M | 47 |     24601 | 127.81K | 36 |  232.0 | 15.0 |  62.3 | 415.5 |   0:01:00 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  98.87% |    125041 |  2.1M | 28 |      1433 |  12.7K | 48 |  232.0 | 13.0 |  64.3 | 406.5 |   0:00:35 |
| 8_mr_spades  |  99.43% |    156585 | 2.44M | 26 |      1433 |  11.5K | 47 |  418.0 | 21.0 | 118.3 | 721.5 |   0:00:37 |
| 8_megahit    |  98.76% |    115906 | 2.44M | 47 |      1433 | 12.97K | 92 |  232.0 | 14.0 |  63.3 | 411.0 |   0:00:31 |
| 8_mr_megahit |  99.63% |    168224 | 2.45M | 25 |      1433 | 13.73K | 46 |  418.0 | 21.0 | 118.3 | 721.5 |   0:00:36 |
| 8_platanus   |  98.26% |    130922 | 2.44M | 31 |      2331 |  8.38K | 54 |  232.0 | 10.0 |  67.3 | 393.0 |   0:00:31 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 2488635 | 2488635 |   1 |
| Paralogs                 |    5627 |   56034 |  18 |
| 7_merge_anchors.anchors  |  107902 | 2436009 |  47 |
| 7_merge_anchors.others   |   59625 |  369844 |  81 |
| glue_anchors             |  107902 | 2434960 |  45 |
| fill_anchors             |  165171 | 2437855 |  24 |
| spades.contig            |  310298 | 2456501 |  29 |
| spades.scaffold          |  310298 | 2456701 |  27 |
| spades.non-contained     |  310298 | 2453467 |  21 |
| mr_spades.contig         |  179012 | 2456591 |  24 |
| mr_spades.scaffold       |  179012 | 2456691 |  23 |
| mr_spades.non-contained  |  179012 | 2454591 |  21 |
| megahit.contig           |  115949 | 2456223 |  57 |
| megahit.non-contained    |  115949 | 2451363 |  45 |
| mr_megahit.contig        |  309887 | 2467470 |  39 |
| mr_megahit.non-contained |  309887 | 2460195 |  21 |
| platanus.contig          |   97966 | 2470000 | 175 |
| platanus.scaffold        |  177051 | 2464880 | 124 |
| platanus.non-contained   |  177051 | 2446544 |  23 |

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
| R.genome.bbtools  | 547.1 |    365 | 2596.0 |                         97.41% |
| R.tadpole.bbtools | 374.0 |    363 |  122.4 |                         95.26% |
| R.genome.picard   | 375.5 |    365 |  108.9 |                             FR |
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

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 1892775 | 1892775 |        1 |
| Paralogs   |   33912 |   93528 |       10 |
| Illumina.R |     101 |   2.14G | 21230270 |
| trim.R     |     100 | 549.68M |  5517578 |
| Q25L60     |     100 | 538.62M |  5415291 |
| Q30L60     |     100 | 516.06M |  5232491 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   2.12G | 21017828 |
| highpass | 101 |   2.12G | 20974022 |
| sample   | 101 | 567.83M |  5622104 |
| trim     | 100 | 549.75M |  5518288 |
| filter   | 100 | 549.68M |  5517578 |
| R1       | 100 | 275.11M |  2758789 |
| R2       | 100 | 274.58M |  2758789 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3585	0.06377%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	710	0.01287%
#Name	Reads	ReadsPct
NC_001422.1 Coliphage phiX174, complete genome	710	0.01287%
```

```text
#R.peaks
#k	31
#unique_kmers	7555499
#error_kmers	5762260
#genomic_kmers	1793239
#main_peak	203
#genome_size_in_peaks	1869798
#genome_size	1883343
#haploid_genome_size	1883343
#fold_coverage	203
#haploid_fold_coverage	203
#ploidy	1
#percent_repeat_in_peaks	4.184
#percent_repeat	4.315
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 549.68M | 5517546 |
| ecco          | 100 | 549.68M | 5517546 |
| eccc          | 100 | 549.68M | 5517546 |
| ecct          | 100 | 547.32M | 5493570 |
| extended      | 140 |  766.6M | 5493570 |
| merged.raw    | 361 | 552.55M | 1617716 |
| unmerged.raw  | 140 | 314.47M | 2258138 |
| unmerged.trim | 140 | 314.46M | 2258134 |
| M1            | 361 | 551.83M | 1615635 |
| U1            | 140 | 157.52M | 1129067 |
| U2            | 140 | 156.95M | 1129067 |
| Us            |   0 |       0 |       0 |
| M.cor         | 311 | 867.91M | 5489404 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 154.6 |    159 |  23.5 |          2.73% |
| M.ihist.merge.txt  | 341.6 |    350 |  59.3 |         58.90% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG | EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|-----:|---------:|----------:|
| Q0L0.R   | 290.4 |  275.8 |    5.05% | "71" | 1.89M | 1.8M |     0.95 | 0:01'06'' |
| Q25L60.R | 284.6 |  272.8 |    4.15% | "71" | 1.89M | 1.8M |     0.95 | 0:01'06'' |
| Q30L60.R | 272.7 |  263.6 |    3.35% | "71" | 1.89M | 1.8M |     0.95 | 0:01'03'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.15% |     27884 | 1.76M | 85 |      6334 | 60.68K | 414 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:26 |
| Q0L0X40P001   |   40.0 |  98.06% |     31596 | 1.76M | 88 |     12504 | 64.43K | 421 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  97.93% |     31601 | 1.76M | 92 |     15111 | 70.17K | 425 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  97.18% |     32151 | 1.76M | 76 |     12504 | 45.34K | 175 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q0L0X80P001   |   80.0 |  97.23% |     32349 | 1.76M | 74 |     26730 | 45.08K | 177 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  97.07% |     32678 | 1.76M | 74 |      8299 | 44.62K | 170 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q25L60X40P000 |   40.0 |  98.20% |     27244 | 1.76M | 93 |     26730 |  81.8K | 414 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  97.93% |     28853 | 1.76M | 91 |     26730 | 84.53K | 405 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:23 |
| Q25L60X40P002 |   40.0 |  98.01% |     31604 | 1.76M | 87 |     26730 | 77.81K | 397 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  97.17% |     32338 | 1.76M | 75 |     26730 | 45.17K | 176 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q25L60X80P001 |   80.0 |  97.12% |     32152 | 1.76M | 79 |      9103 | 48.16K | 176 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q25L60X80P002 |   80.0 |  97.09% |     32696 | 1.76M | 70 |      6334 | 44.46K | 163 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q30L60X40P000 |   40.0 |  98.10% |     24671 | 1.76M | 96 |     26730 | 82.17K | 453 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:25 |
| Q30L60X40P001 |   40.0 |  98.08% |     27545 | 1.76M | 88 |     26730 | 78.92K | 432 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.10% |     26457 | 1.76M | 95 |     26730 | 78.18K | 418 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q30L60X80P000 |   80.0 |  97.23% |     32758 | 1.76M | 73 |     26730 | 66.31K | 172 |   81.0 | 1.0 |  26.0 | 126.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| Q30L60X80P001 |   80.0 |  97.32% |     32692 | 1.76M | 72 |      6334 | 45.16K | 171 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q30L60X80P002 |   80.0 |  97.24% |     32697 | 1.76M | 74 |     26730 | 66.01K | 175 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:21 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.66% |     32251 | 1.75M | 80 |     26730 | 73.39K | 174 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| MRX40P001 |   40.0 |  96.71% |     32598 | 1.75M | 81 |     26730 | 72.57K | 173 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| MRX40P002 |   40.0 |  96.72% |     32562 | 1.75M | 79 |     26730 | 72.18K | 170 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:22 |   0:00:21 |
| MRX80P000 |   80.0 |  96.69% |     32679 | 1.75M | 70 |     26730 | 71.83K | 162 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:23 |
| MRX80P001 |   80.0 |  96.67% |     32680 | 1.75M | 69 |     26730 | 71.34K | 160 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:22 |
| MRX80P002 |   80.0 |  96.66% |     32644 | 1.75M | 73 |     26730 | 71.93K | 168 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:36 |   0:00:22 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.12% |     27884 | 1.75M | 84 |      2719 | 66.51K | 477 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  97.96% |     31555 | 1.76M | 93 |     23818 | 92.39K | 506 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:24 |
| Q0L0X40P002   |   40.0 |  98.02% |     31369 | 1.75M | 94 |      5358 | 62.64K | 473 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  97.73% |     32749 | 1.76M | 76 |     27554 | 48.47K | 356 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  97.55% |     32743 | 1.76M | 71 |     27554 | 45.79K | 340 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:24 |
| Q0L0X80P002   |   80.0 |  97.38% |     35011 | 1.76M | 73 |     27554 | 46.43K | 333 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.14% |     28213 | 1.75M | 93 |      7246 | 65.16K | 503 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  97.83% |     28853 | 1.75M | 91 |      3421 | 66.03K | 455 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  97.96% |     30214 | 1.76M | 88 |      5358 |  59.3K | 478 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:41 |   0:00:25 |
| Q25L60X80P000 |   80.0 |  97.77% |     32756 | 1.76M | 72 |     27554 | 48.65K | 370 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  97.55% |     32748 | 1.76M | 72 |     27554 | 47.07K | 358 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  97.56% |     32743 | 1.76M | 72 |     27554 | 46.19K | 352 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:25 |
| Q30L60X40P000 |   40.0 |  98.05% |     24642 | 1.75M | 98 |      5358 | 64.74K | 498 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q30L60X40P001 |   40.0 |  98.09% |     27538 | 1.76M | 91 |      8029 | 70.39K | 505 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q30L60X40P002 |   40.0 |  97.98% |     23887 | 1.76M | 98 |      5358 | 58.98K | 495 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| Q30L60X80P000 |   80.0 |  97.79% |     32747 | 1.76M | 74 |     27554 | 49.65K | 376 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  97.68% |     32750 | 1.76M | 74 |     27554 | 49.79K | 387 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  97.91% |     32756 | 1.76M | 72 |     27554 | 47.83K | 387 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.32% |     32639 | 1.75M | 78 |     27554 | 45.02K | 149 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:20 |
| MRX40P001 |   40.0 |  96.37% |     32598 | 1.75M | 81 |     27554 | 44.58K | 153 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:43 |   0:00:21 |
| MRX40P002 |   40.0 |  96.34% |     32562 | 1.75M | 78 |     27554 | 44.24K | 149 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:20 |
| MRX80P000 |   80.0 |  96.35% |     32679 | 1.75M | 69 |     27554 | 43.66K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:20 |
| MRX80P001 |   80.0 |  96.34% |     32680 | 1.75M | 69 |     27554 | 43.52K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:21 |
| MRX80P002 |   80.0 |  96.35% |     32691 | 1.75M | 69 |     27554 | 43.36K | 141 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:21 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.86% |     29477 | 1.76M | 84 |     27554 | 52.81K | 357 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q0L0X40P001   |   40.0 |  97.89% |     31564 | 1.76M | 92 |     27554 | 52.04K | 392 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X40P002   |   40.0 |  97.86% |     30656 | 1.76M | 96 |     27091 | 83.59K | 384 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  97.38% |     32761 | 1.76M | 73 |     27554 | 40.82K | 177 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| Q0L0X80P001   |   80.0 |  97.41% |     32349 | 1.76M | 71 |     27554 | 39.91K | 173 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  97.16% |     32678 | 1.76M | 70 |     27554 | 39.37K | 162 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  97.85% |     25229 | 1.67M | 95 |      5358 | 57.22K | 394 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| Q25L60X40P001 |   40.0 |  97.49% |     30147 | 1.76M | 93 |      5358 | 56.39K | 321 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| Q25L60X40P002 |   40.0 |  97.85% |     31576 | 1.76M | 87 |     12974 | 85.51K | 370 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  97.46% |     32757 | 1.76M | 69 |     27554 | 40.15K | 181 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:23 |
| Q25L60X80P001 |   80.0 |  97.10% |     32703 | 1.76M | 71 |     27554 | 39.46K | 162 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:22 |
| Q25L60X80P002 |   80.0 |  97.46% |     32684 | 1.76M | 72 |     27554 |  40.1K | 182 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  97.97% |     24671 | 1.75M | 97 |     20503 | 81.54K | 398 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| Q30L60X40P001 |   40.0 |  97.85% |     27535 | 1.76M | 91 |     27554 | 53.18K | 342 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  97.89% |     25162 | 1.76M | 97 |     27554 | 49.65K | 385 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X80P000 |   80.0 |  97.50% |     32693 | 1.76M | 75 |     27554 | 40.57K | 197 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  97.52% |     32768 | 1.76M | 69 |     27554 |  40.1K | 187 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  97.75% |     32766 | 1.76M | 72 |     27554 | 40.89K | 210 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:24 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.32% |     32251 | 1.75M | 79 |     27554 | 45.12K | 151 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:21 |
| MRX40P001 |   40.0 |  96.36% |     32598 | 1.75M | 81 |     27554 | 44.47K | 153 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:21 |
| MRX40P002 |   40.0 |  96.33% |     32562 | 1.75M | 78 |     27554 | 44.18K | 149 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:20 |
| MRX80P000 |   80.0 |  96.30% |     32679 | 1.75M | 69 |     27554 | 43.71K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:21 |
| MRX80P001 |   80.0 |  96.34% |     32680 | 1.75M | 69 |     27554 | 43.52K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:21 |
| MRX80P002 |   80.0 |  96.30% |     32691 | 1.75M | 69 |     27554 |  43.4K | 141 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:20 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  94.01% |     35101 | 1.75M | 68 |     21288 | 505.04K | 138 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:28 |
| 7_merge_mr_unitigs_bcalm      |  97.70% |     35083 | 1.75M | 68 |     27554 | 160.72K |  18 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:47 |
| 7_merge_mr_unitigs_superreads |  97.64% |     32674 | 1.75M | 69 |     26730 | 105.97K |  16 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:45 |
| 7_merge_mr_unitigs_tadpole    |  97.72% |     35083 | 1.75M | 68 |     27554 | 139.44K |  17 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:48 |
| 7_merge_unitigs_bcalm         |  96.91% |     35148 | 1.76M | 69 |      1144 | 143.52K |  94 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:38 |
| 7_merge_unitigs_superreads    |  97.80% |     35161 | 1.76M | 68 |     16174 | 223.38K |  58 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:49 |
| 7_merge_unitigs_tadpole       |  97.68% |     35150 | 1.76M | 68 |     26471 |  250.6K |  78 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:46 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|----------:|
| 8_spades     |  97.52% |     36658 | 1.76M | 67 |     27548 | 40.34K | 133 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:27 |
| 8_mr_spades  |  97.25% |     37825 | 1.76M | 64 |     27661 | 43.51K | 128 |  462.0 | 7.0 | 147.0 | 724.5 |   0:00:28 |
| 8_megahit    |  97.20% |     36661 | 1.76M | 67 |     27556 | 39.45K | 137 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:27 |
| 8_mr_megahit |  97.41% |     36753 | 1.76M | 65 |     27844 | 53.66K | 134 |  462.0 | 8.0 | 146.0 | 729.0 |   0:00:28 |
| 8_platanus   |  96.96% |     36629 | 1.76M | 66 |     27530 | 39.38K | 131 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:26 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1892775 | 1892775 |   1 |
| Paralogs                 |   33912 |   93528 |  10 |
| 7_merge_anchors.anchors  |   35101 | 1753124 |  68 |
| 7_merge_anchors.others   |   21288 |  505035 | 138 |
| glue_anchors             |   35101 | 1753124 |  68 |
| fill_anchors             |   38074 | 1754726 |  63 |
| spades.contig            |   37811 | 1803326 |  82 |
| spades.scaffold          |   37811 | 1803380 |  80 |
| spades.non-contained     |   37811 | 1799490 |  66 |
| mr_spades.contig         |   37904 | 1811595 |  93 |
| mr_spades.scaffold       |   37904 | 1811595 |  93 |
| mr_spades.non-contained  |   37904 | 1805657 |  64 |
| megahit.contig           |   35250 | 1802645 |  77 |
| megahit.non-contained    |   35250 | 1798804 |  70 |
| mr_megahit.contig        |   35781 | 1821737 |  86 |
| mr_megahit.non-contained |   35781 | 1814792 |  69 |
| platanus.contig          |   35268 | 1808500 | 121 |
| platanus.scaffold        |   37808 | 1805556 |  99 |
| platanus.non-contained   |   37808 | 1798706 |  65 |


# All strains

```shell script
cd ~/data/anchr/fda_argos

cat paralogs/cover.csv |
    grep -v "^#" |
    grep -v "^name" | head -n 10 | tail -n 5 |
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

