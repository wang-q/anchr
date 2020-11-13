# Assemble genomes from FDA-ARGOS data sets

[TOC levels=1-3]: # ""

- [Assemble genomes from FDA-ARGOS data sets](#assemble-genomes-from-fda-argos-data-sets)
- [Download](#download)
  - [Good assemblies before 2010](#good-assemblies-before-2010)
  - [Reference Assemblies](#reference-assemblies)
  - [Illumina Reads](#illumina-reads)
  - [Reference Genomes](#reference-genomes)
  - [Paralogs](#paralogs)
  - [Copy/link files](#copylink-files)
- [Borr_bur_B31](#borr_bur_b31)
- [Bord_pert_Tohama_I](#bord_pert_tohama_i)
- [Ca_jej_jejuni_NCTC_11168_ATCC_700819](#ca_jej_jejuni_nctc_11168_atcc_700819)
- [Clostridio_dif_630](#clostridio_dif_630)
- [Co_dip_NCTC_13129](#co_dip_nctc_13129)
- [Fr_tul_tularensis_SCHU_S4](#fr_tul_tularensis_schu_s4)
- [Ha_inf_Rd_KW20](#ha_inf_rd_kw20)
- [Leg_pneumop_pneumophila_Philadelphia_1](#leg_pneumop_pneumophila_philadelphia_1)
- [N_gon_FA_1090](#n_gon_fa_1090)


* Rsync to hpcc

```shell script
for D in Ftul Hinf Cjej Lpne Cdip Cdif; do
    rsync -avP \
        ~/data/anchr/${D}/ \
        wangq@202.119.37.251:data/anchr/${D}
done

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
SAMN16357376,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
#SAMN16357375,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
#SAMN06173357,Par_dis_ATCC_8503,Parabacteroides distasonis ATCC 8503
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


# Borr_bur_B31


```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Borr_bur_B31

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1467551 \
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
| R.genome.bbtools  | 406.5 |    340 | 1473.4 |                         97.28% |
| R.tadpole.bbtools | 343.4 |    336 |  101.1 |                         81.89% |
| R.genome.picard   | 347.6 |    340 |   99.0 |                             FR |
| R.tadpole.picard  | 343.4 |    336 |   98.9 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 192       | 77069           | 0.0000       | 32.25   |
| R.31 | 185       | 103962          | 0.0000       | 31.29   |
| R.41 | 249       | 0               | 0.0000       | 30.80   |
| R.51 | 250       | 0               | 0.0000       | 30.47   |
| R.61 | 411       | 981665          | 0.1085       | 30.23   |
| R.71 | 303       | 1103645         | 0.0904       | 30.05   |
| R.81 | 213       | 1223735         | 0.0905       | 29.87   |


Table: statReads

| Name       |    N50 |     Sum |        # |
|:-----------|-------:|--------:|---------:|
| Genome     | 910724 | 1521208 |       22 |
| Paralogs   |   5170 |  472496 |      127 |
| Illumina.R |    101 |   1.47G | 14568214 |
| trim.R     |    100 | 423.21M |  4258016 |
| Q25L60     |    100 | 404.87M |  4088324 |
| Q30L60     |    100 | 372.18M |  3805951 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.46G | 14452890 |
| highpass | 101 |   1.45G | 14377628 |
| sample   | 101 | 440.27M |  4359062 |
| trim     | 100 | 423.21M |  4258016 |
| filter   | 100 | 423.21M |  4258016 |
| R1       | 100 | 211.58M |  2129008 |
| R2       | 100 | 211.63M |  2129008 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3878	0.08896%
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
#unique_kmers	9031771
#error_kmers	7801506
#genomic_kmers	1230265
#main_peak	222
#genome_size_in_peaks	3832377
#genome_size	3755242
#haploid_genome_size	1251747
#fold_coverage	77
#haploid_fold_coverage	222
#ploidy	3
#het_rate	0.00099
#percent_repeat_in_peaks	8.071
#percent_repeat	67.038
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 423.21M | 4257974 |
| ecco          | 100 | 423.21M | 4257974 |
| eccc          | 100 | 423.21M | 4257974 |
| ecct          | 100 | 419.38M | 4219142 |
| extended      | 140 | 585.41M | 4219142 |
| merged.raw    | 355 | 438.21M | 1306280 |
| unmerged.raw  | 140 | 221.36M | 1606582 |
| unmerged.trim | 140 | 221.36M | 1606524 |
| M1            | 355 | 437.53M | 1304298 |
| U1            | 140 | 110.78M |  803262 |
| U2            | 140 | 110.57M |  803262 |
| Us            |   0 |       0 |       0 |
| M.cor         | 312 | 660.19M | 4215120 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 152.9 |    157 |  23.0 |          3.47% |
| M.ihist.merge.txt  | 335.5 |    344 |  60.8 |         61.92% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 288.4 |  263.9 |    8.50% | "71" | 1.47M | 1.27M |     0.87 | 0:01'00'' |
| Q25L60.R | 276.0 |  256.5 |    7.04% | "71" | 1.47M | 1.26M |     0.86 | 0:00'53'' |
| Q30L60.R | 253.7 |  238.6 |    5.99% | "71" | 1.47M | 1.25M |     0.85 | 0:00'50'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  89.01% |     12329 | 1.09M | 152 |      1438 | 53.66K | 419 |   44.0 |  5.0 |   9.7 |  88.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:20 |
| Q0L0X40P001   |   40.0 |  89.08% |     11897 | 1.09M | 159 |      1175 | 56.25K | 443 |   44.0 |  4.0 |  10.7 |  84.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:20 |
| Q0L0X40P002   |   40.0 |  88.69% |     10594 |  1.1M | 170 |      1155 |  57.2K | 467 |   44.0 |  5.0 |   9.7 |  88.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:20 |
| Q0L0X80P000   |   80.0 |  86.31% |      8831 | 1.09M | 193 |      1310 | 44.07K | 412 |   89.0 |  9.0 |  20.7 | 174.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:21 |
| Q0L0X80P001   |   80.0 |  86.01% |      8442 | 1.08M | 191 |      1488 | 52.97K | 411 |   89.0 |  8.0 |  21.7 | 169.5 | "31,41,51,61,71,81" |   0:00:20 |   0:00:21 |
| Q0L0X80P002   |   80.0 |  85.40% |      8036 | 1.08M | 199 |      1423 | 43.46K | 423 |   89.0 |  8.0 |  21.7 | 169.5 | "31,41,51,61,71,81" |   0:00:20 |   0:00:20 |
| Q25L60X40P000 |   40.0 |  89.11% |     13123 | 1.09M | 154 |      1575 | 61.34K | 430 |   44.0 |  4.0 |  10.7 |  84.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:20 |
| Q25L60X40P001 |   40.0 |  88.37% |     10838 |  1.1M | 160 |      1139 | 52.75K | 422 |   44.0 |  5.0 |   9.7 |  88.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:19 |
| Q25L60X40P002 |   40.0 |  88.93% |     13110 |  1.1M | 147 |      1102 | 49.44K | 405 |   44.0 |  5.0 |   9.7 |  88.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:19 |
| Q25L60X80P000 |   80.0 |  86.16% |      9040 | 1.09M | 181 |      1514 | 44.41K | 390 |   89.0 |  9.0 |  20.7 | 174.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:20 |
| Q25L60X80P001 |   80.0 |  85.98% |      9587 | 1.09M | 179 |      1514 | 40.11K | 381 |   89.0 |  9.0 |  20.7 | 174.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:20 |
| Q25L60X80P002 |   80.0 |  86.05% |      8493 | 1.09M | 199 |      1393 | 50.25K | 426 |   89.0 |  8.0 |  21.7 | 169.5 | "31,41,51,61,71,81" |   0:00:20 |   0:00:19 |
| Q30L60X40P000 |   40.0 |  88.85% |     13114 | 1.09M | 152 |      1481 | 56.95K | 424 |   44.0 |  4.5 |  10.2 |  86.2 | "31,41,51,61,71,81" |   0:00:13 |   0:00:21 |
| Q30L60X40P001 |   40.0 |  89.18% |     13107 |  1.1M | 147 |      1248 | 51.74K | 425 |   44.0 |  5.0 |   9.7 |  88.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:20 |
| Q30L60X40P002 |   40.0 |  89.32% |     15163 |  1.1M | 148 |      1430 | 45.44K | 404 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:13 |   0:00:20 |
| Q30L60X80P000 |   80.0 |  86.15% |     10412 | 1.09M | 173 |      1364 |  41.6K | 373 |   89.0 |  9.0 |  20.7 | 174.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:21 |
| Q30L60X80P001 |   80.0 |  86.59% |     10158 |  1.1M | 175 |      1236 | 38.75K | 378 |   89.0 | 10.0 |  19.7 | 178.0 | "31,41,51,61,71,81" |   0:00:20 |   0:00:20 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  84.40% |     13608 | 1.09M | 143 |      1334 | 43.76K | 302 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:19 |
| MRX40P001 |   40.0 |  84.18% |     13138 | 1.09M | 144 |      1310 | 40.04K | 305 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:19 |
| MRX40P002 |   40.0 |  84.12% |     11834 |  1.1M | 155 |      1479 | 37.52K | 326 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:19 |
| MRX80P000 |   80.0 |  83.60% |     10703 | 1.09M | 161 |      1293 | 36.93K | 337 |   90.0 | 11.0 |  19.0 | 180.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:20 |
| MRX80P001 |   80.0 |  83.61% |     10521 | 1.09M | 170 |      1148 | 38.19K | 361 |   90.0 | 11.0 |  19.0 | 180.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:19 |
| MRX80P002 |   80.0 |  83.72% |     10149 | 1.09M | 167 |      1257 | 41.62K | 351 |   90.0 | 10.0 |  20.0 | 180.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:19 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  89.72% |     18827 | 1.09M | 127 |      1310 | 46.84K | 407 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:21 |
| Q0L0X40P001   |   40.0 |  89.46% |     18409 | 1.09M | 126 |      1529 | 46.86K | 411 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:20 |
| Q0L0X40P002   |   40.0 |  90.14% |     15591 |  1.1M | 129 |      1543 |  42.9K | 404 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:20 |
| Q0L0X80P000   |   80.0 |  89.19% |     16078 | 1.11M | 131 |      1326 | 38.13K | 348 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:21 |
| Q0L0X80P001   |   80.0 |  89.01% |     14761 | 1.12M | 144 |      1314 | 44.83K | 393 |   89.0 |  9.0 |  20.7 | 174.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:21 |
| Q0L0X80P002   |   80.0 |  89.08% |     17024 | 1.11M | 131 |      1326 | 40.05K | 368 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:21 |
| Q25L60X40P000 |   40.0 |  89.35% |     16620 | 1.09M | 121 |      1505 | 58.28K | 374 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:20 |
| Q25L60X40P001 |   40.0 |  89.03% |     16653 | 1.09M | 129 |      1330 | 42.33K | 420 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:20 |
| Q25L60X40P002 |   40.0 |  89.62% |     18170 |  1.1M | 129 |      1310 | 44.33K | 405 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:20 |
| Q25L60X80P000 |   80.0 |  88.81% |     15061 | 1.11M | 134 |      1310 | 41.45K | 353 |   89.0 | 11.0 |  18.7 | 178.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:20 |
| Q25L60X80P001 |   80.0 |  89.46% |     15592 | 1.11M | 134 |      1163 | 45.11K | 388 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:21 |
| Q25L60X80P002 |   80.0 |  88.71% |     15589 |  1.1M | 127 |      1326 | 40.46K | 342 |   89.0 | 10.0 |  19.7 | 178.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:21 |
| Q30L60X40P000 |   40.0 |  89.28% |     15335 | 1.09M | 133 |      1514 |  43.1K | 405 |   45.0 |  5.5 |   9.5 |  90.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:20 |
| Q30L60X40P001 |   40.0 |  89.37% |     16669 | 1.09M | 123 |      1326 | 47.78K | 433 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:20 |
| Q30L60X40P002 |   40.0 |  89.90% |     15595 |  1.1M | 140 |      1326 | 47.37K | 427 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:20 |
| Q30L60X80P000 |   80.0 |  89.19% |     15608 | 1.11M | 131 |      1310 | 39.69K | 392 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:21 |
| Q30L60X80P001 |   80.0 |  88.83% |     15465 |  1.1M | 127 |      1163 | 38.65K | 373 |   89.0 | 10.0 |  19.7 | 178.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:22 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  85.66% |     15622 | 1.12M | 122 |      1326 | 24.61K | 266 |   45.0 |  9.0 |   6.0 |  90.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:20 |
| MRX40P001 |   40.0 |  85.58% |     16625 | 1.11M | 118 |      1310 | 28.36K | 266 |   45.0 |  8.0 |   7.0 |  90.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:20 |
| MRX40P002 |   40.0 |  85.51% |     18453 | 1.11M | 111 |      1310 | 31.94K | 249 |   45.0 |  6.5 |   8.5 |  90.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:19 |
| MRX80P000 |   80.0 |  84.52% |     14786 | 1.11M | 141 |      1310 | 34.35K | 298 |   90.0 | 12.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:20 |
| MRX80P001 |   80.0 |  84.78% |     15173 | 1.11M | 140 |      1334 | 35.55K | 297 |   90.0 | 12.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:19 |
| MRX80P002 |   80.0 |  84.92% |     15173 | 1.11M | 139 |      1229 | 35.85K | 292 |   90.0 | 12.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:19 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  89.09% |     15594 |  1.1M | 130 |      1326 | 51.95K | 377 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| Q0L0X40P001   |   40.0 |  89.22% |     15601 | 1.11M | 133 |      1163 | 47.71K | 398 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:20 |
| Q0L0X40P002   |   40.0 |  89.02% |     15425 | 1.11M | 129 |      1163 | 44.52K | 374 |   45.0 |  7.0 |   8.0 |  90.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:20 |
| Q0L0X80P000   |   80.0 |  87.97% |     15341 | 1.11M | 133 |      1326 | 36.63K | 293 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:21 |
| Q0L0X80P001   |   80.0 |  87.47% |     15068 | 1.11M | 143 |      1163 | 32.69K | 308 |   89.0 | 11.0 |  18.7 | 178.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:20 |
| Q0L0X80P002   |   80.0 |  87.56% |     14764 | 1.11M | 136 |      1514 |  40.5K | 301 |   89.0 | 11.5 |  18.2 | 178.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:21 |
| Q25L60X40P000 |   40.0 |  89.07% |     16058 |  1.1M | 128 |      1410 | 41.26K | 350 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:20 |
| Q25L60X40P001 |   40.0 |  89.36% |     15655 | 1.11M | 134 |      1092 |  39.2K | 394 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:20 |
| Q25L60X40P002 |   40.0 |  88.67% |     16622 | 1.11M | 134 |      1310 | 39.55K | 360 |   45.0 |  7.0 |   8.0 |  90.0 | "31,41,51,61,71,81" |   0:00:22 |   0:00:19 |
| Q25L60X80P000 |   80.0 |  87.48% |     15065 | 1.11M | 137 |      1438 | 37.56K | 303 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:21 |
| Q25L60X80P001 |   80.0 |  87.66% |     15338 | 1.11M | 137 |      1514 | 30.85K | 291 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:21 |
| Q25L60X80P002 |   80.0 |  87.44% |     15081 | 1.11M | 130 |      1468 | 31.95K | 288 |   89.0 | 12.5 |  17.2 | 178.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:20 |
| Q30L60X40P000 |   40.0 |  90.14% |     15056 |  1.1M | 138 |      1472 | 52.07K | 394 |   45.0 |  5.0 |  10.0 |  90.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:20 |
| Q30L60X40P001 |   40.0 |  89.21% |     16632 |  1.1M | 122 |      1514 | 43.81K | 377 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:19 |
| Q30L60X40P002 |   40.0 |  89.25% |     15641 |  1.1M | 138 |      1147 | 45.18K | 379 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:20 |
| Q30L60X80P000 |   80.0 |  87.88% |     15353 | 1.11M | 137 |      1352 |  32.9K | 308 |   89.0 | 13.0 |  16.7 | 178.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:20 |
| Q30L60X80P001 |   80.0 |  88.24% |     15160 | 1.11M | 132 |      1462 | 34.42K | 312 |   89.0 | 12.0 |  17.7 | 178.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  85.39% |     15384 | 1.11M | 131 |      1611 | 39.28K | 280 |   45.0 |  7.0 |   8.0 |  90.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:20 |
| MRX40P001 |   40.0 |  85.19% |     15160 | 1.12M | 128 |      1272 | 28.69K | 274 |   45.0 |  7.0 |   8.0 |  90.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:20 |
| MRX40P002 |   40.0 |  84.96% |     15225 | 1.11M | 130 |      1585 | 36.63K | 275 |   45.0 |  6.0 |   9.0 |  90.0 | "31,41,51,61,71,81" |   0:00:24 |   0:00:19 |
| MRX80P000 |   80.0 |  84.25% |     13649 |  1.1M | 146 |      1310 | 36.05K | 304 |   90.0 | 11.5 |  18.5 | 180.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:21 |
| MRX80P001 |   80.0 |  84.59% |     13135 | 1.11M | 151 |      1157 | 33.25K | 316 |   90.0 | 12.0 |  18.0 | 180.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:20 |
| MRX80P002 |   80.0 |  84.71% |     14781 |  1.1M | 145 |      1310 | 39.82K | 304 |   90.0 | 10.0 |  20.0 | 180.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:20 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  91.74% |     44341 | 1.13M |  71 |      1639 |  177.4K | 104 |  295.0 | 42.0 |  56.3 | 590.0 |   0:00:32 |
| 7_merge_mr_unitigs_bcalm      |  94.00% |     40409 | 1.13M |  77 |      1376 |   54.5K |  36 |  294.0 | 58.0 |  40.0 | 588.0 |   0:00:37 |
| 7_merge_mr_unitigs_superreads |  93.41% |     15585 | 1.11M | 122 |      2298 |  71.18K |  39 |  294.0 | 29.5 |  68.5 | 573.8 |   0:00:37 |
| 7_merge_mr_unitigs_tadpole    |  93.61% |     27426 | 1.12M | 100 |      1601 |  62.35K |  38 |  294.0 | 41.5 |  56.5 | 588.0 |   0:00:37 |
| 7_merge_unitigs_bcalm         |  94.79% |     24677 | 1.13M | 105 |      1522 | 153.13K |  89 |  296.0 | 31.0 |  67.7 | 583.5 |   0:00:40 |
| 7_merge_unitigs_superreads    |  94.73% |     21807 | 1.12M | 113 |      1958 | 145.62K |  80 |  295.0 | 35.0 |  63.3 | 590.0 |   0:00:42 |
| 7_merge_unitigs_tadpole       |  94.96% |     24506 | 1.13M | 100 |      1486 | 125.68K |  81 |  296.0 | 34.0 |  64.7 | 592.0 |   0:00:41 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |     Sum |  # | N50Others |    Sum |  # | median |   MAD | lower |  upper | RunTimeAN |
|:-------------|--------:|----------:|--------:|---:|----------:|-------:|---:|-------:|------:|------:|-------:|----------:|
| 8_spades     |  93.50% |    353721 | 623.24K | 45 |      2805 | 38.82K | 89 |  296.0 |  95.5 |   3.2 |  592.0 |   0:00:26 |
| 8_mr_spades  |  95.70% |    428465 |   1.21M | 47 |      3591 | 40.62K | 53 |  508.0 | 194.0 |   3.0 | 1016.0 |   0:00:29 |
| 8_megahit    |  94.98% |    856638 |   1.63M | 50 |      4821 | 50.01K | 84 |  298.0 | 105.5 |   3.0 |  596.0 |   0:00:30 |
| 8_mr_megahit |  96.99% |    428371 |   1.24M | 68 |      3051 |  56.4K | 93 |  507.0 | 175.0 |   3.0 | 1014.0 |   0:00:30 |
| 8_platanus   |  89.88% |    353567 |   1.15M | 42 |      1796 | 19.32K | 56 |  295.0 |  96.0 |   3.0 |  590.0 |   0:00:27 |


Table: statFinal

| Name                     |    N50 |     Sum |    # |
|:-------------------------|-------:|--------:|-----:|
| Genome                   | 910724 | 1521208 |   22 |
| Paralogs                 |   5170 |  472496 |  127 |
| 7_merge_anchors.anchors  |  44341 | 1128009 |   71 |
| 7_merge_anchors.others   |   1639 |  177399 |  104 |
| glue_anchors             |  73385 | 1125013 |   64 |
| fill_anchors             | 212315 | 1129797 |   43 |
| spades.contig            | 353804 | 1279499 |  208 |
| spades.scaffold          | 911053 | 1299893 |  206 |
| spades.non-contained     | 353804 | 1235396 |   54 |
| mr_spades.contig         | 910617 | 1252948 |   62 |
| mr_spades.scaffold       | 910617 | 1253148 |   60 |
| mr_spades.non-contained  | 910617 | 1245918 |   50 |
| megahit.contig           | 938112 | 1799189 |  244 |
| megahit.non-contained    | 938112 | 1723011 |   60 |
| mr_megahit.contig        | 906303 | 1437246 |  441 |
| mr_megahit.non-contained | 906303 | 1294724 |   78 |
| platanus.contig          | 213136 | 1357041 | 1071 |
| platanus.scaffold        | 428230 | 1321199 |  885 |
| platanus.non-contained   | 428230 | 1173913 |   43 |


# Bord_pert_Tohama_I


```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Bord_pert_Tohama_I

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 4086189 \
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


# Ha_inf_Rd_KW20

```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Ha_inf_Rd_KW20

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 1830138 \
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
| R.genome.bbtools  | 393.9 |    266 | 2164.6 |                         96.99% |
| R.tadpole.bbtools | 273.5 |    265 |   78.5 |                         94.86% |
| R.genome.picard   | 274.8 |    266 |   77.4 |                             FR |
| R.tadpole.picard  | 273.5 |    265 |   77.2 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 0         | 0               | 0.0000       | 40.91   |
| R.31 | 450       | 1355576         | 0.0000       | 40.46   |
| R.41 | 386       | 1750118         | 0.0000       | 40.21   |
| R.51 | 325       | 1778118         | 0.0000       | 40.02   |
| R.61 | 257       | 1778905         | 0.0000       | 39.87   |
| R.71 | 194       | 1796512         | 0.0000       | 39.73   |
| R.81 | 130       | 1805531         | 0.0000       | 39.59   |


Table: statReads

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 1830138 | 1830138 |        1 |
| Paralogs   |    5432 |   95244 |       29 |
| Illumina.R |     101 |   1.24G | 12231248 |
| trim.R     |     100 | 525.04M |  5286948 |
| Q25L60     |     100 | 507.59M |  5131277 |
| Q30L60     |     100 | 472.62M |  4856022 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.23G | 12143650 |
| highpass | 101 |   1.22G | 12073026 |
| sample   | 101 | 549.04M |  5436052 |
| trim     | 100 | 525.09M |  5287460 |
| filter   | 100 | 525.04M |  5286948 |
| R1       | 100 | 262.87M |  2643474 |
| R2       | 100 | 262.17M |  2643474 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3959	0.07283%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	512	0.00968%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	7923739
#error_kmers	6135262
#genomic_kmers	1788477
#main_peak	201
#genome_size_in_peaks	1828856
#genome_size	1829868
#haploid_genome_size	1829868
#fold_coverage	201
#haploid_fold_coverage	201
#ploidy	1
#percent_repeat_in_peaks	2.212
#percent_repeat	2.067
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 525.03M | 5286814 |
| ecco          | 100 | 525.02M | 5286814 |
| eccc          | 100 | 525.02M | 5286814 |
| ecct          | 100 | 523.18M | 5267922 |
| extended      | 140 | 733.33M | 5267922 |
| merged.raw    | 311 | 714.21M | 2386584 |
| unmerged.raw  | 140 |  68.16M |  494754 |
| unmerged.trim | 140 |  68.16M |  494742 |
| M1            | 311 | 713.36M | 2383764 |
| U1            | 140 |  34.27M |  247371 |
| U2            | 140 |  33.88M |  247371 |
| Us            |   0 |       0 |       0 |
| M.cor         | 303 |  783.9M | 5262270 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 158.4 |    163 |  21.8 |         10.25% |
| M.ihist.merge.txt  | 299.3 |    297 |  61.2 |         90.61% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 286.9 |  271.6 |    5.33% | "71" | 1.83M |  1.8M |     0.98 | 0:01'07'' |
| Q25L60.R | 277.4 |  266.0 |    4.11% | "71" | 1.83M | 1.79M |     0.98 | 0:01'05'' |
| Q30L60.R | 258.4 |  250.3 |    3.12% | "71" | 1.83M | 1.79M |     0.98 | 0:01'02'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.99% |     48258 | 1.77M | 69 |      1032 | 18.69K | 359 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.91% |     55138 | 1.77M | 61 |      1027 |  19.1K | 328 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  99.01% |     52079 | 1.77M | 64 |      1002 | 23.24K | 351 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  98.53% |     36585 | 1.77M | 81 |      1561 | 15.62K | 188 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  98.57% |     39119 | 1.77M | 70 |      1298 | 15.51K | 182 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X80P002   |   80.0 |  98.54% |     33145 | 1.77M | 83 |      1298 | 16.27K | 209 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  99.23% |     58188 | 1.77M | 57 |      1022 | 22.63K | 389 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  99.07% |     57120 | 1.77M | 57 |      1012 | 22.16K | 326 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q25L60X40P002 |   40.0 |  99.06% |     54012 | 1.77M | 59 |      1022 | 20.41K | 319 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  98.83% |     57114 | 1.77M | 58 |      2561 | 15.17K | 153 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  98.85% |     46443 | 1.77M | 64 |      2561 | 15.45K | 158 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  98.91% |     54631 | 1.77M | 63 |      1634 | 15.92K | 176 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  99.28% |     45124 | 1.77M | 64 |       444 | 23.72K | 394 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.22% |     51674 | 1.77M | 61 |      1009 | 24.26K | 368 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:27 |
| Q30L60X40P002 |   40.0 |  99.26% |     55133 | 1.77M | 64 |      1001 | 22.31K | 371 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  99.06% |     55139 | 1.77M | 52 |      3503 | 15.18K | 161 |   78.0 | 8.0 |  18.0 | 153.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q30L60X80P001 |   80.0 |  99.07% |     58187 | 1.77M | 52 |      3503 | 15.44K | 172 |   78.0 | 8.0 |  18.0 | 153.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q30L60X80P002 |   80.0 |  99.15% |     55128 | 1.77M | 56 |      3503 | 16.98K | 195 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.08% |     52505 | 1.77M | 59 |       970 | 16.48K | 135 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| MRX40P001 |   40.0 |  98.09% |     54455 | 1.76M | 56 |      1164 | 19.19K | 127 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| MRX40P002 |   40.0 |  98.00% |     54519 | 1.76M | 59 |      1164 | 18.11K | 133 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| MRX80P000 |   80.0 |  97.97% |     41797 | 1.77M | 64 |      1022 | 17.04K | 143 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:23 |
| MRX80P001 |   80.0 |  97.84% |     42433 | 1.76M | 70 |      1022 | 19.13K | 156 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:23 |
| MRX80P002 |   80.0 |  97.81% |     44335 | 1.76M | 65 |      1164 | 20.14K | 145 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:23 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.93% |     53976 | 1.76M | 57 |      1065 | 22.55K | 436 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.97% |     51072 | 1.77M | 62 |      1022 |  27.1K | 447 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:28 |
| Q0L0X40P002   |   40.0 |  98.86% |     51042 | 1.76M | 62 |      1011 |    27K | 400 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  98.97% |     58150 | 1.77M | 53 |      1022 | 17.51K | 302 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:26 |
| Q0L0X80P001   |   80.0 |  99.08% |     59992 | 1.77M | 55 |      1022 | 19.58K | 303 |   80.0 | 6.5 |  20.2 | 149.2 | "31,41,51,61,71,81" |   0:00:52 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  99.12% |     60011 | 1.77M | 54 |      1023 | 22.29K | 309 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:27 |
| Q25L60X40P000 |   40.0 |  99.03% |     44085 | 1.77M | 63 |      1016 | 25.76K | 469 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  98.98% |     49409 | 1.77M | 64 |      1009 | 26.18K | 477 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  98.98% |     47795 | 1.76M | 61 |      1025 | 26.25K | 410 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  99.03% |     60002 | 1.77M | 53 |      1022 | 20.13K | 304 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q25L60X80P001 |   80.0 |  99.11% |     55137 | 1.77M | 52 |      1634 | 21.05K | 292 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q25L60X80P002 |   80.0 |  99.13% |     54583 | 1.77M | 54 |      1634 | 20.65K | 329 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:27 |
| Q30L60X40P000 |   40.0 |  98.97% |     43965 | 1.76M | 71 |      1021 | 27.07K | 523 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.04% |     49400 | 1.76M | 67 |        55 | 23.33K | 510 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:48 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  98.94% |     34714 | 1.77M | 79 |      1010 | 29.21K | 518 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  99.15% |     55136 | 1.77M | 53 |      1036 |  21.3K | 353 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q30L60X80P001 |   80.0 |  99.10% |     58143 | 1.77M | 53 |      1027 | 21.33K | 359 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:29 |
| Q30L60X80P002 |   80.0 |  99.15% |     55113 | 1.77M | 53 |      1019 | 24.21K | 370 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.98% |     58198 | 1.77M | 50 |      1634 | 14.68K | 102 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:22 |
| MRX40P001 |   40.0 |  98.22% |     59967 | 1.76M | 50 |      3503 | 17.45K | 103 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| MRX40P002 |   40.0 |  98.02% |     59974 | 1.76M | 49 |      1606 | 16.35K | 100 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:22 |
| MRX80P000 |   80.0 |  97.95% |     58097 | 1.77M | 50 |      1103 | 14.75K | 103 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:23 |
| MRX80P001 |   80.0 |  97.94% |     58116 | 1.77M | 50 |      1022 | 14.01K | 103 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:22 |
| MRX80P002 |   80.0 |  97.90% |     58053 | 1.77M | 49 |      1634 | 16.01K | 102 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:25 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.10% |     55092 | 1.77M | 57 |      1634 | 21.13K | 328 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  99.04% |     53972 | 1.77M | 56 |      1047 | 19.59K | 314 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  99.13% |     54596 | 1.77M | 61 |      1031 | 25.89K | 334 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  98.94% |     58171 | 1.77M | 52 |      1634 | 14.45K | 170 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:25 |
| Q0L0X80P001   |   80.0 |  99.03% |     59995 | 1.77M | 49 |      2561 | 19.06K | 193 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  99.10% |     57108 | 1.77M | 53 |      1634 | 17.29K | 212 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  99.15% |     55088 | 1.77M | 60 |      1022 | 21.59K | 349 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  99.07% |     58133 | 1.77M | 57 |      1036 | 21.34K | 327 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q25L60X40P002 |   40.0 |  99.03% |     51111 | 1.77M | 57 |      1071 | 22.73K | 308 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:25 |
| Q25L60X80P000 |   80.0 |  99.05% |     60002 | 1.77M | 49 |      1634 | 16.01K | 175 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:25 |
| Q25L60X80P001 |   80.0 |  99.13% |     58214 | 1.77M | 49 |      3503 |  16.7K | 198 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:26 |
| Q25L60X80P002 |   80.0 |  99.22% |     60017 | 1.77M | 49 |      3503 | 17.12K | 210 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q30L60X40P000 |   40.0 |  99.20% |     45124 | 1.77M | 66 |      1035 | 27.05K | 382 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.07% |     54035 | 1.77M | 62 |      1634 | 21.29K | 363 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:25 |
| Q30L60X40P002 |   40.0 |  99.10% |     51121 | 1.78M | 64 |      1634 | 20.69K | 378 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q30L60X80P000 |   80.0 |  99.18% |     58188 | 1.77M | 50 |      1634 | 18.25K | 258 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q30L60X80P001 |   80.0 |  99.16% |     68727 | 1.77M | 48 |      3503 | 17.39K | 230 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:28 |
| Q30L60X80P002 |   80.0 |  99.18% |     60010 | 1.77M | 48 |      1634 | 18.56K | 248 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:27 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.99% |     58198 | 1.77M | 50 |      1634 | 14.63K | 102 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| MRX40P001 |   40.0 |  98.00% |     59967 | 1.76M | 50 |      3503 | 17.44K | 101 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| MRX40P002 |   40.0 |  98.02% |     59974 | 1.77M | 49 |      1606 | 16.28K | 100 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX80P000 |   80.0 |  97.98% |     58097 | 1.77M | 51 |      1634 | 14.82K | 104 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:24 |
| MRX80P001 |   80.0 |  98.01% |     58116 | 1.77M | 49 |      1634 | 15.43K | 101 |   79.0 | 7.5 |  18.8 | 152.2 | "31,41,51,61,71,81" |   0:00:31 |   0:00:23 |
| MRX80P002 |   80.0 |  97.93% |     57045 | 1.76M | 51 |      3503 | 17.84K | 106 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:33 |   0:00:23 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.04% |     68651 | 1.77M | 47 |     24768 | 265.89K | 62 |  270.0 | 20.0 |  70.0 | 495.0 |   0:00:31 |
| 7_merge_mr_unitigs_bcalm      |  98.71% |     59937 | 1.76M | 47 |      3503 |  11.73K |  4 |  267.0 | 24.0 |  65.0 | 508.5 |   0:00:46 |
| 7_merge_mr_unitigs_superreads |  98.87% |     59937 | 1.76M | 47 |     25128 |  77.65K |  9 |  269.0 | 24.0 |  65.7 | 511.5 |   0:00:53 |
| 7_merge_mr_unitigs_tadpole    |  98.52% |     59935 | 1.76M | 47 |      3503 |  11.73K |  4 |  267.0 | 22.0 |  67.0 | 499.5 |   0:00:43 |
| 7_merge_unitigs_bcalm         |  98.68% |     60041 | 1.77M | 47 |     24786 | 145.26K | 45 |  268.0 | 23.0 |  66.3 | 505.5 |   0:00:48 |
| 7_merge_unitigs_superreads    |  98.70% |     59980 | 1.77M | 47 |      1634 |  29.47K | 17 |  269.0 | 19.5 |  70.2 | 491.2 |   0:00:50 |
| 7_merge_unitigs_tadpole       |  98.79% |     68688 | 1.77M | 47 |     17479 | 113.81K | 35 |  272.0 | 19.0 |  71.7 | 493.5 |   0:00:50 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.33% |    131396 | 1.78M | 25 |      3298 | 11.28K | 53 |  271.0 | 26.5 |  63.8 | 525.8 |   0:00:29 |
| 8_mr_spades  |  99.31% |    161544 | 1.79M | 19 |      1634 |  9.46K | 39 |  428.0 | 47.0 |  95.7 | 853.5 |   0:00:32 |
| 8_megahit    |  98.93% |     60004 | 1.77M | 47 |      1765 | 14.13K | 99 |  271.0 | 26.5 |  63.8 | 525.8 |   0:00:28 |
| 8_mr_megahit |  99.48% |    113733 | 1.79M | 27 |      1634 |  8.38K | 53 |  428.0 | 45.5 |  97.2 | 846.8 |   0:00:32 |
| 8_platanus   |  99.24% |    131369 | 1.78M | 19 |      3251 |  9.38K | 38 |  272.0 | 27.0 |  63.7 | 529.5 |   0:00:29 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1830138 | 1830138 |   1 |
| Paralogs                 |    5432 |   95244 |  29 |
| 7_merge_anchors.anchors  |   68651 | 1768223 |  47 |
| 7_merge_anchors.others   |   24768 |  265887 |  62 |
| glue_anchors             |   68651 | 1766633 |  45 |
| fill_anchors             |  125667 | 1773775 |  27 |
| spades.contig            |  131566 | 1797765 |  53 |
| spades.scaffold          |  131672 | 1798239 |  47 |
| spades.non-contained     |  131566 | 1791727 |  28 |
| mr_spades.contig         |  161695 | 1799629 |  30 |
| mr_spades.scaffold       |  161695 | 1799785 |  28 |
| mr_spades.non-contained  |  161695 | 1796047 |  20 |
| megahit.contig           |   60052 | 1796189 |  75 |
| megahit.non-contained    |   60052 | 1786744 |  52 |
| mr_megahit.contig        |  131755 | 1803952 |  36 |
| mr_megahit.non-contained |  131755 | 1799110 |  26 |
| platanus.contig          |  107567 | 1806867 | 138 |
| platanus.scaffold        |  131416 | 1799789 |  75 |
| platanus.non-contained   |  131416 | 1791351 |  19 |


# Leg_pneumop_pneumophila_Philadelphia_1


```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=Leg_pneumop_pneumophila_Philadelphia_1

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 3397754 \
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
| R.genome.bbtools  | 447.4 |    346 | 1898.4 |                         98.46% |
| R.tadpole.bbtools | 353.4 |    344 |   99.0 |                         93.25% |
| R.genome.picard   | 355.7 |    346 |   98.8 |                             FR |
| R.tadpole.picard  | 353.4 |    344 |   98.5 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 237       | 3333288         | 0.0000       | 40.56   |
| R.31 | 207       | 3385691         | 0.0000       | 40.15   |
| R.41 | 175       | 3387017         | 0.0000       | 39.90   |
| R.51 | 144       | 3387374         | 0.0000       | 39.73   |
| R.61 | 114       | 3389541         | 0.0000       | 39.58   |
| R.71 | 86        | 3418151         | 0.0000       | 39.45   |
| R.81 | 57        | 3424423         | 0.0000       | 39.31   |


Table: statReads

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 3397754 | 3397754 |        1 |
| Paralogs   |    2793 |  100404 |       43 |
| Illumina.R |     101 |   1.06G | 10498482 |
| trim.R     |     100 |  960.7M |  9702334 |
| Q25L60     |     100 | 900.56M |  9151013 |
| Q30L60     |     100 | 796.54M |  8269502 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.06G | 10456976 |
| highpass | 101 |   1.05G | 10354546 |
| sample   | 101 |   1.02G | 10092338 |
| trim     | 100 |  960.7M |  9702334 |
| filter   | 100 |  960.7M |  9702334 |
| R1       | 100 | 480.63M |  4851167 |
| R2       | 100 | 480.06M |  4851167 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6594	0.06534%
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
#unique_kmers	23041325
#error_kmers	19673427
#genomic_kmers	3367898
#main_peak	193
#genome_size_in_peaks	3391554
#genome_size	3399605
#haploid_genome_size	3399605
#fold_coverage	193
#haploid_fold_coverage	193
#ploidy	1
#percent_repeat_in_peaks	0.710
#percent_repeat	0.907
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 960.68M | 9702206 |
| ecco          | 100 | 960.68M | 9702206 |
| eccc          | 100 | 960.68M | 9702206 |
| ecct          | 100 | 955.45M | 9648558 |
| extended      | 140 |   1.34G | 9648558 |
| merged.raw    | 357 |   1.06G | 3119606 |
| unmerged.raw  | 140 |  470.5M | 3409346 |
| unmerged.trim | 140 | 470.49M | 3409244 |
| M1            | 357 |   1.06G | 3113660 |
| U1            | 140 | 235.79M | 1704622 |
| U2            | 140 |  234.7M | 1704622 |
| Us            |   0 |       0 |       0 |
| M.cor         | 322 |   1.53G | 9636564 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 156.1 |    160 |  22.5 |          2.64% |
| M.ihist.merge.txt  | 340.4 |    348 |  57.5 |         64.67% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 282.7 |  256.4 |    9.32% | "71" |  3.4M | 3.49M |     1.03 | 0:01'56'' |
| Q25L60.R | 265.1 |  245.8 |    7.30% | "71" |  3.4M | 3.41M |     1.00 | 0:01'51'' |
| Q30L60.R | 234.6 |  221.0 |    5.81% | "71" |  3.4M | 3.39M |     1.00 | 0:01'40'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.99% |     31266 | 3.28M | 177 |        50 | 39.45K | 746 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:33 |
| Q0L0X40P001   |   40.0 |  98.96% |     31707 | 3.29M | 179 |        62 | 44.62K | 747 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:32 |
| Q0L0X40P002   |   40.0 |  98.94% |     32451 | 3.29M | 185 |        51 | 41.67K | 772 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:31 |
| Q0L0X80P000   |   80.0 |  98.50% |     25528 | 3.35M | 189 |        39 | 21.42K | 505 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| Q0L0X80P001   |   80.0 |  98.41% |     33790 | 3.35M | 171 |        38 |  18.6K | 452 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:30 |
| Q0L0X80P002   |   80.0 |  98.39% |     27814 | 3.35M | 189 |        39 | 22.34K | 492 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  99.13% |     30816 |  3.1M | 160 |       502 | 41.98K | 562 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:32 |
| Q25L60X40P001 |   40.0 |  99.18% |     43831 | 3.05M | 136 |        83 | 33.45K | 527 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:33 |
| Q25L60X40P002 |   40.0 |  99.20% |     39846 |  3.1M | 151 |       558 | 36.43K | 530 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:33 |
| Q25L60X80P000 |   80.0 |  98.89% |     49053 | 3.37M | 113 |        40 |  15.5K | 348 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:34 |
| Q25L60X80P001 |   80.0 |  99.00% |     50518 | 3.35M | 111 |        43 | 17.03K | 349 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q25L60X80P002 |   80.0 |  98.96% |     47852 | 3.39M | 115 |        43 | 18.11K | 377 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:35 |
| Q30L60X40P000 |   40.0 |  99.28% |     24669 | 2.19M | 155 |       676 | 42.71K | 428 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q30L60X40P001 |   40.0 |  99.34% |     26936 |  2.4M | 167 |       759 |  49.5K | 475 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:35 |
| Q30L60X40P002 |   40.0 |  99.30% |     29404 | 2.43M | 151 |       774 | 45.62K | 451 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  99.25% |     79379 | 3.39M |  83 |        73 | 17.63K | 302 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |
| Q30L60X80P001 |   80.0 |  99.18% |     61235 | 3.36M |  90 |        46 | 17.01K | 334 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:37 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.87% |     52467 | 3.26M |  96 |       100 | 24.12K | 252 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:41 |   0:00:35 |
| MRX40P001 |   40.0 |  98.86% |     50553 | 2.95M |  95 |       125 | 26.74K | 246 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| MRX40P002 |   40.0 |  98.84% |     48194 | 3.39M | 101 |       115 | 27.09K | 257 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:40 |   0:00:33 |
| MRX80P000 |   80.0 |  98.60% |     56220 | 3.34M |  99 |        91 | 26.37K | 281 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:07 |   0:00:32 |
| MRX80P001 |   80.0 |  98.53% |     63787 | 3.34M |  92 |       100 | 26.83K | 260 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:34 |
| MRX80P002 |   80.0 |  98.53% |     58258 | 3.34M |  96 |       102 | 27.54K | 262 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:06 |   0:00:35 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.12% |     45897 | 3.15M | 128 |        80 | 35.36K | 586 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:05 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  99.08% |     48772 | 3.09M | 119 |       460 | 35.32K | 546 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:05 |   0:00:34 |
| Q0L0X40P002   |   40.0 |  99.13% |     43476 | 2.82M | 124 |       700 | 36.92K | 586 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  99.08% |    116071 | 3.36M |  62 |        51 | 17.92K | 317 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:01:12 |   0:00:38 |
| Q0L0X80P001   |   80.0 |  99.19% |    145174 | 3.35M |  53 |        52 | 19.34K | 337 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:39 |
| Q0L0X80P002   |   80.0 |  99.09% |    102769 | 3.35M |  61 |        54 | 18.57K | 328 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:01:12 |   0:00:38 |
| Q25L60X40P000 |   40.0 |  99.21% |     41232 | 2.88M | 134 |       513 | 41.74K | 555 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:36 |
| Q25L60X40P001 |   40.0 |  99.16% |     50920 | 2.73M | 112 |       654 | 40.96K | 550 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:36 |
| Q25L60X40P002 |   40.0 |  99.22% |     44676 | 2.84M | 132 |       558 |  39.5K | 593 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:36 |
| Q25L60X80P000 |   80.0 |  99.17% |    132935 | 3.36M |  59 |        47 | 16.77K | 326 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:37 |
| Q25L60X80P001 |   80.0 |  99.24% |    197503 | 3.35M |  57 |        57 | 19.61K | 339 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:42 |
| Q25L60X80P002 |   80.0 |  99.22% |    132915 | 3.35M |  56 |        45 | 17.44K | 334 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:38 |
| Q30L60X40P000 |   40.0 |  99.19% |     23676 | 2.76M | 190 |       747 | 59.97K | 658 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:35 |
| Q30L60X40P001 |   40.0 |  99.20% |     43226 | 3.35M | 115 |        54 | 33.41K | 626 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:33 |
| Q30L60X40P002 |   40.0 |  99.25% |     66952 | 3.35M | 111 |       763 | 45.85K | 677 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:36 |
| Q30L60X80P000 |   80.0 |  99.34% |     95810 | 3.19M |  72 |       567 | 27.87K | 386 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:12 |   0:00:39 |
| Q30L60X80P001 |   80.0 |  99.40% |    106661 | 3.35M |  64 |        46 | 21.52K | 433 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:10 |   0:00:41 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.88% |     94438 | 2.43M | 48 |      1103 | 15.24K | 104 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:07 |   0:00:33 |
| MRX40P001 |   40.0 |  98.59% |     79303 | 2.04M | 51 |      1103 | 12.48K |  82 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:05 |   0:00:29 |
| MRX40P002 |   40.0 |  98.63% |     93655 | 2.09M | 48 |      1103 | 13.71K |  89 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:06 |   0:00:30 |
| MRX80P000 |   80.0 |  98.67% |    145126 | 3.35M | 45 |      1103 | 14.36K | 107 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:20 |   0:00:32 |
| MRX80P001 |   80.0 |  98.54% |    231521 | 3.35M | 41 |      1103 | 14.06K |  98 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:21 |   0:00:30 |
| MRX80P002 |   80.0 |  98.74% |    142054 | 3.35M | 47 |      1103 | 13.54K | 108 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:01:19 |   0:00:34 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.10% |     50441 | 2.97M | 116 |       355 | 29.56K | 415 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  99.16% |     45828 | 2.86M | 108 |       290 | 28.92K | 441 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q0L0X40P002   |   40.0 |  99.20% |     55385 | 2.68M |  98 |        83 | 25.07K | 397 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q0L0X80P000   |   80.0 |  99.06% |    109781 | 3.35M |  62 |        55 | 14.06K | 244 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:38 |
| Q0L0X80P001   |   80.0 |  99.12% |    133554 | 3.53M |  58 |       112 | 13.73K | 209 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q0L0X80P002   |   80.0 |  99.13% |    132943 | 3.47M |  59 |        67 | 15.49K | 248 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:39 |
| Q25L60X40P000 |   40.0 |  99.23% |     41685 | 2.22M | 108 |       896 | 30.63K | 369 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:33 |
| Q25L60X40P001 |   40.0 |  99.31% |     50924 | 2.46M |  96 |       743 | 32.08K | 375 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:35 |
| Q25L60X40P002 |   40.0 |  99.26% |     39194 | 2.58M | 118 |       606 | 27.07K | 380 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:33 |
| Q25L60X80P000 |   80.0 |  99.18% |    121956 | 3.35M |  55 |        50 | 14.22K | 250 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:36 |   0:00:35 |
| Q25L60X80P001 |   80.0 |  99.24% |    142118 | 3.35M |  57 |        49 | 13.87K | 250 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:39 |
| Q25L60X80P002 |   80.0 |  99.26% |    142118 | 3.35M |  55 |        54 | 14.42K | 256 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:39 |
| Q30L60X40P000 |   40.0 |  99.27% |     25555 | 2.14M | 140 |       805 | 40.77K | 392 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:33 |
| Q30L60X40P001 |   40.0 |  99.36% |     28348 | 2.23M | 151 |       759 | 48.19K | 433 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q30L60X40P002 |   40.0 |  99.40% |     27290 | 2.33M | 145 |       760 | 53.12K | 501 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  99.30% |     95815 | 3.18M |  68 |       508 | 18.61K | 266 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:37 |
| Q30L60X80P001 |   80.0 |  99.31% |     96045 | 3.35M |  67 |        75 |  17.6K | 283 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:37 |   0:00:36 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.78% |     94438 |  2.4M | 47 |      1103 | 14.67K |  94 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:36 |   0:00:32 |
| MRX40P001 |   40.0 |  98.68% |     79339 | 2.26M | 56 |      1103 | 13.02K |  91 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| MRX40P002 |   40.0 |  98.68% |     93330 | 2.58M | 58 |      1103 | 14.55K | 103 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:30 |
| MRX80P000 |   80.0 |  98.52% |    142065 |  3.5M | 55 |      1103 |  15.2K | 125 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:31 |
| MRX80P001 |   80.0 |  98.60% |    155319 | 3.35M | 50 |      1103 |  14.9K | 113 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:41 |   0:00:33 |
| MRX80P002 |   80.0 |  98.56% |    127338 | 3.35M | 53 |      1103 | 13.61K | 114 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:33 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  98.68% |    202282 | 3.35M | 38 |     29318 | 421.42K | 92 |  256.0 | 4.0 |  81.3 | 402.0 |   0:00:49 |
| 7_merge_mr_unitigs_bcalm      |  99.33% |    202291 | 3.35M | 37 |      1103 |  11.65K | 10 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:18 |
| 7_merge_mr_unitigs_superreads |  99.44% |    198577 | 3.35M | 38 |      1179 |  11.96K | 10 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:27 |
| 7_merge_mr_unitigs_tadpole    |  99.40% |    198588 | 3.35M | 38 |     89750 | 105.04K | 12 |  256.0 | 4.0 |  81.3 | 402.0 |   0:01:26 |
| 7_merge_unitigs_bcalm         |  99.26% |    198578 | 3.35M | 38 |      1067 |  68.77K | 59 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:09 |
| 7_merge_unitigs_superreads    |  99.44% |    198589 | 3.35M | 41 |     23585 | 228.72K | 55 |  255.0 | 5.0 |  80.0 | 405.0 |   0:01:30 |
| 7_merge_unitigs_tadpole       |  99.26% |    198582 | 3.35M | 39 |     10591 | 124.73K | 48 |  255.0 | 5.0 |  80.0 | 405.0 |   0:01:10 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.28% |    274686 | 2.93M | 24 |      1763 | 30.93K | 51 |  255.0 | 12.5 |  72.5 | 438.8 |   0:00:40 |
| 8_mr_spades  |  99.47% |    274709 | 3.36M | 24 |      1359 | 24.09K | 50 |  451.0 |  6.5 | 143.8 | 705.8 |   0:00:45 |
| 8_megahit    |  99.15% |    202196 | 3.35M | 37 |      1441 | 27.72K | 81 |  255.0 |  6.0 |  79.0 | 409.5 |   0:00:39 |
| 8_mr_megahit |  99.48% |    274969 | 3.36M | 24 |      1203 | 20.85K | 54 |  451.0 | 12.5 | 137.8 | 732.8 |   0:00:47 |
| 8_platanus   |  98.83% |    274614 | 3.35M | 27 |      1072 | 10.07K | 51 |  255.0 |  5.0 |  80.0 | 405.0 |   0:00:40 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 3397754 | 3397754 |   1 |
| Paralogs                 |    2793 |  100404 |  43 |
| 7_merge_anchors.anchors  |  202282 | 3349797 |  38 |
| 7_merge_anchors.others   |   29318 |  421421 |  92 |
| glue_anchors             |  202282 | 3349028 |  37 |
| fill_anchors             |  274595 | 3551622 |  27 |
| spades.contig            |  431777 | 3398852 |  72 |
| spades.scaffold          |  431777 | 3398962 |  70 |
| spades.non-contained     |  431777 | 3388369 |  28 |
| mr_spades.contig         |  351394 | 3390527 |  40 |
| mr_spades.scaffold       |  351394 | 3390527 |  40 |
| mr_spades.non-contained  |  351394 | 3383835 |  26 |
| megahit.contig           |  248588 | 3390749 |  66 |
| megahit.non-contained    |  248588 | 3379367 |  44 |
| mr_megahit.contig        |  275034 | 3402982 |  68 |
| mr_megahit.non-contained |  275034 | 3383640 |  30 |
| platanus.contig          |  198751 | 3392384 | 197 |
| platanus.scaffold        |  363088 | 3385816 | 148 |
| platanus.non-contained   |  363088 | 3364324 |  24 |


# N_gon_FA_1090


```shell script
WORKING_DIR=${HOME}/data/anchr/fda_argos
BASE_NAME=N_gon_FA_1090

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 2153922 \
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
| R.genome.bbtools  | 459.4 |    327 | 2186.9 |                         97.89% |
| R.tadpole.bbtools | 328.6 |    319 |   92.9 |                         77.41% |
| R.genome.picard   | 337.3 |    327 |   94.3 |                             FR |
| R.tadpole.picard  | 328.7 |    319 |   92.8 |                             FR |


Table: statKAT

| k    | mean_freq | est_genome_size | est_het_rate | mean_gc |
|:-----|:----------|:----------------|:-------------|:--------|
| R.21 | 242       | 0               | 0.0000       | 54.49   |
| R.31 | 416       | 1618265         | 0.0000       | 54.41   |
| R.41 | 356       | 1893318         | 0.0000       | 54.32   |
| R.51 | 291       | 1920356         | 0.0000       | 54.23   |
| R.61 | 234       | 1996598         | 0.0000       | 54.15   |
| R.71 | 174       | 2055792         | 0.0000       | 54.07   |
| R.81 | 123       | 2208326         | 0.0000       | 53.99   |


Table: statReads

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 2153922 | 2153922 |        1 |
| Paralogs   |    4318 |  132685 |       45 |
| Illumina.R |     101 |   1.49G | 14768158 |
| trim.R     |     100 | 563.83M |  5767612 |
| Q25L60     |     100 | 494.57M |  5130182 |
| Q30L60     |     100 | 379.87M |  4105219 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.49G | 14705264 |
| highpass | 101 |   1.46G | 14500526 |
| sample   | 101 | 646.18M |  6397788 |
| trim     | 100 | 563.83M |  5767612 |
| filter   | 100 | 563.83M |  5767612 |
| R1       | 100 | 282.49M |  2883806 |
| R2       | 100 | 281.34M |  2883806 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	5556	0.08684%
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
#unique_kmers	14363434
#error_kmers	12325498
#genomic_kmers	2037936
#main_peak	171
#genome_size_in_peaks	2138176
#genome_size	2221208
#haploid_genome_size	2221208
#fold_coverage	171
#haploid_fold_coverage	171
#ploidy	1
#percent_repeat_in_peaks	4.762
#percent_repeat	6.489
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 563.82M | 5767502 |
| ecco          | 100 | 563.81M | 5767502 |
| eccc          | 100 | 563.81M | 5767502 |
| ecct          | 100 |  560.8M | 5736160 |
| extended      | 140 | 787.23M | 5736160 |
| merged.raw    | 345 | 663.83M | 2014412 |
| unmerged.raw  | 140 | 230.53M | 1707336 |
| unmerged.trim | 140 | 230.52M | 1707232 |
| M1            | 345 | 662.75M | 2011113 |
| U1            | 140 | 115.73M |  853616 |
| U2            | 140 | 114.79M |  853616 |
| Us            |   0 |       0 |       0 |
| M.cor         | 318 | 895.28M | 5729458 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 155.7 |    160 |  23.3 |          3.36% |
| M.ihist.merge.txt  | 329.5 |    334 |  58.3 |         70.24% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 261.8 |  235.4 |   10.06% | "51" | 2.15M | 2.24M |     1.04 | 0:01'11'' |
| Q25L60.R | 229.8 |  212.4 |    7.57% | "51" | 2.15M | 2.07M |     0.96 | 0:01'05'' |
| Q30L60.R | 176.7 |  166.6 |    5.73% | "47" | 2.15M | 2.04M |     0.95 | 0:00'54'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  94.65% |      7210 | 1.99M | 390 |      1032 | 136.48K | 1447 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:26 |
| Q0L0X40P001   |   40.0 |  94.41% |      7174 | 1.96M | 380 |      1026 | 129.66K | 1363 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  94.46% |      7927 | 1.98M | 351 |      1087 | 120.54K | 1378 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:19 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  90.86% |      6807 | 1.94M | 390 |      1053 |  67.03K |  892 |   76.0 | 6.0 |  19.3 | 141.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:31 |
| Q0L0X80P001   |   80.0 |  90.68% |      7527 | 1.94M | 379 |      1075 |  62.68K |  878 |   76.0 | 6.0 |  19.3 | 141.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  96.07% |     10141 |  1.9M | 273 |      1258 | 184.37K | 1176 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  96.02% |     10360 |  1.9M | 278 |      1140 | 186.18K | 1215 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  96.03% |     10165 | 1.91M | 269 |      1142 |    179K | 1161 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  95.06% |     12379 | 1.93M | 250 |      1140 |  88.77K |  742 |   76.0 | 6.0 |  19.3 | 141.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:27 |
| Q25L60X80P001 |   80.0 |  95.02% |     10302 | 1.91M | 270 |      1203 | 139.23K |  830 |   76.0 | 5.0 |  20.3 | 136.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:28 |
| Q30L60X40P000 |   40.0 |  96.48% |      6380 | 1.61M | 368 |      1265 | 378.81K | 1342 |   37.0 | 3.0 |   9.3 |  69.0 | "31,41,51,61,71,81" |   0:00:18 |   0:00:28 |
| Q30L60X40P001 |   40.0 |  96.56% |      7632 | 1.77M | 331 |      1223 |  295.8K | 1337 |   37.0 | 4.0 |   8.3 |  73.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:29 |
| Q30L60X40P002 |   40.0 |  96.45% |      8266 | 1.81M | 297 |      1473 | 286.56K | 1267 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:18 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  96.48% |      8027 | 1.74M | 308 |      1244 | 253.96K | 1077 |   75.0 | 6.0 |  19.0 | 139.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:30 |
| Q30L60X80P001 |   80.0 |  96.52% |      9639 |  1.8M | 274 |      1375 | 235.96K | 1078 |   75.0 | 7.0 |  18.0 | 144.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:31 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  92.30% |     15152 | 1.95M | 189 |      1142 | 65.76K | 396 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| MRX40P001 |   40.0 |  92.34% |     17916 | 1.95M | 179 |      1201 | 59.73K | 384 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| MRX40P002 |   40.0 |  92.19% |     15443 | 1.96M | 208 |      1142 | 53.92K | 431 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| MRX80P000 |   80.0 |  91.70% |     12953 | 1.95M | 220 |      1096 | 61.48K | 472 |   77.0 | 6.0 |  19.7 | 142.5 | "31,41,51,61,71,81" |   0:00:43 |   0:00:24 |
| MRX80P001 |   80.0 |  91.59% |     11872 | 1.95M | 246 |      1095 | 62.88K | 514 |   77.0 | 5.5 |  20.2 | 140.2 | "31,41,51,61,71,81" |   0:00:42 |   0:00:25 |
| MRX80P002 |   80.0 |  91.55% |     12754 | 1.96M | 238 |      1095 | 52.92K | 503 |   77.0 | 7.0 |  18.7 | 147.0 | "31,41,51,61,71,81" |   0:00:42 |   0:00:26 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.54% |     11025 | 1.94M | 263 |      1095 | 130.35K | 1375 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:28 |
| Q0L0X40P001   |   40.0 |  95.46% |     10846 | 1.91M | 264 |      1107 | 143.56K | 1327 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q0L0X40P002   |   40.0 |  95.49% |     12038 | 1.93M | 247 |      1095 | 104.27K | 1361 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q0L0X80P000   |   80.0 |  95.85% |     16787 | 1.94M | 211 |      1107 |  118.8K | 1078 |   76.0 | 5.0 |  20.3 | 136.5 | "31,41,51,61,71,81" |   0:00:55 |   0:00:29 |
| Q0L0X80P001   |   80.0 |  95.96% |     15351 | 1.97M | 207 |      1116 |  93.15K | 1026 |   77.0 | 5.0 |  20.7 | 138.0 | "31,41,51,61,71,81" |   0:00:56 |   0:00:31 |
| Q25L60X40P000 |   40.0 |  95.88% |      9104 | 1.91M | 295 |      1244 | 204.58K | 1485 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  95.85% |      8802 |  1.9M | 305 |      1137 |  203.9K | 1453 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |
| Q25L60X40P002 |   40.0 |  95.63% |     10004 | 1.93M | 282 |      1076 | 118.46K | 1502 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  96.48% |     14445 | 1.93M | 219 |      1115 | 132.65K | 1177 |   77.0 | 6.0 |  19.7 | 142.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:32 |
| Q25L60X80P001 |   80.0 |  96.54% |     14172 | 1.93M | 222 |      1096 | 127.58K | 1203 |   77.0 | 6.0 |  19.7 | 142.5 | "31,41,51,61,71,81" |   0:00:55 |   0:00:32 |
| Q30L60X40P000 |   40.0 |  94.37% |      5669 | 1.83M | 447 |      1456 | 268.51K | 1627 |   38.0 | 5.0 |   7.7 |  76.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:26 |
| Q30L60X40P001 |   40.0 |  94.46% |      5633 |  1.8M | 434 |      1123 | 319.81K | 1682 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  93.97% |      5406 |  1.8M | 443 |      1446 | 328.23K | 1581 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  96.60% |      7895 | 1.86M | 327 |      1307 | 310.12K | 1548 |   75.0 | 7.0 |  18.0 | 144.0 | "31,41,51,61,71,81" |   0:00:57 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  96.50% |      9395 | 1.88M | 293 |      1221 | 246.18K | 1467 |   76.0 | 8.0 |  17.3 | 150.0 | "31,41,51,61,71,81" |   0:00:58 |   0:00:31 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  92.54% |     18956 | 1.94M | 157 |      1265 | 55.09K | 330 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:23 |
| MRX40P001 |   40.0 |  92.51% |     19390 | 1.96M | 153 |      1265 | 56.09K | 328 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:22 |
| MRX40P002 |   40.0 |  92.49% |     18488 | 1.96M | 162 |      1429 | 49.33K | 333 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:53 |   0:00:22 |
| MRX80P000 |   80.0 |  92.26% |     19940 | 1.97M | 153 |      1406 | 48.33K | 316 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:59 |   0:00:24 |
| MRX80P001 |   80.0 |  92.31% |     18529 | 1.96M | 160 |      1265 | 49.87K | 325 |   78.0 | 5.0 |  21.0 | 139.5 | "31,41,51,61,71,81" |   0:01:00 |   0:00:24 |
| MRX80P002 |   80.0 |  92.29% |     18413 | 1.96M | 162 |      1361 | 49.51K | 335 |   77.0 | 5.5 |  20.2 | 140.2 | "31,41,51,61,71,81" |   0:00:59 |   0:00:24 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |    # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|--------:|-----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  95.95% |     12484 | 1.95M | 239 |      1093 | 127.74K | 1093 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q0L0X40P001   |   40.0 |  95.92% |     13567 | 1.93M | 237 |      1166 | 138.72K | 1108 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q0L0X40P002   |   40.0 |  96.04% |     13944 | 1.94M | 233 |      1157 | 125.41K | 1110 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:28 |
| Q0L0X80P000   |   80.0 |  95.47% |     15525 | 1.96M | 206 |      1265 |  89.43K |  790 |   76.0 | 5.0 |  20.3 | 136.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:29 |
| Q0L0X80P001   |   80.0 |  95.35% |     14550 | 1.96M | 201 |      1226 |  77.49K |  786 |   77.0 | 5.0 |  20.7 | 138.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |
| Q25L60X40P000 |   40.0 |  96.55% |     10791 | 1.88M | 248 |      1361 | 205.95K | 1274 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:29 |
| Q25L60X40P001 |   40.0 |  96.35% |     10936 | 1.88M | 255 |      1107 |  185.3K | 1191 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:28 |
| Q25L60X40P002 |   40.0 |  96.37% |     11987 |  1.9M | 264 |      1322 | 188.68K | 1219 |   38.0 | 3.0 |   9.7 |  70.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:28 |
| Q25L60X80P000 |   80.0 |  96.12% |     14233 | 1.91M | 210 |      1207 | 111.59K |  830 |   77.0 | 5.0 |  20.7 | 138.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:29 |
| Q25L60X80P001 |   80.0 |  96.30% |     13876 | 1.92M | 215 |      1107 | 128.44K |  917 |   77.0 | 5.0 |  20.7 | 138.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  96.26% |      7709 | 1.86M | 338 |      1367 |  279.4K | 1429 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q30L60X40P001 |   40.0 |  96.38% |      8157 | 1.81M | 327 |      1142 |  278.4K | 1449 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  96.28% |      7726 | 1.84M | 318 |      1436 | 322.36K | 1463 |   38.0 | 4.0 |   8.7 |  75.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:30 |
| Q30L60X80P000 |   80.0 |  96.75% |      8027 | 1.74M | 311 |      1181 | 307.94K | 1257 |   75.0 | 6.0 |  19.0 | 139.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:31 |
| Q30L60X80P001 |   80.0 |  96.76% |      9399 | 1.83M | 280 |      1389 | 256.71K | 1193 |   75.0 | 7.0 |  18.0 | 144.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:32 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  92.33% |     18956 | 1.94M | 156 |      1406 | 54.58K | 313 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| MRX40P001 |   40.0 |  92.34% |     19390 | 1.96M | 154 |      1265 | 54.94K | 318 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:22 |
| MRX40P002 |   40.0 |  92.18% |     18475 | 1.96M | 167 |      1408 | 48.47K | 327 |   39.0 | 3.0 |  10.0 |  72.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| MRX80P000 |   80.0 |  92.25% |     19402 | 1.96M | 158 |      1265 | 49.61K | 331 |   78.0 | 6.0 |  20.0 | 144.0 | "31,41,51,61,71,81" |   0:00:33 |   0:00:24 |
| MRX80P001 |   80.0 |  92.19% |     18299 | 1.96M | 167 |      1265 | 51.47K | 336 |   78.0 | 5.0 |  21.0 | 139.5 | "31,41,51,61,71,81" |   0:00:34 |   0:00:24 |
| MRX80P002 |   80.0 |  92.25% |     16956 | 1.97M | 170 |      1406 | 50.06K | 353 |   77.0 | 6.0 |  19.7 | 142.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:24 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |   # | N50Others |     Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|----:|----------:|--------:|----:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  92.17% |     20064 | 1.95M | 155 |      2854 |  684.4K | 309 |  226.0 | 12.0 |  63.3 | 393.0 |   0:00:35 |
| 7_merge_mr_unitigs_bcalm      |  95.90% |     20003 | 1.96M | 152 |      4287 | 146.23K |  69 |  225.0 | 14.0 |  61.0 | 400.5 |   0:00:48 |
| 7_merge_mr_unitigs_superreads |  95.85% |     19396 | 1.96M | 152 |      2093 | 109.08K |  65 |  225.0 | 15.0 |  60.0 | 405.0 |   0:00:49 |
| 7_merge_mr_unitigs_tadpole    |  95.88% |     18914 | 1.96M | 154 |     11053 | 178.78K |  66 |  225.0 | 14.0 |  61.0 | 400.5 |   0:00:49 |
| 7_merge_unitigs_bcalm         |  94.79% |     18752 | 1.97M | 170 |      1648 | 401.48K | 241 |  226.0 | 14.0 |  61.3 | 402.0 |   0:00:46 |
| 7_merge_unitigs_superreads    |  96.34% |     20565 | 1.97M | 156 |      2077 | 376.84K | 203 |  224.0 | 14.0 |  60.7 | 399.0 |   0:00:57 |
| 7_merge_unitigs_tadpole       |  95.72% |     19919 | 1.96M | 158 |      1642 | 351.72K | 213 |  225.0 | 14.0 |  61.0 | 400.5 |   0:00:50 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  96.66% |     39021 | 1.84M |  93 |      2090 | 55.15K | 165 |  228.0 | 16.0 |  60.0 | 414.0 |   0:00:32 |
| 8_mr_spades  |  97.38% |     44045 | 1.95M |  84 |      1640 | 42.97K | 156 |  406.0 | 25.5 | 109.8 | 723.8 |   0:00:34 |
| 8_megahit    |  95.32% |     20548 | 1.97M | 153 |      1408 | 59.28K | 301 |  227.0 | 14.0 |  61.7 | 403.5 |   0:00:29 |
| 8_mr_megahit |  98.44% |     42044 | 2.03M |  91 |      1238 | 60.95K | 190 |  407.0 | 30.0 | 105.7 | 745.5 |   0:00:33 |
| 8_platanus   |  95.31% |     36989 | 1.89M | 102 |      1384 |  44.2K | 177 |  227.0 | 14.0 |  61.7 | 403.5 |   0:00:30 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 2153922 | 2153922 |   1 |
| Paralogs                 |    4318 |  132685 |  45 |
| 7_merge_anchors.anchors  |   20064 | 1950306 | 155 |
| 7_merge_anchors.others   |    2854 |  684403 | 309 |
| glue_anchors             |   20064 | 1950306 | 155 |
| fill_anchors             |   40639 | 1975927 |  83 |
| spades.contig            |   55104 | 2090954 | 291 |
| spades.scaffold          |   55104 | 2091154 | 289 |
| spades.non-contained     |   57275 | 2054844 |  78 |
| mr_spades.contig         |   50736 | 2085330 | 124 |
| mr_spades.scaffold       |   50736 | 2085530 | 122 |
| mr_spades.non-contained  |   50736 | 2067812 |  74 |
| megahit.contig           |   23014 | 2074576 | 251 |
| megahit.non-contained    |   23067 | 2032021 | 148 |
| mr_megahit.contig        |   41100 | 2160638 | 281 |
| mr_megahit.non-contained |   42380 | 2093013 |  99 |
| platanus.contig          |   20641 | 2142016 | 826 |
| platanus.scaffold        |   46758 | 2104943 | 517 |
| platanus.non-contained   |   46995 | 2039521 |  87 |

