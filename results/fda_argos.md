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
- [Ca_jej_jejuni_NCTC_11168_ATCC_700819](#ca_jej_jejuni_nctc_11168_atcc_700819)
- [Clostridio_dif_630](#clostridio_dif_630)
- [Co_dip_NCTC_13129](#co_dip_nctc_13129)
- [Fr_tul_tularensis_SCHU_S4](#fr_tul_tularensis_schu_s4)
- [Ha_inf_Rd_KW20](#ha_inf_rd_kw20)
- [Legionella pneumophila subsp. pneumophila ATCC 33152D-5](#legionella-pneumophila-subsp-pneumophila-atcc-33152d-5)
  - [Lpne: reference](#lpne-reference)
  - [Lpne: download](#lpne-download)
  - [Lpne: template](#lpne-template)
  - [Lpne: run](#lpne-run)


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
| R.tadpole.bbtools | 310.1 |    299 |  103.7 |                         96.43% |
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
| trim.R     |     100 | 530.04M |  5325364 |
| Q25L60     |     100 |  516.1M |  5196839 |
| Q30L60     |     100 | 489.69M |  4979409 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.54G | 15283902 |
| highpass | 101 |   1.54G | 15224722 |
| sample   | 101 | 549.04M |  5436054 |
| trim     | 100 | 530.08M |  5325832 |
| filter   | 100 | 530.04M |  5325364 |
| R1       | 100 | 265.32M |  2662682 |
| R2       | 100 | 264.72M |  2662682 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	4095	0.07533%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	468	0.00879%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	8176693
#error_kmers	6560383
#genomic_kmers	1616310
#main_peak	221
#genome_size_in_peaks	1640136
#genome_size	1641667
#haploid_genome_size	1641667
#fold_coverage	221
#haploid_fold_coverage	221
#ploidy	1
#percent_repeat_in_peaks	1.463
#percent_repeat	1.525
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 530.03M | 5325284 |
| ecco          | 100 | 530.03M | 5325284 |
| eccc          | 100 | 530.03M | 5325284 |
| ecct          | 100 | 527.04M | 5294872 |
| extended      | 140 |  738.4M | 5294872 |
| merged.raw    | 333 | 688.75M | 2151939 |
| unmerged.raw  | 140 | 137.52M |  990994 |
| unmerged.trim | 140 | 137.52M |  990988 |
| M1            | 333 |  688.1M | 2149922 |
| U1            | 140 |  68.99M |  495494 |
| U2            | 140 |  68.53M |  495494 |
| Us            |   0 |       0 |       0 |
| M.cor         | 317 | 827.77M | 5290832 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 159.6 |    165 |  20.8 |          4.65% |
| M.ihist.merge.txt  | 320.1 |    321 |  58.6 |         81.28% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 289.6 |  273.0 |    5.75% | "71" | 1.83M | 1.62M |     0.89 | 0:01'07'' |
| Q25L60.R | 282.0 |  268.9 |    4.68% | "71" | 1.83M | 1.62M |     0.89 | 0:01'05'' |
| Q30L60.R | 267.7 |  257.5 |    3.81% | "71" | 1.83M | 1.62M |     0.88 | 0:01'04'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.68% |     66811 | 1.59M | 50 |      2340 | 13.67K | 202 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:20 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.79% |     57909 | 1.59M | 51 |      2340 | 14.81K | 230 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q0L0X40P002   |   40.0 |  98.74% |     66531 |  1.6M | 56 |      2340 | 15.47K | 228 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:23 |
| Q0L0X80P000   |   80.0 |  98.08% |     40407 | 1.59M | 69 |      2340 | 11.62K | 148 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q0L0X80P001   |   80.0 |  98.29% |     43899 | 1.59M | 67 |      2340 | 12.16K | 142 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q0L0X80P002   |   80.0 |  98.22% |     46359 | 1.59M | 64 |      2340 | 12.17K | 136 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.73% |     56877 |  1.6M | 54 |      2340 |  14.3K | 221 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q25L60X40P001 |   40.0 |  98.78% |     70625 | 1.59M | 46 |      2340 | 15.27K | 208 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  98.66% |     79970 | 1.59M | 45 |      2340 | 13.26K | 180 |   45.0 | 3.0 |  12.0 |  81.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:22 |
| Q25L60X80P000 |   80.0 |  98.30% |     57893 | 1.59M | 54 |      6071 | 11.65K | 118 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q25L60X80P001 |   80.0 |  98.27% |     66536 | 1.59M | 57 |      6071 | 11.83K | 121 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q25L60X80P002 |   80.0 |  98.23% |     48276 | 1.59M | 54 |      6071 | 11.44K | 118 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q30L60X40P000 |   40.0 |  98.72% |     70556 | 1.59M | 43 |      2340 | 13.37K | 194 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:23 |
| Q30L60X40P001 |   40.0 |  98.86% |     74607 |  1.6M | 45 |      2340 | 15.11K | 223 |   45.0 | 3.0 |  12.0 |  81.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.90% |     60535 |  1.6M | 48 |      2340 | 14.41K | 201 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:18 |   0:00:23 |
| Q30L60X80P000 |   80.0 |  98.42% |     69808 | 1.59M | 53 |      6071 | 12.02K | 138 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:23 |
| Q30L60X80P001 |   80.0 |  98.37% |     70552 | 1.59M | 51 |      2340 |  12.6K | 124 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |
| Q30L60X80P002 |   80.0 |  98.32% |     66680 | 1.59M | 50 |      2340 | 12.66K | 128 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:22 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.67% |     50164 | 1.59M | 55 |      2340 | 14.43K | 115 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:23 |
| MRX40P001 |   40.0 |  97.79% |     75073 | 1.59M | 47 |      2340 | 13.79K |  99 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:21 |
| MRX40P002 |   40.0 |  97.83% |     79910 | 1.59M | 43 |      2340 | 13.28K |  93 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:22 |
| MRX80P000 |   80.0 |  97.58% |     57869 | 1.59M | 62 |      2340 | 14.97K | 130 |   88.0 | 6.0 |  23.3 | 159.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:23 |
| MRX80P001 |   80.0 |  97.60% |     64662 | 1.59M | 57 |      2340 | 14.09K | 120 |   88.0 | 6.0 |  23.3 | 159.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:23 |
| MRX80P002 |   80.0 |  97.49% |     50162 | 1.59M | 60 |      2340 |  14.7K | 124 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:22 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.88% |     79903 | 1.59M | 40 |      2490 | 14.92K | 234 |   45.0 | 2.0 |  13.0 |  76.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| Q0L0X40P001   |   40.0 |  98.69% |     79871 | 1.59M | 45 |      2340 |  16.2K | 224 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q0L0X40P002   |   40.0 |  98.94% |     68219 |  1.6M | 44 |      2340 | 16.04K | 234 |   44.0 | 2.0 |  12.7 |  75.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  98.75% |     79989 |  1.6M | 43 |      2340 | 13.94K | 200 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:26 |
| Q0L0X80P001   |   80.0 |  98.61% |     79949 | 1.59M | 43 |      2340 |    13K | 172 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  98.79% |     80735 |  1.6M | 43 |      2340 | 14.01K | 202 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:25 |
| Q25L60X40P000 |   40.0 |  98.95% |     55758 | 1.59M | 50 |      1024 | 21.82K | 255 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  98.79% |     70854 | 1.59M | 47 |      1053 | 18.11K | 255 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  98.99% |     71573 | 1.59M | 42 |      1069 | 17.64K | 207 |   44.0 | 2.0 |  12.7 |  75.0 | "31,41,51,61,71,81" |   0:00:45 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  98.82% |     71558 |  1.6M | 42 |      2340 | 14.87K | 213 |   90.0 | 6.5 |  23.5 | 164.2 | "31,41,51,61,71,81" |   0:00:49 |   0:00:25 |
| Q25L60X80P001 |   80.0 |  98.80% |     83938 | 1.59M | 39 |      2340 | 14.24K | 181 |   90.0 | 6.0 |  24.0 | 162.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  98.84% |     79976 | 1.59M | 36 |      2340 | 14.76K | 203 |   90.0 | 6.0 |  24.0 | 162.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.71% |     78883 | 1.59M | 48 |      2340 | 14.83K | 212 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:23 |
| Q30L60X40P001 |   40.0 |  98.85% |     70822 |  1.6M | 47 |      2609 |    16K | 258 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:45 |   0:00:23 |
| Q30L60X40P002 |   40.0 |  98.86% |     71545 | 1.59M | 45 |      2340 | 16.76K | 246 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  98.81% |     71579 | 1.59M | 42 |      2340 | 15.26K | 201 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:25 |
| Q30L60X80P001 |   80.0 |  98.88% |     71570 | 1.59M | 43 |      2340 | 15.39K | 205 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:26 |
| Q30L60X80P002 |   80.0 |  98.91% |     79941 | 1.59M | 43 |      2340 | 15.61K | 209 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:26 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.33% |    104120 |  1.6M | 34 |      2340 | 12.94K | 92 |   45.0 | 4.0 |  11.0 |  85.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:24 |
| MRX40P001 |   40.0 |  98.19% |    104091 | 1.59M | 37 |      2340 | 13.55K | 95 |   45.0 | 3.0 |  12.0 |  81.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:22 |
| MRX40P002 |   40.0 |  98.19% |     90197 | 1.59M | 40 |      2340 | 13.36K | 93 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:22 |
| MRX80P000 |   80.0 |  97.84% |     80643 | 1.59M | 41 |      2340 | 13.22K | 89 |   88.0 | 6.0 |  23.3 | 159.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:23 |
| MRX80P001 |   80.0 |  97.82% |     80629 | 1.59M | 42 |      2340 | 13.05K | 92 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:54 |   0:00:23 |
| MRX80P002 |   80.0 |  97.96% |     83989 | 1.59M | 40 |      2340 | 13.26K | 90 |   89.0 | 5.5 |  24.2 | 158.2 | "31,41,51,61,71,81" |   0:00:53 |   0:00:23 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.87% |     79983 | 1.59M | 42 |      2340 | 13.98K | 199 |   45.0 | 2.0 |  13.0 |  76.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q0L0X40P001   |   40.0 |  98.86% |     79885 |  1.6M | 43 |      2340 | 14.57K | 193 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:25 |
| Q0L0X40P002   |   40.0 |  98.99% |     71561 |  1.6M | 45 |      2340 | 15.38K | 188 |   45.0 | 3.0 |  12.0 |  81.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q0L0X80P000   |   80.0 |  98.62% |     75183 | 1.59M | 43 |      6071 | 11.23K | 106 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  98.61% |     80762 | 1.59M | 42 |      6071 | 11.69K | 110 |   88.0 | 4.0 |  25.3 | 150.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q0L0X80P002   |   80.0 |  98.58% |     80745 | 1.59M | 44 |      6071 |  11.7K | 112 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| Q25L60X40P000 |   40.0 |  98.70% |     55840 | 1.59M | 47 |      2340 | 14.94K | 190 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q25L60X40P001 |   40.0 |  98.83% |     71559 |  1.6M | 45 |      2340 | 14.79K | 198 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X40P002 |   40.0 |  98.96% |     80733 |  1.6M | 42 |      2340 | 15.91K | 185 |   44.0 | 2.0 |  12.7 |  75.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q25L60X80P000 |   80.0 |  98.71% |     75174 | 1.59M | 43 |      6071 | 12.05K | 126 |   88.0 | 5.0 |  24.3 | 154.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:23 |
| Q25L60X80P001 |   80.0 |  98.68% |     79976 | 1.59M | 41 |      6071 | 12.01K | 116 |   88.0 | 4.0 |  25.3 | 150.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  98.54% |     75182 | 1.59M | 42 |      6071 | 11.67K | 116 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  98.72% |     70818 | 1.59M | 44 |      2340 | 14.35K | 191 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:22 |
| Q30L60X40P001 |   40.0 |  98.95% |     51994 |  1.6M | 56 |      1037 | 18.46K | 235 |   44.0 | 2.0 |  12.7 |  75.0 | "31,41,51,61,71,81" |   0:00:23 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  98.93% |     79927 |  1.6M | 45 |      2340 | 16.78K | 219 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:23 |   0:00:25 |
| Q30L60X80P000 |   80.0 |  98.82% |     75185 | 1.59M | 40 |      2340 | 12.91K | 142 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:26 |
| Q30L60X80P001 |   80.0 |  98.66% |     75180 | 1.59M | 41 |      2340 | 12.69K | 130 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  98.69% |     79977 | 1.59M | 40 |      6071 |  11.7K | 124 |   89.0 | 5.0 |  24.7 | 156.0 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.90% |     80640 | 1.59M | 40 |      2340 | 12.86K |  87 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX40P001 |   40.0 |  97.97% |    104091 | 1.59M | 39 |      2340 | 13.03K |  85 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX40P002 |   40.0 |  97.84% |     90197 | 1.59M | 41 |      2340 | 12.87K |  85 |   44.0 | 3.0 |  11.7 |  79.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX80P000 |   80.0 |  97.76% |     79879 | 1.59M | 44 |      2340 | 12.81K |  90 |   89.0 | 7.0 |  22.7 | 165.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:22 |
| MRX80P001 |   80.0 |  97.80% |     79887 | 1.59M | 45 |      2340 | 12.86K |  92 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:23 |
| MRX80P002 |   80.0 |  97.97% |     79884 | 1.59M | 46 |      2340 | 13.73K | 102 |   89.0 | 6.0 |  23.7 | 160.5 | "31,41,51,61,71,81" |   0:00:33 |   0:00:23 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.94% |    112408 | 1.59M | 23 |    177862 | 271.07K | 41 |  301.0 | 17.0 |  83.3 | 528.0 |   0:00:34 |
| 7_merge_mr_unitigs_bcalm      |  98.91% |    104621 | 1.59M | 25 |      6071 |  20.27K |  7 |  301.0 | 16.0 |  84.3 | 523.5 |   0:00:48 |
| 7_merge_mr_unitigs_superreads |  98.85% |     90146 | 1.59M | 33 |      6071 |   8.41K |  2 |  303.0 | 15.0 |  86.0 | 522.0 |   0:00:45 |
| 7_merge_mr_unitigs_tadpole    |  98.85% |    104100 | 1.59M | 28 |      6071 |   8.96K |  3 |  300.0 | 15.0 |  85.0 | 517.5 |   0:00:46 |
| 7_merge_unitigs_bcalm         |  98.91% |    104137 |  1.6M | 28 |    177862 | 227.26K | 26 |  302.0 | 18.0 |  82.7 | 534.0 |   0:00:47 |
| 7_merge_unitigs_superreads    |  99.08% |    104147 | 1.59M | 28 |      2340 |  12.47K |  5 |  303.0 | 16.0 |  85.0 | 526.5 |   0:00:54 |
| 7_merge_unitigs_tadpole       |  99.05% |    104147 |  1.6M | 29 |      6079 |   47.4K | 20 |  301.0 | 17.0 |  83.3 | 528.0 |   0:00:52 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |   Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.35% |    104593 | 1.23M | 19 |      2420 | 15.9K | 35 |  303.0 | 19.0 |  82.0 | 540.0 |   0:00:29 |
| 8_mr_spades  |  99.42% |    153954 | 1.61M | 16 |      2340 |   13K | 34 |  501.0 | 45.5 | 121.5 | 956.2 |   0:00:31 |
| 8_megahit    |  98.68% |    112616 |  1.6M | 24 |      6015 | 10.8K | 48 |  303.0 | 17.5 |  83.5 | 533.2 |   0:00:29 |
| 8_mr_megahit |  99.30% |    153889 | 1.61M | 22 |      2340 | 12.8K | 44 |  501.0 | 35.0 | 132.0 | 909.0 |   0:00:30 |
| 8_platanus   |  98.87% |    104629 | 1.11M | 18 |      5981 | 9.77K | 37 |  303.0 | 17.0 |  84.0 | 531.0 |   0:00:27 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1641481 | 1641481 |   1 |
| Paralogs                 |    6079 |   33258 |  13 |
| 7_merge_anchors.anchors  |  112408 | 1594034 |  23 |
| 7_merge_anchors.others   |  177862 |  271070 |  41 |
| glue_anchors             |  112408 | 1594034 |  23 |
| fill_anchors             |  128683 | 1603036 |  19 |
| spades.contig            |  153957 | 1622964 |  34 |
| spades.scaffold          |  189387 | 1623102 |  32 |
| spades.non-contained     |  153957 | 1616696 |  17 |
| mr_spades.contig         |  189482 | 1624572 |  25 |
| mr_spades.scaffold       |  189482 | 1624572 |  25 |
| mr_spades.non-contained  |  189482 | 1621940 |  18 |
| megahit.contig           |  112661 | 1622829 |  60 |
| megahit.non-contained    |  112661 | 1607338 |  24 |
| mr_megahit.contig        |  174584 | 1631643 |  41 |
| mr_megahit.non-contained |  174584 | 1621872 |  22 |
| platanus.contig          |  112554 | 1628911 | 108 |
| platanus.scaffold        |  153895 | 1622467 |  63 |
| platanus.non-contained   |  153895 | 1612680 |  21 |


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
| R.tadpole.bbtools | 374.0 |    363 |  122.5 |                         95.23% |
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

| Name       |     N50 |     Sum |        # |
|:-----------|--------:|--------:|---------:|
| Genome     | 1892775 | 1892775 |        1 |
| Paralogs   |   33912 |   93528 |       10 |
| Illumina.R |     101 |   2.14G | 21230270 |
| trim.R     |     100 | 549.71M |  5517712 |
| Q25L60     |     100 | 538.67M |  5415699 |
| Q30L60     |     100 | 516.17M |  5233711 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   2.12G | 21017828 |
| highpass | 101 |   2.12G | 20971826 |
| sample   | 101 | 567.83M |  5622104 |
| trim     | 100 | 549.74M |  5518040 |
| filter   | 100 | 549.71M |  5517712 |
| R1       | 100 | 275.12M |  2758856 |
| R2       | 100 | 274.59M |  2758856 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3563	0.06337%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	328	0.00594%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	7552021
#error_kmers	5758816
#genomic_kmers	1793205
#main_peak	203
#genome_size_in_peaks	1882939
#genome_size	1883727
#haploid_genome_size	1883727
#fold_coverage	203
#haploid_fold_coverage	203
#ploidy	1
#percent_repeat_in_peaks	4.766
#percent_repeat	4.323
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 549.71M | 5517684 |
| ecco          | 100 | 549.71M | 5517684 |
| eccc          | 100 | 549.71M | 5517684 |
| ecct          | 100 | 547.35M | 5493746 |
| extended      | 140 | 766.64M | 5493746 |
| merged.raw    | 361 | 551.94M | 1615317 |
| unmerged.raw  | 140 | 315.17M | 2263112 |
| unmerged.trim | 140 | 315.17M | 2263108 |
| M1            | 361 | 551.23M | 1613269 |
| U1            | 140 | 157.86M | 1131554 |
| U2            | 140 | 157.31M | 1131554 |
| Us            |   0 |       0 |       0 |
| M.cor         | 311 | 868.02M | 5489646 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 154.4 |    159 |  23.6 |          2.69% |
| M.ihist.merge.txt  | 341.7 |    351 |  59.2 |         58.81% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 290.4 |  275.8 |    5.03% | "71" | 1.89M | 1.81M |     0.95 | 0:01'11'' |
| Q25L60.R | 284.6 |  272.8 |    4.14% | "71" | 1.89M |  1.8M |     0.95 | 0:01'10'' |
| Q30L60.R | 272.8 |  263.7 |    3.34% | "71" | 1.89M |  1.8M |     0.95 | 0:01'07'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.25% |     27264 | 1.75M | 103 |      6263 | 63.93K | 458 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:20 |   0:00:35 |
| Q0L0X40P001   |   40.0 |  98.11% |     28057 | 1.76M |  99 |     12575 | 52.93K | 434 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  97.89% |     24353 | 1.76M |  95 |     26730 | 79.31K | 416 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  97.09% |     32349 | 1.76M |  74 |      6263 | 45.06K | 170 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  97.11% |     32659 | 1.76M |  77 |     15111 | 53.88K | 178 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  97.07% |     32701 | 1.76M |  75 |     15048 | 53.82K | 169 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  98.18% |     31595 | 1.76M |  85 |      8437 | 62.42K | 437 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  98.02% |     27270 | 1.76M |  92 |     15110 | 63.01K | 409 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  97.87% |     27528 | 1.76M |  98 |     15049 | 71.42K | 398 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  97.02% |     32735 | 1.76M |  71 |     26730 |  46.6K | 162 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P001 |   80.0 |  97.06% |     31607 | 1.76M |  79 |      6334 | 45.34K | 180 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  97.14% |     32696 | 1.76M |  72 |      8437 | 47.08K | 167 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:25 |
| Q30L60X40P000 |   40.0 |  98.17% |     31377 | 1.76M |  87 |     26730 | 79.14K | 415 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:28 |
| Q30L60X40P001 |   40.0 |  98.03% |     28685 | 1.76M |  98 |     27554 | 51.51K | 419 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:16 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  98.08% |     32332 | 1.75M |  87 |     26730 | 81.41K | 416 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  97.28% |     32345 | 1.76M |  75 |     15049 | 53.66K | 177 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X80P001 |   80.0 |  97.12% |     32758 | 1.76M |  69 |      8437 | 46.57K | 161 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q30L60X80P002 |   80.0 |  97.08% |     32760 | 1.76M |  72 |      8299 | 44.47K | 167 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.60% |     30610 | 1.75M | 83 |     26730 |  72.9K | 170 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| MRX40P001 |   40.0 |  96.71% |     31535 | 1.75M | 81 |     26730 | 72.73K | 175 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:21 |   0:00:23 |
| MRX40P002 |   40.0 |  96.71% |     30688 | 1.75M | 80 |     26730 | 72.57K | 173 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:21 |   0:00:23 |
| MRX80P000 |   80.0 |  96.62% |     32292 | 1.75M | 72 |     26730 | 71.91K | 166 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:22 |
| MRX80P001 |   80.0 |  96.56% |     32655 | 1.75M | 69 |      6752 | 45.13K | 161 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:23 |
| MRX80P002 |   80.0 |  96.63% |     32668 | 1.75M | 72 |     26730 | 71.99K | 165 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:23 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.08% |     27549 | 1.64M | 99 |      4106 | 67.21K | 518 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q0L0X40P001   |   40.0 |  97.91% |     28739 | 1.75M | 95 |      5358 | 55.75K | 468 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q0L0X40P002   |   40.0 |  97.94% |     24958 | 1.75M | 93 |      5358 | 62.35K | 518 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  97.58% |     35050 | 1.76M | 73 |     27554 | 46.94K | 355 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q0L0X80P001   |   80.0 |  97.61% |     32755 | 1.76M | 72 |     27554 |  46.9K | 347 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:26 |
| Q0L0X80P002   |   80.0 |  97.65% |     32741 | 1.76M | 73 |     27554 | 47.27K | 357 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  97.93% |     31069 | 1.76M | 88 |      5358 |  61.6K | 468 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:26 |
| Q25L60X40P001 |   40.0 |  98.03% |     27529 | 1.76M | 93 |      5358 | 56.03K | 504 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  98.02% |     27518 | 1.75M | 96 |      5358 | 64.01K | 493 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  97.48% |     32738 | 1.76M | 73 |     27554 | 46.01K | 358 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q25L60X80P001 |   80.0 |  97.62% |     35198 | 1.76M | 70 |     27554 | 49.37K | 363 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  97.73% |     32753 | 1.76M | 71 |     27554 | 47.39K | 336 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:25 |
| Q30L60X40P000 |   40.0 |  98.06% |     31364 | 1.76M | 89 |      5358 | 65.52K | 520 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:26 |
| Q30L60X40P001 |   40.0 |  98.02% |     30223 | 1.76M | 96 |     24037 | 88.35K | 495 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:25 |
| Q30L60X40P002 |   40.0 |  98.08% |     32280 | 1.73M | 88 |      5358 |  64.1K | 490 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  97.68% |     32746 | 1.76M | 74 |     27554 |  46.9K | 356 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:25 |
| Q30L60X80P001 |   80.0 |  97.71% |     35027 | 1.76M | 72 |     27554 | 49.47K | 375 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q30L60X80P002 |   80.0 |  97.72% |     32749 | 1.76M | 74 |     27554 | 47.73K | 383 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.33% |     32265 | 1.75M | 81 |     27554 | 45.01K | 152 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:53 |   0:00:24 |
| MRX40P001 |   40.0 |  96.35% |     31535 | 1.75M | 79 |     27554 |  44.4K | 150 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| MRX40P002 |   40.0 |  96.40% |     31312 | 1.75M | 79 |     27554 |  44.7K | 154 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:22 |
| MRX80P000 |   80.0 |  96.31% |     32682 | 1.75M | 69 |     27554 | 43.55K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:23 |
| MRX80P001 |   80.0 |  96.32% |     32655 | 1.75M | 69 |     27554 | 43.76K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:54 |   0:00:22 |
| MRX80P002 |   80.0 |  96.29% |     32668 | 1.75M | 70 |     27554 | 43.58K | 142 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:55 |   0:00:22 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  97.91% |     29005 | 1.63M | 101 |      5358 | 55.81K | 382 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:27 |
| Q0L0X40P001   |   40.0 |  97.73% |     28773 | 1.76M |  99 |      5358 | 56.78K | 374 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  97.83% |     24814 | 1.76M |  97 |     27554 | 54.48K | 390 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:25 |
| Q0L0X80P000   |   80.0 |  97.23% |     32754 | 1.76M |  70 |     27554 | 39.46K | 165 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  97.06% |     32755 | 1.76M |  69 |     27554 | 39.44K | 165 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| Q0L0X80P002   |   80.0 |  97.28% |     32741 | 1.76M |  70 |     27554 | 39.54K | 166 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  97.88% |     31377 | 1.76M |  87 |     18798 | 81.76K | 357 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q25L60X40P001 |   40.0 |  97.62% |     27218 | 1.76M |  94 |     27554 | 52.79K | 350 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:25 |
| Q25L60X40P002 |   40.0 |  97.86% |     30687 | 1.76M |  96 |      5358 | 60.38K | 392 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  97.18% |     32738 | 1.76M |  69 |     27554 | 39.51K | 169 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:25 |
| Q25L60X80P001 |   80.0 |  97.00% |     35202 | 1.76M |  68 |     27554 | 38.83K | 149 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| Q25L60X80P002 |   80.0 |  97.29% |     32757 | 1.76M |  69 |     27554 | 39.11K | 161 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  97.99% |     31179 | 1.76M |  94 |     27554 | 53.72K | 359 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q30L60X40P001 |   40.0 |  97.78% |     30232 | 1.76M |  95 |      5358 | 58.39K | 369 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X40P002 |   40.0 |  97.96% |     32309 | 1.73M |  85 |     27554 | 53.14K | 386 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |
| Q30L60X80P000 |   80.0 |  97.62% |     32746 | 1.76M |  70 |     27554 |  39.9K | 188 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:26 |
| Q30L60X80P001 |   80.0 |  97.38% |     35198 | 1.76M |  70 |     27554 | 39.84K | 181 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:26 |
| Q30L60X80P002 |   80.0 |  97.31% |     32761 | 1.76M |  70 |     27554 | 39.87K | 183 |   81.0 | 2.0 |  25.0 | 130.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:24 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  96.25% |     32265 | 1.75M | 82 |     27554 |  44.8K | 149 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:23 |
| MRX40P001 |   40.0 |  96.35% |     31535 | 1.75M | 79 |     27554 |  44.4K | 150 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX40P002 |   40.0 |  96.32% |     31312 | 1.75M | 79 |     27554 | 44.39K | 150 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:22 |
| MRX80P000 |   80.0 |  96.31% |     32682 | 1.75M | 69 |     27554 | 43.49K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:23 |
| MRX80P001 |   80.0 |  96.32% |     32655 | 1.75M | 69 |     27554 | 43.76K | 141 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| MRX80P002 |   80.0 |  96.29% |     32668 | 1.75M | 70 |     27554 | 43.58K | 142 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:22 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |   # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|----:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  94.49% |     36574 | 1.75M | 67 |     25655 | 553.52K | 137 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:32 |
| 7_merge_mr_unitigs_bcalm      |  97.82% |     36571 | 1.75M | 67 |     39198 | 215.54K |  14 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:56 |
| 7_merge_mr_unitigs_superreads |  97.80% |     32646 | 1.75M | 69 |     39198 | 242.25K |  15 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:55 |
| 7_merge_mr_unitigs_tadpole    |  97.84% |     36571 | 1.75M | 67 |     27554 | 270.54K |  16 |  277.0 | 3.0 |  89.3 | 429.0 |   0:00:58 |
| 7_merge_unitigs_bcalm         |  97.08% |     35147 | 1.76M | 68 |      1050 | 142.98K | 101 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:43 |
| 7_merge_unitigs_superreads    |  97.86% |     35167 | 1.76M | 68 |     22390 | 239.37K |  59 |  277.0 | 4.0 |  88.3 | 433.5 |   0:00:56 |
| 7_merge_unitigs_tadpole       |  97.59% |     35169 | 1.76M | 68 |      1170 | 118.32K |  80 |  278.0 | 4.0 |  88.7 | 435.0 |   0:00:50 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  97.52% |     36656 | 1.76M | 67 |     27548 | 40.28K | 133 |  278.0 |  4.0 |  88.7 | 435.0 |   0:00:29 |
| 8_mr_spades  |  97.22% |     37819 | 1.76M | 64 |     27661 | 43.14K | 128 |  463.0 |  9.0 | 145.3 | 735.0 |   0:00:33 |
| 8_megahit    |  97.19% |     36660 | 1.76M | 67 |     27556 | 39.46K | 137 |  278.0 |  3.5 |  89.2 | 432.8 |   0:00:32 |
| 8_mr_megahit |  97.45% |     36782 | 1.76M | 65 |     27844 | 53.84K | 134 |  463.0 | 11.0 | 143.3 | 744.0 |   0:00:33 |
| 8_platanus   |  96.95% |     36633 | 1.76M | 66 |     27530 | 39.22K | 131 |  278.0 |  4.0 |  88.7 | 435.0 |   0:00:30 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1892775 | 1892775 |   1 |
| Paralogs                 |   33912 |   93528 |  10 |
| 7_merge_anchors.anchors  |   36574 | 1753208 |  67 |
| 7_merge_anchors.others   |   25655 |  553517 | 137 |
| glue_anchors             |   36574 | 1753208 |  67 |
| fill_anchors             |   38073 | 1754747 |  63 |
| spades.contig            |   37811 | 1803249 |  81 |
| spades.scaffold          |   37811 | 1803303 |  79 |
| spades.non-contained     |   37811 | 1799490 |  66 |
| mr_spades.contig         |   37904 | 1811220 |  90 |
| mr_spades.scaffold       |   37904 | 1811220 |  90 |
| mr_spades.non-contained  |   37904 | 1805350 |  64 |
| megahit.contig           |   35250 | 1802644 |  77 |
| megahit.non-contained    |   35250 | 1798803 |  70 |
| mr_megahit.contig        |   35781 | 1822283 |  87 |
| mr_megahit.non-contained |   35781 | 1815258 |  69 |
| platanus.contig          |   35268 | 1808588 | 122 |
| platanus.scaffold        |   37808 | 1805453 |  98 |
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
| R.genome.bbtools  | 393.9 |    266 | 2164.5 |                         96.99% |
| R.tadpole.bbtools | 273.5 |    265 |   78.5 |                         94.87% |
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
| trim.R     |     100 | 525.06M |  5287012 |
| Q25L60     |     100 | 507.64M |  5131691 |
| Q30L60     |     100 |  472.6M |  4856039 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.23G | 12143650 |
| highpass | 101 |   1.22G | 12073026 |
| sample   | 101 | 549.04M |  5436054 |
| trim     | 100 |  525.1M |  5287474 |
| filter   | 100 | 525.06M |  5287012 |
| R1       | 100 | 262.88M |  2643506 |
| R2       | 100 | 262.18M |  2643506 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	3995	0.07349%
#Name	Reads	ReadsPct
```

```text
#R.filter
#Matched	462	0.00874%
#Name	Reads	ReadsPct
```

```text
#R.peaks
#k	31
#unique_kmers	7931177
#error_kmers	6143016
#genomic_kmers	1788161
#main_peak	200
#genome_size_in_peaks	1827985
#genome_size	1829762
#haploid_genome_size	1829762
#fold_coverage	200
#haploid_fold_coverage	200
#ploidy	1
#percent_repeat_in_peaks	2.186
#percent_repeat	2.074
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 525.05M | 5286878 |
| ecco          | 100 | 525.04M | 5286878 |
| eccc          | 100 | 525.04M | 5286878 |
| ecct          | 100 | 523.16M | 5267602 |
| extended      | 140 |  733.3M | 5267602 |
| merged.raw    | 311 |  714.6M | 2387753 |
| unmerged.raw  | 140 |   67.8M |  492096 |
| unmerged.trim | 140 |   67.8M |  492088 |
| M1            | 311 | 713.78M | 2384998 |
| U1            | 140 |  34.09M |  246044 |
| U2            | 140 |  33.71M |  246044 |
| Us            |   0 |       0 |       0 |
| M.cor         | 303 | 783.96M | 5262084 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 158.4 |    163 |  21.8 |         10.25% |
| M.ihist.merge.txt  | 299.3 |    297 |  61.2 |         90.66% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 286.9 |  271.6 |    5.34% | "71" | 1.83M |  1.8M |     0.98 | 0:01'08'' |
| Q25L60.R | 277.4 |  266.0 |    4.12% | "71" | 1.83M | 1.79M |     0.98 | 0:01'06'' |
| Q30L60.R | 258.3 |  250.3 |    3.12% | "71" | 1.83M | 1.79M |     0.98 | 0:01'02'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.01% |     37487 | 1.77M | 72 |      1022 | 24.38K | 365 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:24 |
| Q0L0X40P001   |   40.0 |  99.06% |     41479 | 1.77M | 73 |      1012 | 24.09K | 406 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  99.06% |     47850 | 1.77M | 64 |       970 | 23.63K | 357 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  98.67% |     34621 | 1.77M | 84 |      1298 | 17.91K | 216 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X80P001   |   80.0 |  98.09% |     39114 | 1.77M | 82 |      1824 | 14.24K | 192 |   79.0 | 6.5 |  19.8 | 147.8 | "31,41,51,61,71,81" |   0:00:25 |   0:00:24 |
| Q0L0X80P002   |   80.0 |  98.51% |     42469 | 1.77M | 79 |      1324 | 15.25K | 190 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| Q25L60X40P000 |   40.0 |  99.05% |     58131 | 1.77M | 58 |      1051 | 21.32K | 304 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  99.09% |     46923 | 1.77M | 63 |      1020 | 22.92K | 338 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:17 |   0:00:26 |
| Q25L60X40P002 |   40.0 |  99.13% |     53707 | 1.77M | 59 |      1022 | 20.14K | 360 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:25 |
| Q25L60X80P000 |   80.0 |  98.75% |     46444 | 1.77M | 69 |      1324 | 14.43K | 179 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |
| Q25L60X80P001 |   80.0 |  98.87% |     52624 | 1.77M | 63 |      2412 | 15.61K | 168 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  98.76% |     46780 | 1.77M | 68 |      1634 |  14.6K | 170 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:24 |
| Q30L60X40P000 |   40.0 |  99.26% |     55108 | 1.77M | 60 |       297 |  20.7K | 379 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.20% |     55130 | 1.77M | 59 |      1057 | 22.32K | 358 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:18 |   0:00:26 |
| Q30L60X40P002 |   40.0 |  99.29% |     57099 | 1.77M | 59 |       533 | 22.09K | 390 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:17 |   0:00:29 |
| Q30L60X80P000 |   80.0 |  99.10% |     57112 | 1.77M | 55 |      3503 | 14.67K | 173 |   80.0 | 7.5 |  19.2 | 153.8 | "31,41,51,61,71,81" |   0:00:26 |   0:00:25 |
| Q30L60X80P001 |   80.0 |  99.09% |     54625 | 1.77M | 57 |      1261 | 16.04K | 195 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q30L60X80P002 |   80.0 |  99.10% |     55133 | 1.77M | 59 |      1634 |  16.1K | 186 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:26 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.18% |     57101 | 1.76M | 54 |      1634 | 18.65K | 125 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:24 |
| MRX40P001 |   40.0 |  98.19% |     53977 | 1.76M | 54 |      1634 | 20.39K | 127 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:24 |
| MRX40P002 |   40.0 |  98.01% |     54995 | 1.76M | 56 |      1634 | 19.29K | 125 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:23 |
| MRX80P000 |   80.0 |  98.03% |     44317 | 1.76M | 64 |      1634 | 20.02K | 146 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:23 |
| MRX80P001 |   80.0 |  97.90% |     54521 | 1.76M | 64 |      1634 | 18.75K | 139 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:38 |   0:00:22 |
| MRX80P002 |   80.0 |  98.06% |     37440 | 1.77M | 64 |      1164 |  16.6K | 142 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:39 |   0:00:25 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.00% |     43886 | 1.77M | 69 |      1017 | 28.13K | 441 |   40.0 | 3.5 |   9.8 |  75.8 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q0L0X40P001   |   40.0 |  98.97% |     51078 | 1.77M | 64 |      1017 | 27.34K | 468 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:26 |
| Q0L0X40P002   |   40.0 |  98.97% |     48956 | 1.76M | 63 |      1016 | 28.09K | 443 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:46 |   0:00:27 |
| Q0L0X80P000   |   80.0 |  99.10% |     54600 | 1.77M | 55 |      1026 | 22.96K | 328 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:51 |   0:00:27 |
| Q0L0X80P001   |   80.0 |  99.06% |     53983 | 1.77M | 58 |      1634 |  19.2K | 306 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:27 |
| Q0L0X80P002   |   80.0 |  99.09% |     59991 | 1.77M | 50 |      1036 | 22.74K | 288 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:27 |
| Q25L60X40P000 |   40.0 |  99.05% |     44373 | 1.77M | 65 |      1015 | 23.79K | 485 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q25L60X40P001 |   40.0 |  98.99% |     43983 | 1.76M | 63 |      1012 | 28.33K | 468 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:48 |   0:00:27 |
| Q25L60X40P002 |   40.0 |  99.05% |     44057 | 1.76M | 64 |      1026 | 28.21K | 468 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:27 |
| Q25L60X80P000 |   80.0 |  99.13% |     54612 | 1.77M | 55 |      1048 | 19.01K | 310 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:50 |   0:00:28 |
| Q25L60X80P001 |   80.0 |  99.11% |     55121 | 1.77M | 54 |      1634 | 20.08K | 301 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:50 |   0:00:27 |
| Q25L60X80P002 |   80.0 |  99.15% |     55115 | 1.77M | 54 |      1634 | 20.18K | 297 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:29 |
| Q30L60X40P000 |   40.0 |  99.06% |     43833 | 1.77M | 74 |      1009 | 28.79K | 533 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:46 |   0:00:28 |
| Q30L60X40P001 |   40.0 |  99.05% |     40887 | 1.76M | 71 |      1018 | 25.71K | 486 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  99.00% |     43876 | 1.76M | 73 |      1006 | 26.47K | 485 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:47 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  99.24% |     54584 | 1.77M | 59 |      1030 | 20.07K | 374 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:52 |   0:00:30 |
| Q30L60X80P001 |   80.0 |  99.17% |     55112 | 1.77M | 56 |      1634 | 22.54K | 360 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:52 |   0:00:28 |
| Q30L60X80P002 |   80.0 |  99.22% |     57085 | 1.77M | 55 |      1016 | 22.12K | 370 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:51 |   0:00:29 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.98% |     59929 | 1.76M | 49 |      1634 |  17.3K | 100 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| MRX40P001 |   40.0 |  97.99% |     59937 | 1.76M | 47 |      3503 | 17.57K |  98 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:49 |   0:00:24 |
| MRX40P002 |   40.0 |  97.99% |     59982 | 1.76M | 47 |      3503 |  16.9K |  98 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:49 |   0:00:23 |
| MRX80P000 |   80.0 |  97.99% |     58176 | 1.77M | 52 |       414 | 12.41K | 104 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:55 |   0:00:24 |
| MRX80P001 |   80.0 |  97.99% |     59965 | 1.77M | 49 |      3503 | 16.01K | 101 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:55 |   0:00:25 |
| MRX80P002 |   80.0 |  98.05% |     59935 | 1.77M | 50 |      1164 | 12.07K | 101 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:56 |   0:00:23 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.02% |     54002 | 1.77M | 63 |      1032 |    23K | 339 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |
| Q0L0X40P001   |   40.0 |  99.11% |     44385 | 1.77M | 62 |      1022 | 23.29K | 341 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:27 |
| Q0L0X40P002   |   40.0 |  99.00% |     49029 | 1.77M | 66 |      1634 | 27.33K | 330 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |
| Q0L0X80P000   |   80.0 |  99.16% |     57100 | 1.77M | 52 |      1634 | 18.06K | 223 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:27 |
| Q0L0X80P001   |   80.0 |  99.03% |     57092 | 1.77M | 54 |      3503 | 16.06K | 172 |   79.0 | 6.0 |  20.3 | 145.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:25 |
| Q0L0X80P002   |   80.0 |  99.05% |     68710 | 1.77M | 50 |      3503 | 14.59K | 169 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:26 |
| Q25L60X40P000 |   40.0 |  99.03% |     55080 | 1.77M | 58 |      1038 | 22.95K | 318 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:24 |   0:00:25 |
| Q25L60X40P001 |   40.0 |  99.04% |     53956 | 1.77M | 60 |      1061 | 21.87K | 344 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:26 |
| Q25L60X40P002 |   40.0 |  99.13% |     57098 | 1.77M | 58 |      1634 | 18.38K | 311 |   39.0 | 4.0 |   9.0 |  76.5 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q25L60X80P000 |   80.0 |  99.13% |     59999 | 1.77M | 51 |      3503 | 14.45K | 173 |   79.0 | 7.5 |  18.8 | 152.2 | "31,41,51,61,71,81" |   0:00:29 |   0:00:27 |
| Q25L60X80P001 |   80.0 |  99.01% |     60025 | 1.77M | 48 |      3503 |  15.7K | 169 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:25 |
| Q25L60X80P002 |   80.0 |  99.14% |     58173 | 1.77M | 49 |      3503 | 17.18K | 196 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:27 |   0:00:26 |
| Q30L60X40P000 |   40.0 |  99.25% |     51121 | 1.77M | 68 |       479 | 24.15K | 394 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:26 |   0:00:27 |
| Q30L60X40P001 |   40.0 |  99.15% |     53967 | 1.78M | 59 |      1634 | 19.79K | 341 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q30L60X40P002 |   40.0 |  99.13% |     51061 | 1.77M | 62 |      1033 | 22.91K | 379 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:25 |   0:00:28 |
| Q30L60X80P000 |   80.0 |  99.21% |     58200 | 1.77M | 52 |      1634 | 17.08K | 247 |   80.0 | 6.0 |  20.7 | 147.0 | "31,41,51,61,71,81" |   0:00:30 |   0:00:27 |
| Q30L60X80P001 |   80.0 |  99.13% |     59993 | 1.77M | 51 |      1634 | 18.54K | 228 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:31 |   0:00:26 |
| Q30L60X80P002 |   80.0 |  99.12% |     60025 | 1.77M | 49 |      3503 | 17.24K | 222 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:28 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  97.98% |     59929 | 1.76M | 49 |      1634 | 17.21K | 100 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:23 |
| MRX40P001 |   40.0 |  97.99% |     59937 | 1.76M | 47 |      3503 | 17.61K |  98 |   40.0 | 3.0 |  10.3 |  73.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:22 |
| MRX40P002 |   40.0 |  98.00% |     59982 | 1.76M | 47 |      3503 | 16.87K |  98 |   40.0 | 4.0 |   9.3 |  78.0 | "31,41,51,61,71,81" |   0:00:28 |   0:00:24 |
| MRX80P000 |   80.0 |  97.97% |     57078 | 1.77M | 54 |       414 |  12.7K | 108 |   80.0 | 7.0 |  19.7 | 151.5 | "31,41,51,61,71,81" |   0:00:34 |   0:00:24 |
| MRX80P001 |   80.0 |  97.98% |     55024 | 1.77M | 51 |      3503 |  16.3K | 105 |   79.0 | 7.0 |  19.3 | 150.0 | "31,41,51,61,71,81" |   0:00:32 |   0:00:24 |
| MRX80P002 |   80.0 |  98.04% |     59935 | 1.77M | 52 |      1164 | 12.32K | 105 |   79.0 | 8.0 |  18.3 | 154.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:23 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|-----:|------:|------:|----------:|
| 7_merge_anchors               |  97.54% |     68651 | 1.77M | 47 |     11375 | 187.19K | 57 |  269.0 | 21.0 |  68.7 | 498.0 |   0:00:33 |
| 7_merge_mr_unitigs_bcalm      |  98.86% |     59943 | 1.77M | 48 |      3503 |  12.31K |  5 |  273.0 | 22.0 |  69.0 | 508.5 |   0:00:49 |
| 7_merge_mr_unitigs_superreads |  98.70% |     59941 | 1.76M | 47 |      3503 |  14.16K |  7 |  270.0 | 20.0 |  70.0 | 495.0 |   0:00:45 |
| 7_merge_mr_unitigs_tadpole    |  98.88% |     59942 | 1.77M | 48 |      3503 |  13.14K |  6 |  273.0 | 22.0 |  69.0 | 508.5 |   0:00:50 |
| 7_merge_unitigs_bcalm         |  98.68% |     59976 | 1.77M | 52 |      1033 |  58.49K | 44 |  268.0 | 26.0 |  63.3 | 519.0 |   0:00:45 |
| 7_merge_unitigs_superreads    |  98.95% |     68687 | 1.77M | 46 |      6374 |  62.76K | 24 |  273.0 | 19.0 |  72.0 | 495.0 |   0:00:53 |
| 7_merge_unitigs_tadpole       |  98.73% |     68683 | 1.77M | 47 |     31576 | 131.18K | 33 |  273.0 | 17.0 |  74.0 | 486.0 |   0:00:46 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median |  MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|-----:|------:|------:|----------:|
| 8_spades     |  99.32% |    131395 | 1.78M | 27 |      1765 | 10.44K | 56 |  272.0 | 28.0 |  62.7 | 534.0 |   0:00:29 |
| 8_mr_spades  |  99.26% |    113545 | 1.79M | 23 |      1634 |  8.97K | 45 |  428.0 | 47.5 |  95.2 | 855.8 |   0:00:32 |
| 8_megahit    |  98.89% |     77774 | 1.77M | 47 |      1765 | 14.42K | 98 |  271.0 | 27.0 |  63.3 | 528.0 |   0:00:29 |
| 8_mr_megahit |  99.52% |    121771 | 1.79M | 25 |      1634 |  8.61K | 51 |  428.0 | 44.5 |  98.2 | 842.2 |   0:00:30 |
| 8_platanus   |  99.26% |    131365 | 1.78M | 21 |      3413 |  6.79K | 40 |  272.0 | 30.0 |  60.7 | 543.0 |   0:00:30 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 1830138 | 1830138 |   1 |
| Paralogs                 |    5432 |   95244 |  29 |
| 7_merge_anchors.anchors  |   68651 | 1767887 |  47 |
| 7_merge_anchors.others   |   11375 |  187187 |  57 |
| glue_anchors             |   68651 | 1767887 |  47 |
| fill_anchors             |  107496 | 1775040 |  28 |
| spades.contig            |  131566 | 1797792 |  56 |
| spades.scaffold          |  131625 | 1798243 |  49 |
| spades.non-contained     |  131566 | 1791558 |  29 |
| mr_spades.contig         |  145889 | 1799572 |  31 |
| mr_spades.scaffold       |  147009 | 1799741 |  29 |
| mr_spades.non-contained  |  145889 | 1795278 |  22 |
| megahit.contig           |   68796 | 1795896 |  73 |
| megahit.non-contained    |   77818 | 1785852 |  51 |
| mr_megahit.contig        |  121871 | 1803617 |  35 |
| mr_megahit.non-contained |  121871 | 1799358 |  26 |
| platanus.contig          |  107526 | 1806417 | 132 |
| platanus.scaffold        |  131416 | 1800142 |  79 |
| platanus.non-contained   |  131416 | 1791177 |  19 |


# Legionella pneumophila subsp. pneumophila ATCC 33152D-5

## Lpne: reference

* Reference genome

```shell script
mkdir -p ~/data/anchr/Lpne/1_genome
cd ~/data/anchr/Lpne/1_genome

cp ~/data/anchr/ref/Lpne/genome.fa .
cp ~/data/anchr/ref/Lpne/paralogs.fa .

```

## Lpne: download

```shell script
cd ~/data/anchr/Lpne

mkdir -p ena
cd ena

cat << EOF > source.csv
SRX2179279,Lpne,HiSeq 2500 PE100
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| name | srx        | platform | layout | ilength | srr        | spot    | base      |
|:-----|:-----------|:---------|:-------|:--------|:-----------|:--------|:----------|
| Lpne | SRX2179279 | ILLUMINA | PAIRED | 570     | SRR4272054 | 5249241 | 1,011.23M |


* Illumina

```shell script
cd ~/data/anchr/Lpne

mkdir -p 2_illumina
cd 2_illumina

ln -s ../ena/SRR4272054_1.fastq.gz R1.fq.gz
ln -s ../ena/SRR4272054_2.fastq.gz R2.fq.gz

```

## Lpne: template

* template

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Lpne

cd ${WORKING_DIR}/${BASE_NAME}

rm *.sh
anchr template \
    --genome 3397754 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --fastqc \
    --insertsize \
    --kat \
    \
    --trim "--dedupe --cutoff 50 --cutk 31" \
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

```

## Lpne: run

```shell script
WORKING_DIR=${HOME}/data/anchr
BASE_NAME=Lpne

cd ${WORKING_DIR}/${BASE_NAME}
# rm -fr 4_*/ 6_*/ 7_*/ 8_*/
# rm -fr 2_illumina/trim 2_illumina/merge statReads.md 

# BASE_NAME=Lpne bash 0_bsub.sh
bsub -q mpi -n 24 -J "${BASE_NAME}-0_master" "bash 0_master.sh"
#bkill -J "${BASE_NAME}-*"

#bash 0_master.sh
#bash 0_cleanup.sh
#rm -fr 4_down_sampling 6_down_sampling

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
| trim.R     |     100 |  960.7M |  9702338 |
| Q25L60     |     100 | 900.56M |  9151030 |
| Q30L60     |     100 | 796.47M |  8268802 |


Table: statTrimReads

| Name     | N50 |     Sum |        # |
|:---------|----:|--------:|---------:|
| clumpify | 101 |   1.06G | 10456976 |
| highpass | 101 |   1.05G | 10351590 |
| sample   | 101 |   1.02G | 10092338 |
| trim     | 100 |  960.7M |  9702338 |
| filter   | 100 |  960.7M |  9702338 |
| R1       | 100 | 480.63M |  4851169 |
| R2       | 100 | 480.06M |  4851169 |
| Rs       |   0 |       0 |        0 |


```text
#R.trim
#Matched	6522	0.06462%
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
#unique_kmers	23034343
#error_kmers	19666489
#genomic_kmers	3367854
#main_peak	193
#genome_size_in_peaks	3390924
#genome_size	3399578
#haploid_genome_size	3399578
#fold_coverage	193
#haploid_fold_coverage	193
#ploidy	1
#percent_repeat_in_peaks	0.700
#percent_repeat	0.907
#start	center	stop	max	volume
```


Table: statMergeReads

| Name          | N50 |     Sum |       # |
|:--------------|----:|--------:|--------:|
| clumped       | 100 | 960.68M | 9702210 |
| ecco          | 100 | 960.68M | 9702210 |
| eccc          | 100 | 960.68M | 9702210 |
| ecct          | 100 | 955.48M | 9648822 |
| extended      | 140 |   1.34G | 9648822 |
| merged.raw    | 357 |   1.06G | 3120247 |
| unmerged.raw  | 140 | 470.38M | 3408328 |
| unmerged.trim | 140 | 470.37M | 3408230 |
| M1            | 357 |   1.06G | 3114277 |
| U1            | 140 | 235.73M | 1704115 |
| U2            | 140 | 234.64M | 1704115 |
| Us            |   0 |       0 |       0 |
| M.cor         | 322 |   1.53G | 9636784 |

| Group              |  Mean | Median | STDev | PercentOfPairs |
|:-------------------|------:|-------:|------:|---------------:|
| M.ihist.merge1.txt | 156.1 |    160 |  22.5 |          2.64% |
| M.ihist.merge.txt  | 340.4 |    348 |  57.5 |         64.68% |


Table: statQuorum

| Name     | CovIn | CovOut | Discard% | Kmer | RealG |  EstG | Est/Real |   RunTime |
|:---------|------:|-------:|---------:|-----:|------:|------:|---------:|----------:|
| Q0L0.R   | 282.7 |  256.4 |    9.33% | "71" |  3.4M | 3.47M |     1.02 | 0:01'52'' |
| Q25L60.R | 265.1 |  245.8 |    7.30% | "71" |  3.4M | 3.39M |     1.00 | 0:01'48'' |
| Q30L60.R | 234.6 |  220.9 |    5.81% | "71" |  3.4M | 3.38M |     0.99 | 0:01'38'' |


Table: statUnitigsSuperreads.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  98.90% |     27277 | 3.36M | 201 |       314 | 54.83K | 831 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:33 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  98.98% |     31817 | 3.35M | 176 |        65 | 43.78K | 732 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:32 |
| Q0L0X40P002   |   40.0 |  99.04% |     33858 | 3.28M | 171 |        69 |  47.2K | 788 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:28 |   0:00:32 |
| Q0L0X80P000   |   80.0 |  98.39% |     32256 | 3.35M | 184 |        38 | 20.66K | 483 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:31 |
| Q0L0X80P001   |   80.0 |  98.46% |     31247 | 3.35M | 173 |        37 | 18.67K | 463 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q0L0X80P002   |   80.0 |  98.44% |     26096 | 3.35M | 197 |        39 | 22.18K | 510 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:44 |   0:00:32 |
| Q25L60X40P000 |   40.0 |  99.20% |     33907 | 3.05M | 164 |        82 | 36.16K | 596 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:33 |
| Q25L60X40P001 |   40.0 |  99.24% |     30439 | 2.93M | 155 |        68 | 31.32K | 520 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:35 |
| Q25L60X40P002 |   40.0 |  99.27% |     43217 | 3.34M | 156 |       464 | 48.11K | 592 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:36 |
| Q25L60X80P000 |   80.0 |  98.97% |     60454 | 3.35M | 113 |        42 | 16.81K | 361 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:34 |
| Q25L60X80P001 |   80.0 |  98.93% |     59499 | 3.35M | 107 |        43 | 16.74K | 328 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:34 |
| Q25L60X80P002 |   80.0 |  98.90% |     43501 | 3.35M | 123 |        41 | 15.67K | 358 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:33 |
| Q30L60X40P000 |   40.0 |  99.26% |     21630 | 2.19M | 143 |       462 | 31.82K | 396 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:33 |
| Q30L60X40P001 |   40.0 |  99.28% |     29915 | 2.44M | 163 |       792 |  48.8K | 486 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q30L60X40P002 |   40.0 |  99.27% |     26361 | 2.72M | 167 |       920 | 48.99K | 500 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  99.25% |     77005 | 3.35M |  88 |        62 |  19.1K | 331 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:38 |
| Q30L60X80P001 |   80.0 |  99.29% |     77067 | 3.35M |  81 |        49 | 17.66K | 334 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:44 |   0:00:39 |


Table: statMRUnitigsSuperreads.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.75% |     55161 | 2.87M |  88 |       107 | 22.54K | 224 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:40 |   0:00:33 |
| MRX40P001 |   40.0 |  98.80% |     51560 |  3.2M |  98 |       103 | 24.19K | 241 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:40 |   0:00:34 |
| MRX40P002 |   40.0 |  98.91% |     62303 | 3.31M |  97 |       110 | 26.22K | 233 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:39 |   0:00:32 |
| MRX80P000 |   80.0 |  98.51% |     52536 | 3.34M | 100 |        93 | 24.88K | 265 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:01:05 |   0:00:32 |
| MRX80P001 |   80.0 |  98.57% |     55429 | 3.34M |  91 |        93 | 25.78K | 255 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:34 |
| MRX80P002 |   80.0 |  98.43% |     49545 | 3.34M | 104 |        93 | 25.14K | 271 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:05 |   0:00:31 |


Table: statUnitigsBcalm.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.07% |     41039 | 2.86M | 128 |       821 | 44.44K | 559 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:36 |
| Q0L0X40P001   |   40.0 |  99.14% |     42958 | 3.11M | 124 |        74 |    31K | 532 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:36 |
| Q0L0X40P002   |   40.0 |  99.23% |     43256 | 3.03M | 114 |       529 | 37.67K | 560 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:36 |
| Q0L0X80P000   |   80.0 |  99.05% |    107096 | 3.35M |  59 |        56 |  18.9K | 319 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:01:09 |   0:00:35 |
| Q0L0X80P001   |   80.0 |  99.13% |    132310 | 3.35M |  54 |        48 | 18.17K | 324 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:09 |   0:00:38 |
| Q0L0X80P002   |   80.0 |  99.15% |     90737 | 3.36M |  63 |        46 | 18.65K | 350 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:01:10 |   0:00:39 |
| Q25L60X40P000 |   40.0 |  99.21% |     34855 | 3.14M | 150 |       809 | 62.65K | 644 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:35 |
| Q25L60X40P001 |   40.0 |  99.14% |     29910 | 2.96M | 150 |       129 | 33.27K | 566 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:33 |
| Q25L60X40P002 |   40.0 |  99.22% |     39850 | 3.22M | 141 |       451 | 42.09K | 554 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:02 |   0:00:34 |
| Q25L60X80P000 |   80.0 |  99.28% |    118493 | 3.36M |  61 |        42 | 17.19K | 362 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:39 |
| Q25L60X80P001 |   80.0 |  99.36% |    102804 | 3.35M |  57 |        48 | 19.27K | 357 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:11 |   0:00:43 |
| Q25L60X80P002 |   80.0 |  99.25% |     92172 | 3.35M |  60 |       402 | 21.54K | 332 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:10 |   0:00:39 |
| Q30L60X40P000 |   40.0 |  99.23% |     21394 | 2.65M | 194 |       686 | 60.48K | 651 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:01 |   0:00:34 |
| Q30L60X40P001 |   40.0 |  99.24% |     56650 | 3.35M | 115 |       181 | 36.35K | 641 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:34 |
| Q30L60X40P002 |   40.0 |  99.25% |     60527 | 3.19M | 102 |        63 |  33.5K | 638 |   40.0 | 2.0 |  11.3 |  69.0 | "31,41,51,61,71,81" |   0:01:02 |   0:00:35 |
| Q30L60X80P000 |   80.0 |  99.40% |     97955 | 3.35M |  69 |        59 | 23.84K | 417 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:09 |   0:00:40 |
| Q30L60X80P001 |   80.0 |  99.32% |     97380 | 3.35M |  62 |        42 | 19.16K | 401 |   80.0 | 3.0 |  23.7 | 133.5 | "31,41,51,61,71,81" |   0:01:11 |   0:00:38 |


Table: statMRUnitigsBcalm.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.86% |     89282 |  2.1M | 47 |      1103 | 12.84K |  99 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:04 |   0:00:32 |
| MRX40P001 |   40.0 |  98.65% |     78339 | 2.52M | 59 |      1103 | 12.74K |  97 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:04 |   0:00:29 |
| MRX40P002 |   40.0 |  98.82% |     93690 | 2.35M | 59 |      1103 | 13.92K | 109 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:01:05 |   0:00:30 |
| MRX80P000 |   80.0 |  98.68% |    204299 | 3.35M | 44 |      1103 | 12.99K |  99 |   80.0 | 1.5 |  25.2 | 126.8 | "31,41,51,61,71,81" |   0:01:17 |   0:00:32 |
| MRX80P001 |   80.0 |  98.88% |    219256 | 3.35M | 39 |      1103 | 14.26K | 101 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:37 |
| MRX80P002 |   80.0 |  98.60% |    142049 | 3.35M | 44 |      1103 | 13.03K | 102 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:01:17 |   0:00:31 |


Table: statUnitigsTadpole.md

| Name          | CovCor | Mapped% | N50Anchor |   Sum |   # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:--------------|-------:|--------:|----------:|------:|----:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| Q0L0X40P000   |   40.0 |  99.17% |     42376 | 2.41M | 116 |       821 | 33.51K | 392 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q0L0X40P001   |   40.0 |  99.20% |     39685 | 2.93M | 117 |       760 | 29.41K | 423 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:35 |
| Q0L0X40P002   |   40.0 |  99.21% |     40968 | 2.78M | 114 |       975 | 38.11K | 416 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:32 |   0:00:33 |
| Q0L0X80P000   |   80.0 |  99.10% |    118902 | 3.26M |  64 |        53 | 14.54K | 249 |   79.0 | 1.0 |  25.3 | 123.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:37 |
| Q0L0X80P001   |   80.0 |  99.18% |    118501 | 3.35M |  58 |        67 | 15.91K | 248 |   79.0 | 2.0 |  24.3 | 127.5 | "31,41,51,61,71,81" |   0:00:34 |   0:00:38 |
| Q0L0X80P002   |   80.0 |  99.10% |     78101 |  3.2M |  74 |        60 | 15.56K | 252 |   79.0 | 1.0 |  25.3 | 123.0 | "31,41,51,61,71,81" |   0:00:34 |   0:00:37 |
| Q25L60X40P000 |   40.0 |  99.29% |     40324 | 2.66M | 116 |       623 | 32.08K | 403 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:34 |
| Q25L60X40P001 |   40.0 |  99.27% |     32738 |  2.5M | 114 |       878 | 29.87K | 355 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:34 |
| Q25L60X40P002 |   40.0 |  99.32% |     40626 | 2.96M | 126 |       342 | 32.75K | 419 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:35 |
| Q25L60X80P000 |   80.0 |  99.27% |    118518 | 3.36M |  62 |        50 | 14.45K | 268 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:36 |
| Q25L60X80P001 |   80.0 |  99.27% |    118524 | 3.35M |  54 |        55 | 14.29K | 253 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q25L60X80P002 |   80.0 |  99.23% |     81605 | 3.35M |  70 |        55 |  13.9K | 245 |   80.0 | 1.0 |  25.7 | 124.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:35 |
| Q30L60X40P000 |   40.0 |  99.37% |     21421 | 2.05M | 138 |       590 | 44.16K | 422 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:29 |   0:00:35 |
| Q30L60X40P001 |   40.0 |  99.38% |     27587 | 2.19M | 135 |       807 | 47.51K | 421 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:30 |   0:00:35 |
| Q30L60X40P002 |   40.0 |  99.38% |     32990 | 2.75M | 152 |       682 | 43.85K | 462 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:31 |   0:00:34 |
| Q30L60X80P000 |   80.0 |  99.38% |     79367 | 3.07M |  71 |        63 | 19.47K | 324 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:35 |   0:00:38 |
| Q30L60X80P001 |   80.0 |  99.28% |     97379 | 3.35M |  67 |        53 |  16.9K | 304 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:36 |   0:00:37 |


Table: statMRUnitigsTadpole.md

| Name      | CovCor | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |   # | median | MAD | lower | upper |                Kmer | RunTimeUT | RunTimeAN |
|:----------|-------:|--------:|----------:|------:|---:|----------:|-------:|----:|-------:|----:|------:|------:|--------------------:|----------:|----------:|
| MRX40P000 |   40.0 |  98.67% |     79287 | 2.01M | 50 |      1103 | 12.29K |  89 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:37 |   0:00:29 |
| MRX40P001 |   40.0 |  98.81% |     79280 | 2.33M | 54 |      1103 | 12.13K |  90 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:35 |   0:00:31 |
| MRX40P002 |   40.0 |  98.82% |     65134 | 1.92M | 57 |      1103 | 13.72K | 102 |   40.0 | 1.0 |  12.3 |  64.5 | "31,41,51,61,71,81" |   0:00:38 |   0:00:31 |
| MRX80P000 |   80.0 |  98.56% |    136963 | 3.35M | 53 |       219 | 14.03K | 115 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:43 |   0:00:30 |
| MRX80P001 |   80.0 |  98.55% |    132867 | 3.35M | 48 |      1103 | 14.62K | 106 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:31 |
| MRX80P002 |   80.0 |  98.63% |    118501 | 3.35M | 59 |       112 | 14.78K | 130 |   80.0 | 2.0 |  24.7 | 129.0 | "31,41,51,61,71,81" |   0:00:39 |   0:00:32 |


Table: statMergeAnchors.md

| Name                          | Mapped% | N50Anchor |   Sum |  # | N50Others |     Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:------------------------------|--------:|----------:|------:|---:|----------:|--------:|---:|-------:|----:|------:|------:|----------:|
| 7_merge_anchors               |  98.52% |    198565 | 3.35M | 39 |     30269 | 383.98K | 92 |  256.0 | 4.0 |  81.3 | 402.0 |   0:00:44 |
| 7_merge_mr_unitigs_bcalm      |  99.45% |    198595 | 3.35M | 38 |      1647 |  18.33K | 12 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:25 |
| 7_merge_mr_unitigs_superreads |  99.41% |    198597 | 3.35M | 38 |     17171 |  28.59K | 11 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:19 |
| 7_merge_mr_unitigs_tadpole    |  99.45% |    198595 | 3.35M | 38 |     17853 |  67.26K | 12 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:24 |
| 7_merge_unitigs_bcalm         |  99.36% |    198577 | 3.35M | 39 |     30269 | 163.74K | 65 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:16 |
| 7_merge_unitigs_superreads    |  99.37% |    198595 | 3.35M | 39 |      1452 |  85.47K | 53 |  256.0 | 5.0 |  80.3 | 406.5 |   0:01:16 |
| 7_merge_unitigs_tadpole       |  99.32% |    198597 | 3.35M | 38 |     87392 | 169.91K | 66 |  255.0 | 4.0 |  81.0 | 400.5 |   0:01:10 |


Table: statOtherAnchors.md

| Name         | Mapped% | N50Anchor |   Sum |  # | N50Others |    Sum |  # | median | MAD | lower | upper | RunTimeAN |
|:-------------|--------:|----------:|------:|---:|----------:|-------:|---:|-------:|----:|------:|------:|----------:|
| 8_spades     |  99.32% |    248211 | 2.93M | 23 |      1689 | 16.55K | 41 |  256.0 | 6.0 |  79.3 | 411.0 |   0:00:39 |
| 8_mr_spades  |  99.48% |    274713 | 3.36M | 24 |      1380 |  17.3K | 45 |  451.0 | 6.0 | 144.3 | 703.5 |   0:00:45 |
| 8_megahit    |  99.16% |    202194 | 3.35M | 37 |      1432 | 14.31K | 73 |  256.0 | 4.5 |  80.8 | 404.2 |   0:00:39 |
| 8_mr_megahit |  99.50% |    274967 | 3.36M | 24 |      1380 | 16.07K | 50 |  451.0 | 7.5 | 142.8 | 710.2 |   0:00:44 |
| 8_platanus   |  98.86% |    274612 | 3.35M | 27 |      1072 | 10.07K | 51 |  256.0 | 5.5 |  79.8 | 408.8 |   0:00:41 |


Table: statFinal

| Name                     |     N50 |     Sum |   # |
|:-------------------------|--------:|--------:|----:|
| Genome                   | 3397754 | 3397754 |   1 |
| Paralogs                 |    2793 |  100404 |  43 |
| 7_merge_anchors.anchors  |  198565 | 3349937 |  39 |
| 7_merge_anchors.others   |   30269 |  383981 |  92 |
| glue_anchors             |  198565 | 3348430 |  37 |
| fill_anchors             |  261682 | 3351256 |  30 |
| spades.contig            |  431777 | 3381613 |  59 |
| spades.scaffold          |  431777 | 3381713 |  58 |
| spades.non-contained     |  431777 | 3375259 |  19 |
| mr_spades.contig         |  351394 | 3379512 |  30 |
| mr_spades.scaffold       |  351394 | 3379512 |  30 |
| mr_spades.non-contained  |  351394 | 3376822 |  21 |
| megahit.contig           |  248588 | 3374114 |  52 |
| megahit.non-contained    |  248588 | 3366478 |  36 |
| mr_megahit.contig        |  275034 | 3390001 |  51 |
| mr_megahit.non-contained |  275034 | 3378696 |  26 |
| platanus.contig          |  198751 | 3392430 | 197 |
| platanus.scaffold        |  363088 | 3385815 | 148 |
| platanus.non-contained   |  363088 | 3364323 |  24 |

