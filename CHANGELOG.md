# Change Log

## Unreleased - ReleaseDate

* Add `anchr dazzname`
* Move `ovlpr` from intspan here
* Use `mlr` instead of manual `printf` to create markdown files
* `anchr template`
    * Place scripts to `0_script/`
    * Place markdowns to `9_markdown/`

## 0.3.18 - 2022-04-21

* Update Github actions
    * Use a container with GLIBC 2.17 to build linux-gnu binary
    * Codecov with cargo-tarpaulin
* Simplify `quorum.sh`
* Add `anchr template` scripts
    * 2_fastk.sh
* Add bifrost as a unitigger

## 0.3.17 - 2022-04-15

* Update check_dep.sh and install_dep.sh
* Add fastk.md. It will replace KAT in the near future.
    * FastK
    * Merqury.FK
    * GeneScope.FK

## 0.3.16 - 2021-01-12

* Add `anchr template` scripts
    * 3_bwa.sh
    * 3_gatk.sh
    * Add dependency on bwa, samtools, and gatk
* Add --ascp to ena_prep.pl

## 0.3.15 - 2020-12-15

* bbtools 37.77
* sort keys of env.json
* `anchr anchors`
    * Use tsv-summarize to replace Perl codes
    * Rewrite --keepedge
    * Non-covered regions should be ignored
* Set `dazz --idt` to a higher level
    * `dazz contained` - 0.9999
    * `dazz orient` - 0.999
    * `dazz merge` - 0.9999
    * `dazz group` - 0.999
* Adjust `anchr template` default options
    * --gluemin 30
    * --fillmax 100
* Add `anchr template` scripts
    * 9_busco.sh
* Perl version of `readlink -f`

## 0.3.6 - 2020-11-18

* `anchr anchors`
    * Proportional read coverages near edges of contigs
    * Rename opt `--scale` to `--mscale`
    * Add opts `--readl`, `--lscale`, and `--uscale`
    * Remove opt `--lower`
* `anchr template`
    * Add opts `--readl`, `--lscale`, and `--uscale`

## 0.3.1 - 2020-11-16

* Adjust `anchr template` default options
    * --gluemin 10
    * --fillmax 500
* Support SAMN and PRJNA in `anchr ena`
* Add `repetitives.fa` to `1_genome`
* Add docs
    * fda_argos.md
    * yeast.md

## 0.3.0 - 2020-11-08

* Add `anchr template` scripts
    * 0_bsub.sh
* Features completed

## 0.2.8 - 2020-11-07

* Update `anchr dep` scripts
* Support multiple unitiggers
    * superreads, tadpole, or bcalm
* Prebuilt resources
* Add docs
    * se.md
    * qlx.md
    * gage_b.md

## 0.2.6 - 2020-10-29

* Add `anchr template` scripts
    * 7_glue_anchors.sh
    * 7_fill_anchors.sh
* Remove `anchr template` scripts
    * 2_kmergenie.sh
    * 2_sga_preqc.sh
* Add `kat hist` to 2_kat.sh
* Fix bugs

## 0.2.4 - 2020-10-28

* Save RUNTIME in `anchr anchors`
* Add func time_format()

* Add `anchr template` scripts
    * 2_insert_size.sh
    * 2_kat.sh
    * 2_sga_preqc.sh
    * 8_spades.sh
    * 8_mr_spades.sh
    * 8_megahit.sh
    * 8_mr_megahit.sh
    * 8_platanus.sh
    * 9_stat_other_anchors.sh
    * 9_quast.sh
    * 9_stat_final.sh

## 0.2.3 - 2020-10-27

* Remove dependency on samtools
* Fix bugs

* Add `anchr template` scripts
    * 7_merge_anchors.sh
    * 9_stat_merge_anchors.sh

## 0.2.2 - 2020-10-27

* Avoid camelCases

* Add `anchr template` scripts
    * 9_stat_anchors.sh
    * 9_stat_mr_anchors.sh

## 0.2.1 - 2020-10-26

* Rename `environment.json` to `env.json`

* Add `anchr anchors`

* Add `anchr template` scripts
    * 4_unitigs.sh
    * 6_unitigs.sh
    * 4_anchors.sh
    * 6_anchors.sh

## 0.2.0 - 2020-10-26

* Binary name `anchr`

* Add `anchr template` scripts
    * 0_cleanup.sh
    * 0_real_clean.sh
    * 0_master.sh
    * 4_down_sampling.sh
    * 6_down_sampling.sh

* Add `anchr unitigs`

## 0.1.5 - 2020-10-24

* Add `anchr ena`
* Add `anchr dep`
* Add `anchr trim`
* Add `anchr quorum`
* Add `anchr merge`

* Add `anchr template` scripts
    * 2_fastqc.sh
    * 2_kmergenie.sh
    * 2_trim.sh
    * 9_statReads.sh
    * 2_quorum.sh
    * 2_merge.sh

* Tests of `anchr template`

* Use `tsv-sample` to replace `shuf`
* Fix quotes

## 0.0.1 - 2020-10-22

* Skeletons, need to be filled
* Setup github actions

