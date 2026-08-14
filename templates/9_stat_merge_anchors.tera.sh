{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_merge_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [FILENAME_MD]"

DIR_PREFIX=${1:-"7_merge"}
FILENAME_MD=${2:-"statMergeAnchors.md"}

tempfile=$(mktemp /tmp/stat_merge_anchor_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "Mapped" \
    "N50Anchor" "Sum" "#" \
    "SumOthers" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$/\n/' \
    > ${tempfile}

for D in $( find . -type d -name "${DIR_PREFIX}*" | sort ); do
	if [ ! -e ${D}/anchor.merge.fasta ]; then
		continue;
	fi

	pushd ${D}/ > /dev/null

    # merged anchors 的覆盖统计（老流程 anchors 化等价物）
    anchr asm anchor \
        anchor.merge.fasta \
        ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz \
        --mincov 5 --mscale 3 \
        --lscale {{ opt.lscale }} \
        --uscale {{ opt.uscale }} \
        -p {{ opt.parallel }} \
        --stats merge.stats.tsv \
        -o /dev/null

    printf "%s\t" \
        $(basename "${D}") \
        $( cut -f1 merge.stats.tsv ) \
        $( stat_format anchor.merge.fasta ) \
        $( cut -f6 merge.stats.tsv ) \
        $( cut -f2 merge.stats.tsv ) \
        $( cut -f3 merge.stats.tsv ) \
        $( cut -f4 merge.stats.tsv ) \
        $( cut -f5 merge.stats.tsv ) |
        sed 's/\t$/\n/'

    rm -f merge.stats.tsv

	popd > /dev/null
done \
>> ${tempfile}

tva to md ${tempfile} --right 2-10 -o ${FILENAME_MD}
echo -e "\nTable: ${FILENAME_MD}\n" >> ${FILENAME_MD}

cat ${FILENAME_MD}
mkdir -p ${BASH_DIR}/../9_markdown
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
