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

echo -e "Table: ${FILENAME_MD}\n" > ${FILENAME_MD}
printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" \
    "Name" "Mapped%" \
    "N50Anchor" "Sum" "#" \
    "N50Others" "Sum" "#" \
    "median" "MAD" "lower" "upper" \
    "RunTimeAN" \
    >> ${FILENAME_MD}
printf "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n" \
    >> ${FILENAME_MD}

for D in $( find . -type d -name "${DIR_PREFIX}*" | sort ); do
	if [ ! -e ${D}/anchor.merge.fasta ]; then
		continue;
	fi

	pushd ${D}/ > /dev/null

	MAPPED_RATIO=$( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber' )
    SECS_AN=$( cat anchor/env.json | jq '.RUNTIME | tonumber' )

	printf "| %s | %s | %s | %s | %s | %s | %s | %s | %.1f | %.1f | %.1f | %.1f | %s |\n" \
		$(basename "${D}") \
		$( perl -e "printf qq(%.2f%%), ${MAPPED_RATIO} * 100;" ) \
		$( stat_format anchor.merge.fasta ) \
		$( stat_format others.non-contained.fasta ) \
		$( cat anchor/env.json | jq '.median | tonumber' ) \
		$( cat anchor/env.json | jq '.MAD | tonumber' ) \
		$( cat anchor/env.json | jq '.lower | tonumber' ) \
		$( cat anchor/env.json | jq '.upper | tonumber' ) \
        $( time_format ${SECS_AN} )

	popd > /dev/null
done \
>> ${FILENAME_MD}

cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/9_markdown
