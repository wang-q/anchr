{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_other_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [FILENAME_MD]"

FILENAME_MD=${1:-"statOtherAnchors.md"}

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

for D in 8_spades 8_mr_spades 8_megahit 8_mr_megahit 8_platanus; do
	if [ ! -e ${D}/anchor/anchor.fasta ]; then
		continue;
	fi

	pushd ${D}/ > /dev/null

	MAPPED_RATIO=$( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber' )
    SECS_AN=$( cat anchor/env.json | jq '.RUNTIME | tonumber' )

	printf "| %s | %s | %s | %s | %s | %s | %s | %s | %.1f | %.1f | %.1f | %.1f | %s |\n" \
		$(basename "${D}") \
		$( perl -e "printf qq(%.2f%%), ${MAPPED_RATIO} * 100;" ) \
		$( stat_format anchor/anchor.fasta ) \
		$( stat_format anchor/pe.others.fa ) \
		$( cat anchor/env.json | jq '.median | tonumber' ) \
		$( cat anchor/env.json | jq '.MAD | tonumber' ) \
		$( cat anchor/env.json | jq '.lower | tonumber' ) \
		$( cat anchor/env.json | jq '.upper | tonumber' ) \
        $( time_format ${SECS_AN} )

	popd > /dev/null
done \
>> ${FILENAME_MD}

cat ${FILENAME_MD}
