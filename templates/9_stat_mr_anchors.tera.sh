{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_mr_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [FILENAME_MD]"

DIR_PREFIX=${1:-"6_unitigs"}
FILENAME_MD=${2:-"statMRAnchors.md"}

echo -e "Table: ${FILENAME_MD}\n" > ${FILENAME_MD}
printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" \
    "Name" "CovCor" "Mapped%" \
    "N50Anchor" "Sum" "#" \
    "N50Others" "Sum" "#" \
    "median" "MAD" "lower" "upper" \
    "Kmer" "RunTimeUT" "RunTimeAN" \
    >> ${FILENAME_MD}
printf "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n" \
    >> ${FILENAME_MD}

for X in {{ opt.cov }}; do
	for P in $(printf "%03d " {0..{{ opt.statp }}}); do
		if [ ! -e ${DIR_PREFIX}/MRX${X}P${P}/anchor/anchor.fasta ]; then
			continue;
		fi

		pushd ${DIR_PREFIX}/MRX${X}P${P}/ > /dev/null

		SUM_COR=$( cat env.json | jq '.SUM_COR | tonumber' )
		MAPPED_RATIO=$( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber' )
		SECS_UT=$( cat env.json | jq '.RUNTIME | tonumber' )
        SECS_AN=$( cat anchor/env.json | jq '.RUNTIME | tonumber' )

		printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %.1f | %.1f | %.1f | %.1f | %s | %s | %s |\n" \
			"MRX${X}P${P}" \
			$( perl -e "printf qq(%.1f), ${SUM_COR} / {{ opt.genome }};" ) \
			$( perl -e "printf qq(%.2f%%), ${MAPPED_RATIO} * 100;" ) \
			$( stat_format anchor/anchor.fasta ) \
			$( stat_format anchor/pe.others.fa ) \
			$( cat anchor/env.json | jq '.median | tonumber' ) \
			$( cat anchor/env.json | jq '.MAD | tonumber' ) \
			$( cat anchor/env.json | jq '.lower | tonumber' ) \
			$( cat anchor/env.json | jq '.upper | tonumber' ) \
			$( cat env.json | jq '.KMER' ) \
            $( time_format ${SECS_UT} ) \
            $( time_format ${SECS_AN} )

		popd > /dev/null
	done
done \
>> ${FILENAME_MD}

cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
