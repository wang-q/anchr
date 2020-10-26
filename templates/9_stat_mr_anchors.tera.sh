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
    "Kmer" "RunTimeKU" "RunTimeAN" \
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
		SECS_KU=$( cat env.json | jq '.RUNTIME | tonumber' )
		SECS_AN=$( expr $(stat -c %Y anchor/anchor.success) - $(stat -c %Y anchor/anchors.sh) )

		printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %.1f | %.1f | %.1f | %.1f | %s | %s | %s |\n" \
			"MRX${X}P${P}" \
			$( perl -e "printf qq(%.1f), ${SUM_COR} / [% opt.genome %];" ) \
			$( perl -e "printf qq(%.2f%%), ${MAPPED_RATIO} * 100;" ) \
			$( stat_format anchor/anchor.fasta ) \
			$( stat_format anchor/pe.others.fa ) \
			$( cat anchor/env.json | jq '.median | tonumber' ) \
			$( cat anchor/env.json | jq '.MAD | tonumber' ) \
			$( cat anchor/env.json | jq '.lower | tonumber' ) \
			$( cat anchor/env.json | jq '.upper | tonumber' ) \
			$( cat env.json | jq '.KMER' ) \
			$( printf "%d:%02d'%02d''\n" $((${SECS_KU}/3600)) $((${SECS_KU}%3600/60)) $((${SECS_KU}%60)) ) \
			$( printf "%d:%02d'%02d''\n" $((${SECS_AN}/3600)) $((${SECS_AN}%3600/60)) $((${SECS_AN}%60)) )

		popd > /dev/null
	done
done \
>> ${FILENAME_MD}

cat ${FILENAME_MD}
