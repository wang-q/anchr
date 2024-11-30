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

tempfile=$(mktemp /tmp/stat_mr_anchor_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "CovCor" "Mapped" \
    "N50Anchor" "Sum" "#" \
    "SumOthers" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$/\n/' \
    > ${tempfile}

for X in {{ opt.cov }}; do
	for P in $(printf "%03d " {0..{{ opt.statp }}}); do
		if [ ! -e ${DIR_PREFIX}/MRX${X}P${P}/anchor/anchor.fasta ]; then
			continue;
		fi

		pushd ${DIR_PREFIX}/MRX${X}P${P}/ > /dev/null

		SUM_COR=$( cat env.json | jq '.SUM_COR | tonumber' )

        printf "%s\t" \
			"MRX${X}P${P}" \
			$( perl -e "printf qq(%.1f), ${SUM_COR} / {{ opt.genome }};" ) \
            $( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber | (. * 1000 | round) / 1000' ) \
			$( stat_format anchor/anchor.fasta ) \
            $( stat_format anchor/pe.others.fa | cut -f 2 ) \
            $( cat anchor/env.json | jq '.median | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.MAD    | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.lower  | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.upper  | tonumber | (. * 10 | round) / 10' ) |
            sed 's/\t$/\n/'

		popd > /dev/null
	done
done \
>> ${tempfile}

rgr md ${tempfile} --right 2-11 -o ${FILENAME_MD}
echo -e "\nTable: ${FILENAME_MD}\n" >> ${FILENAME_MD}

cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
