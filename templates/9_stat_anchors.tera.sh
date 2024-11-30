{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 9_stat_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [FILENAME_MD]"

DIR_PREFIX=${1:-"4_unitigs_superreads"}
FILENAME_MD=${2:-"statAnchors.md"}

tempfile=$(mktemp)

echo -e "Table: ${FILENAME_MD}\n" > ${FILENAME_MD}

printf "%s\t" \
    "Name" "CovCor" "Mapped" \
    "N50Anchor" "Sum" "#" \
    "SumOthers" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$/\n/' \
    > tempfile

for Q in 0 {{ opt.qual }}; do
    for L in 0 {{ opt.len }}; do
	    for X in {{ opt.cov }}; do
		    for P in $(printf "%03d " {0..{{ opt.statp }}}); do
		        if [ ! -e ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}/anchor/anchor.fasta ]; then
			        continue;
			    fi

		        pushd ${DIR_PREFIX}/Q${Q}L${L}X${X}P${P}/ > /dev/null

		        SUM_COR=$( cat env.json | jq '.SUM_COR | tonumber' )

                printf "%s\t" \
		            "Q${Q}L${L}X${X}P${P}" \
		            $( perl -e "printf qq(%.1f), ${SUM_COR} / {{ opt.genome }};" ) \
		            $( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber | (. * 1000 | round) / 1000' ) \
		            $( stat_format anchor/anchor.fasta ) \
		            $( stat_format anchor/pe.others.fa | head -n 2 | tail -n 1) \
		            $( cat anchor/env.json | jq '.median | tonumber | (. * 10 | round) / 10' ) \
		            $( cat anchor/env.json | jq '.MAD    | tonumber | (. * 10 | round) / 10' ) \
		            $( cat anchor/env.json | jq '.lower  | tonumber | (. * 10 | round) / 10' ) \
		            $( cat anchor/env.json | jq '.upper  | tonumber | (. * 10 | round) / 10' ) |
                    sed 's/\t$/\n/'

		        popd > /dev/null
		    done
	    done
    done
done \
>> tempfile

rgr md tempfile --right 2-13 -o ${FILENAME_MD}
cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
