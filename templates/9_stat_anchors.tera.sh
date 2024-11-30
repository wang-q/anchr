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
    "Name" "CovCor" "Mapped%" \
    "N50Anchor" "Sum" "#" \
    "N50Others" "Sum" "#" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$//' \
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
		        MAPPED_RATIO=$( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber' )

                printf "%s\t" \
		            "Q${Q}L${L}X${X}P${P}" \
		            $( perl -e "printf qq(%.1f), ${SUM_COR} / {{ opt.genome }};" ) \
                    $( perl -e "printf qq(%.2f%%), ${MAPPED_RATIO} * 100;" ) \
		            $( stat_format anchor/anchor.fasta ) \
		            $( stat_format anchor/pe.others.fa ) \
		            $( cat anchor/env.json | jq '.median | tonumber | map((. * 10 | round) / 10)' ) \
		            $( cat anchor/env.json | jq '.MAD | tonumber | map((. * 10 | round) / 10)' ) \
		            $( cat anchor/env.json | jq '.lower | tonumber | map((. * 10 | round) / 10)' ) \
		            $( cat anchor/env.json | jq '.upper | tonumber | map((. * 10 | round) / 10)' ) |
                    sed 's/\t$//'

		        popd > /dev/null
		    done
	    done
    done
done \
>> tempfile

rgr md tempfile -o ${FILENAME_MD}
cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
