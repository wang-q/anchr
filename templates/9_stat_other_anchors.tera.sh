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

tempfile=$(mktemp /tmp/stat_merge_anchor_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "Mapped" \
    "N50Anchor" "Sum" "#" \
    "SumOthers" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$/\n/' \
    > ${tempfile}

for D in 8_spades 8_mr_spades 8_megahit 8_mr_megahit; do
	if [ ! -e ${D}/anchor/anchor.fasta ]; then
		continue;
	fi

	pushd ${D}/ > /dev/null

    printf "%s\t" \
		$(basename "${D}") \
        $( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber | (. * 1000 | round) / 1000' ) \
		$( stat_format anchor/anchor.fasta ) \
		$( stat_format anchor/pe.others.fa | cut -f 2  ) \
        $( cat anchor/env.json | jq '.median | tonumber | (. * 10 | round) / 10' ) \
        $( cat anchor/env.json | jq '.MAD    | tonumber | (. * 10 | round) / 10' ) \
        $( cat anchor/env.json | jq '.lower  | tonumber | (. * 10 | round) / 10' ) \
        $( cat anchor/env.json | jq '.upper  | tonumber | (. * 10 | round) / 10' ) |
        sed 's/\t$/\n/'

	popd > /dev/null
done \
>> ${tempfile}

rgr md ${tempfile} --right 2-10 -o ${FILENAME_MD}
echo -e "\nTable: ${FILENAME_MD}\n" > ${FILENAME_MD}

cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
