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

    printf "%s\t" \
        $(basename "${D}") \
        $( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber | (. * 1000 | round) / 1000' ) \
        $( stat_format anchor.merge.fasta ) \
        $( stat_format others.non-contained.fasta | cut -f 2 ) \
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
