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
    "Name" "N50" "Sum" "#" |
    sed 's/\t$/\n/' \
    > ${tempfile}

for D in 8_spades 8_mr_spades 8_megahit 8_mr_megahit; do
    if [ ! -e ${D}/anchor.fasta ]; then
        continue;
    fi

    printf "%s\t%s\n" \
        $(basename "${D}") \
        $( stat_format ${D}/anchor.fasta ) \
        >> ${tempfile}
done

rgr md ${tempfile} --right 2-4 -o ${FILENAME_MD}
echo -e "\nTable: ${FILENAME_MD}\n" >> ${FILENAME_MD}

cat ${FILENAME_MD}
mkdir -p ${BASH_DIR}/../9_markdown
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown
