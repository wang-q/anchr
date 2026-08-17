{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 7_merge_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [DIR_MERGE] [MIN_CONTIG_LEN]"

DIR_PREFIX=${1:-"4_unitigs_multik"}
DIR_MERGE=${2:-"7_merge_anchors"}
MIN_CL=${3:-1000}

if [ -e ${DIR_MERGE}/anchor.merge.fasta ]; then
    echo >&2 "${DIR_MERGE}/anchor.merge.fasta presents"
    exit;
fi

#----------------------------#
# merge anchors（现代：asm olc --unitigs）
#----------------------------#
log_info "merge anchors with anchr asm olc --unitigs"

mkdir -p ${DIR_MERGE}

# reversely sorted files, so that Q30L60X80 will be first
find . -path "*${DIR_PREFIX}*" \
    \( -name "anchor.fasta" -o -name "anchor.merge.fasta" \) |
    sort -r \
    > ${DIR_MERGE}/anchors.list

anchr asm olc --unitigs --cross-validate --list-files ${DIR_MERGE}/anchors.list \
    --min-overlap 1000 --min-contig-len ${MIN_CL} \
    -o ${DIR_MERGE}/anchor.merge.fasta

log_info Done.

exit 0
