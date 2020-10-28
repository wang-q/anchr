{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 7_merge_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [DIR_MERGE]"

DIR_PREFIX=${1:-"4_unitigs"}
DIR_MERGE=${2:-"7_merge_anchors"}

if [ -e ${DIR_MERGE}/anchor.merge.fasta ]; then
    echo >&2 "${DIR_MERGE}/anchor.merge.fasta presents"
    exit;
fi

#----------------------------#
# merge anchors
#----------------------------#
log_info "anchor.non-contained"

mkdir -p ${DIR_MERGE}

# reverse sorted files, so that Q30L60X80 will be infile_0
dazz contained \
    $( find . -path "*${DIR_PREFIX}*" -name "anchor.fasta" -or -path "*${DIR_PREFIX}*" -name "anchor.merge.fasta" | sort -r ) \
    --len 1000 --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
    -o stdout |
    faops filter -a 1000 -l 0 stdin ${DIR_MERGE}/anchor.non-contained.fasta

{% if opt.redo == "0" -%}
dazz orient \
    ${DIR_MERGE}/anchor.non-contained.fasta \
    --len 1000 --idt 0.98 --parallel {{ opt.parallel }} \
    -o ${DIR_MERGE}/anchor.intermediate_0.fasta
dazz merge \
    ${DIR_MERGE}/anchor.intermediate_0.fasta \
    --len 1000 --idt 0.999 --parallel {{ opt.parallel }} \
    -o ${DIR_MERGE}/anchor.intermediate_1.fasta
dazz contained \
    ${DIR_MERGE}/anchor.intermediate_1.fasta \
    --len 1000 --idt 0.98 --proportion 0.99 --parallel {{ opt.parallel }} \
    -o stdout \
    | faops filter -a 1000 -l 0 stdin ${DIR_MERGE}/anchor.merge.fasta
{% else -%}
#----------------------------#
# anchors with Q0L0 reads
#----------------------------#
log_info "anchors with Q0L0 reads"

mkdir -p ${DIR_MERGE}/anchor
cd ${DIR_MERGE}/anchor

anchr anchors \
    ../anchor.non-contained.fasta \
    ${BASH_DIR}/2_illumina/trim/pe.cor.fa.gz \
    -p {{ opt.parallel }} \
    --keepedge \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

mv anchor.fasta ../anchor.merge.fasta
{% endif -%}
{# Keep a blank line #}

#----------------------------#
# others
#----------------------------#
log_info "others"

cd ${BASH_DIR}

dazz contained \
    $( find . -path "*${DIR_PREFIX}*" -name "pe.others.fa" -or -path "*${DIR_PREFIX}*" -name "others.non-contained.fasta" | sort -r ) \
{% if opt.redo == "1" -%}
    ${DIR_MERGE}/anchor/pe.others.fa \
{% endif -%}
    --len 500 --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
    -o stdout |
    faops filter -a 500 -l 0 stdin ${DIR_MERGE}/others.intermediate_0.fasta

dazz contained \
    ${DIR_MERGE}/anchor.merge.fasta \
    ${DIR_MERGE}/others.intermediate_0.fasta \
    --len 500 --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
    -o stdout \
    | faops filter -a 500 -l 0 stdin ${DIR_MERGE}/others.intermediate_1.fasta

cat ${DIR_MERGE}/others.intermediate_1.fasta |
    grep '>infile_1/' |
    sed 's/>//' \
    > ${DIR_MERGE}/others.txt

faops some -l 0 \
    ${DIR_MERGE}/others.intermediate_1.fasta \
    ${DIR_MERGE}/others.txt \
    ${DIR_MERGE}/others.non-contained.fasta

{% if opt.redo == "1" -%}
find ${DIR_MERGE}/anchor -name "*.fasta" -or -name "*.fa" | parallel --no-run-if-empty -j 1 rm
{% endif -%}
find ${DIR_MERGE} -name "anchor.intermediate*" | parallel --no-run-if-empty -j 1 rm
find ${DIR_MERGE} -name "others.intermediate*" | parallel --no-run-if-empty -j 1 rm

log_info Done.

exit 0
