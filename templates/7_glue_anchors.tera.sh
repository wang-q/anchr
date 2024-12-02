{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
log_warn 7_glue_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 FILE_ANCHOR FILE_LONG GAP_COV"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

FILE_ANCHOR=$1
FILE_LONG=$2
GAP_COV=${3:-3}

if [ -e 7_glue_anchors/contig.fasta ]; then
    echo >&2 "7_glue_anchors/contig.fasta presents"
    exit;
fi

#----------------------------#
# glue anchors
#----------------------------#
mkdir -p 7_glue_anchors

log_info "overlap: between anchor-long"

anchr overlap2 \
    --parallel {{ opt.parallel }} \
    ${FILE_ANCHOR} \
    ${FILE_LONG} \
    -d 7_glue_anchors \
    -b 50 --len 1000 --idt 0.999 --all

cd 7_glue_anchors

log_info "overlap: within anhcors"
anchr overlap \
    anchor.fasta \
    --serial --len {{ opt.gluemin }} --idt 0.9999 \
    -o stdout |
    perl -nla -e '
        BEGIN {
            our %seen;
            our %count_of;
        }

        @F == 13 or next;
        $F[3] > 0.9999 or next;

        my $pair = join( "-", sort { $a <=> $b } ( $F[0], $F[1], ) );
        next if $seen{$pair};
        $seen{$pair} = $_;

        $count_of{ $F[0] }++;
        $count_of{ $F[1] }++;

        END {
            for my $pair ( keys %seen ) {
                my ($f_id, $g_id) = split "-", $pair;
                next if $count_of{$f_id} > 2;
                next if $count_of{$g_id} > 2;
                print $seen{$pair};
            }
        }
    ' |
    sort -k 1n,1n -k 2n,2n \
    > anchor.ovlp.tsv

log_info "group: anchor-long"
rm -fr group
dazz group \
    anchorLong.db \
    anchorLong.ovlp.tsv \
    --oa anchor.ovlp.tsv \
    --parallel {{ opt.parallel }} \
    --range "1-$(hnsm n50 -H -N 0 -C anchor.fasta)" \
    --len 1000 --idt 0.999 --max "-{{ opt.gluemin }}" -c ${GAP_COV}

log_info "Processing each groups"
{% set parallel2 = opt.parallel | int / 2 -%}
{% set parallel2 = parallel2 | round(method="floor") -%}
{% if parallel2 < 2 %}{% set parallel2 = 2 %}{% endif -%}
cat group/groups.txt |
    parallel --no-run-if-empty --linebuffer -k -j {{ parallel2 }} '
        echo {};
        anchr orient \
            --len 1000 --idt 0.999 \
            group/{}.anchor.fasta \
            group/{}.long.fasta \
            -r group/{}.restrict.tsv \
            -o group/{}.strand.fasta;

        anchr overlap --len 1000 --idt 0.9999 \
            group/{}.strand.fasta \
            -o stdout |
            anchr restrict \
                stdin group/{}.restrict.tsv \
                -o group/{}.ovlp.tsv;

        anchr overlap --len {{ opt.gluemin }} --idt 0.9999 \
            group/{}.strand.fasta \
            -o stdout |
            perl -nla -e '\''
                @F == 13 or next;
                $F[3] > 0.9999 or next;
                $F[9] == 0 or next;
                $F[5] > 0 and $F[6] == $F[7] or next;
                /anchor.+anchor/ or next;
                print;
            '\'' \
            > group/{}.anchor.ovlp.tsv

        dazz layout \
            group/{}.strand.fasta \
            group/{}.ovlp.tsv \
            group/{}.relation.tsv \
            --oa group/{}.anchor.ovlp.tsv \
            -o group/{}.contig.fasta
    '

log_info "Build contigs"
cat \
   group/non_grouped.fasta \
   group/*.contig.fasta |
   hnsm filter -a 1000 stdin -o contig.fasta

log_info Done.

exit 0
