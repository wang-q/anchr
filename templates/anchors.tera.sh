{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Prepare UT
#----------------------------#
START_TIME=$(date +%s)

log_info Symlink input files

if [ ! -e UT.fasta ]; then
    ln -s {{ args.0 }} UT.fasta
fi

log_debug "UT sizes"
faops size UT.fasta > ut.chr.sizes
spanr genome ut.chr.sizes -o ut.chr.json

#----------------------------#
# Mapping reads
#----------------------------#
log_info "Mapping reads"

log_debug "bbmap via bbwrap"
rm -f *.sam
rm -f basecov.txt
rm -f basecov.hist.tsv
bbwrap.sh \
    maxindel=0 strictmaxindel perfectmode \
    threads={{ opt.parallel }} \
    ambiguous=all \
    nodisk append \
    ref=UT.fasta \
    in={{ args | slice(start=1) | join(sep=",") }} \
    outm=mapped.sam outu=unmapped.sam \
    basecov=basecov.txt \
    1>bbmap.err 2>&1

COUNT_MAPPED=$( cat mapped.sam | wc -l | sed 's/ //g' )
COUNT_UNMAPPED=$( cat unmapped.sam | wc -l | sed 's/ //g' )
MAPPED_RATIO=$(
    echo "print ${COUNT_MAPPED}/(${COUNT_MAPPED} + ${COUNT_UNMAPPED})" | perl
)
save MAPPED_RATIO
find . -type f -name "mapped.sam"   | parallel --no-run-if-empty -j 1 rm
find . -type f -name "unmapped.sam" | parallel --no-run-if-empty -j 1 rm

#----------------------------#
# basecov
#----------------------------#
log_info "basecov"

# How to best eliminate values in a list that are outliers
# http://exploringdatablog.blogspot.com/2013/02/finding-outliers-in-numerical-data.html
#
# basecov.txt
# Pos is 0-based
#RefName	Pos	Coverage
cat basecov.txt |
    grep -v '^#' |
    tsv-filter --ne 3:0 | # Non-covered regions should be ignored
    tsv-summarize --median 3 --mad 3 --quantile 3:0.25,0.75 |
    perl -MPath::Tiny -MJSON::PP -e '
        my $json = JSON::PP::decode_json( Path::Tiny::path( q(env.json) )->slurp );

        my $line = <>;
        my @fields = split qq(\t), $line;

        $json->{median} = $fields[0];
        $json->{MAD} = $fields[1];
        $json->{IQR} = $fields[3] - $fields[2];

        my $lower = ( $json->{median} - {{ opt.mscale }} * $json->{MAD} ) / {{ opt.lscale }};
        $lower = {{ opt.mincov }} if $lower < {{ opt.mincov }};
        my $upper = ( $json->{median} + {{ opt.mscale }} * $json->{MAD} ) * {{ opt.uscale }};

        $json->{mscale} = {{ opt.mscale }};
        $json->{lscale} = {{ opt.lscale }};
        $json->{uscale} = {{ opt.uscale }};
        $json->{lower} = $lower;
        $json->{upper} = $upper;

        Path::Tiny::path(q(env.json))->spew( JSON::PP::encode_json($json) );
    '

#----------------------------#
# Properly covered regions by reads
#----------------------------#
# at least some reads covered
log_debug "covered"
cat basecov.txt |
    grep -v '^#' |
    tsv-filter --ne 3:0 | # Non-covered regions should be ignored
    perl -nla -MPath::Tiny -MJSON::PP -e '
        BEGIN {
            our $name;
            our @list;
            our @edges;
            our $limit = JSON::PP::decode_json(
                Path::Tiny::path( q(env.json) )->slurp
            );
            our $length_of;
            for my $line ( Path::Tiny::path( q(ut.chr.sizes) )->lines({ chomp => 1 }) ) {
                my ( $key, $value ) = split /\t/, $line;
                $length_of{$key} = $value;
            }
        }

        sub list_to_ranges {
            my @ranges;
            my @list = sort { $a <=> $b } @_;
            my $count = scalar @list;
            my $pos   = 0;
            while ( $pos < $count ) {
                my $end = $pos + 1;
                $end++ while $end < $count && $list[$end] <= $list[ $end - 1 ] + 1;
                push @ranges, ( $list[$pos], $list[ $end - 1 ] );
                $pos = $end;
            }

            return @ranges;
        }

        if ( !defined $name ) {
            $name = $F[0];
            @list = ();
            @edges = ();
        }

{% if opt.keepedge == "1" -%}
        # proportionally decreases the limits
        if ( $F[1] < {{ opt.readl }} / 2 ) { # left edges
            if ( $F[2] >= {{ opt.mincov }} ) {
                my $lower = $limit->{lower} * $F[1] * 2 / {{ opt.readl }};
                my $upper = $limit->{upper} * $F[1] * 2 / {{ opt.readl }};
                if ( $F[2] >= $lower and $F[2] <= $upper ) {
                    push @edges, $F[1];
                }
            }
        }
        elsif ( $F[1] >= $length_of->{$name} - {{ opt.readl }} / 2 ) { # right edges
            if ( $F[2] >= {{ opt.mincov }} ) {
                my $lower = $limit->{lower} * ($length_of->{$name} - $F[1]) * 2 / {{ opt.readl }};
                my $upper = $limit->{upper} * ($length_of->{$name} - $F[1]) * 2 / {{ opt.readl }};
                if ( $F[2] >= $lower and $F[2] <= $upper ) {
                    push @edges, $F[1];
                }
            }
        }
{% endif -%}

        if ( $F[2] < $limit->{lower} or $F[2] > $limit->{upper} ) {
            next;
        }

        if ( $name eq $F[0] ) {
            push @list, $F[1];
        }
        else {
            my @ranges = list_to_ranges(@list{% if opt.keepedge == "1" %}, @edges{% endif %});
            for ( my $i = 0; $i < $#ranges; $i += 2 ) {
                if ( $ranges[$i] == $ranges[ $i + 1 ] ) {
                    printf qq(%s:%s\n), $name, $ranges[$i] + 1;
                }
                else {
                    printf qq(%s:%s-%s\n), $name, $ranges[$i] + 1, $ranges[ $i + 1 ] + 1;
                }
            }

            $name = $F[0];
            @list = ( $F[1] );
{% if opt.keepedge == "1" -%}
            if ( $F[2] >= {{ opt.mincov }} ) {
                my $lower = $limit->{lower} * $F[1] * 2 / {{ opt.readl }};
                my $upper = $limit->{upper} * $F[1] * 2 / {{ opt.readl }};
                if ( $F[2] >= $lower and $F[2] <= $upper ) {
                    @edges = ( $F[1] );
                }
                else {
                    @edges = ();
                }
            }
            else {
                @edges = ();
            }
{% endif -%}

        }

        END {
            my @ranges = list_to_ranges(@list{% if opt.keepedge == "1" %}, @edges{% endif %});
            for ( my $i = 0; $i < $#ranges; $i += 2 ) {
                if ( $ranges[$i] == $ranges[ $i + 1 ] ) {
                    printf qq(%s:%s\n), $name, $ranges[$i] + 1;
                }
                else {
                    printf qq(%s:%s-%s\n), $name, $ranges[$i] + 1, $ranges[ $i + 1 ] + 1;
                }
            }
        }
    ' \
    > contig.covered.txt
find . -type f -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm

#----------------------------#
# anchor
#----------------------------#
log_info "anchor - proper covered regions"

OPT_FILL={{ opt.fill }}
OPT_MIN={{ opt.min }}
save OPT_FILL
save OPT_MIN

# covered region
spanr cover contig.covered.txt -o contig.covered.json
spanr stat ut.chr.sizes contig.covered.json -o contig.covered.csv

# fill all holes of {{ opt.ratio }} covered contigs
cat contig.covered.csv |
    perl -nla -F"," -e '
        $F[0] eq q{chr} and next;
        $F[0] eq q{all} and next;
        $F[2] < {{ opt.min }} and next;
        $F[3] < {{ opt.ratio }} and next;
        print $F[0];
    ' |
    sort -n \
    > fill_all.txt

spanr some contig.covered.json fill_all.txt -o contig.fill_all.json
rm -f contig.fill_all.temp.json; ln -s contig.fill_all.json contig.fill_all.temp.json

if [ -s fill_all.txt ]; then
    spanr span contig.fill_all.temp.json --op fill -n {{ opt.fill | int * 10 }} -o contig.fill_all.1.json
    rm contig.fill_all.temp.json; ln -s contig.fill_all.1.json contig.fill_all.temp.json

    spanr span contig.fill_all.temp.json --op excise -n {{ opt.min }} -o contig.fill_all.2.json
    rm contig.fill_all.temp.json; ln -s contig.fill_all.2.json contig.fill_all.temp.json
fi

# fill small holes
cat ut.chr.sizes |
    cut -f 1 |
    grep -Fx -f fill_all.txt -v \
    > fill_hole.txt

spanr some contig.covered.json fill_hole.txt -o contig.fill_hole.json
rm -f contig.fill_hole.temp.json; ln -s contig.fill_hole.json contig.fill_hole.temp.json

if [ -s fill_hole.txt ]; then
    spanr span contig.fill_hole.json --op fill -n {{ opt.fill }} -o contig.fill_hole.1.json
    rm contig.fill_hole.temp.json; ln -s contig.fill_hole.1.json contig.fill_hole.temp.json

    spanr span contig.fill_hole.1.json --op excise -n {{ opt.min }} -o contig.fill_hole.2.json
    rm contig.fill_hole.temp.json; ln -s contig.fill_hole.2.json contig.fill_hole.temp.json
fi

# get proper regions
spanr compare --op union contig.fill_all.temp.json contig.fill_hole.temp.json -o contig.proper.json
spanr stat ut.chr.sizes contig.proper.json -o contig.proper.csv

{% if opt.longest == "1" -%}
spanr convert contig.proper.json --longest |
    spanr cover stdin -o anchor.json
{% else -%}
ln -s contig.proper.json anchor.json
{% endif -%}
spanr convert anchor.json -o anchor.regions.txt
{# Keep a blank line #}

#----------------------------#
# others
#----------------------------#
spanr compare ut.chr.json anchor.json --op diff -o others.json
spanr convert others.json -o others.regions.txt

#----------------------------#
# Split UT.fasta to anchor and others
#----------------------------#
log_info "pe.anchor.fa & pe.others.fa"

hnsm range UT.fasta -r anchor.regions.txt -o pe.anchor.fa
hnsm range UT.fasta -r others.regions.txt -o pe.others.fa

#----------------------------#
# Merging anchors
#----------------------------#
log_info "Merging anchors"
anchr contained \
    pe.anchor.fa \
    --len {{ opt.min }} --idt 0.9999 --ratio 0.99999 --parallel {{ opt.parallel }} \
    -o anchor.non-contained.fasta
anchr orient \
    anchor.non-contained.fasta \
    --len {{ opt.min }} --idt 0.999 --parallel {{ opt.parallel }} \
    -o anchor.orient.fasta
anchr merge \
    anchor.orient.fasta --len {{ opt.min }} --idt 0.9999 --parallel {{ opt.parallel }} \
    -o anchor.merge0.fasta

# loss idt to remove duplicates
anchr contained \
    anchor.merge0.fasta \
    --len {{ opt.min }} --idt 0.98 --ratio 0.99 --parallel {{ opt.parallel }} \
    -o anchor.fasta

#----------------------------#
# Done.
#----------------------------#
rm -f "pe.anchor.fa"
rm -f "anchor.*.fasta"
rm -f "contig.fill*.json"

save START_TIME

END_TIME=$(date +%s)
save END_TIME

RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info "Done."

exit 0
