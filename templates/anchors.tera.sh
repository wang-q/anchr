{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Prepare SR
#----------------------------#
START_TIME=$(date +%s)

log_info Symlink input files

if [ ! -e SR.fasta ]; then
    ln -s {{ args.0 }} SR.fasta
fi

log_debug "SR sizes"
faops size SR.fasta > sr.chr.sizes
spanr genome sr.chr.sizes -o sr.chr.yml

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
    ref=SR.fasta \
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
cat basecov.txt |
    grep -v '^#' |
    perl -nla -MPath::Tiny -MJSON -MApp::Fasops::Common -MApp::Dazz::Common -e '
        BEGIN {
            our $name;
            our @list;
            our %hist_of;
            our $json = JSON->new->decode( Path::Tiny::path( q{env.json} )->slurp );
        }

        if ( !defined $name ) {
            $name = $F[0];
            @list = ( $F[2] );
        }
        elsif ( $name eq $F[0] ) {
            push @list, $F[2];
        }
        else {
            my $mean_cov = App::Fasops::Common::mean(@list);
            printf qq(%s\t%d\n), $name, int $mean_cov;

            $name = $F[0];
            @list = ( $F[2] );
        }

        $hist_of{$F[2]}++;

        END {
            my $mean_cov = App::Fasops::Common::mean(@list);
            printf qq(%s\t%d\n), $name, int $mean_cov;

            my $content;
            for my $key ( sort {$a <=> $b} keys %hist_of ) {
                $content .= sprintf qq(%s\t%d\n), $key, $hist_of{$key};
            }
            Path::Tiny::path(q{basecov.hist.tsv})->spew( $content );

            # Non-covered regions should be ignored
            delete $hist_of{0};

            my $median = App::Dazz::Common::histogram_percentile(\%hist_of, 0.5);
            my $q25 = App::Dazz::Common::histogram_percentile(\%hist_of, 0.25);
            my $q75 = App::Dazz::Common::histogram_percentile(\%hist_of, 0.75);
            my $IQR = $q75 - $q25;

            $json->{median} = $median;
            $json->{IQR} = $IQR;
            Path::Tiny::path(q{env.json})->spew( JSON->new->encode($json) );
        }
    ' \
    > contigs.coverage.tsv

# How to best eliminate values in a list that are outliers
# http://www.perlmonks.org/?node_id=1147296
# http://exploringdatablog.blogspot.com/2013/02/finding-outliers-in-numerical-data.html
cat contigs.coverage.tsv |
    perl -nla -MPath::Tiny -MJSON -MStatistics::Descriptive -e '
        BEGIN {
            our $stat   = Statistics::Descriptive::Full->new();
            our %cov_of = ();
            our $json = JSON->new->decode( Path::Tiny::path( q{env.json} )->slurp );
        }

        next if $F[1] < {{ opt.mincov }};

        $cov_of{ $F[0] } = $F[1];
        $stat->add_data( $F[1] );

        END {
            my $contig_median = $stat->median();

            my $median = $json->{median};
            my @abs_res      = map { abs( $median - $_ ) } $stat->get_data();
            my $abs_res_stat = Statistics::Descriptive::Full->new();
            $abs_res_stat->add_data(@abs_res);
            my $MAD = $abs_res_stat->median();

            my $lower = ( $median - {{ opt.mscale }} * $MAD ) / {{ opt.lscale }};
            $lower = {{ opt.mincov }} if $lower < {{ opt.mincov }};
            my $upper = ( $median + {{ opt.mscale }} * $MAD ) * {{ opt.uscale }};

            $json->{contig_median} = $contig_median;
            $json->{MAD} = $MAD;
            $json->{mscale} = {{ opt.mscale }};
            $json->{lscale} = {{ opt.lscale }};
            $json->{uscale} = {{ opt.uscale }};
            $json->{lower} = $lower;
            $json->{upper} = $upper;

            Path::Tiny::path(q{env.json})->spew( JSON->new->encode($json) );
        }
    '

#----------------------------#
# Properly covered regions by reads
#----------------------------#
# at least some reads covered
# basecov.txt
# Pos is 0-based
#RefName	Pos	Coverage
log_debug "covered"
cat basecov.txt |
    grep -v '^#' |
    perl -nla -MPath::Tiny -MJSON -MApp::RL::Common -e '
        BEGIN {
            our $name;
            our @list;
            our $limit = JSON->new->decode(
                Path::Tiny::path( q(env.json) )->slurp
            );
            our $length_of = App::RL::Common::read_sizes( q(sr.chr.sizes) );
        }

        sub list_to_ranges {
            my @ranges;
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
            @list = ( $F[1] );
        }

        if ( $F[1] < {{ opt.readl }} ) { # left edges
            if ( $F[2] < {{ opt.mincov }} ) {
                next;
            }

            my $lower = $limit->{lower} * $F[1] / {{ opt.readl }};
            my $upper = $limit->{upper} * $F[1] / {{ opt.readl }};
            if ( $F[2] < $lower or $F[2] > $upper ) {
                next;
            }
        }
        elsif ( $F[1] >= $length_of->{$name} - {{ opt.readl }} ) { # right edges
            if ( $F[2] < {{ opt.mincov }} ) {
                next;
            }

            my $lower = $limit->{lower} * ($length_of->{$name} - $F[1]) / {{ opt.readl }};
            my $upper = $limit->{upper} * ($length_of->{$name} - $F[1]) / {{ opt.readl }};
            if ( $F[2] < $lower or $F[2] > $upper ) {
                next;
            }
        }
        else {
            if ( $F[2] < $limit->{lower} or $F[2] > $limit->{upper} ) {
                next;
            }
        }

        if ( $name eq $F[0] ) {
            push @list, $F[1];
        }
        else {
            my @ranges = list_to_ranges();
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
        }

        END {
            my @ranges = list_to_ranges();
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
spanr cover contig.covered.txt -o contig.covered.yml
spanr stat sr.chr.sizes contig.covered.yml -o contig.covered.csv

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

spanr some contig.covered.yml fill_all.txt -o contig.fill_all.yml
rm -f contig.fill_all.temp.yml; ln -s contig.fill_all.yml contig.fill_all.temp.yml

if [ -s fill_all.txt ]; then
    spanr span contig.fill_all.temp.yml --op fill -n {{ opt.fill | int * 10 }} -o contig.fill_all.1.yml
    rm contig.fill_all.temp.yml; ln -s contig.fill_all.1.yml contig.fill_all.temp.yml

    spanr span contig.fill_all.temp.yml --op excise -n {{ opt.min }} -o contig.fill_all.2.yml
    rm contig.fill_all.temp.yml; ln -s contig.fill_all.2.yml contig.fill_all.temp.yml
fi

# fill small holes
cat sr.chr.sizes |
    cut -f 1 |
    grep -Fx -f fill_all.txt -v \
    > fill_hole.txt

spanr some contig.covered.yml fill_hole.txt -o contig.fill_hole.yml
rm -f contig.fill_hole.temp.yml; ln -s contig.fill_hole.yml contig.fill_hole.temp.yml

if [ -s fill_hole.txt ]; then
    spanr span contig.fill_hole.yml --op fill -n {{ opt.fill }} -o contig.fill_hole.1.yml
    rm contig.fill_hole.temp.yml; ln -s contig.fill_hole.1.yml contig.fill_hole.temp.yml

    spanr span contig.fill_hole.1.yml --op excise -n {{ opt.min }} -o contig.fill_hole.2.yml
    rm contig.fill_hole.temp.yml; ln -s contig.fill_hole.2.yml contig.fill_hole.temp.yml
fi

# get proper regions
perl -MYAML::Syck -MAlignDB::IntSpan -e '
    my $yml1 = YAML::Syck::LoadFile( q{contig.fill_all.temp.yml} );
    my $yml2 = YAML::Syck::LoadFile( q{contig.fill_hole.temp.yml} );

    my $yml_chr = YAML::Syck::LoadFile( q{sr.chr.yml} );

    my $yml = {};

    for my $key ( sort keys %{$yml1} ) {
        $yml->{$key} = $yml1->{$key};
    }
    for my $key ( sort keys %{$yml2} ) {
        $yml->{$key} = $yml2->{$key};
    }

{% if opt.keepedge == "1" -%}
    # repeats on the edge hurt assembling, so fill edges near big island by {{ opt.fill }}
    for my $key ( sort keys %{$yml} ) {
        my $length = AlignDB::IntSpan->new( $yml_chr->{$key} )->max;

        my $set = AlignDB::IntSpan->new( $yml->{$key} );
        next if $set->is_empty;

        my @sets = $set->sets;
        my $begin = $set->min;
        my $end = $set->max;

        if ( $sets[0]->size >= {{ opt.min }} ) {
            my $new_begin = $begin - {{ opt.fill }} * 10;
            $new_begin = 1 if $new_begin < 1;

            $set->add_pair($new_begin, $begin);
        }

        if ( $sets[-1]->size >= {{ opt.min }} ) {
            my $new_end = $end + {{ opt.fill }} * 10;
            $new_end = $length if $new_end > $length;

            $set->add_pair($end, $new_end);
        }

        $yml->{$key} = $set->runlist;
    }
{% endif -%}

    YAML::Syck::DumpFile( q{contig.proper.yml}, $yml );
    '
spanr stat sr.chr.sizes contig.proper.yml -o contig.proper.csv

{% if opt.longest == "1" -%}
perl -MYAML::Syck -MAlignDB::IntSpan -e '
    my $yml = YAML::Syck::LoadFile( q{contig.proper.yml} );

    for my $key ( sort keys %{$yml} ) {
        my $runlist = $yml->{$key};
        my $intspan = AlignDB::IntSpan->new($runlist);
        my @sets = sort { $b->size <=> $a->size } $intspan->sets;
        printf "%s:%s\n", $key, $sets[0]->runlist;
    }
    ' \
    > longest.regions.txt

spanr cover longest.regions.txt -o anchor.yml
{% else -%}
perl -MYAML::Syck -MAlignDB::IntSpan -e '
    my $yml = YAML::Syck::LoadFile( q{contig.proper.yml} );

    for my $key ( sort keys %{$yml} ) {
        my $runlist = $yml->{$key};
        printf "%s:%s\n", $key, $runlist;
    }
    ' \
    > anchor.regions.txt

ln -s contig.proper.yml anchor.yml
{% endif -%}
{# Keep a blank line #}

#----------------------------#
# others
#----------------------------#
spanr compare sr.chr.yml anchor.yml --op diff -o others.yml

perl -MYAML::Syck -MAlignDB::IntSpan -e '
    my $yml = YAML::Syck::LoadFile( q{others.yml} );

    for my $key ( sort keys %{$yml} ) {
        my $runlist = $yml->{$key};
        printf "%s:%s\n", $key, $runlist;
    }
    ' \
    > others.regions.txt

#----------------------------#
# Split SR.fasta to anchor and others
#----------------------------#
log_info "pe.anchor.fa & pe.others.fa"

faops region -l 0 SR.fasta anchor.regions.txt pe.anchor.fa
faops region -l 0 SR.fasta others.regions.txt pe.others.fa

#----------------------------#
# Merging anchors
#----------------------------#
log_info "Merging anchors"
dazz contained \
    pe.anchor.fa \
    --len {{ opt.min }} --idt 0.98 --proportion 0.99999 --parallel {{ opt.parallel }} \
    -o anchor.non-contained.fasta
dazz orient \
    anchor.non-contained.fasta \
    --len {{ opt.min }} --idt 0.98 --parallel {{ opt.parallel }} \
    -o anchor.orient.fasta
dazz merge \
    anchor.orient.fasta --len {{ opt.min }} --idt 0.999 --parallel {{ opt.parallel }} \
    -o anchor.merge0.fasta
dazz contained \
    anchor.merge0.fasta \
    --len {{ opt.min }} --idt 0.98 --proportion 0.99 --parallel {{ opt.parallel }} \
    -o anchor.fasta

#----------------------------#
# Done.
#----------------------------#
find . -type f -name "pe.anchor.fa"   | parallel --no-run-if-empty -j 1 rm
find . -type f -name "anchor.*.fasta" | parallel --no-run-if-empty -j 1 rm

find . -type f -name "contig.fill_hole.*.yml" | parallel --no-run-if-empty -j 1 rm
find . -type l -name "contig.fill_hole.*.yml" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "contig.fill_all.*.yml"  | parallel --no-run-if-empty -j 1 rm
find . -type l -name "contig.fill_all.*.yml"  | parallel --no-run-if-empty -j 1 rm


save START_TIME

END_TIME=$(date +%s)
save END_TIME

RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info "Done."

exit 0
