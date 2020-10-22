{%- include "header" -%}
{# Keep a blank line #}
#----------------------------#
# Run
#----------------------------#
START_TIME=$(date +%s)

# Add masurca to $PATH
export PATH="$(readlink -f $(which masurca) | xargs dirname):$PATH"

#----------------------------#
# Renaming reads
#----------------------------#
log_info 'Processing pe and/or se library reads'

{% if args | length == 2 -%}
rename_filter_fastq \
    'pe' \
    <(exec expand_fastq '{{ args.0 }}' ) \
    <(exec expand_fastq '{{ args.1 }}' ) \
    > 'pe.renamed.fastq'
{% elif args | length == 3 -%}
rename_filter_fastq \
    'pe' \
    <(exec expand_fastq '{{ args.0 }}' ) \
    <(exec expand_fastq '{{ args.1 }}' ) \
    > 'pe.renamed.fastq'

rename_filter_fastq \
    'se' \
    <(exec expand_fastq '{{ args.2 }}' ) \
    '' \
    > 'se.renamed.fastq'
{% else -%}
rename_filter_fastq \
    'pe' \
    <(exec expand_fastq '{{ args.0 }}' ) \
    '' \
    > 'pe.renamed.fastq'
{% endif -%}
{# Keep a blank line #}
#----------------------------#
# Stats of reads
#----------------------------#
head -n 80000 pe.renamed.fastq > pe_data.tmp

KMER=$(
    tail -n 40000 pe_data.tmp |
        perl -e '
            my @lines;
            while ( my $line = <> ) {
                $line = <>;
                chomp($line);
                push( @lines, $line );
                $line = <>;
                $line = <>;
            }
            my @legnths;
            my $min_len    = 100000;
            my $base_count = 0;
            for my $l (@lines) {
                $base_count += length($l);
                push( @lengths, length($l) );
                for $base ( split( "", $l ) ) {
                    if ( uc($base) eq "G" or uc($base) eq "C" ) { $gc_count++; }
                }
            }
            @lengths  = sort { $b <=> $a } @lengths;
            $min_len  = $lengths[ int( $#lengths * .75 ) ];
            $gc_ratio = $gc_count / $base_count;
            $kmer     = 0;
            if ( $gc_ratio < 0.5 ) {
                $kmer = int( $min_len * .7 );
            }
            elsif ( $gc_ratio >= 0.5 && $gc_ratio < 0.6 ) {
                $kmer = int( $min_len * .5 );
            }
            else {
                $kmer = int( $min_len * .33 );
            }
            $kmer++ if ( $kmer % 2 == 0 );
            $kmer = 31  if ( $kmer < 31 );
            $kmer = 127 if ( $kmer > 127 );
            print $kmer;
    ' )
save KMER
log_debug "Choosing kmer size of $KMER"

MIN_Q_CHAR=$(
    head -n 40000 pe_data.tmp |
        awk 'BEGIN{flag=0}{if($0 ~ /^\+/){flag=1}else if(flag==1){print $0;flag=0}}' |
        perl -ne '
            BEGIN { $q0_char = "@"; }

            chomp;
            for $v ( split "" ) {
                if ( ord($v) < ord($q0_char) ) { $q0_char = $v; }
            }

            END {
                $ans = ord($q0_char);
                if   ( $ans < 64 ) { print "33\n" }
                else               { print "64\n" }
            }
    ')
save MIN_Q_CHAR
log_debug "MIN_Q_CHAR: $MIN_Q_CHAR"

#----------------------------#
# Error correct reads
#----------------------------#
JF_SIZE=$(
    ls -l *.fastq |
        awk '{n+=$5} END{s=int(n / 50 * 1.1); if(s>{{ opt.jf }})print s;else print "{{ opt.jf }}";}'
)
perl -e '
    if(int('$JF_SIZE') > {{ opt.jf }}) {
        print "WARNING: JF_SIZE set too low, increasing JF_SIZE to '$JF_SIZE'.\n";
    }
    '

log_info Creating mer database for Quorum.
quorum_create_database \
    -t {{ opt.parallel }} \
    -s $JF_SIZE -b 7 -m 24 -q $((MIN_Q_CHAR + 5)) \
    -o quorum_mer_db.jf.tmp \
    pe.renamed.fastq {% if args | length == 3 %}se.renamed.fastq {% endif %}\
    && mv quorum_mer_db.jf.tmp quorum_mer_db.jf
if [ $? -ne 0 ]; then
    log_warn Increase JF_SIZE by --jf, the recommendation value is genome_size*coverage/2
    exit 1
fi

# -m Minimum count for a k-mer to be considered "good" (1)
# -g Number of good k-mer in a row for anchor (2)
# -a Minimum count for an anchor k-mer (3)
# -w Size of window (10)
# -e Maximum number of error in a window (3)
# As we have trimmed reads with sickle, we lower `-e` to 1 from original value of 3,
# remove `--no-discard`.
# And we only want most reliable parts of the genome other than the whole genome, so dropping rare
# k-mers is totally OK for us. Raise `-m` from 1 to 3, `-g` from 1 to 3, and `-a` from 1 to 4.
log_info Error correct reads.
quorum_error_correct_reads \
    -q $((MIN_Q_CHAR + 40)) \
    -m 3 -s 1 -g 3 -a 4 -t {{ opt.parallel }} -w 10 -e 1 \
    quorum_mer_db.jf \
    pe.renamed.fastq {% if args | length == 3 %}se.renamed.fastq {% endif %}\
    -o {{ opt.prefix }}.cor --verbose 1>quorum.err 2>&1 \
|| {
    mv {{ opt.prefix }}.cor.fa {{ opt.prefix }}.cor.fa.failed;
    log_warn Error correction of reads failed.;
    exit 1;
}

log_debug "Discard any reads with subs"
mv {{ opt.prefix }}.cor.fa {{ opt.prefix }}.cor.sub.fa
cat {{ opt.prefix }}.cor.sub.fa |
    grep -E '^>\w+\s*$' -A 1 |
    sed '/^--$/d' |
    perl -nl -e '
        BEGIN {
            our $prev_name, $prev_idx, $prev_seq, $cur_name, $cur_idx, $cur_seq;
        }

        if ( substr( $_, 0, 1 ) eq q{>} ) {
            ($cur_name) = split /\s+/, $_;
            $cur_name =~ s/^>//;
            $cur_idx = substr $cur_name, 2;
        }
        else {
            $cur_seq = $_;

            if ( $cur_idx & 1 ) { # odd
                if ( defined $prev_name ) {
                    if ( $cur_idx - $prev_idx == 1 ) {
                        print qq{>$prev_name/1};
                        print qq{$prev_seq};
                        print qq{>$prev_name/2};
                        print qq{$cur_seq};
                    }
                    else {
                        print qq{>$prev_name/1};
                        print qq{$prev_seq};
                        print qq{>$prev_name/2};
                        print qq{N};

                        print qq{>$cur_name/1};
                        print qq{N};
                        print qq{>$cur_name/2};
                        print qq{$cur_seq};
                    }

                    undef $prev_name;
                    undef $prev_idx;
                    undef $prev_seq;
                }
                else {
                    print qq{>$cur_name/1};
                    print qq{N};
                    print qq{>$cur_name/2};
                    print qq{$cur_seq};
                }
            }
            else { # even
                if ( defined $prev_name ) {
                    print qq{>$prev_name/1};
                    print qq{$prev_seq};
                    print qq{>$prev_name/2};
                    print qq{N};
                }

                $prev_name = $cur_name;
                $prev_idx = $cur_idx;
                $prev_seq = $cur_seq;
            }
        }
    ' \
    > {{ opt.prefix }}.cor.fa

rm {{ opt.prefix }}.cor.sub.fa

#----------------------------#
# Estimating genome size.
#----------------------------#
log_info Estimating genome size.

{% if opt.estsize == 'auto' -%}
jellyfish count -m 31 -t {{ opt.parallel }} -C -s $JF_SIZE -o k_u_hash_0 {{ opt.prefix }}.cor.fa
ESTIMATED_GENOME_SIZE=$(
    jellyfish histo -t {{ opt.parallel }} -h 1 k_u_hash_0 |
        tail -n 1 |
        awk '{print $2}'
)
save ESTIMATED_GENOME_SIZE
log_debug "Estimated genome size: $ESTIMATED_GENOME_SIZE"
{% else -%}
ESTIMATED_GENOME_SIZE={{ opt.estsize }}
save ESTIMATED_GENOME_SIZE
log_debug "You set ESTIMATED_GENOME_SIZE of $ESTIMATED_GENOME_SIZE"
{% endif -%}
{# Keep a blank line #}
log_debug "Reads stats with faops"
SUM_IN=$( faops n50 -H -N 0 -S pe.renamed.fastq {% if args | length == 3 %}se.renamed.fastq {% endif %})
save SUM_IN
SUM_OUT=$( faops n50 -H -N 0 -S {{ opt.prefix }}.cor.fa )
save SUM_OUT

#----------------------------#
# Shuffle interleaved reads.
#----------------------------#
log_info Shuffle interleaved reads.
mv {{ opt.prefix }}.cor.fa {{ opt.prefix }}.interleave.fa
cat {{ opt.prefix }}.interleave.fa |
    awk '{
        OFS="\t"; \
        getline seq; \
        getline name2; \
        getline seq2; \
        print $0,seq,name2,seq2}' |
    shuf |
    awk '{OFS="\n"; print $1,$2,$3,$4}' \
    > {{ opt.prefix }}.cor.fa
rm {{ opt.prefix }}.interleave.fa
pigz -p {{ opt.parallel }} {{ opt.prefix }}.cor.fa

#----------------------------#
# Done.
#----------------------------#
find . -type f -name "quorum_mer_db.jf" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "k_u_hash_0"       | parallel --no-run-if-empty -j 1 rm
find . -type f -name "*.tmp"            | parallel --no-run-if-empty -j 1 rm
find . -type f -name "pe.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "se.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "pe.cor.sub.fa"    | parallel --no-run-if-empty -j 1 rm
find . -type f -name "*.cor.log"        | parallel --no-run-if-empty -j 1 rm

save START_TIME

END_TIME=$(date +%s)
save END_TIME

RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info Done.

exit 0
