#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use YAML::Syck qw();

use LWP::Simple qw();
use JSON::PP qw();
use Number::Format;
use List::MoreUtils::PP;
use Path::Tiny qw();
use Text::CSV_XS;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    <<'MARKDOWN',
Grab information from ENA.

Usage: perl %c [options] <infile.csv> > <outfile.json>

* <infile> == stdin means read from STDIN
* <infile> format
    * first column is one SRA object ID, /[DES]R\w\d+/, SRP or SRX
    * second column is the name of one group
    * other columns are optional
MARKDOWN

    [ 'sra', "download sra instead of fastq", ],
    [],
    [ 'help|h',    'display this message' ],
    [ 'verbose|v', 'verbose mode' ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

if ( @ARGV != 1 ) {
    my $message = "This command need one filename.\n\tIt found";
    $message .= sprintf " [%s]", $_ for @ARGV;
    $message .= ".\n";
    $usage->die( { pre_text => $message } );
}
for (@ARGV) {
    next if lc $_ eq "stdin";
    if ( !Path::Tiny::path($_)->is_file ) {
        $usage->die( { pre_text => "The input file [$_] doesn't exist.\n" } );
    }
}

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $master = {};

my $csv_fh;
if ( lc $ARGV[0] eq 'stdin' ) {
    $csv_fh = *STDIN;
}
else {
    open $csv_fh, "<", $ARGV[0];
}

my $csv = Text::CSV_XS->new( { binary => 1 } )
    or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
while ( my $row = $csv->getline($csv_fh) ) {
    next if $row->[0]     =~ /^#/;
    next unless $row->[0] =~ /[DES]R\w\d+/;

    my ( $key, $name ) = ( $row->[0], $row->[1] );
    if ( !defined $name ) {
        $name = $key;
    }
    warn "key: [$key]\tname: [$name]\n";

    my @srx = erp_worker( $key, $opt->{verbose} );
    warn "@srx\n";

    my $sample
        = exists $master->{$name}
        ? $master->{$name}
        : {};
    for (@srx) {
        $sample->{$_} = erx_worker( $_, $opt->{sra}, $opt->{verbose} );
    }
    $master->{$name} = $sample;
    warn "\n";
}
close $csv_fh;

print YAML::Syck::Dump($master);
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub erp_worker {
    my $term    = shift;
    my $verbose = shift;

    my $url_part1 = "http://www.ebi.ac.uk/ena/portal/api/filereport?accession=";
    my $url_part2
        = "&result=read_run&fields=secondary_study_accession,experiment_accession" . "&format=tsv";
    my $url = $url_part1 . $term . $url_part2;
    warn "$url\n" if $verbose;

    my @lines = split /\n/, LWP::Simple::get($url);

    my @srx;
    for (@lines) {
        if (/([DES]RX\d+)/) {
            push @srx, $1;
        }
    }
    @srx = List::MoreUtils::PP::uniq(@srx);

    return @srx;
}

sub erx_worker {
    my $term    = shift;
    my $use_sra = shift;
    my $verbose = shift;

    my $url_part1 = "http://www.ebi.ac.uk/ena/portal/api/filereport?accession=";
    my $url_part2
        = "&result=read_run&fields=secondary_study_accession,secondary_sample_accession,"
        . "experiment_accession,run_accession,scientific_name,"
        . "instrument_platform,instrument_model,"
        . "library_name,nominal_length,library_layout,library_source,library_selection,"
        . "read_count,base_count,"
        . ( $use_sra ? "sra_md5,sra_ftp" : "fastq_md5,fastq_ftp" )
        . "&format=json";
    my $url = $url_part1 . $term . $url_part2;
    warn "$url\n" if $verbose;

    my $content = LWP::Simple::get($url);
    if ( !scalar $content ) {
        warn "Can't get any SRR, please check.\n";
        return;
    }
    my $json = JSON::PP::decode_json($content);

    my $info = {
        srp                 => $json->[0]{secondary_study_accession},
        srs                 => $json->[0]{secondary_sample_accession},
        srx                 => $json->[0]{experiment_accession},
        scientific_name     => $json->[0]{scientific_name},
        instrument_platform => $json->[0]{instrument_platform},
        instrument_model    => $json->[0]{instrument_model},
        library_name        => $json->[0]{library_name},
        nominal_length      => $json->[0]{nominal_length},
        library_layout      => $json->[0]{library_layout},
        library_source      => $json->[0]{library_source},
        library_selection   => $json->[0]{library_selection},
        srr_info            => {},
    };
    my ( @srrs, @downloads, @md5s );
    for my $elem ( @{$json} ) {
        my $srr = $elem->{run_accession};
        push @srrs, $srr;

        # ftp path and md5
        my @parts_ftp = map { "ftp://" . $_ } grep {defined} split ";",
            $elem->{ $use_sra ? "sra_ftp" : "fastq_ftp" };
        push @downloads, @parts_ftp;

        my @basenames = map  { ( split "/", $_ )[-1] } @parts_ftp;
        my @parts_md5 = grep {defined} split ";", $elem->{ $use_sra ? "sra_md5" : "fastq_md5" };
        for my $i ( 0 .. $#basenames ) {
            push @md5s, ( sprintf "%s %s", $parts_md5[$i], $basenames[$i] );
        }

        $info->{srr_info}{$srr} = {
            read_count => $elem->{read_count},
            base_count => Number::Format::format_bytes( $elem->{base_count} ),
        };
    }
    $info->{srr}       = \@srrs;
    $info->{downloads} = \@downloads;
    $info->{md5s}      = \@md5s;

    return $info;
}

__END__
