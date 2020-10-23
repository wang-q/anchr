#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use YAML::Syck qw();

use Path::Tiny qw();
use Text::CSV_XS qw();

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

(    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    <<'MARKDOWN',
Create downloading files for ENA.

Usage: perl %c [options] <infile.yml>

* Three files will be generated
    * .csv for information
    * .ftp.txt for aria2c
    * .md5.txt for checkup
* Example:
    * `sra_prep.pl sra_info.yml`
    * `aria2c -x 9 -s 3 -c -i sra_info.ftp.txt`
    * `md5sum --check sra_info.md5.txt`
MARKDOWN

    [ 'platform|p=s', 'illumina, 454 or pacbio', ],
    [ 'layout|l=s',   'pair or single', ],
    [],
    [ 'help|h', 'display this message' ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

if ( @ARGV != 1 ) {
    my $message = "This command need one input file.\n\tIt found";
    $message .= sprintf " [%s]", $_ for @ARGV;
    $message .= ".\n";
    $usage->die( { pre_text => $message } );
}
for (@ARGV) {
    if ( !Path::Tiny::path($_)->is_file ) {
        $usage->die( { pre_text => "The input file [$_] doesn't exist.\n" } );
    }
}

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $yml = YAML::Syck::LoadFile( $ARGV[0] );
my $basename = Path::Tiny::path( $ARGV[0] )->basename( ".yml", ".yaml" );
$basename .= "." . $opt->{platform} if $opt->{platform};
$basename .= "." . $opt->{layout}   if $opt->{layout};

my $csv = Text::CSV_XS->new( { binary => 1 } )
    or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
$csv->eol("\n");
open my $csv_fh, ">", "$basename.csv";

my $ftp_fn = "$basename.ftp.txt";
my $md5_fn = "$basename.md5.txt";
Path::Tiny::path($ftp_fn)->remove;
Path::Tiny::path($md5_fn)->remove;

$csv->print( $csv_fh, [qw{ name srx platform layout ilength srr spot base }] );
for my $name ( sort keys %{$yml} ) {
    print "$name\n";

    for my $srx ( sort keys %{ $yml->{$name} } ) {
        my $info = $yml->{$name}{$srx};
        print " " x 4, "$srx\n";
        if ( !defined $yml->{$name}{$srx} ) {
            print " " x 8, "Empty record\n";
            next;
        }

        my $platform = $info->{instrument_platform};
        my $layout   = $info->{library_layout};
        my $ilength  = $info->{nominal_length};
        print " " x 8, $platform, " " x 8, $layout, "\n";

        if ( $opt->{platform} ) {
            next unless $platform =~ qr/$opt->{platform}/i;
        }
        if ( $opt->{layout} ) {
            next unless $layout =~ qr/$opt->{layout}/i;
        }

        for my $i ( 0 .. scalar @{ $info->{srr} } - 1 ) {
            my $srr = $info->{srr}[$i];

            my $spot = $info->{srr_info}{$srr}{read_count};
            my $base = $info->{srr_info}{$srr}{base_count};

            $csv->print( $csv_fh,
                [ $name, $srx, $platform, $layout, $ilength, $srr, $spot, $base, ] );
        }

        Path::Tiny::path($ftp_fn)->append( map {"$_\n"} @{ $info->{downloads} } );
        Path::Tiny::path($md5_fn)->append( map {"$_\n"} @{ $info->{md5s} } );
    }
}

close $csv_fh;

__END__
