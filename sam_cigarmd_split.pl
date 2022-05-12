#' This Source Code Form is subject to the terms of the Mozilla Public
#' License, v. 2.0. If a copy of the MPL was not distributed with this
#' file, You can obtain one at http://mozilla.org/MPL/2.0/.
#'
#' Youtao Lu@Kim Lab, 2016-2020

use warnings;
use strict;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

our $VERSION = '0.8';
our $LOGGER  = get_logger(__PACKAGE__);
my ( $inFile, $inFh ) = ( undef, undef );
my $splice = 0;
my $reference;
my $nameSorted = 0;
my $ncores     = 1;
my ( $outFile, $outFh ) = ( '', undef );
my $endType = "PE";
my ( $help, $version ) = ( 0, 0 );
my $debug = "info";

sub usage {
    print STDERR <<DOC;
Summary:
    Splits query sequence and phred scores into softclip (S), insertion (I), match (=) and mismatch (X) accoriding to CIGAR and tag MD.

Usage:
    perl $0 --inFile|-i <input.bam|sam> [--outFile|-o <output.txt>] [--endType PE] [--debug info] [--version|-v] [--help|-h] [--ncores 1]

Options:
    --inFile, -i  BAM or SAM input, has to be name sorted if endType == PE. If multiple mapping, only primary alignment will be analyzed;
    --endType     PE or SE (default PE);
    --nameSorted  whether or not the input is name sorted. If not, it will be name sorted and saved to another file with '.nameSorted.bam' suffix;
    --ncores      how many threads for sorting;
    --splice      whether or not output the sliced bases, requires --reference to be present. It can be very slow if enabled as it calls `samtools faidx` for every splice junction ('N'); 
    --reference   path of the reference genome, required if --splice is on;
    --debug       debug level, choose from fatal, error, warn, info, debug, trace;
    --outFile, -o (optional, if omitted, write to STDOUT; otherwise, if ending with '.gz', will be GZ compressed)

    The output is a table with 13 (SE) or 25 (PE) columns. For SE, they are
        1) read ID,
        2) aligned query seq., 3) aligned ref. seq., 4) aligned phred scores,
        5) softclipped query splits, 6) softclipped phred splits,
        7) inserted query splits, 8) inserted phred splits,
        9) matched query splits, 10) matched phred splits,
        11) mismatched query splits, 12) mismatched ref. seq, 13) mismatched phred splits, 
        For PE, the ouput replicates column 2-13 for R2. 
         

Examples:
    Suppose a SE read has the following record:
    query seq.   = GTCTATTCTCTTCTCTTCTATTCTCTTCTGTACTCTTCTCTTATGTTCTCTTCTC
    phred scores = 6EE//EE/A/EEEEAE//A/AAEE/EEEE/E//AE6/E///6/E<//<E/6//E6
    CIGAR        = 5S39M2I1M1D5M3S
    MD tag       = 26t10c2^5

    below are colunm 2-13 from the output (~ denotes a one-base deletion in query or in ref.; * denotes any single base in the ref.)
    2)  GTCTATTCTCTTCTCTTCTATTCTCTTCTGTACTCTTCTCTTATGTT~CTCTTCTC 
    3)  *****==========================t==========c=~~=G=====***
    4)  6EE//EE/A/EEEEAE//A/AAEE/EEEE/E//AE6/E///6/E<//~<E/6//E6
    5)  GTCTA,CTC
    6)  6EE//,/E6
    7)  GT
    8)  </
    9)  TTCTCTTCTCTTCTATTCTCTTCTGT,CTCTTCTCTT,T,T,CTCTT
    10) EE/A/EEEEAE//A/AAEE/EEEE/E,/AE6/E///6,E,/,<E/6/
    11) A,A
    12) t,c
    13) /,/
    
    Another example shows the effect of --splice. Suppose a SE read that has undergone splicing has the following record:
    query seq.   = CATCCCGCCTCCGTCCCCGTTGCTGCCGCCATACACGCTCGCAGTTCTTAGCTCTTCTGTCGGAAACTGGTGTCTTTCCCCTTTCTGTTCT
    phred scores = F,,F:,,F:FFFF:F,F,:F::F,,:F,FF:::FF,,F:::,,,:,F,FFF,F:,,:FFFF:FFFF,F:,FF,:F,:F:,,:,,FFFFF,F
    CIGAR        = 51M10N40M
    MD tag       = 1T43G37G7

    when --splice is off (by default), column 2-4 will be:
    2) CATCCCGCCTCCGTCCCCGTTGCTGCCGCCATACACGCTCGCAGTTCTTAGCTCTTCTGTCGGAAACTGGTGTCTTTCCCCTTTCTGTTCT
    3) =T===========================================G=====================================G=======
    4) F,,F:,,F:FFFF:F,F,:F::F,,:F,FF:::FF,,F:::,,,:,F,FFF,F:,,:FFFF:FFFF,F:,FF,:F,:F:,,:,,FFFFF,F
    
    when --splice is on, column 2-4 will be (~ denotes the missing bases in query seq. due to splicing):
    2) CATCCCGCCTCCGTCCCCGTTGCTGCCGCCATACACGCTCGCAGTTCTTAG~~~~~~~~~~CTCTTCTGTCGGAAACTGGTGTCTTTCCCCTTTCTGTTCT
    3) =T===========================================G=====GTAAGCTTTG================================G=======
    4) F,,F:,,F:FFFF:F,F,:F::F,,:F,FF:::FF,,F:::,,,:,F,FFF~~~~~~~~~~,F:,,:FFFF:FFFF,F:,FF,:F,:F:,,:,,FFFFF,F

Extension: 
    Given the splits, we can easily compute stats such as average phred score in matches and mismatches, or base composition difference between soft-clipping and insertion. For example, 

    *) to extract the triplet of aligned query, ref. and phred:
    perl sam_cigarmd_split.pl -i test.bam 2> err.log | awk 'NR > 1 {print \$2 "\\n" \$3 "\\n" \$4 "\\n"}' > test_align.txt

    *) to calculate the average phred scores in matches and mismatches:
    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ use List::Util "sum" }; \@F = split("\\t", \$_, -1); \$phred = \$F[9] . \$F[21]; \$phred =~ s/,//g; \@p = map{ ord(\$_) - 33 } split(\"\", \$phred); print \$F[0], "\\t", \$#p > -1 ? sum(\@p) / (\$#p+1) : ""' > test_avgphred_match.txt
perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ use List::Util "sum" }; \@F = split("\\t", \$_, -1); \$phred = \$F[12] . \$F[24]; \$phred =~ s/,//g; \@p = map{ ord(\$_) - 33 } split(\"\", \$phred); print \$F[0], "\\t", \$#p > -1 ? sum(\@p) / (\$#p+1) : ""' > test_avgphred_mismatch.txt

    *) to calculate the base composition of mismatched ref.:
    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ %comp }; \@F = split("\\t", \$_, -1); \$x = \$F[11] . \$F[23]; \$x =~ s/,//g; \@b = split("", \$x); for \$n (\@b) { \$comp{uc(\$n)}++ }; END{ for \$k (sort {\$comp{\$b} <=> \$comp{\$a}} keys(%comp)) { print \$k, "\\t", \$comp{\$k}} }'

    *) to calculate the mutation transition probability of the mismatches:
    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -ne 'BEGIN{ %comp; print "\\t" . join("\\t", qw(A C G T)), "\n" }; \@F = split("\\t", \$_, -1); \$x = \$F[11] . \$F[23]; \$x =~ s/,//g; \$y = \$F[10] . \$F[22]; \$y =~ s/,//g; \@b = split("", \$x); \@m = split("", \$y); for \$i (0..\$#b) { \$comp{uc(\$b[\$i])}{uc(\$m[\$i])}++ }; END{ for \$r (qw(A C G T)) { print \$r; for \$c (qw(A C G T)) { print "\\t", \$comp{\$r}{\$c} // 0 }; print "\n" } }'
DOC
}

GetOptions(
    "inFile|i=s"  => \$inFile,
    "outFile|o:s" => \$outFile,
    "endType:s"   => \$endType,
    "splice"      => \$splice,
    "reference:s" => \$reference,
    "nameSorted"  => \$nameSorted,
    "ncores:i"    => \$ncores,
    "help|h"      => \$help,
    "version|v"   => \$version,
    "debug:s"     => \$debug,
) or &usage() && exit(-1);

&usage() && exit(-1) if $help or !defined($inFile);
print "$0 v$VERSION\n" && exit(0) if $version;
die("--endType can only choose from PE and SE!\n")
  if $endType ne "PE" && $endType ne "SE";
die("--inFile does not exist!\n") if !defined($inFile) || !-e $inFile;
die("--inFile does not look like SAM or BAM!\n") if $inFile !~ /\.(s|b)am$/i;
die("--reference is missing or file not found!\n")
  if $splice && ( !defined($reference) || !-e $reference );
$reference = '' if !defined($reference);

if ( $debug eq "fatal" ) {
    $LOGGER->level($FATAL);
}
elsif ( $debug eq "error" ) {
    $LOGGER->level($ERROR);
}
elsif ( $debug eq "warn" ) {
    $LOGGER->level($WARN);
}
elsif ( $debug eq "info" ) {
    $LOGGER->level($INFO);
}
elsif ( $debug eq "debug" ) {
    $LOGGER->level($DEBUG);
}
elsif ( $debug eq "trace" ) {
    $LOGGER->level($TRACE);
}
my $appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::Screen");
my $layout   = Log::Log4perl::Layout::PatternLayout->new(
    "[%d{yyyy-MM-dd HH:mm:ss.SSS Z}] %m");    #%F %L
$appender->layout($layout);
$LOGGER->add_appender($appender);

$LOGGER->info(
"{ VERSION = $VERSION, inFile = $inFile, outFile = $outFile, endType == $endType, splice = $splice, reference = $reference, nameSorted = $nameSorted, ncores = $ncores, help = $help, debug = $debug, version = $version }\n"
);

{

    package Tag_md;

    sub new {
        my $class = shift;
        my $tag   = shift;
        my $self  = bless {
            tag      => $tag,
            oplen    => 0,
            opchr    => '',
            mmstr    => '',
            mmout    => '',
            consumed => 0,
            dels     => [],
        }, $class;
        return $self;
    }

    sub shift_op
    {    # each call will take one string of digits or one string of bases
        my $self = shift;
        return 0 if ( length( $self->{tag} ) == 0 );
        if ( ord( substr( $self->{tag}, 0, 1 ) ) == 48 )
        { # in 10A3T0T10, 0 indicates neighboring SNPs (https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files); we simply skip it.
            $self->{tag} = substr( $self->{tag}, 1 );
            $LOGGER->trace(
"(shift_op) \$self->{tag} = $self->{tag}, \$self->{oplen} = $self->{oplen}, \$self->{opchr} = $self->{opchr}, \$self->{mmstr} = $self->{mmstr}, \$self->{mmout} = $self->{mmout}, \$self->{consumed} = $self->{consumed}, \@{\$self->{dels}} = @{$self->{dels}}\n"
            );
        }
        elsif ( ord( substr( $self->{tag}, 0, 1 ) ) < 58 ) {    # is digit
            $self->{tag} =~ s/^(\d+)//;
            $self->{oplen} = $1;
            $self->{opchr} = '=';
            $LOGGER->trace(
"(shift_op) \$self->{tag} = $self->{tag}, \$self->{oplen} = $self->{oplen}, \$self->{opchr} = $self->{opchr}, \$self->{mmstr} = $self->{mmstr}, \$self->{mmout} = $self->{mmout}, \$self->{consumed} = $self->{consumed}, \@{\$self->{dels}} = @{$self->{dels}}\n"
            );
        }
        elsif ( substr( $self->{tag}, 0, 1 ) eq '^' )
        { # deletion ^[acgtnACGTN]+, e.g. in 56^AGG45, ^AGG is a string of 3 deletions in the read, corresponding to '3D' in cigar, AGG are deleted letters in the reference genome.
            $self->{tag} =~ s/\^(\D+)//;
            push @{ $self->{dels} }, $1;
            $LOGGER->trace(
"(shift_op) \$self->{tag} = $self->{tag}, \$self->{oplen} = $self->{oplen}, \$self->{opchr} = $self->{opchr}, \$self->{mmstr} = $self->{mmstr}, \$self->{mmout} = $self->{mmout}, \$self->{consumed} = $self->{consumed}, \@{\$self->{dels}} = @{$self->{dels}}\n"
            );
        }
        else {
            $self->{tag} =~ s/^(\D+)//;
            $self->{mmstr} = $1;
            $self->{oplen} = length( $self->{mmstr} );
            $self->{opchr} = 'X';
            $LOGGER->trace(
"(shift_op) \$self->{tag} = $self->{tag}, \$self->{oplen} = $self->{oplen}, \$self->{opchr} = $self->{opchr}, \$self->{mmstr} = $self->{mmstr}, \$self->{mmout} = $self->{mmout}, \$self->{consumed} = $self->{consumed}, \@{\$self->{dels}} = @{$self->{dels}}\n"
            );
        }
        return $self->{oplen};
    }

    sub consume_mapped {
        my $self  = shift;
        my $m_len = shift;
        if ( $m_len < $self->{oplen} )
        {    # consume part of current op and cache the leftover
            $self->{oplen} -= $m_len;
            $self->{consumed} = $m_len;
            if ( $self->{opchr} eq 'X' )
            {    # (5M)... : (AGCATGTA)6... -> ()... : (GTA)6...
                $self->{mmout} = substr( $self->{mmstr}, 0, $m_len );
                $self->{mmstr} = substr( $self->{mmstr}, $m_len );
            }
            else {    # (5M)... : (15)AGT... -> ()... : (10)AGT...
                $self->{mmout} = $self->{mmstr} = '';
            }
        }
        else {        # consume all current op
            $self->{consumed} = $self->{oplen};
            if ( $self->{opchr} eq 'X' )
            {         # (15M)... : (ATC)3... -> (12M)... : (3)...
                $self->{mmout} = $self->{mmstr};
                $self->{mmstr} = '';
                $self->{oplen} = 0;
            }
            else {    # (15M)... : 4ATG... -> (11M)... : (ATG)
                $self->{mmout} = $self->{mmstr} = '';
                $self->{oplen} = 0;
            }
        }
        return $self->{consumed};
    }
}

{

    package Read;

    sub new {
        my $class = shift;
        my ( $readID, $whichInPair, $chr, $pos, $cigar, $seq, $phred, $md ) =
          @_;
        my $self = bless {
            readID      => $readID,
            whichInPair => $whichInPair,
            chr         => $chr,
            pos         => $pos,
            seq         => $seq,
            phred       => $phred,
            cigar       => $cigar,
            md          => $md,
            seq_s       => { S => [], I => [], '=' => [], X => [] },
            ref_s       => { X => [] },
            phred_s     => { S => [], I => [], '=' => [], X => [] },
            seq_a       => '',
            ref_a       => '',
            phred_a     => '',
        }, $class;
    }

    sub parse {
        my $self   = shift;
        my $tag_md = Tag_md->new( $self->{md} );
        while ( $self->{cigar} =~ s/(\d+)(M|I|D|N|S|H|P|=|X)// ) {
            my $oplen = $1;
            my $opchr = $2;
            if ( $opchr eq 'D' )
            { #  || $opchr eq 'N' || $opchr eq 'H' || $opchr eq 'P'; does not consume query, only consumes reference; we don't print the reference bases since we don't have data in SAM;
                $LOGGER->fatal("\$tag_md->{oplen} (= $tag_md->{oplen}) must be 0!") && die() if $tag_md->{oplen};
                $tag_md->shift_op(); #if !$tag_md->{oplen};
                $self->{seq_a}   .= '~' x $oplen;
                $self->{phred_a} .= '~' x $oplen;
                my $del = shift @{ $tag_md->{dels} };
                $self->{ref_a} .= $del;
                $self->{pos} += $oplen;
                next;
            }
            if ( $opchr eq 'S' || $opchr eq 'I' )
            { # soft clipped, remove the part from $seq and $phred, keep $tag_md unchanged
                push @{ $self->{seq_s}->{$opchr} },
                  substr( $self->{seq}, 0, $oplen );
                push @{ $self->{phred_s}->{$opchr} },
                  substr( $self->{phred}, 0, $oplen );
                $self->{seq_a} .= substr( $self->{seq}, 0, $oplen );
                if ( $opchr eq 'S' ) {
                    $self->{ref_a} .= '*' x $oplen;
                    $self->{pos} += $oplen;
                }
                else {
                    $self->{ref_a} .= '~' x $oplen;
                }
                $self->{phred_a} .= substr( $self->{phred}, 0, $oplen );
                $self->{seq}   = substr( $self->{seq},   $oplen );
                $self->{phred} = substr( $self->{phred}, $oplen );
                next;
            }
            if ( $opchr eq '=' || $opchr eq 'X' )
            { # match, remove the part from $seq and $phred, remove the part from $tag_md
                $self->{pos} += $oplen;
                while ( $oplen > 0 ) {
                    $tag_md->consume_mapped($oplen);
                    push @{ $self->{seq_s}->{$opchr} },
                      substr( $self->{seq}, 0, $tag_md->{consumed} );
                    push @{ $self->{phred_s}->{$opchr} },
                      substr( $self->{phred}, 0, $tag_md->{consumed} );
                    $self->{seq_a}   .= substr( $self->{seq},   0, $oplen );
                    $self->{phred_a} .= substr( $self->{phred}, 0, $oplen );
                    if ( $opchr eq 'X' ) {
                        push @{ $self->{ref_s}->{$opchr} }, $tag_md->{mmout};
                        $self->{ref_a} .= $tag_md->{mmout};
                    }
                    else
                    { # substr($self->{seq}, 0, $tag_md->{consumed}); # output original bases
                        $self->{ref_a} .=
                          "=" x $tag_md->{consumed};    # output = instead
                    }
                    $self->{seq} = substr( $self->{seq}, $tag_md->{consumed} );
                    $self->{phred} =
                      substr( $self->{phred}, $tag_md->{consumed} );
                    $oplen -= $tag_md->{consumed};
                    $tag_md->shift_op() if !$tag_md->{oplen};
                } # if oplen is greater than $md_tag->{oplen}, there is mismatch in $md_tag->oplen
                next;
            }
            if ( $opchr eq 'M' ) {    # mapped, can be match or mismatch
                $self->{pos} += $oplen;
                while ( $oplen > 0 ) {
                    $tag_md->consume_mapped($oplen);
                    push @{ $self->{seq_s}->{ $tag_md->{opchr} } },
                      substr( $self->{seq}, 0, $tag_md->{consumed} );
                    push @{ $self->{phred_s}->{ $tag_md->{opchr} } },
                      substr( $self->{phred}, 0, $tag_md->{consumed} );
                    $self->{seq_a} .=
                      substr( $self->{seq}, 0, $tag_md->{consumed} );
                    $self->{phred_a} .=
                      substr( $self->{phred}, 0, $tag_md->{consumed} );
                    if ( $tag_md->{opchr} eq 'X' ) {    # mismatch
                        push @{ $self->{ref_s}->{ $tag_md->{opchr} } },
                          $tag_md->{mmout};
                        $self->{ref_a} .= $tag_md->{mmout};
                    }
                    else {    # substr($self->{seq}, 0, $tag_md->{consumed});
                        $self->{ref_a} .= '=' x $tag_md->{consumed};
                    }
                    $self->{seq} = substr( $self->{seq}, $tag_md->{consumed} );
                    $self->{phred} =
                      substr( $self->{phred}, $tag_md->{consumed} );
                    $oplen -= $tag_md->{consumed};
                    $tag_md->shift_op() if $tag_md->{oplen} == 0;
                } # if oplen is greater than $md_tag->{oplen}, there is mismatch in $md_tag->oplen
                next;
            }
            if ( $splice && $opchr eq 'N' ) {
                my $l = $self->{pos};
                my $r = $l + $oplen - 1;
                my $fh;
                open( $fh, "samtools faidx $reference $self->{chr}:$l-$r |" )
                  or $LOGGER->fatal(
                    "Cannot run samtools faidx $reference $self->{chr}:$l-$r\n")
                  && die
                  ();
                <$fh>;
                my $ref = "";
                while (<$fh>) {
                    chomp;
                    $ref .= $_;
                }
                close $fh;
                $self->{ref_a}   .= $ref;
                $self->{seq_a}   .= '~' x $oplen;
                $self->{phred_a} .= '~' x $oplen;
            }
        }
    }
}

{

    package ReadPair;

    sub new {
        my $class = shift;
        my $pair  = shift;
        my $self  = bless {
            R1 => $pair->{R1},
            R2 => $pair->{R2}
        }, $class;
    }

    sub parse {
        my $self = shift;
        $self->{R1}->parse() if defined( $self->{R1} );
        $self->{R2}->parse() if defined( $self->{R2} );
    }

    sub output {
        my $self  = shift;
        my $outFh = shift;
        my $readID =
          defined( $self->{R1} )
          ? $self->{R1}->{readID}
          : $self->{R2}->{readID};
        my $l1 =
          defined( $self->{R1} )
          ? "\t"
          . $self->{R1}->{seq_a} . "\t"
          . $self->{R1}->{ref_a} . "\t"
          . $self->{R1}->{phred_a} . "\t"
          . join( ",", @{ $self->{R1}->{seq_s}->{S} } ) . "\t"
          . join( ",", @{ $self->{R1}->{phred_s}->{S} } ) . "\t"
          . join( ",", @{ $self->{R1}->{seq_s}->{I} } ) . "\t"
          . join( ",", @{ $self->{R1}->{phred_s}->{I} } ) . "\t"
          . join( ",", @{ $self->{R1}->{seq_s}->{'='} } ) . "\t"
          . join( ",", @{ $self->{R1}->{phred_s}->{'='} } ) . "\t"
          . join( ",", @{ $self->{R1}->{seq_s}->{X} } ) . "\t"
          . join( ",", @{ $self->{R1}->{ref_s}->{X} } ) . "\t"
          . join( ",", @{ $self->{R1}->{phred_s}->{X} } )
          : "\t" x 12;
        my $l2 =
          defined( $self->{R2} )
          ? "\t"
          . $self->{R2}->{seq_a} . "\t"
          . $self->{R2}->{ref_a} . "\t"
          . $self->{R2}->{phred_a} . "\t"
          . join( ",", @{ $self->{R2}->{seq_s}->{S} } ) . "\t"
          . join( ",", @{ $self->{R2}->{phred_s}->{S} } ) . "\t"
          . join( ",", @{ $self->{R2}->{seq_s}->{I} } ) . "\t"
          . join( ",", @{ $self->{R2}->{phred_s}->{I} } ) . "\t"
          . join( ",", @{ $self->{R2}->{seq_s}->{'='} } ) . "\t"
          . join( ",", @{ $self->{R2}->{phred_s}->{'='} } ) . "\t"
          . join( ",", @{ $self->{R2}->{seq_s}->{X} } ) . "\t"
          . join( ",", @{ $self->{R2}->{ref_s}->{X} } ) . "\t"
          . join( ",", @{ $self->{R2}->{phred_s}->{X} } )
          : "\t" x 12;
        my $l = $readID . $l1 . $l2 . "\n";
        print $outFh $l;
    }
}

{

    package SamReader;
    use IO::Zlib;

    sub new {
        my $class = shift;
        my ( $inFile, $outFile, $endType, $nameSorted ) = @_;
        my $self = bless {
            inFile     => $inFile,
            inFh       => undef,
            outFile    => $outFile,
            outFh      => undef,
            endType    => $endType,
            nameSorted => $nameSorted,
            ncores     => $ncores,
        }, $class;
        return $self;
    }

    sub init_fh {
        my $self = shift;
        if ( $self->{endType} eq "PE" && !$self->{nameSorted} ) {
            my $tmpFile = $self->{inFile};
            $tmpFile =~ s/\.(bam|sam)$//i;
            my $suffix = $1;
            $tmpFile .= ".nameSorted" . ".$suffix";
            $LOGGER->warn(
"$self->{inFile} is not name sorted. Name sorting it and saving to $tmpFile...\n"
            );
            $LOGGER->warn("$tmpFile exists already! It will be overwritten.\n")
              if -e $tmpFile;
            my $exit = system(
"samtools sort -n -\@ $self->{ncores} -o $tmpFile $self->{inFile}"
            );
            $LOGGER->fatal("Error happended when sorting $self->{inFile}!")
              && die()
              if $exit;
            $self->{inFile} = $tmpFile;
        }
        $LOGGER->info("Opening file handle for $self->{inFile}...\n");
        open( $self->{inFh}, "samtools view -F 0x100 $self->{inFile} |" )
          or $LOGGER->fatal("Cannot open $self->{inFile} for read!\n") && die();
        if ( $self->{outFile} eq '' ) {
            $LOGGER->info("Opening file handle for STDOUT...\n");
            $self->{outFh} = *STDOUT;
        }
        else {
            $LOGGER->info("Opening file handle for $self->{outFile}...\n");
            if ( $self->{outFile} =~ /\.gz$/i ) {
                ( $self->{outFh} = IO::Zlib->new( $self->{outFile}, 'w' ) )
                  or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
                  && die();
            }
            else {
                open( $self->{outFh}, ">", $self->{outFile} )
                  or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
                  && die();
            }
        }
    }

    sub fin_fh {
        my $self = shift;
        $LOGGER->info("Closing file handle for $self->{outFile}...\n");
        close $self->{outFh}
          or $LOGGER->warn("Cannot close $self->{outFile}!\n");
        $LOGGER->info("Closing file handle for $self->{inFile}...\n");
        close $self->{inFh} or $LOGGER->warn("Cannot close $self->{inFile}!\n");
    }

    sub get_md {
        my $self = shift;
        my $tags = shift;
        my @md   = grep { /^MD/ } @$tags;    # should be one and only one MD tag
        return undef if $#md != 0;
        my $md = substr( $md[0], 5 );
        return $md;
    }

    sub iterate {
        my $self = shift;
        $self->{endType} eq "PE" ? $self->iterate_PE() : $self->iterate_SE();
    }

    sub iterate_PE {
        my $self = shift;
        my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
        my ( $readID, $flag, $chr, $pos, $cigar, $seq, $phred, $md, @F, @tags );
        my ( $read, $prev_read, $readPair );
        my $whichInPair;
        print $outFh (
            join( "\t",
                qw( ReadID R1_align_seq R1_align_ref R1_align_phred R1_split_S_seq R1_split_S_phred R1_split_I_seq R1_split_I_phred R1_split_=_seq R1_split_=_phred R1_split_X_seq R1_split_X_ref R1_split_X_phred )
            ),
            "\t",
            join( "\t",
                qw( R2_align_seq R2_align_ref R2_align_phred R2_split_S_seq R2_split_S_phred R2_split_I_seq R2_split_I_phred R2_split_=_seq R2_split_=_phred R2_split_X_seq R2_split_X_ref R2_split_X_phred )
            ),
            "\n"
        );
        while ( readline($inFh) ) {
            chomp;
            @F           = split( "\t", $_ );
            $readID      = $F[0];
            $whichInPair = $F[1] & 0x40 ? "R1" : "R2";
            $chr         = $F[2];
            $pos         = $F[3];
            $cigar       = $F[5];
            $seq         = $F[9];
            $phred       = $F[10];

   # SAM spec cannot guarantee a tag is always at the same column in each record
   # see htslib-1.9/sam.c sam_format1
            @tags = @F[ 11 .. $#F ];
            $md   = $self->get_md( \@tags );
            if ( !defined($md) ) {
                $LOGGER->warn("$readID has no MD tag! Skipping it...\n");
                next;
            }
            $read =
              Read->new( $readID, $whichInPair, $chr, $pos, $cigar, $seq,
                $phred, $md );
            $LOGGER->debug(
"(iterate_PE) $read->{readID}, $read->{whichInPair}, $read->{chr}, $read->{pos}, $read->{cigar}, $read->{seq}, $read->{phred}, $read->{md}\n"
            );
            if ( !defined( $prev_read->{readID} ) ) {
                $prev_read = $read;
                next;
            }
            if ( $prev_read->{readID} eq $read->{readID} ) {
                $readPair = ReadPair->new(
                    {
                        $prev_read->{whichInPair} => $prev_read,
                        $read->{whichInPair}      => $read
                    }
                );
                $readPair->parse();
                $readPair->output( $self->{outFh} );
                $prev_read = undef;
            }
            else {
                $readPair =
                  ReadPair->new( { $prev_read->{whichInPair} => $prev_read } );
                $readPair->parse();
                $readPair->output( $self->{outFh} );
                $prev_read = $read;
            }
        }
        if ( defined( $prev_read->{readID} ) ) {
            $readPair =
              ReadPair->new( { $prev_read->{whichInPair} => $prev_read } );
            $readPair->parse();
            $readPair->output( $self->{outFh} );
        }
    }

    sub iterate_SE {    #
        my $self = shift;
        my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
        my ( $readID, $flag, $chr, $pos, $cigar, $seq, $phred, $md, @F, @tags );
        my ( $read, $readPair );
        print $outFh (
            join( "\t",
                qw( ReadID R1_align_seq R1_align_ref R1_align_phred R1_split_S_seq R1_split_S_phred R1_split_I_seq R1_split_I_phred R1_split_=_seq R1_split_=_phred R1_split_X_seq R1_split_X_ref R1_split_X_phred )
            ),
            "\n"
        );
        while ( readline($inFh) ) {
            chomp;
            @F      = split( "\t", $_ );
            $readID = $F[0];
            $chr    = $F[2];
            $pos    = $F[3];
            $cigar  = $F[5];
            $seq    = $F[9];
            $phred  = $F[10];

   # SAM spec cannot guarantee a tag is always at the same column in each record
   # see htslib-1.9/sam.c sam_format1
            @tags = @F[ 11 .. $#F ];
            $md   = $self->get_md( \@tags );
            $LOGGER->debug(
"(iterate_SE) $read->{readID}, $read->{whichInPair}, $read->{chr}, $read->{pos}, $read->{cigar}, $read->{seq}, $read->{phred}, $read->{md}\n"
            );
            $read =
              Read->new( $readID, "R1", $chr, $pos, $cigar, $seq, $phred, $md );
            $readPair = ReadPair->new( { R1 => $read } );
            $readPair->parse();
            $readPair->output( $self->{outFh} );
        }
    }
}

my $samReader = SamReader->new( $inFile, $outFile, $endType, $nameSorted );
$samReader->init_fh();
$samReader->iterate();
$samReader->fin_fh();
$LOGGER->info("All done.\n");
