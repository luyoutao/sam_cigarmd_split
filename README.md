## Summary:
Split query sequence and phred scores into softclip (S), insertion (I), match (=) and mismatch (X) accoriding to CIGAR and tag MD.

## Usage:

    perl Source/functions/sam_cigarmd_split.pl --inFile|-i <input.bam|sam> [--outFile|-o <output.txt>] [--endType PE] [--debug info] [--version|-v] [--help|-h] [--ncores 1]

## Options:

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
         

## Examples:
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

## Extension: 
Given the splits, we can easily compute stats such as average phred score in matches and mismatches, or base composition difference between soft-clipping and insertion. For example, 

\*) to extract the triplet of aligned query, ref. and phred:

    perl sam_cigarmd_split.pl -i test.bam 2> err.log | awk 'NR > 1 {print $2 "\n" $3 "\n" $4 "\n"}' > test_align.txt

\*) to calculate the average phred scores in matches and mismatches:

    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ use List::Util "sum" }; @F = split("\t", $_, -1); $phred = $F[9] . $F[21]; $phred =~ s/,//g; @p = map{ ord($_) - 33 } split("", $phred); print $F[0], "\t", $#p > -1 ? sum(@p) / ($#p+1) : ""' > test_avgphred_match.txt
    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ use List::Util "sum" }; @F = split("\t", $_, -1); $phred = $F[12] . $F[24]; $phred =~ s/,//g; @p = map{ ord($_) - 33 } split("", $phred); print $F[0], "\t", $#p > -1 ? sum(@p) / ($#p+1) : ""' > test_avgphred_mismatch.txt

\*) to calculate the base composition of mismatched ref.:

    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -lne 'BEGIN{ %comp }; @F = split("\t", $_, -1); $x = $F[11] . $F[23]; $x =~ s/,//g; @b = split("", $x); for $n (@b) { $comp{uc($n)}++ }; END{ for $k (sort {$comp{$b} <=> $comp{$a}} keys(%comp)) { print $k, "\t", $comp{$k}} }'

\*) to calculate the mutation transition probability of the mismatches:

    perl sam_cigarmd_split.pl -i test.bam 2> err.log | sed 1d | perl -ne 'BEGIN{ %comp; print "\t" . join("\t", qw(A C G T)), "\n" }; @F = split("\t", $_, -1); $x = $F[11] . $F[23]; $x =~ s/,//g; $y = $F[10] . $F[22]; $y =~ s/,//g; @b = split("", $x); @m = split("", $y); for $i (0..$#b) { $comp{uc($b[$i])}{uc($m[$i])}++ }; END{ for $r (qw(A C G T)) { print $r; for $c (qw(A C G T)) { print "\t", $comp{$r}{$c} // 0 }; print "\n" } }'
