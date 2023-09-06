#!/usr/bin/env perl

# Recording the position where the 3' end of each read is mapped as a WIG file
# 
# usage:
# perl countReads3Edge_calibrated.readscount.pl yeast.Ribo.seq_calibrated.bam > inforam_trans.bed positive_strand_info_trans.wig negative_strand_info_trans.wig

use strict;

my $bam_file = $ARGV[0];
my %read_3pos_plus;
my %read_3pos_minus;
my %max_pos;
my $offset = 0;
my $strand;
my $readscounts= 0;

#print "$offset";
open(my $fp, "-|", "samtools view $bam_file") || die;
while(<$fp>) {
    chomp;
    my @line = split(/\t/, $_);
    my $refname = $line[2];
    my $ciger = $line[5];
    my $flags = $line[1];
    my $pos = $line[3];
    my $seq = $line[9];
    my $seqlen = length($seq);

    next if ($flags & 0x4); # unmapped

    my $dir = 1;
    $dir = -1 if ($flags & 0x10);

# The original codes only use reads with exact match, which maybe too stringent
#    next if ($ciger ne "${seqlen}M"); # partial match

# This codes allow reads with soft clipping at 3' end
    my $match_len = 0;
    my $softclip_len = 0;
    if ($ciger eq "${seqlen}M") {
#	print STDERR "read with exact match\n";
	$match_len = $seqlen;
	$readscounts ++;
	next if ($match_len < $offset);
    }
    elsif ($dir == 1 && $ciger =~ /(\d+)M(\d+)S/) {
#	print STDERR "read on plus strand with 3' soft-clipping\n";
	$match_len = $1;
	$softclip_len = $2;
	$readscounts ++;
	next if ($match_len + $softclip_len != $seqlen || $match_len< $offset);
    }
    elsif ($dir == -1 && $ciger =~ /^(\d+)S(\d+)M$/) {
#	print STDERR "read on minus strand with 3' soft-clipping\n";
	$match_len = $2;
	$softclip_len = $1;
	$readscounts ++;
	next if ($match_len + $softclip_len != $seqlen || $match_len< $offset);
    }
    else {
	next;
    }

    my $rib_last_pos;
    if ($dir == 1) {
	$rib_last_pos = $pos + $seqlen - 1 - $softclip_len;
	#print "$pos\n";
	#print "$seqlen\n";
	#print "$softclip_len\n";
        #print "$rib_last_pos\n";
	$rib_last_pos = $rib_last_pos - $offset;
	#print "$rib_last_pos\n";
	$read_3pos_plus{$refname}->[$rib_last_pos]++;
    } else {
	$rib_last_pos = $pos + $softclip_len;
	$rib_last_pos = $rib_last_pos + $offset;
	$read_3pos_minus{$refname}->[$rib_last_pos]++;
    }
    #$read_3pos{$refname}->[$rib_last_pos]++;
    $max_pos{$refname} = $rib_last_pos if ($max_pos{$refname} < $rib_last_pos);
}

print STDERR "loading finished\n";
#print STDERR "max position: ", $max_pos{"chr"}, "\n";
print STDERR "total reads counts:",$readscounts,"\n";

# BED format 
# chr1  start  stop  some_value   strand ....
my %read_3pos = (%read_3pos_plus, %read_3pos_minus);
open(my $ibf, '>inforam_trans.bed');
foreach my $refname (sort keys %read_3pos) {
  for (my $i = 1; $i < $max_pos{$refname}; $i++) {
    if ($read_3pos_plus{$refname}->[$i] > 0) {
        print $ibf join("\t", $refname, $i, $i+1, "SomeDummyName", $read_3pos_plus{$refname}->[$i], "+") , "\n";
    }
    if ($read_3pos_minus{$refname}->[$i] > 0) {
        print $ibf join("\t", $refname, $i, $i+1, "SomeDummyName", $read_3pos_minus{$refname}->[$i], "-") , "\n";
    }
}
}
close $ibf;


# --> inforam.bed

# WIG format
# positive
open(my $pwig, '>positive_strand_info_trans.wig');
foreach my $refname (sort keys %read_3pos_plus) {
    for (my $i = 1; $i < $max_pos{$refname}; $i++) {
	if ($read_3pos_plus{$refname}->[$i] > 0) {
	    print $pwig join("\t", $refname, $i, $read_3pos_plus{$refname}->[$i]) , "\n";
	}
    }
}
close $pwig;

# --> positive_strand_info

# WIG format
# minus
open(my $nwig, '>negative_strand_info_trans.wig');
foreach my $refname (sort keys %read_3pos_minus) {
    for (my $i = 1; $i < $max_pos{$refname}; $i++) {
	if ($read_3pos_minus{$refname}->[$i] > 0) {
	    print $nwig join("\t", $refname, $i, $read_3pos_minus{$refname}->[$i]) , "\n";
	}
    }
}
close $nwig
# --> negative_strand_info
