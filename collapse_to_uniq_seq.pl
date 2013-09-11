#! /usr/bin/perl

#Tristan Lefebure, Licenced under the GPL v3
#tristan.lefebure@gmail.com

#2008-01-23, probably the first version
#5May2011: also make a table with name conversions
#1Sep2011: add noamb option, export haplo with the least ambiguity, and tble reports the ref haplotype

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use Bio::SeqIO;

my $format = 'fasta';
my $help= '';
my $noamb;

GetOptions(
	'help|?' => \$help,
	'format=s'  => \$format,
	'noamb' => \$noamb,
);

my $printhelp = "\nUsage: $0 [options] <alignment> <seq output> <convertion table output>\n\nWill collapse identical sequences. Ns and leading and ending gaps are not considered as differences. BioPerl is required.

Options:
	-h
		print this help
	-format: fasta|phylip|...
		give output format (default=fasta)
	-noamb, will convert all the ambiguity code into missing data (only for DNA data)\n";

if ($#ARGV<2) {
	print $printhelp;
	exit;
}

my $in  = Bio::AlignIO->new(-file   => $ARGV[0] ,
                             -format => $format);

# my $out = Bio::AlignIO->new(-file   => ">$ARGV[1].old" ,
#                              -format => $format);

my $out2 = Bio::SeqIO->new(-file   => ">$ARGV[1]" ,
                             -format => 'fasta');


open TAB, ">$ARGV[2]";
print TAB "ST\tSequence\tChosenHaplotype\n";

while ( my $aln = $in->next_aln() ) {
	#transform the ambiguity code, the leading and ending gaps to ?
	#make this into a new aligment
	my $tmpAln = Bio::SimpleAlign->new();
	if ($noamb) {
	    foreach my $seq ( $aln->each_seq() ) {
		my $str = $seq->seq();
		$str =~ s/[nkmbvswdyrh]/\?/gi;
		my $new = Bio::LocatableSeq->new(
			    -id      => $seq->id(),
			    -alphabet=> $seq->alphabet,
			    -seq     => $str,
			    -start   => $seq->start,
			    -end     => $seq->end
			);
		$tmpAln->add_seq($new);
	    }
	}
	else {
	    $tmpAln = $aln;
	}

	my ($red_aln, $st) = $tmpAln->uniq_seq;
#         $out->write_aln($red_aln);



	my %selectSeq;
	#export the sequences that have the least ambiguity code
	foreach my $stn (sort { $a <=> $b } keys %$st) {
	    my @tax = @{$st->{$stn}};
	    my $bestseq;
	    my $leastamb = 999999;
	    #find the best sequence
	    foreach my $t (@tax) {
		my $seq = $aln->get_seq_by_id($t);
		my $str = $seq->seq;
		my $count = 0;
		#count ambiguous code
		while($str =~ /[nkmbvswdyrh\?]/ig) { ++$count };
		#add the leading and ending gaps
		if($str =~ /^(-+)/ig) { $count += length($1) };
		if($str =~ /(-+)$/ig) { $count += length($1) };
		if($count < $leastamb) {
		    $bestseq = $seq;
		    $leastamb = $count;
		}
	    }
	    #print out the best sequence
	    print "The best seq for group $stn is ", $bestseq->id, "\n";
	    $selectSeq{$stn} = $bestseq->id;
	    $out2->write_seq($bestseq);
	}

	foreach my $stn (sort { $a <=> $b } keys %$st) {
	    my @tax = @{$st->{$stn}};
	    foreach (@tax) {
		print TAB "$stn\t$_\t$selectSeq{$stn}\n";
	    }
	}
}


