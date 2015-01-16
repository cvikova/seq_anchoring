#!/usr/bin/perl -w

use warnings;
use Getopt::Long;
use File::Basename;
use vars qw($help $mtp_file $fpc_file $match $pool_list $cov_folder $sequence_list $prefix); 


&GetOptions ( "h"      => \$help,
	      "mtp:s"  => \$mtp_file,
	      "fpc:s"  => \$fpc_file,
	      "match:s"  => \$match_file,	      
	      "pool:s"  => \$pool_list,
	      "cov:s"  => \$cov_folder,
	      "seq:s"  => \$sequence_list,	
	      "alp:i"    => \$aln_perc,
	      "p:s" => \$prefix );

$help and &help;
$mtp_file or print "\n!!!MISSING FILE WITH MTP CLONES!!!\n" and  &help;
$fpc_file or print "\n!!!MISSING FPC FILE!!!\n" and  &help;
$pool_list or print "\n!!!MISSING FILE WITH LIST OF POOLS!!!\n" and  &help;
$cov_folder or print "\n!!!MISSING FOLDER WITH MosaikCoverage RESULTS!!!\n" and  &help;
$sequence_list or print "\n!!!MISSING FILE WITH LIST OF SEQUENCES!!!\n" and  &help;

my $alignment_percentage_treshold;
my $alignment_coverage_treshold;

$aln_perc and $alignment_percentage_treshold = $aln_perc or print "\nsetting alignment length treshold to 80%\n" and $alignment_percentage_treshold = 80;

$prefix or print "\nsetting prefix of output files to analysis\n" and $prefix = 'analysis';



my %positive_pools;				# hash made from marker names and arrays of corresponding positive pools


open (POOL, $pool_list);
while (my $pool = <POOL>) {
	chomp $pool;
	-e "$cov_folder/${pool}_coverage.txt" or die "File $cov_folder/${pool}_coverage.txt does not exist!";
	open (COV, "$cov_folder/${pool}_coverage.txt");

	while (my $entry = <COV>) {
		if ($entry =~ /coverage statistics/) {
			$entry =~ /coverage statistics for (\S+) \((\d+) bp\): (\d+).+mean/;
			my $seq_name = $1;
			my $seq_length = $2;
			my $align_lenght = $3;

			my $alignment_percentage = 100*$align_lenght/$seq_length;
			if ($alignment_percentage >= $alignment_percentage_treshold){

				push (@{$positive_pools{$seq_name}}, "$pool");

			}
		} 				# end of porcessing of one line in coverage file

		else {next;}

	} 					# end of processing of complete coverage file for one pool 

} 						# end of positive pool deteciton
close POOL;
# end feeding of hash %positive_pools

open (OUT, ">$prefix.results.out");

print OUT "\"Status\" column refer to how the marker was assigned to clone(s)\n\"anchor_type_1\" - single positive pool in each dimmension (plate, row, column)\n\"anchor_type_2\" - multiple positive pools, two or more clones among candidate share same contig of the physical map\n\"anchor_type_3\" - multiple positive pools, two or more clones have a hit in End2End analysis in FPC, if used\n\"no_clone_overlap\" - CANDIDATE instead of POSITIVE clones are listed in column 7\n\n";
print OUT "marker\tstatus\tpositive_plate_pools\tpositive_row_pools\tpositive_column_pools\tpositive_pool_list\tpositive_clones";


open (FPC, $fpc_file);				# load FPC file
my @fpc = <FPC>;				#
close (FPC);					#
my $f_map = join ("", @fpc);			#



open (SEQ, $sequence_list);			# open sequence list

while (my $seq = <SEQ>) {
	chomp $seq;

	if (exists($positive_pools{"$seq"})) {


		($pos_plates_count, $pos_rows_count, $pos_columns_count, @candidate_clones) = BAC_candidates(@{$positive_pools{$seq}});


		if ($pos_plates_count == 0 or $pos_rows_count == 0 or $pos_columns_count == 0){
			print OUT "\n$seq\tno_positive_BAC\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count\t@{$positive_pools{$seq}}";
		}

		elsif ($pos_plate_count == 1 and $pos_row_count == 1 and $pos_column_count == 1){
			print OUT "\n$seq\tanchor_type_1\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count\t@{$positive_pools{$seq}}\t$candidate_clones[0]";
		}


		elsif ($pos_plate_count > 5 or $pos_row_count > 5 or $pos_column_count > 5){
			print OUT "\n$seq\tto_many_positive_pools\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count";
		}

		else { 

		
			($match, @positive_BACs) = contig_BACs ($f_map, @candidate_clones); # sub - BAC identification - same contig

			if ($match == 1) {
			print OUT "\n$seq\tanchor_type_2\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count\t@{$positive_pools{$seq}}\t@positive_BACs";
			}
			
			else { 
 				$match_file or next;

				
				($pos, @positives) = matching_BACs ($match_file, @candidate_clones); # sub - BAC identification - matching_clones

				if ($pos == 1) {
					print OUT "\n$seq\tanchor_type_3\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count\t@{$positive_pools{$seq}}\t@positives";
				}

				else {
					print OUT "\n$seq\tno_clone_overlap\t$pos_plates_count\t$pos_rows_count\t$pos_columns_count\t@{$positive_pools{$seq}}\t\t@candidate_clones";
				}

			}
		}

	}

	else {
		print OUT "\n$seq\tno_positive_pool";
	}

}

exit(0);


###############
###  BAC candidate selection   ###
###############
sub BAC_candidates {

	my (@positive_pools) = @_;

	my @plate;
	my @row;
	my @column;
	foreach $pool (@positive_pools) {

		if ($pool =~ /p/) {
			push (@plate, $pool);
		}

		if ($pool =~ /r([A-P])/) {push (@row, $1);}

		if ($pool =~ /c(\d{2})/) {push (@column, $1);}
	}

	$pos_plate_count = scalar@plate;
	$pos_row_count = scalar@row;
	$pos_column_count = scalar@column;

	my $p;
	my $r;
	my $c;

	open (MTP, $mtp_file);				# load MTP adressess of clones
	my @mtp_clones = <MTP>;				#
	close (MTP);					#
	my $addresses = join ("", @mtp_clones);		#

	my @candidate_BACs;

	for ($p = 0; $p < $pos_plate_count; $p++) {
		for ($r = 0; $r < $pos_row_count; $r++) {
			for ($c = 0; $c < $pos_column_count; $c++) {
				my $clone = $plate[$p] . $row[$r] . $column[$c];
				chomp $clone;

				if ($addresses =~ /$clone\t(\w+)/i) {
					my $fpc_name = $1;
					push (@candidate_BACs, $fpc_name);
				}
			}
		}
	}

return ($pos_plate_count, $pos_row_count, $pos_column_count, @candidate_BACs);

}


###############
###  BAC identification - same contig   ###
###############

sub contig_BACs {
	my ($fpc_data, @candidates) = @_;

	my $number_of_clones = scalar(@candidates);
	my $n;
	my @clones;
	my @contigs;

	for ($n = 0; $n < $number_of_clones; $n++) {

		if ($fpc_data =~ /BAC\s+:\s+"$candidates[$n]"\s+Map\s+"(ctg\d+)"/){

			my $contig = $1;
			push (@clones, $candidates[$n]);
			push (@contigs, $contig);
		}

		else {print "contig for $candidates[$n] not found\n";}


	}

	my $i;
	my $j;
	my $match = 0;
	my @positive_clones;

	for ($i = 0; $i < $number_of_clones-1; $i++) {

		for ($j = $i+1; $j < $number_of_clones; $j++) {

			if ($contigs[$i] eq $contigs[$j]) {
				$match = 1;
				push (@positive_clones, $clones[$i], $clones[$j]);
			}
		}
	}

	my @positive_BACs;

	if ($match == 1){
		sort @positive_clones;

		my $i;
		my $j;

		push (@positive_BACs, $positive_clones[0]);
		for ($i = 0; $i < scalar(@positive_clones)-1; $i++) {

			unless ($positive_clones[$i] eq $positive_clones[$i+1]) {
				push (@positive_BACs,$positive_clones[$i+1]);
			}
		}


	}
	
return ($match, @positive_BACs);

} #end of sub



###############
###  BAC identification - matching_clones   ###
###############

sub matching_BACs {
	my ($match_file, @candidates) = @_;

	my $number_of_clones = scalar(@candidates);
	my $match = 0;
	my @positive_clones;

open (MATCH, $match_file);

	while (my $line = <MATCH>) {
		my $i;
		my $j;
		for ($i = 0; $i < $number_of_clones-1; $i++) {

			if ($line =~ /$candidates[$i]/) {

				for ($j = $i+1; $j < $number_of_clones; $j++) {
					if ($line =~ /$candidates[$j]/) {

						$match = 1;
						push (@positive_clones, $candidates[$i], $candidates[$j]);
						print $line;

					}
				}
			}
		}
	}


	my @positive_BACs;

	if ($match == 1){
		sort @positive_clones;

		my $i;
		my $j;

		push (@positive_BACs, $positive_clones[0]);
		for ($i = 0; $i < scalar(@positive_clones)-1; $i++) {

			unless ($positive_clones[$i] eq $positive_clones[$i+1]) {
				push (@positive_BACs,$positive_clones[$i+1]);
			}
		}


	}

close MATCH;
return ($match, @positive_BACs);

} #end of sub



###############
###  HELP   ###
###############
sub help() {

my $prog = basename($0);
print STDERR <<EOF ;

$prog 
--------------------------------------
This script is used for high-througput anchoring of sequences in silico to physical map based on sequencing of 3-dimensional pools and Read mapping to sequences.  

USAGE: 
       $prog  [OPTIONS]
			
EXAMPLE:
       $prog  -mtp  MTP_clones.txt  -fpc 3DS.fpc  -match End2End.txt -pool Pool_list.txt -cov Zipper_coverage/ -seq Sequence_list.txt -alp 75 -alc 3 -p Zipper_pct_75_cov_3

OPTIONS:

       -h       print this help

       -mtp   file with MTP address of each clone and corresponding address in BAC library

			[MTP address]	[BAC library address(in format used by FPC)]
			------------------------------------------------------------
			p01A01		TaaCsp3DS001A20
			p01C11		TaaCsp3DS001B13
			....		...............	
			....		...............
			p10N20		TaaCsp3DS096P10
	
       -fpc   physical map *.fpc file

       -match	log file with matching clones from fpc

       -pool	file with pool list, each pool on separate line (pool names must correspond to MTP addresses in "mtp" file; plate pools must start with p, the number must reflect numbering of plates in MTP addresses; row pools must start with r followed by letter of row; column pools must start with c followed by two digits (if the mtp address if p03A12, corresponding pools should be p03, rA and c12) 

       -cov	folder with coverage files, file name must be in format [pool_name]_coverage.txt (i.e. c05_coverage.txt for column pool c05; name must be same with pool name in "pool" file) 	 

       -seq	file with marker list (same names as in coverage files)

       -alp	alignment percentage threshold (default 80)

       -p	prefix of output file

EOF

exit(1);
}
