# seq_anchoring
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
