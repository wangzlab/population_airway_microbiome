#!/usr/bin/perl -w
# This scripts calculate AMHI by adopting a slightly modified scheme from that of Gupta and colleagues to calculate GMHI. 
# Both bacterial and fungal profiles are used and the thresholds for their fold-changes were simutaneously optimized.

open (IN, "bacteria_name_mapping");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$name{$a[0]}="B_".$a[1];
}
open (IN, "../data/bacteria_L6.txt");
$dump=<IN>;
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next unless (exists $name{$a[0]});
	for my $i (1..$#a) {
		$hash{$name{$a[0]}}{$headers[$i]}=$a[$i];
		$transhash{$headers[$i]}{$name{$a[0]}}=$a[$i];
		$bac{$headers[$i]}=1;
	}
}
open (IN, "fungi_name_mapping");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$name{$a[0]}="F_".$a[1];
}
open (IN, "../data/fungi_L6.txt");
$dump=<IN>;
$header=<IN>;
chop $header;
@headers=split("\t",$header);
while (<IN>) {
	chop;
	@a=split("\t",$_);
	next unless (exists $name{$a[0]});
	for my $i (1..$#a) {
		$hash{$name{$a[0]}}{$headers[$i]}=$a[$i];
		$transhash{$headers[$i]}{$name{$a[0]}}=$a[$i];
		$fun{$headers[$i]}=1;
	}
}

for my $key (keys %bac) {
	if (exists $fun{$key}) {
		$select{$key}=1; ## get bac fungi overlap samples
	}
}

open (IN, "../meta_airwayHealth_overallHealth.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	if ($a[1] == 1) {
		$rd{$a[0]}="Y";
	}
	if ($a[1] == 0) {
		$rd{$a[0]}="N";
	}
	if ($a[2] == 1) {
		$ad{$a[0]}="Y";
	}
	if ($a[2] == 0) {
		$ad{$a[0]}="N";
	}
}
for my $key (keys %hash) {
	$totalY=0;
	$numY=0;
	$totalN=0;
	$numN=0;
	$totalY2=0;
	$numY2=0;
	$totalN2=0;
	$numN2=0;	
	for my $key2 (keys %select) {
		if ($rd{$key2} eq 'Y') {
			$totalY += $hash{$key}{$key2};
			$numY++;
			if ($hash{$key}{$key2}>0) {
				#print $hash{$key}{$key2}."\n";
				$numY2++;
				$totalY2 += abs($hash{$key}{$key2}*log($hash{$key}{$key2}));
			}
		}
		if ($rd{$key2} eq 'N') {
			$totalN += $hash{$key}{$key2};
			$numN++;
			if ($hash{$key}{$key2}>0) {
				$numN2++;
				$totalN2 += abs($hash{$key}{$key2}*log($hash{$key}{$key2}));
			}
		}
	}
	$averY=$totalY/$numY;
	$averN=$totalN/$numN;
	$averY2=$totalY2/$numY * ($numY2/$numY);
	$averN2=$totalN2/$numN * ($numN2/$numN);
	$fc=log(($averY+0.000001)/($averN+0.000001))/log(2);
	$fc2=log(($averY2+0.000001)/($averN2+0.000001))/log(2);
	print $key."\t".$fc."\t".$fc2."\n";
	$fc{$key}=$fc;
	$fc2{$key}=$fc2;
}

### the above steps calculate fold change of each taxa in disease vs healthy ###

LOOP:for (my $i=0;$i<=5;$i+=0.05) {  # threshold for selecting cutoff
	LOOP:for (my $j=0;$j<=5;$j+=0.05) {
		%disease=();
		%health=();
		$m=0;
		$n=0;
		for my $key (keys %fc) { # obtain disease-enriched and healthy enriched taxa based on cutoff
			next if ($key =~ /F_/);
			if ($fc{$key}<=-1*$i) {
				$health{$key}=1;
				$m++;
				print STDERR $i."\t".$j."\t".$key."\tN\n";
			}	
			if ($fc{$key}>=$i) {
				$disease{$key}=1;
				$n++;
				print STDERR $i."\t".$j."\t".$key."\tY\n";
			}
		}
		for my $key (keys %fc) { # obtain disease-enriched and healthy enriched taxa based on cutoff
			next if ($key =~ /B_/);
			if ($fc{$key}<=-1*$j) {
				$health{$key}=1;
				$m++;
				print STDERR $i."\t".$j."\t".$key."\tN\n";
			}	
			if ($fc{$key}>=$j) {
				#next if ($key =~ /F_/);
				$disease{$key}=1;
				$n++;
				print STDERR $i."\t".$j."\t".$key."\tY\n";
			}
		}
		if ($m==0 or $n==0) {
			next LOOP;
		}
		$h=0;
		$hh=0;
		$d=0;
		$dd=0;
		#print "SampleID\tAirway_disease\tAll_disease\tGMHI\n";
		for my $key (keys %select) { # for each selected sample calculate GMHI
			next unless (exists $ad{$key});
			$totalY=0;
			$posY=0;
			$sumY=0;
			for my $key2 (keys %disease) { # disease-enriched taxa
				$totalY ++;
				if ($transhash{$key}{$key2}>0) {
					$posY++;
					$tmp=$transhash{$key}{$key2};
					$sumY+=abs($tmp*log($tmp));	
				}
			}
			$phiY=$posY/$totalY*$sumY; # calculate phi score of all disease-enriched taxa in each sample
			$totalN=0;
			$posN=0;
			$sumN=0;
			for my $key2 (keys %health) { # healthy-enriched taxa
				$totalN ++;
				if ($transhash{$key}{$key2}>0) {
					$posN++;
					$tmp=$transhash{$key}{$key2};
					$sumN+=abs($tmp*log($tmp));	
				}
			}
			$phiN=$posN/$totalN*$sumN; # calculate phi score of all healthy-enriched taxa in each sample
			$diff=log(($phiN+0.000001)/($phiY+0.000001))/log(10); # get GMHI (health over disease)
			print $key."\t".$ad{$key}."\t".$phiN."\t".$phiY."\t".$diff."\n";
			if ($ad{$key} eq 'Y') { # if a sample is a disease one
				$d++;
				if ($diff<0) { # if its GMHI < 0
					$dd++;
				}
			}
			if ($ad{$key} eq 'N') { # if a sample is a healthy one 
				$h++;
				if ($diff>0) { # if its GMHI > 0
					$hh++;
				}
			}
		}
		$accuracy=($dd/$d+$hh/$h)*0.5;
		#print $i."\t".$j."\t".$accuracy."\n";
	}
}
