#!/usr/bin/env perl

use warnings;
use strict;

if ($#ARGV !=0)
{
	print "use: karplus_couplings.pl angaver.xvg\n\nwhere 'angaver.xvg' contains the average angle per time, as produced by 'gmx angle'\n";
	exit;
}
my @time=();
my @angle=();
my @j_coupling=();

my $averJcoup=0;
my $stderrJcoup=0;

open XVG, $ARGV[0];

while (<XVG>)
{
	if ($_=~/^(\d+\.\d+)\s+(\S+)/)
	{
		push @time, $1;
		push @angle, $2;
	}
}

close XVG;

for (my $i=0; $i<=$#angle; $i+=1)
{
	my $Jcoup=7.49*cos($angle[$i])**2 -0.96*cos($angle[$i])+0.15;
	push @j_coupling, $Jcoup;	
	$averJcoup+=$Jcoup;
}

$averJcoup=$averJcoup/($#angle+1);

for (my $i=0; $i<=$#j_coupling; $i+=1)
{
	$stderrJcoup+=($averJcoup-$j_coupling[$i])**2
}
$stderrJcoup=sqrt($stderrJcoup/(($#j_coupling+1)*$#j_coupling));

print "Aver: $averJcoup +/- $stderrJcoup Hz\n";
