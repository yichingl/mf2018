#!/usr/bin/perl -w


if (@ARGV != 4) {
	die "Usage: <dir> <dim> <mix> <covtyp>\n";
}

$dir	= shift @ARGV;
$dim	= shift @ARGV;
$mix	= shift @ARGV;
$covtyp	= shift @ARGV;

if ($covtyp == 1) {
	$file = "$dir/g.d$dim.m$mix.f.pro";
}
else {
	$file = "$dir/g.d$dim.m$mix.d.pro";
}

open OUT, ">$file" || die $!;
print OUT "<NAME> gmmproto USER\n";
if ($covtyp == 1) {
	print OUT "<VECSIZE> $dim <FULL>\n";
}
else {
	print OUT "<VECSIZE> $dim <DIAG>\n";
}

print OUT "<NUMMIXES> $mix\n";
$r = 1.0/$mix;
for ($i = 0; $i < $mix; $i++) {
	printf OUT "\n<MIXTURE> $i %.3f\n", $r;
	print OUT "<MEAN>\n";
	for ($j = 0; $j < $dim; $j++) {
		print OUT "0.0\t";
	}
	print OUT "\n";
	print OUT "<VARIANCE>\n";
	for ($j = 0; $j < $dim; $j++) {
		for ($k = $j; $k < $dim; $k++) {
			if ($k == $j) {
				print OUT "1.0\t";
			}
			else {
				print OUT "0.0\t";
			}
		}
		print OUT "\n";
	}
	print OUT "\n";
}
print OUT "<ENDGMM>\n";
close OUT;
