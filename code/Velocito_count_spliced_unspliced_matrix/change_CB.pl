#! /usr/bin/env perl
#
# Short description for change_CB.pl
#
# use strict;
use Getopt::Long;
use warnings;
my %hash;
my @a;
my $aa;
my @b;
my @c;
my @d;
my ($config,$bam);

GetOptions(
        "i=s" => \$config,
        "b=s" => \$bam,
);

open IN,$config;
while(<IN>){
	chmod;
	@a = split(/,/,$_);
	$aa = $a[0].'_'.$a[1];
	$hash{$aa} = $a[2];
         
}
#$bam=shift;
open BAM,"samtools view -h $bam|";
while(<BAM>){
	#chmod;
        if(/^@/){
                print $_;
        }
        else{
		#chmod();
		@d =split(/\n/,$_);
		@b = split(/\t/,$d[0]);
		@c = split(/CB:Z:/,$b[15]);
		if(exists($hash{$c[1]})){
                #$_=~s/\tCB:Z:([^\t]*)\t/\tCB:Z:$hash{$c[1]}\t/;
                print"$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b[5]\t$b[6]\t$b[7]\t$b[8]\t$b[9]\t$b[10]\t$b[11]\t$b[12]\t$b[13]\t$b[14]\t$b[16]\t$b[17]\tCB:Z:$hash{$c[1]}"};
          
        }
}
