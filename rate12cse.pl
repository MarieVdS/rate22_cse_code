#! /usr/bin/perl 
# C-rich CSE ODE writer for Rate12 (Feb 2012)
# Works only with with McElroy colon-separated rate12 file format (including multipl temperature ranges)
# Based on original code by Andrew J Markwick June 2000
# Writes ODEs to filename specified after -o (eg. -o odes.f)
# Writes species list to <ratefile_name>.specs, including conserved species and parents
# Parent abundances (with respect to H2) are specified at the bottom of this script

$version = "rate10odecsT.pl";
$date="Feb 2012";
$verbose=1;
$showspecs=0;
$showelems=0;
$cons{"e-"}=1;
$cons{"H2"}=2;
@csp=("e-","H2");
@dsp=();
$nc=2;
$rfile="";
$output="";

# parse options

while ($_ = shift(@ARGV)) {
    if (/\.rates$/) {
	$rfile=$_;
    }
    elsif (/^-o$/) {
	$output=shift;
    }
    elsif (/^-c=([a-zA-Z0-9,+-]+)$/) {
	@t=split(/,/,substr($_,3));
	foreach $f (@t) {
	    push(@csp,$f);
	    $cons{$f}=$nc;
	    $nc++;
	}
    }
    elsif (/^-q$/ || /^--quiet$/) {
	$verbose=0;
    }
    elsif (/^-s$/ || /^--species$/) {
	$showspecs=1;
    }
    elsif (/^-e$/ || /^--elements$/) {
	$showelems=1;
    }
    elsif (/^--help$/) { 
	&usage && die "\n"; 
    }
    elsif (/^--version$/) { 
	die "$0: rate version $version ($date)\n"; 
    }
    else {
	print "$0: [safe] unrecognized switch: $_\n";
    }
}

$rfile =~ m/\.[^.]+$/;
$filename = $`;

if ($rfile eq "") { &usage && die "$0: [serious] no ratefile specified\n"; }
if ($output eq "") { die "$0: [serious] no output file specified\n"; }

if ($verbose) { print "$version ($date)\nA. J. Markwick and M. A. Cordiner\nFor use with Rate10 carbon-rich CSE model\n\n"; }

# try to read ratefile

@lines=();
open(RATE,"<$rfile") || die "$0: [serious] file not found: $rfile\n";
if ($verbose) { print "$0: reading ratefile $rfile\n"; } 
while (<RATE>) {
    push(@lines,$_);
}
close(RATE);

$nreacs=$#lines+1;

$ct=1;
$elems="";
$ec=0;
$nc--;

# extract species

foreach (@lines) {
    @data=split(/\:/);
    @temp=($data[2], $data[3], $data[4], $data[5], $data[6], $data[7] );
    foreach $item (@temp) {
	if (exists $specs{$item}) {
	} 
	else {
	    unless ($item eq "" || $item eq "PHOTON" || $item eq "CRPHOT" || 
        $item eq "INPHOTON" || $item eq "ACPHOTON" ||
		    $item eq "M" || $item eq "CRP" || exists $cons{$item}) {
		$specs{$item}=$ct;
		$ct++;
		push(@dsp,$item);
	    }
	}
    }
}


$nsp=$ct-1;

# mass array

%mel= (
       "H" => 1,
       "He" => 4,
       "C" => 12,
       "T" => 13,
       "N" => 14,
       "O" => 16,
       "X" => 18,
       "F" => 19,
       "Na" => 23,
       "Mg" => 24,
       "Al" => 27,
       "Si" => 28,
       "P" => 31,
       "S" => 32,
       "Cl" => 35,
       "Ca" => 40,
       "Ar" => 40,
       "Ti" => 48,
       "Fe" => 56
       );

# extract elements

foreach $key (keys %cons) {
    if ($key =~ m/^[A-Z]$/) {
	$ec++;
	$elems=$elems.$key." :";
	$els{$key}=0;
    }
    elsif ($key =~ m/^[A-Z][a-z]$/) {
	$ec++;
	$elems=$elems.$key.":";
	$els{$key}=0;
    }
}
foreach $key (keys %specs) {
    if ($key =~ m/^[A-Z]$/) {
	$ec++;
	$elems=$elems.$key." :";
	$els{$key}=0;
    }
    elsif ($key =~ m/^[A-Z][a-z]$/) {
	$ec++;
	$elems=$elems.$key.":";
	$els{$key}=0;
    }
} 

foreach $key (keys %els) {
    unless (exists $mel{$key}) { die "$0: [serious] mass of $key unknown\n"; }
}

# check conserved species are elements or H2

foreach $c (keys %cons) {
    if (length($c) < 3) {
	if ((length($c) == 1 && index($elems,$c." :") eq -1) || 
	    (length($c) == 2 && index($elems,$c.":") eq -1)) {
	    unless ($c eq "H2" || $c eq "e-" ) {
		die "$0: [serious] can only conserve elements or H2\n";
	    }
	}
    }
    else {
	print $c;
	die "$0: [serious] can only conserve elements or H2\n";
    }
}

# report statistics

if ($verbose) { 
    print "$0: ratefile $rfile - $nreacs reactions, ";
    print $nsp+$nc." species, $ec elements\n"; 
}
if ($showspecs && $verbose) { 
    print "$0: species  - "; 
    foreach $key (keys %specs) {
	print $key.":";
    } 
    print "\n"; 
}
if ($showelems && $verbose) { print "$0: elements - $elems\n"; }
$elems=$elems."+ :- :";
$ec=$ec+2;

# decompose species into element array + work out masses

if ($verbose) { print "$0: decomposing species ... \n"; }

foreach $sc (keys %specs) {
    @dc=(0) x $ec;
    $i=0;
    while ( $i < length($sc) ) {
	$te=substr($sc,$i,1);
	$tt=substr($sc,$i+1,1);
	$tu=substr($sc,$i+2,1);
	if (index($elems,$te.$tt.":") eq -1) {
# not a 2 char elem
	    $p=index($elems,$te." :")/3;
	    if ($tt =~ m/[0-9]/) {
# second char is number
		if ($tu =~ m/[0-9]/) {
# char 3 is number
		    $dc[$p]=$dc[$p]+($tt*10)+$tu;
		    $i=$i+3;
		    $masses{$sc}=$masses{$sc}+($mel{$te}*(($tt*10)+$tu));
		}
		else {
		    $dc[$p]=$dc[$p]+$tt;
		    $i=$i+2;
		    $masses{$sc}=$masses{$sc}+($mel{$te}*$tt);
		}
	    }
	    else {
# second char is not number
		$dc[$p]++;
		$i++;
		$masses{$sc}=$masses{$sc}+$mel{$te};
	    }
	}
       else {
# 2 char element
	   $p=index($elems,$te.$tt.":")/3;
	   $tu=substr($sc,$i+2,1);
	   if ($tu =~ m/[0-9]/) {
# char 3 number
	       $dc[$p]=$dc[$p]+$tu;
	       $i=$i+3;
	       $masses{$sc}=$masses{$sc}+($mel{$te.$tt}*$tu);
	   }
	   else {
# char 3 not a number
	       $dc[$p]++;
	       $i=$i+2;
	       $masses{$sc}=$masses{$sc}+$mel{$te.$tt};
	   }
       }
    }
    $decomp{$sc} = [ @dc ];
    $nat{$sc}=&total(@dc);
    if ( $sc =~ m/[+-]/ ) { $nat{$sc}--; }
}

@sorted = sort { $masses{$a} <=> $masses{$b} } (keys %masses);
@sortat = sort { $nat{$a} <=> $nat{$b} } (keys %nat);
$ix=1;
$tmx=0;
foreach $it ( @sorted ) {
    $specs{$it}=$ix;
    $ix++;
}

if ($verbose) { print "$0: biggest species is ".$sortat[$#sortat]." (".
		    $nat{$sortat[$#sortat]}.")\n"; }
if ($verbose) { print "$0: heaviest species is ".$sorted[$#sorted]." (".
		    $masses{$sorted[$#sorted]}.")\n"; }

# evaluate algebraic terms

if ($verbose) { print "$0: preparing conservation terms ... \n"; }

foreach $cn (keys %cons) {
    if ($cn eq "e-") {
	$cterm{$cn}="+TOTAL(".$cons{$cn}.")+(";
	$p=index($elems,"+ :")/3;
    } 
    elsif ($cn eq "H2") {
	$cterm{$cn}="+TOTAL(".$cons{$cn}.")-0.5*(";
	$p=index($elems,"H :")/3;
    }
    else {
	$cterm{$cn}="+TOTAL(".$cons{$cn}.")-(";
	if (length($cn) == 1) {
	    $p=index($elems,$cn." :")/3;
	}
	else {
	    $p=index($elems,$cn.":")/3;
	}
    }
    foreach $sp ( @sorted ) {
	if ( $sp ne $cn ) {
	    $djp = $decomp{$sp}[$p];
	    if ($cn eq "e-" && $djp == 0 && $decomp{$sp}[$p+1] > 0) {
		$xi="-";
		$djp=$decomp{$sp}[$p+1];
	    }
	    else {
		$xi="+";
	    }
	    if ( $djp > 0 ) {
		if ( $djp > 1 ) {
		    $cterm{$cn}=$cterm{$cn}.$xi.$djp."*Y(".$specs{$sp}.")";
		}
		else {
		    $cterm{$cn}=$cterm{$cn}.$xi."Y(".$specs{$sp}.")";
		}
	    }
	}
    }
    $cterm{$cn}=$cterm{$cn}.")";
    
}

$ct=1;
$noch=0;

if ($verbose) { print "$0: evaluating formation/destruction terms ... \n"; }

# evaluate differential terms

foreach (@lines) { 
    @data=split(/\:/);
    @temp=( $data[2], $data[3], $data[4], $data[5], $data[6], $data[7] );
    $flag = $data[1];
#     print "$flag\n";
    for $i ( 0 .. 1 ){
	$te=$temp[$i];
	$tt=$temp[1-$i];
	unless ($specs{$te} eq "" || exists $cons{$te}) {
	    $ex="+K(".$ct.")";
	    unless ($specs{$tt} eq "" && not exists $cons{$tt}) {
		if (exists $cons{$tt}) {
		    $ex=$ex."*X(".$cons{$tt}.")*HNR";
		}
		else {
		    $ex=$ex."*Y(".$specs{$tt}.")*HNR";
		}
	    }
	    
	    if ($specs{$tt} eq "" && $flag eq "FO") {
		    $ex=$ex."*HNR";
		    
		}
	    if ($specs{$tt} eq "" && $flag eq "TD") {
		    $ex=$ex."*HNR";
		}
# 		
	
	    $dest{$te}=$dest{$te}.$ex;
	    $noch=$noch+length($ex);
	    if ((split(/ /,$temp[2]))[0] eq "M" || (split(/ /,$temp[1]))[0] eq "M") {
		$dest{$te}=$dest{$te}."*DN";
		$noch=$noch+3;
	    }
	}
    }
    for $i ( 2 .. 5 ){
	$te=$temp[0];
	$tt=$temp[1];
	$tu=$temp[$i];
	unless ($tu eq "M" || $specs{$tu} eq "" || exists $cons{$tu}) {
	    $ex="+K(".$ct.")";
	    if (exists $cons{$te}) {
		$ex=$ex."*X(".$cons{$te}.")";
	    }
	    else {
		$ex=$ex."*Y(".$specs{$te}.")";
	    }
	    unless ($specs{$tt} eq "" && not exists $cons{$tt}) {
		if (exists $cons{$tt}) {
		    $ex=$ex."*X(".$cons{$tt}.")*HNR";
		}
		else {
		    $ex=$ex."*Y(".$specs{$tt}.")*HNR";
		}
	    }
	    if ($specs{$tt} eq "" && $flag eq "FO") {
		    $ex=$ex."*HNR";
# 		    print "doing it";
		    
		}
	    if ($specs{$tt} eq "" && $flag eq "TD") {
		    $ex=$ex."*HNR";
# 		    print "doing it";
		}
# 	    if ($specs{$tt} eq "" && $flag eq "PD") {
# 		    $ex=$ex."*HNR";
# 		}		    
		
	    $form{$tu}=$form{$tu}.$ex;
	    $noch=$noch+length($ex);
	    if ($temp[2] eq "M" || $temp[1] eq "M") {
		$form{$tu}=$form{$tu}."*DN";
		$noch=$noch+3;
	    }
	}
    }
    $ct++;
}


$noc=int($noch/30000)+1;

# check ratefile is closed

if ($verbose) {
    print "$0: checking terms ... \n";
}
foreach $it (keys %dest) {
    unless ($it eq "" || exists $form{$it} || exists $cons{$it}) { 
	if ($verbose) { print "$0: [safe] species $it is never formed\n"; }
    }
}
foreach $it (keys %form) {
    unless ($it eq "CH2CO+" || exists $dest{$it} || exists $cons{$it}) { 
	die "$0: [serious] species $it is never destroyed\n";
    }
}


###- CH2CO+ is only formed by companion photons.
###- Still compile the ODE file (see exception above),
###- but print this out
foreach $it (keys %form) {
    if ($it eq "CH2CO+" ) { 
	print "$0: [serious] species $it is never destroyed\n";
    }
}


# output section

if ($verbose) { print "$0: writing formatted output ...\n "; }

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime;
$made=(Sun,Mon,Tue,Wed,Thu,Fri,Sat,Sun)[$wday];
$made=$made.", ".(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)[$mon];
$made=$made." ".$mday." ".($year+1900)." at ";
$h=sprintf("%02d", $hour);
$m=sprintf("%02d", $min);
$s=sprintf("%02d", $sec);
$made=$made.$h.":".$m.":".$s;
$args="(N,R,Y,YDOT)";
$spc="      ";
$splitter=$spc."END\nC \nC \n";
$subname="DIFFUN";
$chunk=65;
open(FORT,">$output") || die "$0: [serious] unable to write to file $output\n";
print FORT $spc."SUBROUTINE $subname".$args."\n";
&header;
&comments;
print FORT "C  Begin Conservation equations\n\n";

foreach $it (@csp) {
    print FORT "C $it\n";
    $th="X(".$cons{$it}.")";
    &doit($cterm{$it},0,$th);
}

print FORT "C\nC  Begin differential equations\nC\n";

for $i ( 1 .. $noc+1 ) {
    print FORT $spc."CALL ODES".$i."(N,Y,YDOT)\n";
}

$couted=0;
$ncot=1;
foreach $it ( @sorted ) {
    if ($couted > 30000 || $ncot == 1) {
	if ($ncot > 1) {
	    print FORT $spc."RETURN\n";
	}
	print FORT $splitter;
	print FORT $spc."SUBROUTINE ODES".$ncot."(N,Y,YDOT)\n";
	&header;
	$couted=0;
	$ncot++;
    }
    print FORT "C ".$it." \n";
    &doit($form{$it},1,"F");
    &doit($dest{$it},1,"D");
    $couted=$couted+length($form{$it})+length($dest{$it});
    if ($it eq "H") {
	print FORT $spc."YDOT(".$specs{$it}.")=(F+HLOSS)-(D*Y(".$specs{$it}."))\n";
    } else {
	print FORT $spc."YDOT(".$specs{$it}.")=F-(D*Y(".$specs{$it}."))\n";
    }
}
print FORT $spc."RETURN\n".$spc."END\n";

while ($ncot <= $noc) {
    print FORT "C \nC \n";
    print FORT $spc."SUBROUTINE ODES".$ncot."(N,Y,YDOT)\n";
    &header;
    print FORT $spc."RETURN\n".$spc."END\n";
    $ncot++;
}

    print FORT "C \nC \n";
    print FORT $spc."SUBROUTINE ODES".$ncot."(N,Y,YDOT)\n";
    &headerPlusFooter;
    $ncot++;

close(FORT);

# species file

$sfile=$filename.".specs";
open(SPEC,">$sfile") || die "$0: [serious] unable to write to file $sfile\n";
print SPEC "Index Species      Mass\n";
foreach $it (@sorted) {
    printf SPEC " %4d %-12s %4d\n", ($specs{$it},$it,$masses{$it});
}

$i=0;

print SPEC  " 9999 Conserved:\n";
foreach $it (@csp) {
    $i++;
    if($it eq "e-"){$initval=0.0}
#     If abundances are expressed relative to H2
#     if($it eq "H2"){$initval=1.0}
#     If abundances are expressed relative to H
    if($it eq "H2"){$initval=0.5}
    printf SPEC " %4d %-12s %8.2E\n", ($i,$it,$initval);
}

&parents;

close(SPEC);

if ($verbose) { print "$0: completed successfully\n"; }

sub doit {
#
# line=0, split=1, thing=2
#
    $ln=$_[0];
    print FORT $spc.$_[2]."=0.\n";
    $ps=0;
    while ($ps < length($ln)) {
	$f=substr($ln,$ps,$chunk);
	$pt="     * ".$f;
	if (length($pt) > 7) {
	    if (substr($f,0,1) eq "+" && $ps > $chunk && $_[1]) {
		print FORT $spc.$_[2]."=".$_[2]."\n";
		print FORT "     * ".$f."\n";
	    }
	    else {
		print FORT "     * ".$f."\n";
	    }
	}
	$ps=$ps+$chunk;
    }
}

sub comments {
    print FORT <<_EOC_

C  CALL SUBROUTINE TO CALCULATE RATE COEFFICIENTS
      CALL RATES(R)

C  Made by: $version on $made
C  using ratefile: $rfile

_EOC_
}

sub total {
    $tot=0;
    foreach $t ( @_ ) {
	$tot=$tot+$t;
    }
    return $tot;
}

sub header {
    print FORT <<_EOT_
      DOUBLE PRECISION R,F,D,K,Y,YDOT,X,TOTAL,ACCR,HNR,HLOSS
      INTEGER N
      DIMENSION Y(N),YDOT(N),K(10000),X(10),TOTAL(10)
      COMMON/BL1/ K,X,TOTAL,ACCR,HNR
      HLOSS=-ACCR*Y($specs{"H"})
_EOT_
}

sub headerPlusFooter {
    print FORT <<_EOT_
      DOUBLE PRECISION Y,YDOT,MLOSS,V
      INTEGER N
      DIMENSION Y(N),YDOT(N)
      COMMON/BL2/ MLOSS,V
      DO 20 I=1,N
      YDOT(I)=(1.0/V)*YDOT(I)
 20   CONTINUE
      RETURN
      END
_EOT_
}

sub usage {
    print <<_EOM_
Usage:> rate10odecsT.pl [options]... <ratefile_name> 
Reads colon-separated reaction rate file and writes FORTRAN ODE source code.

  -c=<species>    comma separated list of species to conserve
  -o <file>       place the output ODEs into <file> 
  -q, --quiet     refrain from printing messages   
  -e, --elements  show elements in ratefile
  -s, --species   list species in ratefile
      --help      display this help and exit
      --version   output version information and exit

Report bugs to andrew.markwick\@manchester.ac.uk
_EOM_
}


sub parents {
print SPEC <<_EOA_
 9999 Parents:
He    0.85E-1
CO    1.50E-04
N2    2.00E-05
H2O   1.075E-04
HCN   1.295E-07
CO2   1.50E-07
NH3   3.125E-07
SO    1.53E-06
CS    2.785E-08
H2S   0.875E-05
SO2   1.86E-06
SiO   1.355E-05
SiS   4.765E-07
PO    3.875E-08
PN    0.75E-08
Cl    0.50E-08
F     0.50E-08
_EOA_
}
