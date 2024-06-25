#!/usr/bin/perl  -w

# extract-pl:  perl version of extract.f
#

$help =  "\n".
"use: extract [outfile] [-h]\n".
"   'extract' with no file writes to output file 'c.tmp'\n".
"   '-h' flag gives this help message\n".
"\n".
"        output file contains (in columns):\n".
"        1. TOA day number in 1966 days (i.e., 1-Jan-1966 = 0)\n".
"        2. postfit residual (in P)\n".
"        3. postfit residual (in microsec)\n".
"        4. timing error used by Tempo (in microsec)\n".
"        5. prefit residual (in microsec)\n".
"        6. orbital phase\n".
"        7. mjd\n".
"        8. ptype (use as point type in mongo/sm to distinguish points by freq.)\n".
"        9. serial number\n".
"        10. day of year (day number mod 365.25)\n".
"        11. observing frequency (in MHz)\n".
"        12. time of day in seconds\n".
"        13. serial number, except that this counter skips ahead 10\n".
"            each time there is a gap of more than 0.25 days or a\n".
"            negative jump in the MJDs\n".
"        14. MJD-50000   (less quantization problems in mongo....)\n".
"        15. weight of point in tempo fit\n".
"        16. red noise model (in microsec) or ddm*1e6\n".
"        17. residual minus red noise model (in microsec)\n".
"\n";


$outfile = 'c.tmp';
@parlist = ();

foreach $i (0..$#ARGV) {
  $a = $ARGV[$i];
  if ($a=~/^-/) {
    $a =~ s/^.//;  # remove leading character
    while (length($a)>0) {
      if ($a=~/^h/) {
        die $help; 
      } else {
        # could process other flags here....
        # default, can't understand flag ==> die with help message
        die $help;   
      }
      $a =~ s/^.//;  # remove leading character
    }
  } else {
    push @parlist, $a;
  }
}
die $help if ($#parlist>0);
$outfile = $parlist[0] if ($#parlist==0);
open (B,">$outfile");

open (A,"resid2.tmp\n");
binmode A;
$p = 0;
$mjd0 = 0 ;
$i = 0;
$i2 = -10;
while (read(A,$buf,80)) { 
  ($mjd,$res,$ressec,$phase,$freq,$wgt,$terr,$dt1,$resred) = 
    unpack("x4d9",$buf);
  $mjd+= 39126 if ($mjd<20000);  # convert 1966 day number to MJD if needed
  $ct = $mjd-39126;              # now convert MJD to 1966 day number
  $p = $ressec/$res if ($p==0 && $res!=0);
  $ptype = 52;
  $ptype = 30 if ($freq>320);
  $ptype = 31 if ($freq>380);
  $ptype = 52 if ($freq>500);
  $ptype = 32 if ($freq>607);
  $ptype = 41 if ($freq>609);
  $ptype = 800 if ($freq>750);
  $ptype = 40 if ($freq>1000);
  $ptype = 61 if ($freq>1500);
  $i++;
  $i2++;
  $i2 += 10 if ($mjd>$mjd0+0.25 || $mjd<$mjd0);
  $mjd0 = $mjd;
  $dayyr = $ct - 365.25*int($ct/365.25);
  $fmt = "%17.11f %14.8f %13.2f %10.3f %12.3f %9.7f %16.10f %6.1f ".
    "%5d %8.4f %9.3f %12.6f %5d %16.10f %10.5f %13.2f %13.2f\n";
  printf B $fmt,
    $ct, $res, $ressec*1.e6,$terr,$dt1*1.e6*$p, $phase, $mjd, $ptype, 
    $i, $dayyr, $freq, ($ct-int($ct))*86400., $i2, $mjd-50000, $wgt, $resred*1.e6,  
    ($ressec-$resred)*1.e6;
  
}

close A;
close B;

