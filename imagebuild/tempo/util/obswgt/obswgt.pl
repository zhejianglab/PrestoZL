#!/usr/bin/perl -w

# obswgt -- calculate relative weights of summed over
#           all points from each data set in a tempo run,
#           where "data set" is defined as a group of
#           points with the same entry in info.tmp
#
#           also calculate the relative chi^2s for each
#           data set.  these are simply
#              rchi2 = sum [(residual/est.err)**2]/npts_in_set
#           and a slightly re-scaled version, to take
#           into account degrees of freedom of the whole set:
#              rchi2x = chi2*(total_#_TOAs)/(total_#_degress of freedom)
#
# use: obswgt

$useage =  "\n".
           "use: obswgt [flags]\n".
           "flags:\n".
           "-e: print N_eff=(sum of wgt)^2/(sum of wgt^2) in addition to everything else\n".
           "-f: segregate by frequency bin number [int(freq/100)] as well as info.tmp\n".
           "-m: print start and end MJDs for each set\n".
           "\n";

$infofile = "info.tmp";
$resfile = "resid2.tmp";

# the usual parsing code is here in case it is needed for future development

$freqbinflag = 0;
$neffflag = 0;
$mjdflag = 0;

foreach $par (@ARGV) {
  if ($par=~/^-/) {           # flags
    $par = substr($par,1);      #   strip off the leading "-"
    while ($par ne "") {      #   loop through the string
      $f = substr($par,0,1);  #     strip off the current flag
      $par = substr($par,1);    #     process flags:
      if ($f eq "e") {
        $neffflag = 1;
      } elsif ($f eq "f") {
        $freqbinflag = 1;
      } elsif ($f eq "m") {
        $mjdflag     = 1;
      } elsif ($f eq "h") {
        die $useage   
      } else {
	print "\nError: Invalid flag -$f\n\n";
        die $useage;
      }
    }
  } else {                  # parameters
    push @param, $par;        #   push onto parameter list
  }
}

die $useage if ($#param>-1);
($npts) = @param if ($#param==0);

$mtxtime = (stat("matrix.tmp"))[9];
$tlistime = (stat("tempo.lis"))[9];

if($mtxtime =~ $tlistime) {
  open(C,"matrix.tmp");
  if(read(C,$buf,8)) {
    $nparam = (unpack("x4i1",$buf))[0];
  } else {
    print "Warning: can't read matrix.tmp;  chi^2 values will not\n";
    print "  reflect true number of degrees of freedom\n";
    $nparam = 0;
  }
  close C;
}
else {
  print "Matrix.tmp was not generated with tempo.lis;\n";
  print "assuming only 1 degree of freedom.";
}
$nparam++;  # accounts for fit for phase

open(D,"pwd|");
$dir = <D>;
chomp $dir;
close D;

open(E,"tempo.lis");
while (<E>) {
  last if (/Input data from/);
}
chomp;
$file = $_;
$file =~ s/Input data from *//;

die "Cannot open residual file $resfile\n" if (! -e $resfile);
die "Cannot open info file $infofile\n" if (! -e $infofile);

open(A,$resfile);
binmode A;
open(B,$infofile);

$ninfo = -1;
$infostr = "ALL";  # default code -- all observatories combined
@obsname = ();
$ntoa=0;

while(read(A,$buf,80)) {
  # (undef,undef,$res,undef,$freq,$tpowgt,$terr,undef)=
  ($mjd,undef,$res,undef,$freq,$tpowgt,$terr,undef)=
    unpack("x4d8",$buf);
  $mjd += 39126 if ($mjd<20000); # convert 1966 day number to MJD if needed
  $infostr = <B>;
  next if ($tpowgt==0);
  $res = $res * 1.e6;
  chomp $infostr;
  if ($freqbinflag) {
    $freqbin = int($freq/100);
    $infostr = $infostr."_".$freqbin;
  }
  if (!exists($idxobs{$infostr})) {
    $idxobs{$infostr} = ++$ninfo;
    $obsname[$ninfo] = $infostr;
    $npts{$infostr} = 0;
    $chi2{$infostr} = 0;
    $res2{$infostr} = 0;
    $res2w{$infostr} = 0;
    $sumwgt{$infostr} = 0;
    $sumwgt2{$infostr} = 0;
    $mjdstart{$infostr} = 999999.;
    $mjdend{$infostr} = -999999.;
  }
  $npts{$infostr}++;  
  # $sumwgt{$infostr} += 1/($terr*$terr);
  # $sumwgt2{$infostr} += 1/($terr*$terr);
  $sumwgt{$infostr} += $tpowgt;
  $sumwgt2{$infostr} += $tpowgt*$tpowgt;
  $chi2{$infostr} += ($res*$res)/($terr*$terr);
  # $chi2{$infostr} += ($res*$res)*$tpowgt;
  $res2{$infostr} += $res*$res;
  $res2w{$infostr} += $res*$res*$tpowgt;
  $ntoa++;
  $mjdstart{$infostr} = $mjd if ($mjd<$mjdstart{$infostr});
  $mjdend{$infostr} = $mjd if ($mjd>$mjdend{$infostr});
}
close A;
close B;

$ndof = $ntoa - $nparam;

$wgttot = 0;
foreach $i (keys %sumwgt) {
  $wgttot += $sumwgt{$i};
}

print "\n";
print "Working Directory:    $dir\n";
print "Tempo Input File:     $file\n";
print "Number of TOAs:       $ntoa\n";
print "Number of parameters: $nparam\n";
print "Number of DOF:        $ndof\n";

print "\n";

print "total  avg weight  number    chi^2   adjusted chi^2   rms         rms     ";
print "    neff    " if ($neffflag==1);
print " mjd range   nyear     " if ($mjdflag==1);
print "info  \n"; 

print "weight   per TOA   of TOAs  per TOA    (per DOF)    unweighted  weighted  ";
print "            " if ($neffflag==1);
print "                    " if ($mjdflag==1);
print "      \n"; 

print "\n";

$swgtnorm = 0;
$ssumwgt  = 0;
$ssumwgt2 = 0;
$savgwgt  = 0;
$snpts    = 0;
$schi2    = 0;
$sres2    = 0;
$sres2w   = 0;
$smjdstart = 999999.;
$smjdend = -999999.;

foreach $i (@obsname) {
  $wgtnorm = $sumwgt{$i}/$wgttot;
  $avgwgt = $wgtnorm/$npts{$i};
  $rchi2 = $chi2{$i}/$npts{$i};
  $rchi2x = $rchi2*($ntoa/$ndof);
  $resrms = sqrt($res2{$i}/$npts{$i});
  $resrmsw = sqrt($res2w{$i}/$sumwgt{$i});
  $neff = ($sumwgt{$i})**2/$sumwgt2{$i};

  printf "%7.5f  %7.5f  %5d  %10.4f  %10.4f  %10.4f  %10.4f  ",
    $wgtnorm,$avgwgt,$npts{$i},$rchi2,$rchi2x,$resrms,$resrmsw;
  printf "  %8.2f  ",$neff if ($neffflag==1);
  printf " %5d-%5d %5.2f      ", 
    int($mjdstart{$i}), int($mjdend{$i}), ($mjdend{$i}-$mjdstart{$i})/365.25
         if ($mjdflag==1);
  printf " %s\n", $i;

  $swgtnorm += $wgtnorm;           # this had better sum to 1!
  $ssumwgt  += $sumwgt{$i};
  $ssumwgt2 += $sumwgt2{$i};
  $savgwgt  += $avgwgt*$npts{$i};  # this had better sum to 1/$ntoa after noralization!
  $snpts    += $npts{$i};          # this had better sum to $ntoa!
  $schi2    += $chi2{$i};           
  $sres2    += $res2{$i};           
  $sres2w   += $res2w{$i};
  $smjdstart = $mjdstart{$i} if ($mjdstart{$i}<$smjdstart);
  $smjdend   = $mjdend{$i}   if ($mjdend{$i}  >$smjdend);
}

$savgwgt = $savgwgt/$ntoa;
$srchi2 = $schi2/$ntoa;
$srchi2x = $schi2/($ntoa-$nparam);
$sresrms = sqrt($sres2/$ntoa);
$sresrmsw = sqrt($sres2w/$ssumwgt);
$neff = ($ssumwgt)**2/$ssumwgt2;

print "\n";
printf "%7.5f  %7.5f  %5d  %10.4f  %10.4f  %10.4f  %10.4f  ",
    $swgtnorm,$savgwgt,$snpts,$srchi2,$srchi2x,$sresrms,$sresrmsw;
  printf "  %8.2f  ",$neff if ($neffflag==1);
printf " %5d-%5d %5.2f      ", 
  int($smjdstart), int($smjdend), ($smjdend-$smjdstart)/365.25
       if ($mjdflag==1);
printf " %s\n", "overall";
print "\n";

