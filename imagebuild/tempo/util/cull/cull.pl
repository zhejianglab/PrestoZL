#!/usr/bin/perl -w

# cull -- removes high-residual points from a tempo file
#         based on 20 Jan 99 version of cull.f
#         which was based on 13 Feb 93 version of avtime
#         innumerable modifications since then
#         DJN is to blame for all code written before 18 March 2010
#         After that date, check the CVS log 
#
#

$usage =
  "Use:\n" .
  "  cull infile outfile [flags]\n" .
  "Parameters:\n".
  "  infile:  input TOA file\n" .
  "  outfile: output TOA file(*)\n" .
  "           (*) exception:  when using -o or -O flag, outfile parameter is an extension\n".
  "                    added to input file names to create output files names\n".
  "Flags:\n" .
  "  -r: -rx.xxx = maximum residaual as a fraction of pulse period\n" .
  "      -rxxxs  = maximum residual, seconds\n" .
  "      -rxxxm  = maximum residual, milliseconds\n" .
  "      -rxxxu  = maximum residual, microseconds\n" .
  "         in all cases, if 'xxx' starts with '+' (or '-'), remove points\n".
  "         only with high positive (or negative) residuals.\n".
  "  -sxxx: use xxx-day smoothed timescale for residual test\n".
  "  -axxx: exclude any TOAs before MJD xxx\n".
  "  -bxxx: exclude any TOAs after MJD xxx\n".
  "  -A: \"and\" -a and -b options (default is \"or\").\n".
  "  -pxxx -qyyy: exclude any TOAs between orbital phases xxx and yyy\n".
  "  -c: remove comment lines starting with 'C' in TOA section\n".
  "  -C: Add a copy of the cull arguments as a comment at start of outfile\n".
  "  -d: convert any 1966 day numbers (day<20000) to MJD\n".
  "  -exxx: remove points for which DMX delta-DM error is greater than xxx\n".
  "  -fxxx: keep only sets of TOAs with multiple freqs within xxx days\n" .
  "  -g: remove points for which DMX delta-DM ranges have <25% freq. variation\n".
  "  -h: help message\n".
  "  -i: incorporate tempo's final erros (after efac, etc) into TOA lines\n".
  "  -j: replace JUMP statements with TIME offsets (after running tempo -j)\n".
  "  -k: remove SKIP'd sections\n".
  "  -l: keep only one orbit of TOAs, from phase zero to one\n".
  "  -m: monitor progress of program with xxx/xxx display\n".
  "  -n: convert Parkes/Jodrell format TOAs to Princeton format\n".
  "  -o: preserve the INCLUDE file structure, add 'outext' parameter to each file name\n".
  "  -O: preserve the INCLUDE file structure, substitute 'outext' parameter for each file name extension\n".
  "  -t: incorporate TIME offsets directly into TOAs\n".
  "  -uxxx -vyyy: include only TOAs with frequency between xxx and yyy MHz\n".
  "  -w: remove points with zero weight\n".
  "Notes:\n" . 
  "  + residuals are read from tempo output file resid2.tmp\n".
  "  + extra output (ASCII MJD, residuals, flags) to e.tmp\n".
  "  + options -d and -t work only with Princeton-style TOAs\n";

$cflag = 0;
$Cflag = 0;
$dflag = 0;
$eflag = 0;
$fflag = 0;
$gflag = 0;
$iflag = 0;
$jflag = 0;
$kflag = 0;
$lflag = 0;
$mflag = 0;
$nflag = 0;
$oflag = 0;
$Oflag = 0;
$pflag = 0;  # residual cut defined as fraciton of a pulse period
$rflag = 0;
$rnflag = 1;
$rpflag = 1;
$sflag = 0;
$tflag = 0;
$uflag = 0;
$vflag = 0;
$gdiff = 0;
$lflag1 = 0;
$lflag2 = 0;
$wflag = 0;
$hdrflag = 0;

$slop = 600./86400;  # slop in days allowed when checking for membership
                     # in a DMX range.  the issue is that DMX ranges are
                     # defined in observatory coordinates whereas the
                     # TOAs in the residual file are in SSB coordinates
                     # solution is to allow a little slop when interpreting
                     # DMX ranges, though this could be a problem if
                     # they are intentionally defined with high precision


$mjdmin = 0;
$mjdmax = 9999999;
$mjdand = 0;

$pha1in = -1;
$pha2in = -1;

$freqmin = 0.;
$freqmax = 1.e20;   # a big number


foreach $par (@ARGV) {
  if ($par=~/^-/) {           # flags
    push @flags, $par;        # Keep a list of these for reference
    $par = substr($par,1);      #   strip off the leading "-"
    while ($par ne "") {      #   loop through the string
      $f = substr($par,0,1);  #     strip off the current flag
      $par = substr($par,1);    #     process flags:
      if ($f eq "a")      {#       -a
        $mjdmin = $par;
        last;
      } elsif ($f eq "b") {#       -b
        $mjdmax = $par;
        last;
      } elsif ($f eq "C") {#       -C
        $Cflag = 1;
      } elsif ($f eq "c") {#       -c
        $cflag = 1;
      } elsif ($f eq "d") {#       -d
        $dflag = 1;
      } elsif ($f eq "e") {#       -e
        $eflag = 1;
        $dmxmax = $par;
        last;
      } elsif ($f eq "f") {#       -f
        $fflag = 1;
        $nday = $par;
        $fdiff = 0.25;
        last;
      } elsif ($f eq "g") {#       -g
        $gflag = 1;
        $gdiff = 0.25;
        last;        
      } elsif ($f eq "h") {#       -h
        die $usage;
      } elsif ($f eq "i") {#       -i
        $iflag = 1;
      } elsif ($f eq "j") {#       -j
        $jflag = 1;
      } elsif ($f eq "k") {#       -k
        $kflag = 1;
      } elsif ($f eq "l") {#       -l
        $lflag = 1;
      } elsif ($f eq "m") {#       -m
        $mflag = 1;
      } elsif ($f eq "n") {#       -n
        $nflag = 1;
      } elsif ($f eq "o") {#       -o
        $oflag = 1;
      } elsif ($f eq "O") {#       -o
        $Oflag = 1;
      } elsif ($f eq "r") {#       -r
        $rflag = 1;
        if ($par=~/s$/) {             #   parameter ends in 's'
          $r = substr($par,0,-1);                  
        } elsif ($par=~/m$/) {        #   parameter ends in 'm'
          $r = substr($par,0,-1)*1.e-3;
        } elsif ($par=~/u$/) {        #   parameter ends in 'u'
          $r = substr($par,0,-1)*1.e-6;
        } elsif ($par=~/\d$/) {       #   parameter ends in a digit
          $r = $par;
          $pflag = 1;
        } else {
          print "\nError: Cannot parse residual parameter $r\n\n";
          die $usage;
        }
        if (substr($par,0,1) eq "-") {
          $rpflag = 0; # no cut of positive residuals
          $r = -$r;
        } elsif (substr($par,0,1) eq "+") {
          $rnflag = 0; # no cut of negative residuals
        }
        last;
      } elsif ($f eq "p") {
        $pha1in = $par;
        last;
      } elsif ($f eq "q") {
        $pha2in = $par;
        last;
      } elsif ($f eq "s") {#       -s
        $sflag = 1; 
        $nsmooth = $par;
        last;
      } elsif ($f eq "t") {#       -t
        $tflag = 1;
      } elsif ($f eq "u") {#       -u
        $uflag = 1;
        $freqmin = $par;
        last;
      } elsif ($f eq "v") {#       -v
        $vflag = 1;
        $freqmax = $par;
        last;
      } elsif ($f eq "w") {#       -w
        $wflag = 1;
        last;
      } elsif ($f eq "A") {#       -A
        $mjdand = 1;
        last;
      }  else {            #       invalid flag
	print "\nERROR: Invalid flag -$f\n\n";
        die $usage;
      }
    }
  } else {                  # parameters
    push @param, $par;        #   push onto parameter list
  }
}

die $usage if ($#param!=1);

$oOflag = $oflag || $Oflag;

if ($pha1in>=0 && $pha1in<=1 && $pha2in>=0 && $pha2in<=1) {
  if ($pha1in<$pha2in) {
    $pha1 = $pha1in;
    $pha2 = $pha2in;
    $pha3 = 2;
    $pha4 = -1;
  } else {
    $pha1 = 0;
    $pha2 = $pha2in;
    $pha3 = $pha1in;
    $pha4 = 1;
  }
} elsif ($pha1in==-1 && $pha2in==-1) {
  $pha1 = 2;  # phase will never be more than 2 and less than -1
  $pha2 = -1; # so these numbers will allow phase testing without
  $pha3 = 2;  # danger of zapping points
  $pha4 = -1;
} else {
  die "ERROR: For orbital phase cuts, both -p and -q must be used\n".
    "       Each must have a parameter in the range 0<=x<=1\n".
      "$usage";
}

($infile,$outfile) = @param;

$outext = $outfile if ($oOflag);

# read in the residuals
open (A,"resid2.tmp") or die "Error opening resid2.tmp";
binmode A;
$n = 0;
$p = 0;
# while ($q=read(A,$buf,80)) {
while (read(A,$buf,80)) {
  # undef'd unpacked variable 7 (zero indexed) is dt1:
  ($mjd[$n],$res[$n],$dt2sec,$phase[$n],$freq[$n],$wgt[$n],$terr[$n],undef) = 
    unpack("x4d8",$buf);
  $mjd[$n]+= 39126 if ($mjd[$n]<20000);  # convert futs to MJDs
  $p = $dt2sec/$res[$n] if ($p==0 && $res[$n]!=0);
  $n++;
}
close A;


if ($eflag || $gflag) {
  @eok = &dmxcheck(\@mjd,\@freq,$dmxmax,$eflag,$gflag,$gdiff);
} else {
  foreach $i (0..$n-1) {
    $eok[$i] = 1;
  }
}

if ($fflag) {
  @fok = &dualfreq(\@mjd,\@freq,$nday,$fdiff);
} else {
  foreach $i (0..$n-1) {
    $fok[$i] = 1;
  }
}

@jump = &getjump if ($jflag);

if ($rflag) {  # preparatory work for cull-by-residual
  die "All residuals exactly zero ?!\n" if ($p==0);
  $r = $r/$p if (!$pflag);  # max residual, $r$, always fraction of period
  if (!$sflag) {
    $iaold = -1; # force re-calculation from scratch ... 
    $ibold = -1; # ... for efficiency reasons only
    $avall=&getav(0,$n-1);    
  }
}

# now we are ready to process TOAs


# check to see if it has a header with parameters at the top of the file
# (rather than just being a .tim file)
# search for TOA at start of line.  Don't need to consider INCLUDEd files

open (A,$infile);
while (<A>) {
  $_ =~ s/^\s+|\s+$//g ; # remove leading or trailing whitespace
  if (length>1 && uc((split)[0]) eq "TOA") {
    $hdrflag = 1;
    last;
  }
}
close A;

# now open the input file for real

open (A,$infile);

# open the output files

if ($oflag) {
  $outfile = "$infile.$outext";
  print "oflag section set outfile to $outfile\n";
}
if ($Oflag) {
  @infilex = split(/\./,$infile);
  $outfile = join(".",@infilex[0..$#infilex-1]).".".$outext;
}
open (Z,">$outfile");
open (Y,">e.tmp");


# copy the header if one is present

if ($hdrflag) {
  while (<A>) {
    print Z $_;
    chomp;
    next if (/^ *$/);  # avoid warnings from split function in next line
    last if (uc((split)[0]) eq "TOA");
  }
}

# Print a copy of the cull arguments as a comment
if ($Cflag) { print Z "C cull " . join(' ',@flags) . "\n"; }
  
# now loop through lines, allowing switches into and out of INCLUDE'd files

$skip = 0;

$fhin[0] = *A;
$fhinx = $fhin[0];

$fhout[0] = *Z;
$fhoutx = $fhout[0];

$i=0;
$| = 1 if ($mflag);
print "     " if ($mflag);  

$iaold = -1; 
$ibold = -1;

$toff = 0.;

$jcounter = 0;

$emax = 1.e99;

$tempo2 = 0;

for(;;) {

  printf "\b\b\b\b\b\b\b\b\b\b\b%5d/%5d", $i,$n if ($mflag);
  if (eof($fhinx)) {     # file processing and input section
    close pop @fhin;
    close pop @fhout if ($oOflag);
    last if ($#fhin==-1);
    $fhinx = $fhin[$#fhin];
    $fhoutx = $fhout[$#fhout] if ($oOflag);
    next;
  }
  $a = <$fhinx>;
  if ($a=~/^ *$/) {
    $a1 = "";
  } else {
    $a1 = uc((split(' ',$a))[0]);
  }
  if ($a1=~/^INCLUDE/) {
    if ($oOflag) {  # print INCLUDE line in new file if preserving INCLUDE structure
      $a2 = $a;
      chomp $a2;
      $a2 =~ s/ *$//;
      if ($oflag) {
        print $fhoutx "$a2.$outext\n";  
      }
      if ($Oflag) {
        @a2x = split(/\./,$a2);
        print $fhoutx join(".",@a2x[0..$#a2x-1]).".$outext\n";
      }
    }
    local *A;
    local $infile;
    $infile = (split(" ",$a))[1];
    open (A,$infile);
    $fhin[$#fhin+1] = *A;
    $fhinx = $fhin[$#fhin];
    if ($oOflag) {
      local *Z;
      local $outfile;
      if ($oflag) {
        $outfile = "$infile.$outext";
      }
      if ($Oflag) { 
        @infilex = split(/\./,$infile);
        $outfile = join(".",@infilex[0..$#infilex-1]).".".$outext;
      }
      open (Z,">$outfile");
      $fhout[$#fhout+1] = *Z;
      $fhoutx = $fhout[$#fhout];
    }
    next;
  } 
  #                      process commands
  last if ($a1=~/^END/);
  $skip = 1 if ($a1=~/^SKIP/);
  $skip = 0 if ($a1=~/^NOSK/);

  next if ($cflag && $a=~/^C /); # skip comments if $cflag
  next if ($kflag && ($skip || $a1=~/^NOSK/)); # skip SKIP'd sections if $kflag

  if ($jflag && $a1=~/^JUMP/) {
    $jval = $jump[int($jcounter/2)];
    $jval = -$jval if (($jcounter%2)==1);
    $a = "TIME $jval\n";
    $jcounter++;
  }


  #                      process TIME cards  
  if ($tflag && $a1=~/^TIME/) {
    if ($skip) {
      printf "\b\b\b\b\b\b\b\b\b\b\b" if ($mflag);
      printf "Warning: ignoring TIME card in skipped section:\n$a";
    } else {
      $toff += (split(' ',$a))[1];
      next;
    }
  }

  #                      skip EFACS, etc., in iflag=1 mode
  if ($iflag && ($a1=~/^EFAC/ || $a1=~/^EQUAD/ || $a1=~/^EMIN/)) {
    next;
  }

  if ($a1=~/^EMAX/) { 
    $emax = (split(' ',$a))[1];
  }

  if ($a1=~/^FORMAT/) {
    $tempo2 = 1;
    print "Warning: FORMAT line found, assuming tempo2 format.\n";
    print "  Note that some features of cull may not work correctly.\n";
    print "  But actually they usually do, so don't worry too much.\n";
  }

  #                      process TOAs
  if (($a=~/^[0-9a-z@ ]/ || ($tempo2 && &t2_istoa($a)))
          && !($a=~/^ *$/) && !$skip) {
    die "Error: $n TOAs in resid2.tmp, but more in $infile\n" if ($i>$n);

    #                    optionally convert Parkes/Jodrell fmt to Princeton fmt
    if ($nflag) {
      $nfmt = &getformat($a); 
      if ($nfmt==1) {
        $jbfreq = substr($a,25,9);
	$jbmjd = substr($a,36,19);
	$jbphs = substr($a,55,8);
        die "nonzero PHS in TOA at $jbmjd.\n" if ($jbphs != 0);
	$jbterr = substr($a,63,8);
	$jbsite = substr($a,79,1);
        $a = sprintf("%1s              %9s %19s %8s                     \n",
                   $jbsite, $jbfreq, $jbmjd, $jbterr);
      }
    }

    if ($tempo2) {
      $err = (split(' ',$a))[3];
    } else {
      $err = substr($a,44,9);
    }

    if ($err<$emax) {

      #                    in iflag==1 mode, substitute tempo error
      if ($iflag) {       
          substr($a,44,9) = sprintf (" %8.2f",$terr[$i]);
      }

      $cull = 0;
      # check for zero weight
      if ($wflag && $wgt[$i]==0) {
        $cull = 1;
      }
      #                    check date range
      if ($mjdand==0 && ($mjd[$i]<$mjdmin || $mjd[$i]>$mjdmax)) {
        substr($a,0,1) = "C";
        $cull = 1;
      }
      if ($mjdand==1 && ($mjd[$i]<$mjdmin && $mjd[$i]>$mjdmax)) {
        substr($a,0,1) = "C";
        $cull = 1;
      }
      #                    check dmx flag, and dual-frequency flag 
      if (!$fok[$i] || !$eok[$i]) { 
        substr($a,0,1) = "C";
        $cull = 1;
      }
      #                    check orbital phase
      if (($phase[$i]>$pha1 && $phase[$i]<$pha2) || 
          ($phase[$i]>$pha3 && $phase[$i]<$pha4))  {
        substr($a,0,1) = "C";
        $cull = 1;
      }
      if ($vflag>0 && $freq[$i]<$freqmin || 
          $vflag>0 && $freq[$i]>$freqmax )        {
        substr($a,0,1) = "C";
        $cull = 1;
      }
      if ($lflag) {
        $lflag1++  if ($i>0 && $phase[$i]<$phase[$i-1]);
        $lflag2 = 1 if ($lflag1>1 && $phase[$i]<$phase[$i-1]);
        $cull = 1 if ($lflag1==0 || $lflag2==1);
      }
      if ($rflag && (!$cull)) {# check residual, but skip if date is out of range
        if ($sflag) { # use only a subset of residuals around current point
          for ($ia=$i-1; $ia>=0; $ia--) {
            last if ($mjd[$ia]<$mjd[$i]-$nsmooth || 
                     abs(($freq[$ia]-$freq[$i])/$freq[$i])>0.10 ||
                     $mjd[$ia]>$mjd[$ia+1]);
          }
          $ia++;
          for ($ib=$i+1; $ia<=$n; $ib++) {
            last if ($mjd[$ib]>$mjd[$i]+$nsmooth || 
                     abs(($freq[$ib]-$freq[$i])/$freq[$i])>0.10 ||
                     $mjd[$ib]<$mjd[$ib-1]);
          }
          $ib--;
          $av = &getav($ia,$ib);
        } else {
          $av = $avall;
        }    
        $r1 =  ($res[$i]-$av);
        $r1 = $r1 - int($r1);
        $r1++ if ($r1<0.);
        $cull = 1 if ($rpflag && $r1>$r && $r1<=0.5);
        $cull = 1 if ($rnflag && $r1>0.5 && $r1<=1-$r);
        $r2 =  ($res[$i]-$avall);  # calculate res relative to average ...
        $r2 = $r2 - int($r1);      # ...of all points;  used only for printout
        $r2++ if ($r2<-0.5);
        $r2-- if ($r2>=0.5);
        $r1-- if ($r1>0.5);   # go from 0<$r1<1 to -0.5<$r1<0.5
        printf Y "%17.11f %13.11f %13.11f %13.2f %13.2f %1d\n",
        $mjd[$i], $r1, $r2, $r1*$p*1.e6, $r2*$p*1.e6, $cull;
      }
  
      # Append a C rather than erasing 1st char
      if ($cull) { $a = "C " . $a; }
  
      if ($dflag || ($tflag && $toff!=0)) { # adjustments to TOA
        $mjdidx = 0;
        $mjdidx = 24 if (substr($a,29,1) eq ".");
        $mjdidx = 25 if (substr($a,30,1) eq ".");
        if ($mjdidx==0) {
          print "Warning: can't parse TOA for offset or MJD adjustment for:\n$a";
        } else {
          $nmjd = substr($a,$mjdidx,5);
          $fmjd = substr($a,$mjdidx+5,15);
          $fmjd =~ s/ *$//;               # chop trailing spaces from $fmjd
          $lfmjd = length($fmjd);         # only write out this many digits
          if ($nmjd<20000 || $toff!=0) {  # skip next part if no adjustment needed
            $fmjd += $toff/86400. if ($tflag);
            $nmjd += 39126 if ($dflag && $nmjd<20000);
            if ($fmjd>=1) {
              $nmjd = $nmjd + int($fmjd);
              $fmjd = $fmjd - int($fmjd);
            } elsif ($fmjd<0) {
              $nmjd = $nmjd + int($fmjd) - 1;
              $fmjd = $fmjd - int($fmjd) + 1;
            }
            substr($a,$mjdidx,5) = sprintf("%5d",$nmjd);
            $tmp = sprintf("%16.14f",$fmjd);
            if ($tmp =~ /^0/) {
              substr($a,$mjdidx+5,$lfmjd) = substr($tmp,1,$lfmjd);
            } else {
              die "Internal error:  fmjd=$tmp not in range 0 to 1 on card:\n$a\n";
            }
                                          #  avoid collision of freq, MJD fields:
            substr($a,23,1) = " " if (substr($a,23,1) eq "0");
          }

        }
      }
      $i++;
    }
  }
  
  print $fhoutx $a if (!($cflag && $cull));
    
}

print "\n" if ($mflag);


# dualfreq:  check that for a given TOA there is at least one other
#    TOA in the data set within nday days that is differs from it
#    in frequency by fractional amount fdiff (i.e., 
#      2*abs(f1-f2)/(f1+f2)>=fdiff).  This is useful to find epochs
#    on which dual-frequency data were taken.
# inputs:
#    \@mjd[0..n]
#    \@freq[0..n]
#    $nday
#    $fdiff
# returns:
#    @freqflag[0..n] = 1 for TOAs from dual-freq epochs; 0 otherwise
# note: code is not very efficient.  should sort by MJD and then
#       compare nearby array elements.  

sub dualfreq {
  my ($mjd,$freq,$nday,$fdiff) = @_;
  my ($i,$j);
  my (@ok);

  foreach $i (0..$#{$mjd}) {
    $ok[$i] = 0; 
  }
  foreach $i (0..$#{$mjd}-1) {
    foreach $j ($i+1..$#{$mjd}) {
      if (abs(${$mjd}[$i]-${$mjd}[$j])<$nday  &&
        2*abs(${$freq}[$i]-${$freq}[$j])/(${$freq}[$i]+${$freq}[$j])>$fdiff) {
        $ok[$i] = 1;
        $ok[$j] = 1;
      }
    }
  }
@ok;
}






#
# getav: find "average" residual point of an array 
# (to get a subarray, pass only part of the array to this routine)
#

sub getav {
  my ($ia,$ib) = @_;
  my ($i,$i1,$i2);
  my ($r1,$r2);
  my ($n2);
  # note: s, n1, iaold, ibold not local, so that they can be re-read

  if ($ib<$iaold || $ia>$ibold) {  # start completely anew
    foreach $i (0..15) {
      $s[$i] = 0;
    $n1[$i] = 0;
    }
    foreach $i ($ia..$ib) {
      &addpt($i);
    }
  } else {
    foreach $i ($ia .. $iaold-1, $ibold+1 .. $ib) {
      &addpt($i);
    }
    foreach $i ($iaold .. $ia-1, $ib+1 .. $ibold) {
      &subpt($i);
    }
  }
  $iaold = $ia;
  $ibold = $ib;
  $n2=$n1[0]+$n1[15];
  $i2 = 0;
  foreach $i (1..15) {
    if ($n1[$i-1]+$n1[$i] > $n2) {
      $n2 = $n1[$i-1]+$n1[$i];
      $i2 = $i;
    }
  }
  $i1 = $i2 - 1;
  if ($i1==-1) {
      # note: "-$n1[15]" in next line takes residuals like 0.999
      # and makes them -0.001, for proper addition to the s[0] resids.
    return ($s[0]+$s[15]-$n1[15])/($n1[0]+$n1[15]);
  } else {
    return ($s[$i1]+$s[$i2])/($n1[$i1]+$n1[$i2]); 
  }
}


sub addpt {
  my ($i) = @_;
  my ($r1,$r1idx);
  $r1 = (10.+$res[$i])-int(10.+$res[$i]);   # ensure 0<=$a<1
  $r1idx = int(16*$r1);
  $s[$r1idx] += $r1;
  $n1[$r1idx]++;
}

sub subpt {
  my ($i) = @_;
  my ($r1,$r1idx);
  $r1 = (10.+$res[$i])-int(10.+$res[$i]);   # ensure 0<=$a<1
  $r1idx = int(16*$r1);
  $s[$r1idx] -= $r1;
  $n1[$r1idx]--;
}


# getdmx: --get ranges and error bars for DMX values
#         --use whatever xxx.par file corresponds to the
#           pulsar name in tempo.lis
#         --code liberally borrowed from dmxoff
sub dmxcheck {

  my ($mjd,$freq,$dmxmax,$eflag,$gflag,$gdiff) = @_;
  my ($i, $j);
  my (@ok);
  my (@dm, @dmx, @mjd1, @mjd2);
  my (@minfreq, @maxfreq);
  my (@dmxrange);

  foreach $i (0..$#{mjd}) {
    $ok[$i] = 1;
  }

# generate .par name from tempo.lis
  my($psr) = "";
  open (GETDMX_A,"tempo.lis");
  while (<GETDMX_A>) {
    if (/^Assumed parameters/) {
      $psr = (split)[4];
      last;
    }
  }
  close GETDMX_A;
  die "Error: Can't find pulsar name in tempo.lis to ".
    "construct input file name\n" if ($psr eq "");

# read in dmx values
  die "Error, no parameter file $psr.par\n" if (!(-e "$psr.par"));
  open (GETDMX_B,"$psr.par");
  $dmxinuse = 0;
  while (<GETDMX_B>) {
    exit if (/^ *TOA/);
    $dmxinuse = 1 if (/^ *DMX/);
    my($a,$b,$c) = (split)[0,1,3];
    $b = uc($b);
    $b =~ s/D/E/g;
    $c = uc($c);
    $c =~ s/D/E/g;
    my($d,$e) = (split("_",$a))[0,1];
    $e--;  # fortran->perl array index convention  
    if ($d eq "DMX" && $e > -1) {
      $dm[$e] = $b;
      $dme[$e] = $c;
    } elsif ($d eq "DMXR1") {
      $mjd1[$e] = $b-$slop;  # expand range by $slop; see notes above
    } elsif ($d eq "DMXR2") {
      $mjd2[$e] = $b+$slop;
    }
  }
  close GETDMX_B;
  die "Error: DMX not used in parameter file $psr.par\n" if (!$dmxinuse);

  if ($eflag || $gflag) {
    foreach $i (0..$#{$mjd}) {
      $dmxrange[$i] = -1;
      foreach $j (0..$#mjd1) {
        $dmxrange[$i] = $j if (${$mjd}[$i]>$mjd1[$j] && ${$mjd}[$i]<$mjd2[$j]);        
      }
      die "Can't find DMX range for TOA at MJD ${$mjd}[$i]\n" if ($dmxrange[$i]==-1);
    }

  }

  if ($gflag) {
    foreach $j (0..$#mjd1) {
      $minfreq[$j]=9.e50;
      $maxfreq[$j]=0;
    }
    # find minimum and maximum frequency in each dmx range
    foreach $i (0..$#{$mjd}) {
      $maxfreq[$dmxrange[$i]]=${freq}[$i] if (${$freq}[$i]>$maxfreq[$dmxrange[$i]]);
      $minfreq[$dmxrange[$i]]=${freq}[$i] if (${$freq}[$i]<$minfreq[$dmxrange[$i]]);
    }
    foreach $j (0..$#mjd1) {
      $usedmx[$j] = 1;
      $usedmx[$j] = 0 if (($maxfreq[$j]-$minfreq[$j])/$maxfreq[$j] < $gdiff);
    }
    foreach $i (0..$#{$mjd}) {   # this used to loop through $#{$mjd}-1.  Why???
      $ok[$i] = 0 if ($usedmx[$dmxrange[$i]]==0);
    }
  }

  # check MJD of each toa to ensure that the relevant dme is less than dmxmax
  if ($eflag) {
    foreach $i (0..$#{$mjd}) {
      $ok[$i] = 1 if ($dme[$dmxrange[$i]]<$dmxmax);
    }
  }

@ok;  
}

# getformat: Figure out the format of a TOA line,
#   0 = Princeton format
#   1 = Parkes/Jodrell format
#   2 = ITOA format
# These are equivalent to the nfmt variable in $TEMPO/src/arrtim.f,
# and the same algorithm is used as arrtim.f

sub getformat {  
  my($a) = $_[0];
  my($nfmt) = 0;  # default, 0 princeton format
  my ($i);
  my ($b);

  if ($a =~ /^ /) {  # first char blank, 1 pks format
    $nfmt = 1 
  } else {
    if ($a=~/^. /) {  # first char not blank, second blank
      $i = 1;
      while (($b=substr($a,$i++,1)) ne '\0') {
        if ($b eq " ") {
          $nfmt = 0;
          last;
        } elsif ($b eq "+" || $b eq "-") {
          $nfmt = 2;
          last; 
        }
      }
    }
  }
  $nfmt;
}

# t2_istoa: Check whether a line appears to be a tempo2 TOA line
sub t2_istoa {
  my @stuff = split(' ',$_[0]);

  # Needs at least 5 fields
  if (scalar(@stuff) < 5) { return 0; }
  # 2nd thru 4th should all be numeric (freq, mjd, err)
  if ($stuff[1] !~ /^[0-9.]*$/) { return 0; }
  if ($stuff[2] !~ /^[0-9.]*$/) { return 0; }
  if ($stuff[3] !~ /^[0-9.]*$/) { return 0; }

  # Assume that's good enough?
  return 1;
}


#
#
# getjump: --read xxx.par file, get JUMP values
#
sub getjump{

  my(@jump);

  # generate .par name from tempo.lis
  my($psr) = "";
  open (GETDMX_A,"tempo.lis");
  while (<GETDMX_A>) {
    if (/^Assumed parameters/) {
      $psr = (split)[4];
      last;
    }
  }
  close GETDMX_A;
  die "Error: Can't find pulsar name in tempo.lis to ".
    "construct input file name\n" if ($psr eq "");

# read in jump values.  
  die "Error, no parameter file $psr.par\n" if (!(-e "$psr.par"));
  open (GETDMX_B,"$psr.par");
  my($jumpsinuse) = 0;
  while (<GETDMX_B>) {
    exit if (/^ *TOA/);
    if (/^JUMP/) {
      $jumpsinuse = 1;
      my($a,$b) = (split)[0,1];
      my($c) = (split('_',$a))[1];
      $jump[$c-1] = $b;  # so jump_1 goes into $jump[0], etc.
    }
  }
  die "Error, no JUMP values in $psr.par.\n".
      "Run 'tempo -j' on the input file before running cull.\n" if (!$jumpsinuse);
  @jump;
}
