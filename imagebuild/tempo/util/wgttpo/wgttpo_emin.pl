#!/usr/bin/perl
($parfile,$toafile) = @ARGV;
if(!(-e mode.dummy)) {
  open(OUT,">mode.dummy");
  printf OUT "MODE 1\n\n";
  close OUT;
}
open(IN,$toafile) || die "cannot open $toafile for reading: $!";
open(OUT,">wgttpo.out") || die "cannot open wgttpo.out for writing: $!";
printf OUT "C %s %s\n\n", $parfile, $toafile;
printf OUT "MODE 1\n\n";
while (<IN>) {
  ($order,$toas) = split;
  print $order," ",$toas,"\n";
  if($order =~ /INCLUDE/) {
    system "cat mode.dummy $toas > $toas.tmp";
    system "tempo -f $parfile $toas.tmp";
    system "tail -2 tempo.lis > tempo.lis.end";
    open(IN2,"tempo.lis.end");
    $_ = <IN2>;
    ($jnq,$jnq2,$jnq3,$jnq4,$wresp,$jnq5,$jnq6,$jnq7,$wres,$jnq8) = split;
    $_ = <IN2>;
    $goodstring = substr($_,33);
    ($chisq,$jnq3,$pp,$jnq4,$wmax) = split(" ",$goodstring);
    print "Initial: ",$wres," ",$wmax," ",$chisq," ",$emin,"\n";
    $emin = $wres*$chisq/20.;
    print $wres," ",$wmax," ",$chisq," ",$emin,"\n";
    if ($wmax <= 2 || $chisq <= 1.05) {
      $emin = 0.;
    }
    while ($wmax > 2 && $chisq > 1.05) {
      open(OUT2,">emin.dummy") || die "cannot open emin.dummy for writing: $!";
      printf OUT2 "EMIN %8.4f\n", $emin;
      close(OUT2);
      system "cat emin.dummy";
      system "cat mode.dummy emin.dummy $toas > $toas.tmp"; 
      system "tempo -f $parfile $toas.tmp";
      system "tail -2 tempo.lis > tempo.lis.end";
      open(IN2,"tempo.lis.end");
      $_ = <IN2>;
      ($jnq,$jnq2,$jnq3,$jnq4,$wresp,$jnq5,$jnq6,$jnq7,$wres,$jnq8) = split;
      $_ = <IN2>;
      $goodstring = substr($_,33);
      ($chisq,$jnq3,$pp,$jnq4,$wmax) = split(" ",$goodstring);
      print "Before reweighting: ",$wres," ",$wmax," ",$chisq," ",$emin,"\n";
      $emin += $wres/40.;
      print $wres," ",$wmax," ",$chisq," ",$emin,"\n";
    }
    open(OUT2,">emin.dummy") || die "cannot open emin.dummy for writing: $!";
    printf OUT2 "EMIN %8.4f\n", $emin;
    close(OUT2);
    system "cat emin.dummy";
    system "cat mode.dummy emin.dummy $toas > $toas.tmp";
    system "tempo -f $parfile $toas.tmp";
    system "tail -2 tempo.lis > tempo.lis.end";
    open(IN2,"tempo.lis.end");
    $_ = <IN2>;
    ($jnq,$jnq2,$jnq3,$jnq4,$wresp,$jnq5,$jnq6,$jnq7,$wres,$jnq8) = split;
    $_ = <IN2>;
    $goodstring = substr($_,33);
    ($chisq,$jnq3,$pp,$jnq4,$wmax) = split(" ",$goodstring);
    print "Final weighting: ",$wres," ",$wmax," ",$chisq," ",$emin,"\n";
    $efac = sqrt($chisq);
    printf OUT "EMIN %8.4f\n", $emin; 
    printf OUT "EFAC %8.4f\n", $efac; 
    printf OUT "INCLUDE %s\n",$toas;
    system "rm $toas.tmp";
  }
  elsif(!($order =~ /MODE/) && !($order =~ /EMIN/) && !($order =~ /EQUAD/) && !($order =~ /EFAC/)) {
    print OUT
  }
}
