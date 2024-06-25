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
    $equad = $wres*$chisq/20.;
    print $wres," ",$wmax," ",$chisq," ",$equad,"\n";
    if ($wmax < 1 || $chisq <= 1.00) {
      $equad = 0.;
    }
    while ($wmax >= 1 && $chisq > 1.01) {
      open(OUT2,">equad.dummy") || die "cannot open equad.dummy for writing: $!";
      printf OUT2 "EQUAD %8.4f\n", $equad;
      close(OUT2);
      system "cat equad.dummy";
      system "cat mode.dummy equad.dummy $toas > $toas.tmp"; 
      system "tempo -f $parfile $toas.tmp";
      system "tail -2 tempo.lis > tempo.lis.end";
      open(IN2,"tempo.lis.end");
      $_ = <IN2>;
      ($jnq,$jnq2,$jnq3,$jnq4,$wresp,$jnq5,$jnq6,$jnq7,$wres,$jnq8) = split;
      $_ = <IN2>;
      $goodstring = substr($_,33);
      ($chisq,$jnq3,$pp,$jnq4,$wmax) = split(" ",$goodstring);
      print $wres," ",$wmax," ",$chisq," ",$equad,"\n";
      $equad += $wres/40.;
      print $wres," ",$wmax," ",$chisq," ",$equad,"\n";
    }
    open(OUT2,">equad.dummy") || die "cannot open equad.dummy for writing: $!";
    printf OUT2 "EQUAD %8.4f\n", $equad;
    close(OUT2);
    system "cat equad.dummy";
    system "cat mode.dummy equad.dummy $toas > $toas.tmp";
    system "tempo -f $parfile $toas.tmp";
    system "tail -2 tempo.lis > tempo.lis.end";
    open(IN2,"tempo.lis.end");
    $_ = <IN2>;
    ($jnq,$jnq2,$jnq3,$jnq4,$wresp,$jnq5,$jnq6,$jnq7,$wres,$jnq8) = split;
    $_ = <IN2>;
    $goodstring = substr($_,33);
    ($chisq,$jnq3,$pp,$jnq4,$wmax) = split(" ",$goodstring);
    print $wres," ",$wmax," ",$chisq," ",$equad,"\n";
    $efac = sqrt($chisq);
    printf OUT "EQUAD %8.4f\n", $equad; 
#    printf OUT "EFAC %8.4f\n", $efac;
    printf OUT "EFAC 1.0\n";
    printf OUT "INCLUDE %s\n",$toas;
    system "rm $toas.tmp";
  }
  elsif(!($order =~ /MODE/) && !($order =~ /EQUAD/) && !($order =~ /EFAC/) && !($order =~ /EMIN/)) {
    print OUT
  }
}
