#!/usr/bin/perl

open ($fh, ">", "logs/log") or die $!;

for my $p (1..5) {
  for my $n (500) {
    for my $m (1..30) {
      print $fh $p . ' ' . $n . ' ' . $m . "\n";
      system ("mpirun -n $p ./bin/a.out $n $m");
    }
  }
}
