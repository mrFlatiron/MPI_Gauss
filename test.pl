#!/usr/bin/perl

for my $p (1..5) {
  for my $n (500) {
    for my $m (1..30) {
      print $p . ' ' . $n . ' ' . $m . "\n";
      system ("mpirun -n $p ./app_exec $n $m");
    }
  }
}
