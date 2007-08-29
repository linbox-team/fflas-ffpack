#!/usr/bin/perl

$fichier1 = $ARGV[0];
$fichier2 = $ARGV[1];
open(FIC1, $fichier1) or die "Impossible d'ouvrir $fichier...\n";
open(FIC2, $fichier2) or die "Impossible d'ouvrir $fichier...\n";

while(($l1 = <FIC1>) && ($l2 = <FIC2>)) {
  @t1 = split(/\s+/, $l1); $v1 = $t1[5]; $v0 = $t1[1];
  @t2 = split(/\s+/, $l2); $v2 = $t2[5];
  print "$v0 $v1 $v2\n";
}

close(FIC1);
close(FIC2);
