#!/usr/bin/perl

# Copyright (c) FFLAS-FFPACK
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/

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
