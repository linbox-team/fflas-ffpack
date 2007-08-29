#!/bin/csh -f

set block = $1

#foreach i (10 12 14 16 18 20 22 24 26 28 30 50 75 100 125 150)
foreach i (175 200 225)
	@ N = $i * $block
	tblockmat-$block $N 1 
end
	
