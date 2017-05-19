#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK CharPoly Autotuning =========
echo =================================================
echo 
./charpoly-LUK-ArithProg 2> charpoly-LUK-ArithProg-threshold.h  | tee charpoly-LUK-ArithProg-autotune.log
./charpoly-Danilevskii-LUK 2> charpoly-Danilevskii-LUK-threshold.h  | tee charpoly-Danilevskii-LUK-autotune.log
./arithprog 8 160 4 4 1000 3 2> arithprog-blocksize.h  | tee arithprog-blocksize-autotune.log
