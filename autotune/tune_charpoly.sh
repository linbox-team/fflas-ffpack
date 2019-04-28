#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK CharPoly Autotuning =========
echo =================================================
echo 
(./arithprog 16 64 2 4 800 3 > arithprog-blocksize.h) 2>&1 | tee arithprog-blocksize-autotune.log

(./charpoly-LUK-ArithProg > charpoly-LUK-ArithProg-threshold.h) 2>&1 | tee charpoly-LUK-ArithProg-autotune.log

(./charpoly-Danilevskii-LUK > charpoly-Danilevskii-LUK-threshold.h) 2>&1 | tee charpoly-Danilevskii-LUK-autotune.log
