#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fgemm Autotuning =========
echo =================================================
echo 
echo "== Tuning fgemm over Modular<double> =="
./winograd-modular-double 2> fgemm-thresholds.h  | tee fgemm-autotune.log
echo 
echo "== Tuning fgemm over Modular<float> =="
./winograd-modular-float 2>>  fgemm-thresholds.h  | tee -a fgemm-autotune.log
echo 
echo "== Tuning fgemm over ModularBalanced<double> =="
./winograd-modularbalanced-double 2>> fgemm-thresholds.h  | tee -a fgemm-autotune.log
echo 
echo "== Tuning fgemm over ModularBalanced<float> =="
./winograd-modularbalanced-float 2>>  fgemm-thresholds.h  | tee -a fgemm-autotune.log
