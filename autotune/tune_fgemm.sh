#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fgemm Autotuning =========
echo =================================================
echo 
echo "== Tuning fgemm over Modular<double> =="
(./winograd-modular-double > fgemm-thresholds.h) 2>&1 | tee fgemm-autotune.log
echo 
echo "== Tuning fgemm over Modular<float> =="
(./winograd-modular-float >>  fgemm-thresholds.h) 2>&1 | tee -a fgemm-autotune.log
echo 
echo "== Tuning fgemm over ModularBalanced<double> =="
(./winograd-modularbalanced-double >> fgemm-thresholds.h) 2>&1 | tee -a fgemm-autotune.log
echo 
echo "== Tuning fgemm over ModularBalanced<float> =="
(./winograd-modularbalanced-float >>  fgemm-thresholds.h) 2>&1 | tee -a fgemm-autotune.log
