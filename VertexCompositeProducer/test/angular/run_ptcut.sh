#!/bin/bash

echo "==== Start ===="

# 第一个任务
(
  cd /afs/cern.ch/user/j/jianjie/scratch1/CMSSW_15_1_0_patch3/src/VertexCompositeAnalysis/VertexCompositeProducer/test/angular || exit
  echo "Running fourpair_rho_x_ptcut.C"
  root -l -b -q fourpair_rho_x_ptcut.C
)

# 第二个任务
(
  cd /afs/cern.ch/user/j/jianjie/scratch1/CMSSW_15_1_0_patch3/src/VertexCompositeAnalysis/VertexCompositeProducer/test/angular/fourpair_rho_x_ptcut_PNG || exit
  echo "Running plot_pt_vs_fitparams.C"
  root -l -b -q plot_pt_vs_fitparams.C
)

echo "==== Done ===="