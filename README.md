# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_13_2_11

cd CMSSW_13_2_11/src

cmsenv

git clone -b 13_2_X https://github.com/stahlleiton/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2023_DiMuContBothGammaGamma_mc_cfg.py
