# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_14_1_9_patch2

cd CMSSW_14_1_9_patch2/src

cmsenv

git clone -b ParticleFitter_14_1_X https://github.com/stahlleiton/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

cmsRun VCTree_PbPb2023_UPCDiKa_UPCReco_cfg.py
