# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_10_patch1
cd CMSSW_15_0_10_patch1/src

cmsenv

git cms-addpkg DataFormats/PatCandidates ; git fetch cmssw git@github.com:stahlleiton/cmssw.git ParticleAnalyzer_CMSSW_15_1_X ; git cherry-pick 198311ccf2a5943c3086868838a49cd77d6bc5da

git clone git@github.com:stahlleiton/VertexCompositeAnalysis.git -b ParticleFitter_15_0_X

scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2025_SingleObject_MC_ParticleAnalyzer_cfg.py
