# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_10_patch1

cd CMSSW_15_0_10_patch1/src

cmsenv

git cms-addpkg DataFormats/PatCandidates ; git fetch git@github.com:stahlleiton/cmssw.git ParticleAnalyzer_CMSSW_15_1_X ; git cherry-pick ce0e4ae41f60f84dc814aa94ff39e89419da18ef

git clone git@github.com:stahlleiton/VertexCompositeAnalysis.git -b ParticleFitter_15_0_X

scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2025_SingleObject_MC_ParticleAnalyzer_cfg.py
