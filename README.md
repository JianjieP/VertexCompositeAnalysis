# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_15_0_10_patch1
cd CMSSW_15_0_10_patch1/src

cmsenv

git cms-addpkg DataFormats/PatCandidates

git remote add cmssw git@github.com:stahlleiton/cmssw.git

git fetch cmssw --no-tags

git cherry-pick 91562810cee10bac976afe057a878addc088493e

git clone git@github.com:stahlleiton/VertexCompositeAnalysis.git -b ParticleFitter_15_0_X

scram b -j8

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2025_SingleObject_MC_ParticleAnalyzer_cfg.py
