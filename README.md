# Pi0Tuplizer
Simple tuplizer for low energy photon clusters (3x3) and diphoton events (pi0/eta)
=============================
-----------------------------
setup
-----------------------------
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
if [ -r CMSSW_8_0_25/src ] ; then
 echo release CMSSW_8_0_25 already exists
else
scram p CMSSW CMSSW_8_0_25
fi
cd CMSSW_8_0_25/src
eval `scram runtime -sh`
git clone git@github.com:RazorCMS/Pi0Tuplizer.git PiZero/Pi0Tuplizer
scram b
```
-----------------------------
cmsRun
-----------------------------
```
cd python
cmsRun Pi0Tuplizer_80X_AlCaP0_RAW.py
```

-----------------------------
submit crab jobs
-----------------------------
```
crab submit -c crab.py
```


