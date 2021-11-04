# Xray-analysis
Analysis of X-ray data

Script originally developed to run the analysis of the first extended X-ray campaign (January-March 2021): python fitVfb.py
Then it was modified to be able to manage also the annealing part, but its main purpose is still to analyse the irradiation data.

### Get the code

```
export SCRAM_ARCH="slc7_amd64_gcc900"  # for bash
cmsrel CMSSW_11_3_1
cd cmsrel CMSSW_11_3_1/src/
cmsenv
```
Set this environment variable with your github use (or set it by hand if it doesn't work)
Currently the code is in the branch named __mciprian_devel__ from the repository of cippy
```
YOUR_GITHUB_REPOSITORY=$(git config user.github)
git remote add cippy git@github.com:cippy/Xray-analysis.git
git clone -o cippy git@github.com:cippy/Xray-analysis.git -b mciprian_devel xRayAnalysis
cd xRayAnalysis
git remote add origin git@github.com:$YOUR_GITHUB_REPOSITORY/Xray-analysis.git

### Examples

Testing some options to customize output
```
python3 fitVfb.py -o /some/output/folder/on/eos/ --outfiles /local/folder/for/output/rootfiles/ --skip MOS --samples 1012_UL
```

For annealing studies, use a command like the following
```
python3 fitVfb.py -o /some/output/folder/on/eos/MOShalf/ --skip GCD --samples 1011_LR -a --annealing-path-regexp "1011_LR_ann_60min_60C_MOShalf_20C.*"
```


