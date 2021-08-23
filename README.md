# Xray-analysis
Analysis of X-ray data

to run the analysis of the first extended X-ray campaign (January-March 2021): python fitVfb.py

the script for the new setup (fitNewSetup.py) is still in preparation

### Examples

Testing some options to customize output
```
python3 fitVfb.py -o /eos/user/m/mciprian/www/xRayAnalysis/testCode/test_11August2021/ --outfiles toremove2 --skip MOS --samples 1012_UL
```

For annealing studies, use a command like the following
```
python3 fitVfb.py -o /eos/user/m/mciprian/www/xRayAnalysis/testCode/test_17August2021/MOShalf/ --skip GCD --samples 1011_LR -a --annealing-path-regexp "1011_LR_ann_60min_60C_MOShalf_20C.*"
```