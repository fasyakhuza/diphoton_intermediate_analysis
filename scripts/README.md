#
To produce parquet file using Coffea Framework, you can follow below steps.
##
1. Put the `photon_to_parquet_EB_EE.py` to your folder that coffea framework installed in.

2. Activate higgs-dna environment.
```
micromamba activate higgs-dna
```

3. Set your proxy.
```
voms-proxy-init --rfc --voms cms -valid 192:00
```

4. Run the script.
```
python photon_to_parquet_EB_EE.py -y XXXX -model AAA -isTest True
```
-y : the year of the analysis. Choices for XXX: 2016preVFP, 2016postVFP, 2017, or 2018

-model : Choices for AAA: GGH, RSG

-isTest : to only scan the signal MC samples with M130, M400, and M800 for each year and model.
