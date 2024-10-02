put the photon_to_parquet_EB_EE.py to your folder that coffea framework installed in
activate higgs-dna
micromamba activate higgs-dna
set your proxy
run the script
python photon_to_parquet_EB_EE.py -y XXXX -model AAA -isTest True
-y : the year of the analysis. Choices: 2016preVFP, 2016postVFP, 2017, or 2018
-model : Choices: GGH, RSG
-isTest : to only scan the signal MC samples with M130, M400, and M800 for each year and model
