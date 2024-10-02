import numpy as np
import awkward as ak
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
import argparse
import os
import json
import pandas as pd


def binning_width_center_1d(nbins, xlow, xhigh):
    binning_width_center = {}
    binning = np.linspace(xlow, xhigh, nbins + 1)
    width = binning[1] - binning[0]
    center = (binning[:-1] + binning[1:]) / 2
    binning_width_center["binning"]=binning
    binning_width_center["width"]=width
    binning_width_center["center"]=center
    return binning_width_center

def weighted_kinematics(var, width):
    kinematics = {}
    norm_weights = 1.0 / (sum(var) * width)
    kinematics["weighted"] = var * norm_weights
    
    error = np.sqrt(var) * norm_weights
    kinematics["ydn"] = np.array([kinematics["weighted"][i] - x for i, x in enumerate(error)])
    kinematics["yup"] = np.array([kinematics["weighted"][i] + x for i, x in enumerate(error)])
    
    return kinematics



###add parser###
parser = argparse.ArgumentParser()

#add command
#parser.add_argument("-cat", dest='category', help='Add category of the analysis.')
parser.add_argument("-y", dest='year', help='Add year of the analysis.')
parser.add_argument("-model", dest='model', help='Choose model (GGH or RSG).')
parser.add_argument("-isTest", help='The script is run for testing or not.')

#Get arguments from the user
args = parser.parse_args()


yearinput = str(args.year)
model = args.model

if args.isTest:
    if yearinput == "2016preVFP":
        if model == "GGH":
            samplefile = open('../configs/GGH_M130_M400_M800_MC_UL2016preVFP.json')
        if model == "RSG":
            samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2016preVFP.json')
    if yearinput == "2016postVFP":
        if model == "GGH":
            samplefile = open('../configs/GGH_M130_M400_M800_MC_UL2016postVFP.json')
        if model == "RSG":
            samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2016postVFP.json')
    if yearinput == "2017":
        if model == "GGH":
            samplefile = open('../configs/GGH_M130_M400_M800_MC_UL2017.json')
            #samplefile = open('../configs/M130_GluGlu_MC_UL2017.json')
        if model == "RSG":
            samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2017.json')
        #samplefile = open('../configs/M130_M400_M800_MC_UL2017.json')
    if yearinput == "2018":
        if model == "GGH":
            samplefile = open('../configs/GGH_M130_M400_M800_MC_UL2018.json')
        if model == "RSG":
            samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2018.json')
else:
    if yearinput == "2016preVFP":
        if model == "GGH":
            samplefile = open('../configs/GGHSpin0_MC_UL2016preVFP.json')
        if model == "RSG":
            samplefile = open('../configs/RSGraviton_MC_UL2016preVFP.json')
    if yearinput == "2016postVFP":
        if model == "GGH":
            samplefile = open('../configs/GGHSpin0_MC_UL2016postVFP.json')
        if model == "RSG":
            samplefile = open('../configs/RSGraviton_MC_UL2016postVFP.json')
    if yearinput == "2017":
        if model == "GGH":
            samplefile = open('../configs/GGHSpin0_MC_UL2017.json')
        if model == "RSG":
            samplefile = open('../configs/RSGraviton_MC_UL2017.json')
    if yearinput == "2018":
        if model == "GGH":
            samplefile = open('../configs/GGHSpin0_MC_UL2018.json')
        if model == "RSG":
            samplefile = open('../configs/RSGraviton_MC_UL2018.json')


# check directory
#folder="photon_kinematics"
directory="../parquet/%s/%s" %(yearinput, model)
print("Directory: ", directory)
if os.path.exists(directory):
    print("Directory '%s' exists" %(directory))
else:
    os.makedirs(directory)
    print("Directory '%s' just created" %(directory))


#samplefile = open('signallist_2017.json')

bin_dict = {
    "mvaID": {"nbins":50, "xlow":-1, "xhigh":1},
    "hoe": {"nbins":50, "xlow":0, "xhigh":1},
    "mass": {"nbins":50, "xlow":-0.1, "xhigh":0.1},
    "pt": {"nbins":50, "xlow":0, "xhigh":1000},
    "eta": {"nbins":30, "xlow":-3., "xhigh":3.},
    "phi": {"nbins":35, "xlow":-3.5, "xhigh":3.5},
    "pfRelIso03_all": {"nbins":25, "xlow":0, "xhigh":5},
    "pfRelIso03_chg": {"nbins":25, "xlow":0, "xhigh":5},
    "r9": {"nbins":40, "xlow":0, "xhigh":1.6},
    "sieie": {"nbins":50, "xlow":0, "xhigh":0.1},
    "mvaID_Fall17V1p1": {"nbins":50, "xlow":-1, "xhigh":1},
    "mvaID_WP80": {"nbins":50, "xlow":0, "xhigh":1},
    "mvaID_WP90": {"nbins":50, "xlow":0, "xhigh":1}
}

xfig, yfig = 14, 14

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass


    def process(self, events):
        dataset = events.metadata["dataset"]

        # define the histogram
        results={}
        results[dataset]={
            "count": len(events)
        }

        # Declare Hist for ECAL Barrel

        h_mvaID_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID", label="Primary dataset")
            .Reg(bin_dict["mvaID"]["nbins"], bin_dict["mvaID"]["xlow"], bin_dict["mvaID"]["xhigh"], overflow=False, underflow=False, name="x_mvaID", label = r"mvaID")
            .Weight()
        )

        h_eta_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_eta", label="Primary dataset")
            .Reg(bin_dict["eta"]["nbins"], bin_dict["eta"]["xlow"], bin_dict["eta"]["xhigh"], overflow=False, underflow=False, name="x_eta", label = r"$\eta$")
            .Weight()
        )

        h_hoe_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_hoe", label="Primary dataset")
            .Reg(bin_dict["hoe"]["nbins"], bin_dict["hoe"]["xlow"], bin_dict["hoe"]["xhigh"], overflow=False, underflow=False, name="x_hoe", label = r"H/E")
            .Weight()
        )

        h_mass_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_mass", label="Primary dataset")
            .Reg(bin_dict["mass"]["nbins"], bin_dict["mass"]["xlow"], bin_dict["mass"]["xhigh"], overflow=False, underflow=False, name="x_mass", label = r"m_$\gamma$ (GeV)")
            .Weight()
        )
        
        h_pfRelIso03_all_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_all", label="Primary dataset")
            .Reg(bin_dict["pfRelIso03_all"]["nbins"], bin_dict["pfRelIso03_all"]["xlow"], bin_dict["pfRelIso03_all"]["xhigh"], overflow=False, underflow=False, name="x_pfRelIso03_all", label = r"pfRelIso03_all")
            .Weight()
        )

        h_pfRelIso03_chg_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_chg", label="Primary dataset")
            .Reg(bin_dict["pfRelIso03_chg"]["nbins"], bin_dict["pfRelIso03_chg"]["xlow"], bin_dict["pfRelIso03_chg"]["xhigh"], overflow=False, underflow=False, name="x_pfRelIso03_chg", label = r"pfRelIso03_chg")
            .Weight()
        )

        h_phi_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_phi", label="Primary dataset")
            .Reg(bin_dict["phi"]["nbins"], bin_dict["phi"]["xlow"], bin_dict["phi"]["xhigh"], overflow=False, underflow=False, name="x_phi", label = r"$\phi$")
            .Weight()
        )

        h_pt_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_pt", label="Primary dataset")
            .Reg(bin_dict["pt"]["nbins"], bin_dict["pt"]["xlow"], bin_dict["pt"]["xhigh"], overflow=False, underflow=False, name="x_pt", label = r"$\pt$")
            .Weight()
        )

        h_r9_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_r9", label="Primary dataset")
            .Reg(bin_dict["r9"]["nbins"], bin_dict["r9"]["xlow"], bin_dict["r9"]["xhigh"], overflow=False, underflow=False, name="x_r9", label = r"$\r9$")
            .Weight()
        )

        h_sieie_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_sieie", label="Primary dataset")
            .Reg(bin_dict["sieie"]["nbins"], bin_dict["sieie"]["xlow"], bin_dict["sieie"]["xhigh"], overflow=False, underflow=False, name="x_sieie", label = r"$\sigma_{i\etai\eta}$")
            .Weight()
        )
        
        h_mvaID_Fall17V1p1_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_Fall17V1p1", label="Primary dataset")
            .Reg(bin_dict["mvaID_Fall17V1p1"]["nbins"], bin_dict["mvaID_Fall17V1p1"]["xlow"], bin_dict["mvaID_Fall17V1p1"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_Fall17V1p1", label = r"mvaID_Fall17V1p1")
            .Weight()
        )

        h_mvaID_WP80_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP80", label="Primary dataset")
            .Reg(bin_dict["mvaID_WP80"]["nbins"], bin_dict["mvaID_WP80"]["xlow"], bin_dict["mvaID_WP80"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_WP80", label = r"mvaID_WP80")
            .Weight()
        )

        h_mvaID_WP90_EB = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP90", label="Primary dataset")
            .Reg(bin_dict["mvaID_WP90"]["nbins"], bin_dict["mvaID_WP90"]["xlow"], bin_dict["mvaID_WP90"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_WP90", label = r"mvaID_WP90")
            .Weight()
        )
        

        # Declare Hist for ECAL Endcap

        h_mvaID_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID", label="Primary dataset")
            .Reg(bin_dict["mvaID"]["nbins"], bin_dict["mvaID"]["xlow"], bin_dict["mvaID"]["xhigh"], overflow=False, underflow=False, name="x_mvaID", label = r"mvaID")
            .Weight()
        )

        h_eta_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_eta", label="Primary dataset")
            .Reg(bin_dict["eta"]["nbins"], bin_dict["eta"]["xlow"], bin_dict["eta"]["xhigh"], overflow=False, underflow=False, name="x_eta", label = r"$\eta$")
            .Weight()
        )

        h_hoe_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_hoe", label="Primary dataset")
            .Reg(bin_dict["hoe"]["nbins"], bin_dict["hoe"]["xlow"], bin_dict["hoe"]["xhigh"], overflow=False, underflow=False, name="x_hoe", label = r"H/E")
            .Weight()
        )

        h_mass_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_mass", label="Primary dataset")
            .Reg(bin_dict["mass"]["nbins"], bin_dict["mass"]["xlow"], bin_dict["mass"]["xhigh"], overflow=False, underflow=False, name="x_mass", label = r"m_$\gamma$ (GeV)")
            .Weight()
        )
        
        h_pfRelIso03_all_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_all", label="Primary dataset")
            .Reg(bin_dict["pfRelIso03_all"]["nbins"], bin_dict["pfRelIso03_all"]["xlow"], bin_dict["pfRelIso03_all"]["xhigh"], overflow=False, underflow=False, name="x_pfRelIso03_all", label = r"pfRelIso03_all")
            .Weight()
        )

        h_pfRelIso03_chg_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_chg", label="Primary dataset")
            .Reg(bin_dict["pfRelIso03_chg"]["nbins"], bin_dict["pfRelIso03_chg"]["xlow"], bin_dict["pfRelIso03_chg"]["xhigh"], overflow=False, underflow=False, name="x_pfRelIso03_chg", label = r"pfRelIso03_chg")
            .Weight()
        )

        h_phi_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_phi", label="Primary dataset")
            .Reg(bin_dict["phi"]["nbins"], bin_dict["phi"]["xlow"], bin_dict["phi"]["xhigh"], overflow=False, underflow=False, name="x_phi", label = r"$\phi$")
            .Weight()
        )

        h_pt_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_pt", label="Primary dataset")
            .Reg(bin_dict["pt"]["nbins"], bin_dict["pt"]["xlow"], bin_dict["pt"]["xhigh"], overflow=False, underflow=False, name="x_pt", label = r"$\pt$")
            .Weight()
        )

        h_r9_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_r9", label="Primary dataset")
            .Reg(bin_dict["r9"]["nbins"], bin_dict["r9"]["xlow"], bin_dict["r9"]["xhigh"], overflow=False, underflow=False, name="x_r9", label = r"$\r9$")
            .Weight()
        )

        h_sieie_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_sieie", label="Primary dataset")
            .Reg(bin_dict["sieie"]["nbins"], bin_dict["sieie"]["xlow"], bin_dict["sieie"]["xhigh"], overflow=False, underflow=False, name="x_sieie", label = r"$\sigma_{i\etai\eta}$")
            .Weight()
        )
        
        h_mvaID_Fall17V1p1_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_Fall17V1p1", label="Primary dataset")
            .Reg(bin_dict["mvaID_Fall17V1p1"]["nbins"], bin_dict["mvaID_Fall17V1p1"]["xlow"], bin_dict["mvaID_Fall17V1p1"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_Fall17V1p1", label = r"mvaID_Fall17V1p1")
            .Weight()
        )

        h_mvaID_WP80_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP80", label="Primary dataset")
            .Reg(bin_dict["mvaID_WP80"]["nbins"], bin_dict["mvaID_WP80"]["xlow"], bin_dict["mvaID_WP80"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_WP80", label = r"mvaID_WP80")
            .Weight()
        )

        h_mvaID_WP90_EE = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP90", label="Primary dataset")
            .Reg(bin_dict["mvaID_WP90"]["nbins"], bin_dict["mvaID_WP90"]["xlow"], bin_dict["mvaID_WP90"]["xhigh"], overflow=False, underflow=False, name="x_mvaID_WP90", label = r"mvaID_WP90")
            .Weight()
        )

        
        # get photons
        photons = events.Photon    
       
        # add selections
        # for ECAL Barrel
        #if args.category == "EB":
        photons_cat_EB = photons[(photons.isScEtaEB)]
        
        count_number_of_photons_EB = ak.num(photons_cat_EB, axis=1)

        photons_cat_selected_EB = photons_cat_EB[count_number_of_photons_EB > 0]


        # fill histogram ECAL Barrel
        h_mvaID_EB.fill(dataset_mvaID=dataset,x_mvaID=photons_cat_selected_EB.mvaID[:,0]) 
        results["mvaID_EB"] = h_mvaID_EB

        h_eta_EB.fill(dataset_eta=dataset,x_eta=photons_cat_selected_EB.eta[:,0])
        results["eta_EB"] = h_eta_EB

        h_hoe_EB.fill(dataset_hoe=dataset,x_hoe=photons_cat_selected_EB.hoe[:,0])
        results["hoe_EB"] = h_hoe_EB

        h_mass_EB.fill(dataset_mass=dataset,x_mass=photons_cat_selected_EB.mass[:,0])
        results["mass_EB"] = h_mass_EB

        h_pfRelIso03_all_EB.fill(dataset_pfRelIso03_all=dataset,x_pfRelIso03_all=photons_cat_selected_EB.pfRelIso03_all[:,0])
        results["pfRelIso03_all_EB"] = h_pfRelIso03_all_EB
        
        h_pfRelIso03_chg_EB.fill(dataset_pfRelIso03_chg=dataset,x_pfRelIso03_chg=photons_cat_selected_EB.pfRelIso03_chg[:,0])
        results["pfRelIso03_chg_EB"] = h_pfRelIso03_chg_EB

        h_phi_EB.fill(dataset_phi=dataset,x_phi=photons_cat_selected_EB.phi[:,0])
        results["phi_EB"] = h_phi_EB

        h_pt_EB.fill(dataset_pt=dataset,x_pt=photons_cat_selected_EB.pt[:,0])
        results["pt_EB"] = h_pt_EB
        
        h_r9_EB.fill(dataset_r9=dataset,x_r9=photons_cat_selected_EB.r9[:,0])
        results["r9_EB"] = h_r9_EB

        h_sieie_EB.fill(dataset_sieie=dataset,x_sieie=photons_cat_selected_EB.sieie[:,0])
        results["sieie_EB"] = h_sieie_EB

        h_mvaID_Fall17V1p1_EB.fill(dataset_mvaID_Fall17V1p1=dataset,x_mvaID_Fall17V1p1=photons_cat_selected_EB.mvaID_Fall17V1p1[:,0])
        results["mvaID_Fall17V1p1_EB"] = h_mvaID_Fall17V1p1_EB
        
        h_mvaID_WP80_EB.fill(dataset_mvaID_WP80=dataset,x_mvaID_WP80=photons_cat_selected_EB.mvaID_WP80[:,0])
        results["mvaID_WP80_EB"] = h_mvaID_WP80_EB

        h_mvaID_WP90_EB.fill(dataset_mvaID_WP90=dataset,x_mvaID_WP90=photons_cat_selected_EB.mvaID_WP90[:,0])
        results["mvaID_WP90_EB"] = h_mvaID_WP90_EB


        # for ECAL Endcap
        #if args.category == "EE":
        photons_cat_EE = photons[(photons.isScEtaEE)]
    
        count_number_of_photons_EE = ak.num(photons_cat_EE, axis=1)

        photons_cat_selected_EE = photons_cat_EE[count_number_of_photons_EE > 0]


        # fill histogram ECAL Endcap
        h_mvaID_EE.fill(dataset_mvaID=dataset,x_mvaID=photons_cat_selected_EE.mvaID[:,0]) 
        results["mvaID_EE"] = h_mvaID_EE

        h_eta_EE.fill(dataset_eta=dataset,x_eta=photons_cat_selected_EE.eta[:,0])
        results["eta_EE"] = h_eta_EE

        h_hoe_EE.fill(dataset_hoe=dataset,x_hoe=photons_cat_selected_EE.hoe[:,0])
        results["hoe_EE"] = h_hoe_EE

        h_mass_EE.fill(dataset_mass=dataset,x_mass=photons_cat_selected_EE.mass[:,0])
        results["mass_EE"] = h_mass_EE

        h_pfRelIso03_all_EE.fill(dataset_pfRelIso03_all=dataset,x_pfRelIso03_all=photons_cat_selected_EE.pfRelIso03_all[:,0])
        results["pfRelIso03_all_EE"] = h_pfRelIso03_all_EE
        
        h_pfRelIso03_chg_EE.fill(dataset_pfRelIso03_chg=dataset,x_pfRelIso03_chg=photons_cat_selected_EE.pfRelIso03_chg[:,0])
        results["pfRelIso03_chg_EE"] = h_pfRelIso03_chg_EE

        h_phi_EE.fill(dataset_phi=dataset,x_phi=photons_cat_selected_EE.phi[:,0])
        results["phi_EE"] = h_phi_EE

        h_pt_EE.fill(dataset_pt=dataset,x_pt=photons_cat_selected_EE.pt[:,0])
        results["pt_EE"] = h_pt_EE
        
        h_r9_EE.fill(dataset_r9=dataset,x_r9=photons_cat_selected_EE.r9[:,0])
        results["r9_EE"] = h_r9_EE

        h_sieie_EE.fill(dataset_sieie=dataset,x_sieie=photons_cat_selected_EE.sieie[:,0])
        results["sieie_EE"] = h_sieie_EE

        h_mvaID_Fall17V1p1_EE.fill(dataset_mvaID_Fall17V1p1=dataset,x_mvaID_Fall17V1p1=photons_cat_selected_EE.mvaID_Fall17V1p1[:,0])
        results["mvaID_Fall17V1p1_EE"] = h_mvaID_Fall17V1p1_EE
        
        h_mvaID_WP80_EE.fill(dataset_mvaID_WP80=dataset,x_mvaID_WP80=photons_cat_selected_EE.mvaID_WP80[:,0])
        results["mvaID_WP80_EE"] = h_mvaID_WP80_EE

        h_mvaID_WP90_EE.fill(dataset_mvaID_WP90=dataset,x_mvaID_WP90=photons_cat_selected_EE.mvaID_WP90[:,0])
        results["mvaID_WP90_EE"] = h_mvaID_WP90_EE

        return results

    def postprocess(self, accumulant):
        pass

#########
# get input samples

sample_dict = json.load(samplefile)

# run processor

run = processor.Runner(
    # executor=processor.IterativeExecutor(),
    executor=processor.FuturesExecutor(workers=4), # user 4 cores
    schema=NanoAODSchema
)

results = run(
    sample_dict,
    treename="Events",
    processor_instance=MyProcessor(),
)


# results
#print(results)


#############################

variables = ["mvaID", "eta", "hoe", "mass", "pfRelIso03_all", "pfRelIso03_chg", "phi", "pt", "r9", "sieie", "mvaID_Fall17V1p1"]

for key in sample_dict.keys():
    print("Sample: " + key)
    key_ntuple={}

    # for ECAL Barrel
    for variable in variables:
        var_string = "EB_"+str(variable)
        var_results = "%s_EB" %(variable)
        dataset_var = "dataset_%s" %(variable) #(f'dataset_{variable}')

        key_var_hist = results[var_results][{dataset_var:key}]
        key_var_value = key_var_hist.values()
        key_var_edges = key_var_hist.axes.edges
        
        key_ntuple[var_string]=[{
            'bin_center':[(key_var_edges[0][i] + key_var_edges[0][i+1]) / 2 for i in range(len(key_var_edges[0])-1)],
            'value':key_var_value,
            'nbins': [bin_dict[variable]["nbins"]],
            'xlow': [bin_dict[variable]["xlow"]],
            'xhigh': [bin_dict[variable]["xhigh"]]

        #key_ntuple[var_string]['bin_center']=
        }] #if include '[]' it will create bin_center, nbins, xlow, and xhigh as subfields. But, if not, those will being as index at 0 level

    # for ECAL Endcap
    for variable in variables:
        var_string = "EE_"+str(variable)
        var_results = "%s_EE" %(variable)
        dataset_var = "dataset_%s" %(variable) #(f'dataset_{variable}')
        
        key_var_hist = results[var_results][{dataset_var:key}]
        key_var_value = key_var_hist.values()
        key_var_edges = key_var_hist.axes.edges
        
        key_ntuple[var_string]=[{
            'bin_center':[(key_var_edges[0][i] + key_var_edges[0][i+1]) / 2 for i in range(len(key_var_edges[0])-1)],
            'value':key_var_value,
            'nbins': [bin_dict[variable]["nbins"]],
            'xlow': [bin_dict[variable]["xlow"]],
            'xhigh': [bin_dict[variable]["xhigh"]]

        #key_ntuple[var_string]['bin_center']=
        }] #if include '[]' it will create bin_center, nbins, xlow, and xhigh as subfields. But, if not, those will being as index at 0 level

    df_key = pd.DataFrame(data=key_ntuple)
    df_key.to_parquet("%s/%s.parquet" %(directory, key))
    # if args.category == "EE":
    #     df_key.to_parquet(f'../parquet/{key}_EE.parquet')
    # if args.category == "EB":
    #     df_key.to_parquet(f'../parquet/{key}_EB.parquet')