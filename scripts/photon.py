import numpy as np
import awkward as ak
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
import argparse
import os
import json


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
parser.add_argument("-cat", dest='category', help='Add category of the analysis.')
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
        if model == "RSG":
            samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2017.json')
        #samplefile = open('../configs/M130_M400_M800_MC_UL2017.json')
        #samplefile = open('../configs/M130_GluGlu_MC_UL2017.json')
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
folder=".../plots"
subfolder="photon_kinematics"
directory = os.path.join(folder, yearinput, subfolder, model, "/")
print("Directory: ", directory)
if os.path.exists(directory):
    print("Directory '%s' exists" %directory)
else:
    os.makedirs(directory)


#samplefile = open('signallist_2017.json')


#variables = ["mvaID", "eta", "hoe", "mass", "pfRelIso03_all", "pfRelIso03_chg", "phi", "pt", "r9", "sieie", "mvaID_Fall17V1p1", "mvaID_WP80", "mvaID_WP90"]
variables = ["mvaID", "eta"]

nbins_mvaID, nbins_hoe, nbins_mass, nbins_pt, nbins_mvaID_Fall17V1p1, nbins_mvaID_WP80, nbins_mvaID_WP90 = 50, 50, 50, 50, 50, 50, 50
nbins_eta = 30
nbins_phi = 35
nbins_pfRelIso03_all, nbins_pfRelIso03_chg = 40, 40
nbins_r9 = 40
nbins_sieie = 50

xlow_hoe, xlow_pfRelIso03_all, xlow_pfRelIso03_chg, xlow_pt, xlow_r9, xlow_sieie, xlow_mvaID_WP80, xlow_mvaID_WP90 = 0, 0, 0, 0, 0, 0, 0, 0
xlow_mvaID, xlow_mvaID_Fall17V1p1 = -1., -1.
xlow_mass = -0.1
xlow_eta = -3.
xlow_phi = -3.5
xhigh_mvaID, xhigh_hoe, xhigh_mvaID_Fall17V1p1, xhigh_mvaID_WP80, xhigh_mvaID_WP90 = 1, 1, 1, 1, 1
xhigh_eta=3.
xhigh_pfRelIso03_all, xhigh_pfRelIso03_chg = 4, 4
xhigh_phi = 3.5
xhigh_pt = 1000
xhigh_mass = 0.1
xhigh_r9 = 1.6
xhigh_sieie = 0.1

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
            
        h_mvaID = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID", label="Primary dataset")
            .Reg(nbins_mvaID, xlow_mvaID, xhigh_mvaID, overflow=False, underflow=False, name="x_mvaID", label = r"mvaID")
            .Weight()
        )

        h_eta = (
            Hist.new.StrCat([], growth=True, name="dataset_eta", label="Primary dataset")
            .Reg(nbins_eta, xlow_eta, xhigh_eta, overflow=False, underflow=False, name="x_eta", label = r"$\eta$")
            .Weight()
        )

        h_hoe = (
            Hist.new.StrCat([], growth=True, name="dataset_hoe", label="Primary dataset")
            .Reg(nbins_hoe, xlow_hoe, xhigh_hoe, overflow=False, underflow=False, name="x_hoe", label = r"H/E")
            .Weight()
        )

        h_mass = (
            Hist.new.StrCat([], growth=True, name="dataset_mass", label="Primary dataset")
            .Reg(nbins_mass, xlow_mass, xhigh_mass, overflow=False, underflow=False, name="x_mass", label = r"m_$\gamma$ (GeV)")
            .Weight()
        )
        
        h_pfRelIso03_all = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_all", label="Primary dataset")
            .Reg(nbins_pfRelIso03_all, xlow_pfRelIso03_all, xhigh_pfRelIso03_all, overflow=False, underflow=False, name="x_pfRelIso03_all", label = r"pfRelIso03_all")
            .Weight()
        )

        h_pfRelIso03_chg = (
            Hist.new.StrCat([], growth=True, name="dataset_pfRelIso03_chg", label="Primary dataset")
            .Reg(nbins_pfRelIso03_chg, xlow_pfRelIso03_chg, xhigh_pfRelIso03_chg, overflow=False, underflow=False, name="x_pfRelIso03_chg", label = r"pfRelIso03_chg")
            .Weight()
        )

        h_phi = (
            Hist.new.StrCat([], growth=True, name="dataset_phi", label="Primary dataset")
            .Reg(nbins_phi, xlow_phi, xhigh_phi, overflow=False, underflow=False, name="x_phi", label = r"$\phi$")
            .Weight()
        )

        h_pt = (
            Hist.new.StrCat([], growth=True, name="dataset_pt", label="Primary dataset")
            .Reg(nbins_pt, xlow_pt, xhigh_pt, overflow=False, underflow=False, name="x_pt", label = r"$\pt$")
            .Weight()
        )

        h_r9 = (
            Hist.new.StrCat([], growth=True, name="dataset_r9", label="Primary dataset")
            .Reg(nbins_r9, xlow_r9, xhigh_r9, overflow=False, underflow=False, name="x_r9", label = r"$\r9$")
            .Weight()
        )

        h_sieie = (
            Hist.new.StrCat([], growth=True, name="dataset_sieie", label="Primary dataset")
            .Reg(nbins_sieie, xlow_sieie, xhigh_sieie, overflow=False, underflow=False, name="x_sieie", label = r"$\sigma_{i\etai\eta}$")
            .Weight()
        )
        
        h_mvaID_Fall17V1p1 = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_Fall17V1p1", label="Primary dataset")
            .Reg(nbins_mvaID_Fall17V1p1, xlow_mvaID_Fall17V1p1, xhigh_mvaID_Fall17V1p1, overflow=False, underflow=False, name="x_mvaID_Fall17V1p1", label = r"mvaID_Fall17V1p1")
            .Weight()
        )

        h_mvaID_WP80 = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP80", label="Primary dataset")
            .Reg(nbins_mvaID_WP80, xlow_mvaID_WP80, xhigh_mvaID_WP80, overflow=False, underflow=False, name="x_mvaID_WP80", label = r"mvaID_WP80")
            .Weight()
        )

        h_mvaID_WP90 = (
            Hist.new.StrCat([], growth=True, name="dataset_mvaID_WP90", label="Primary dataset")
            .Reg(nbins_mvaID_WP90, xlow_mvaID_WP90, xhigh_mvaID_WP90, overflow=False, underflow=False, name="x_mvaID_WP90", label = r"mvaID_WP90")
            .Weight()
        )

        
        # get photons
        photons = events.Photon    
       
        # add selections
        if args.category == "EE":
            photons_cat = photons[(photons.isScEtaEE)]
        
        if args.category == "EB":
            photons_cat = photons[(photons.isScEtaEB)]
        
        count_number_of_photons = ak.num(photons_cat, axis=1)

        photons_cat_selected = photons_cat[count_number_of_photons > 0]


        # fill histogram
        h_mvaID.fill(dataset_mvaID=dataset,x_mvaID=photons_cat_selected.mvaID[:,0]) 
        results["mvaID"] = h_mvaID

        h_eta.fill(dataset_eta=dataset,x_eta=photons_cat_selected.eta[:,0])
        results["eta"] = h_eta

        h_hoe.fill(dataset_hoe=dataset,x_hoe=photons_cat_selected.hoe[:,0])
        results["hoe"] = h_hoe

        h_mass.fill(dataset_mass=dataset,x_mass=photons_cat_selected.mass[:,0])
        results["mass"] = h_mass

        h_pfRelIso03_all.fill(dataset_pfRelIso03_all=dataset,x_pfRelIso03_all=photons_cat_selected.pfRelIso03_all[:,0])
        results["pfRelIso03_all"] = h_pfRelIso03_all
        
        h_pfRelIso03_chg.fill(dataset_pfRelIso03_chg=dataset,x_pfRelIso03_chg=photons_cat_selected.pfRelIso03_chg[:,0])
        results["pfRelIso03_chg"] = h_pfRelIso03_chg

        h_phi.fill(dataset_phi=dataset,x_phi=photons_cat_selected.phi[:,0])
        results["phi"] = h_phi

        h_pt.fill(dataset_pt=dataset,x_pt=photons_cat_selected.pt[:,0])
        results["pt"] = h_pt
        
        h_r9.fill(dataset_r9=dataset,x_r9=photons_cat_selected.r9[:,0])
        results["r9"] = h_r9

        h_sieie.fill(dataset_sieie=dataset,x_sieie=photons_cat_selected.sieie[:,0])
        results["sieie"] = h_sieie

        h_mvaID_Fall17V1p1.fill(dataset_mvaID_Fall17V1p1=dataset,x_mvaID_Fall17V1p1=photons_cat_selected.mvaID_Fall17V1p1[:,0])
        results["mvaID_Fall17V1p1"] = h_mvaID_Fall17V1p1
        
        h_mvaID_WP80.fill(dataset_mvaID_WP80=dataset,x_mvaID_WP80=photons_cat_selected.mvaID_WP80[:,0])
        results["mvaID_WP80"] = h_mvaID_WP80

        h_mvaID_WP90.fill(dataset_mvaID_WP90=dataset,x_mvaID_WP90=photons_cat_selected.mvaID_WP90[:,0])
        results["mvaID_WP90"] = h_mvaID_WP90

        return results

    def postprocess(self, accumulant):
        pass

#########
# cell 20

sample_dict = json.load(samplefile)

#sample_dict = {
#    "glugluSpin0toGG_W-0p014_M-650_2017":[
#        "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluSpin0ToGammaGamma_W-0p014_M-650_TuneCP2_13TeV_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2560000/7C5BA93A-8B0D-1F44-A2CA-6E9577FBAE99.root",
#        "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluSpin0ToGammaGamma_W-0p014_M-650_TuneCP2_13TeV_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2560000/93E2ACF4-DC11-EE4F-908C-9D9B5882C3AC.root",
#        "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluSpin0ToGammaGamma_W-0p014_M-650_TuneCP2_13TeV_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2560000/A21BCAD7-126F-CF4E-B61C-36836A09E754.root",
#        "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluSpin0ToGammaGamma_W-0p014_M-650_TuneCP2_13TeV_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2560000/E2959584-3BF8-D943-84A1-E548CF67452C.root",
#        "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluSpin0ToGammaGamma_W-0p014_M-650_TuneCP2_13TeV_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2560000/F4CDCC58-B38F-494B-A977-9AA7DBDBA5AF.root"
#    ]#,
#    #"DY": [
#    #    "/eos/user/f/fkhuzaim/higgsdna_finalfits_tutorial_24/01_columnar_introduction/DY/44449913-A842-E811-9863-0CC47A545060.root"
#    #]
#}

# cell 21

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

# cell 22

print(results)

# cell 23 
# get kinematics from hist and get the density
# plot histogram

import mplhep as hep
import matplotlib.pyplot as plt

hep.style.use(hep.style.CMS)

'''for key in sample_dict.keys():
    print("Sample: " + key)

    f_mvaID, ax_mvaID = plt.subplots(figsize=(10,10))
    ax_mvaID.set_ylabel("Counts")
    #results["mvaID"][{"dataset_mvaID":"glugluSpin0toGG_W-0p014_M-650_2017"}].plot(ax=ax,label="Spin0toGG")
    results["mvaID"][{"dataset_mvaID":key}].plot(ax=ax_mvaID,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    plt.legend()

    f_eta, ax_eta = plt.subplots(figsize=(10,10))
    ax_eta.set_ylabel("Counts")
    #results["eta"][{"dataset_eta":"glugluSpin0toGG_W-0p014_M-650_2017"}].plot(ax=ax_eta,label="Spin0toGG")
    results["eta"][{"dataset_eta":key}].plot(ax=ax_eta,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_hoe, ax_hoe = plt.subplots(figsize=(10,10))
    ax_hoe.set_ylabel("Counts")
    results["hoe"][{"dataset_hoe":key}].plot(ax=ax_hoe,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_mass, ax_mass = plt.subplots(figsize=(10,10))
    ax_mass.set_ylabel("Counts")
    results["mass"][{"dataset_mass":key}].plot(ax=ax_mass,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_pfRelIso03_all, ax_pfRelIso03_all = plt.subplots(figsize=(10,10))
    ax_pfRelIso03_all.set_ylabel("Counts")
    results["pfRelIso03_all"][{"dataset_pfRelIso03_all":key}].plot(ax=ax_pfRelIso03_all,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_pfRelIso03_chg, ax_pfRelIso03_chg = plt.subplots(figsize=(10,10))
    ax_pfRelIso03_chg.set_ylabel("Counts")
    results["pfRelIso03_chg"][{"dataset_pfRelIso03_chg":key}].plot(ax=ax_pfRelIso03_chg,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_phi, ax_phi = plt.subplots(figsize=(10,10))
    ax_phi.set_ylabel("Counts")
    results["phi"][{"dataset_phi":key}].plot(ax=ax_phi,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_pt, ax_pt = plt.subplots(figsize=(10,10))
    ax_pt.set_ylabel("Counts")
    results["pt"][{"dataset_pt":key}].plot(ax=ax_pt,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_r9, ax_r9 = plt.subplots(figsize=(10,10))
    ax_r9.set_ylabel("Counts")
    results["r9"][{"dataset_r9":key}].plot(ax=ax_r9,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_sieie, ax_sieie = plt.subplots(figsize=(10,10))
    ax_sieie.set_ylabel("Counts")
    results["sieie"][{"dataset_sieie":key}].plot(ax=ax_sieie,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_mvaID_Fall17V1p1, ax_mvaID_Fall17V1p1 = plt.subplots(figsize=(10,10))
    ax_mvaID_Fall17V1p1.set_ylabel("Counts")
    results["mvaID_Fall17V1p1"][{"dataset_mvaID_Fall17V1p1":key}].plot(ax=ax_mvaID_Fall17V1p1,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_mvaID_WP80, ax_mvaID_WP80 = plt.subplots(figsize=(10,10))
    ax_mvaID_WP80.set_ylabel("Counts")
    results["mvaID_WP80"][{"dataset_mvaID_WP80":key}].plot(ax=ax_mvaID_WP80,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()
    
    f_mvaID_WP90, ax_mvaID_WP90 = plt.subplots(figsize=(10,10))
    ax_mvaID_WP90.set_ylabel("Counts")
    results["mvaID_WP90"][{"dataset_mvaID_WP90":key}].plot(ax=ax_mvaID_WP90,label=key)
    hep.cms.label("Preliminary",loc=0,com=13)
    #ax.set_yscale("log")
    plt.legend()


    #for var in variables:
    #    ax_name = "ax_"+var
    #    print(ax_name)
    if args.category == "EE":
        ax_mvaID.figure.savefig("%s/%s_EE_mvaID.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID)
        ax_eta.figure.savefig("%s/%s_EE_eta.png" %(directory,key), bbox_inches="tight")
        plt.close(f_eta)
        ax_hoe.figure.savefig("%s/%s_EE_hoe.png" %(directory,key), bbox_inches="tight")
        plt.close(f_hoe)
        ax_mass.figure.savefig("%s/%s_EE_mass.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mass)
        ax_pfRelIso03_all.figure.savefig("%s/%s_EE_pfRelIso03_all.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pfRelIso03_all)
        ax_pfRelIso03_chg.figure.savefig("%s/%s_EE_pfRelIso03_chg.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pfRelIso03_chg)
        ax_phi.figure.savefig("%s/%s_EE_phi.png" %(directory,key), bbox_inches="tight")
        plt.close(f_phi)
        ax_pt.figure.savefig("%s/%s_EE_pt.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pt)
        ax_r9.figure.savefig("%s/%s_EE_r9.png" %(directory,key), bbox_inches="tight")
        plt.close(f_r9)
        ax_sieie.figure.savefig("%s/%s_EE_sieie.png" %(directory,key), bbox_inches="tight")
        plt.close(f_sieie)
        ax_mvaID_Fall17V1p1.figure.savefig("%s/%s_EE_mvaID_Fall17V1p1.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_Fall17V1p1)
        ax_mvaID_WP80.figure.savefig("%s/%s_EE_mvaID_WP80.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_WP80)
        ax_mvaID_WP90.figure.savefig("%s/%s_EE_mvaID_WP90.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_WP80)

    if args.category == "EB":
        ax_mvaID.figure.savefig("%s/%s_EB_mvaID.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID)
        ax_eta.figure.savefig("%s/%s_EB_eta.png" %(directory,key), bbox_inches="tight")
        plt.close(f_eta)
        ax_hoe.figure.savefig("%s/%s_EB_hoe.png" %(directory,key), bbox_inches="tight")
        plt.close(f_hoe)
        ax_mass.figure.savefig("%s/%s_EB_mass.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mass)
        ax_pfRelIso03_all.figure.savefig("%s/%s_EB_pfRelIso03_all.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pfRelIso03_all)
        ax_pfRelIso03_chg.figure.savefig("%s/%s_EB_pfRelIso03_chg.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pfRelIso03_chg)
        ax_phi.figure.savefig("%s/%s_EB_phi.png" %(directory,key), bbox_inches="tight")
        plt.close(f_phi)
        ax_pt.figure.savefig("%s/%s_EB_pt.png" %(directory,key), bbox_inches="tight")
        plt.close(f_pt)
        ax_r9.figure.savefig("%s/%s_EB_r9.png" %(directory,key), bbox_inches="tight")
        plt.close(f_r9)
        ax_sieie.figure.savefig("%s/%s_EB_sieie.png" %(directory,key), bbox_inches="tight")
        plt.close(f_sieie)
        ax_mvaID_Fall17V1p1.figure.savefig("%s/%s_EB_mvaID_Fall17V1p1.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_Fall17V1p1)
        ax_mvaID_WP80.figure.savefig("%s/%s_EB_mvaID_WP80.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_WP80)
        ax_mvaID_WP90.figure.savefig("%s/%s_EB_mvaID_WP90.png" %(directory,key), bbox_inches="tight")
        plt.close(f_mvaID_WP80)
'''


################################

for key in sample_dict.keys():
    print("Sample: " + key)

    # mvaID
    mvaID = results["mvaID"][{"dataset_mvaID":key}].values()
    mvaID_bwc = binning_width_center_1d(nbins=nbins_mvaID, xlow=xlow_mvaID, xhigh=xhigh_mvaID)
    mvaID_stat = weighted_kinematics(var=mvaID, width=mvaID_bwc["width"])

    f_mvaID, ax_mvaID = plt.subplots(figsize=(xfig,yfig))

    # draw error bar for MC
    '''for i, x in enumerate(mvaID_bwc["center"]):
        if i == 0:
            ax_mvaID.fill_between([x - mvaID_bwc["width"] / 2., x + mvaID_bwc["width"] / 2.], mvaID_stat["ydn"][i], mvaID_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_mvaID.fill_between([x - mvaID_bwc["width"] / 2., x + mvaID_bwc["width"] / 2.],  mvaID_stat["ydn"][i], mvaID_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(mvaID_bwc["center"], bins=mvaID_bwc["binning"], weights=mvaID_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('mvaID')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_mvaID, xhigh_mvaID)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_mvaID.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_mvaID.png" %(directory,key))
    #plt.close(f_mvaID)
    plt.close()

    # eta
    eta = results["eta"][{"dataset_eta":key}].values()
    eta_bwc = binning_width_center_1d(nbins=nbins_eta, xlow=xlow_eta, xhigh=xhigh_eta)
    eta_stat = weighted_kinematics(var=eta, width=eta_bwc["width"])

    f_eta, ax_eta = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(eta_bwc["center"]):
        if i == 0:
            ax_eta.fill_between([x - eta_bwc["width"] / 2., x + eta_bwc["width"] / 2.], eta_stat["ydn"][i], eta_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_eta.fill_between([x - eta_bwc["width"] / 2., x + eta_bwc["width"] / 2.],  eta_stat["ydn"][i], eta_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(eta_bwc["center"], bins=eta_bwc["binning"], weights=eta_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$\eta$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_eta, xhigh_eta)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_eta.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_eta.png" %(directory,key))
    #plt.close(f_eta)
    plt.close()

    # H over E
    hoe = results["hoe"][{"dataset_hoe":key}].values()
    hoe_bwc = binning_width_center_1d(nbins=nbins_hoe, xlow=xlow_hoe, xhigh=xhigh_hoe)
    hoe_stat = weighted_kinematics(var=hoe, width=hoe_bwc["width"])

    f_hoe, ax_hoe = plt.subplots(figsize=(xfig,yfig))

    # draw error bar for MC
    '''for i, x in enumerate(hoe_bwc["center"]):
        if i == 0:
            ax_hoe.fill_between([x - hoe_bwc["width"] / 2., x + hoe_bwc["width"] / 2.], hoe_stat["ydn"][i], hoe_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_hoe.fill_between([x - hoe_bwc["width"] / 2., x + hoe_bwc["width"] / 2.],  hoe_stat["ydn"][i], hoe_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(hoe_bwc["center"], bins=hoe_bwc["binning"], weights=hoe_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('H/E')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_hoe, xhigh_hoe)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_hoe.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_hoe.png" %(directory,key))
    #plt.close(f_hoe)
    plt.close()

    # mass
    mass = results["mass"][{"dataset_mass":key}].values()
    mass_bwc = binning_width_center_1d(nbins=nbins_mass, xlow=xlow_mass, xhigh=xhigh_mass)
    mass_stat = weighted_kinematics(var=mass, width=mass_bwc["width"])

    f_mass, ax_mass = plt.subplots(figsize=(xfig,yfig))

    # draw error bar for MC
    '''for i, x in enumerate(mass_bwc["center"]):
        if i == 0:
            ax_mass.fill_between([x - mass_bwc["width"] / 2., x + mass_bwc["width"] / 2.], mass_stat["ydn"][i], mass_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_mass.fill_between([x - mass_bwc["width"] / 2., x + mass_bwc["width"] / 2.],  mass_stat["ydn"][i], mass_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(mass_bwc["center"], bins=mass_bwc["binning"], weights=mass_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$m_{\gamma}$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_mass, xhigh_mass)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_mass.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_mass.png" %(directory,key))
    #plt.close(f_mass)
    plt.close()

    # pfRelIso03_all
    pfRelIso03_all = results["pfRelIso03_all"][{"dataset_pfRelIso03_all":key}].values()
    pfRelIso03_all_bwc = binning_width_center_1d(nbins=nbins_pfRelIso03_all, xlow=xlow_pfRelIso03_all, xhigh=xhigh_pfRelIso03_all)
    pfRelIso03_all_stat = weighted_kinematics(var=pfRelIso03_all, width=pfRelIso03_all_bwc["width"])
    
    f_pfRelIso03_all, ax_pfRelIso03_all = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(pfRelIso03_all_bwc["center"]):
        if i == 0:
            ax_pfRelIso03_all.fill_between([x - pfRelIso03_all_bwc["width"] / 2., x + pfRelIso03_all_bwc["width"] / 2.], pfRelIso03_all_stat["ydn"][i], pfRelIso03_all_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_pfRelIso03_all.fill_between([x - pfRelIso03_all_bwc["width"] / 2., x + pfRelIso03_all_bwc["width"] / 2.],  pfRelIso03_all_stat["ydn"][i], pfRelIso03_all_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(pfRelIso03_all_bwc["center"], bins=pfRelIso03_all_bwc["binning"], weights=pfRelIso03_all_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('pfRelIso03_all')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_pfRelIso03_all, xhigh_pfRelIso03_all)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_pfRelIso03_all.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_pfRelIso03_all.png" %(directory,key))
    #plt.close(f_pfRelIso03_all)
    plt.close()

    # pfRelIso03_chg
    pfRelIso03_chg = results["pfRelIso03_chg"][{"dataset_pfRelIso03_chg":key}].values()
    pfRelIso03_chg_bwc = binning_width_center_1d(nbins=nbins_pfRelIso03_chg, xlow=xlow_pfRelIso03_chg, xhigh=xhigh_pfRelIso03_chg)
    pfRelIso03_chg_stat = weighted_kinematics(var=pfRelIso03_chg, width=pfRelIso03_chg_bwc["width"])

    f_pfRelIso03_chg, ax_pfRelIso03_chg = plt.subplots(figsize=(xfig,yfig))

    # draw error bar for MC
    '''for i, x in enumerate(pfRelIso03_chg_bwc["center"]):
        if i == 0:
            ax_pfRelIso03_chg.fill_between([x - pfRelIso03_chg_bwc["width"] / 2., x + pfRelIso03_chg_bwc["width"] / 2.], pfRelIso03_chg_stat["ydn"][i], pfRelIso03_chg_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_pfRelIso03_chg.fill_between([x - pfRelIso03_chg_bwc["width"] / 2., x + pfRelIso03_chg_bwc["width"] / 2.],  pfRelIso03_chg_stat["ydn"][i], pfRelIso03_chg_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(pfRelIso03_chg_bwc["center"], bins=pfRelIso03_chg_bwc["binning"], weights=pfRelIso03_chg_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('pfRelIso03_chg')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_pfRelIso03_chg, xhigh_pfRelIso03_chg)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_pfRelIso03_chg.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_pfRelIso03_chg.png" %(directory,key))
    #plt.close(f_pfRelIso03_chg)
    plt.close()

    # phi
    phi = results["phi"][{"dataset_phi":key}].values()
    phi_bwc = binning_width_center_1d(nbins=nbins_phi, xlow=xlow_phi, xhigh=xhigh_phi)
    phi_stat = weighted_kinematics(var=phi, width=phi_bwc["width"])
    
    f_phi, ax_phi = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(phi_bwc["center"]):
        if i == 0:
            ax_phi.fill_between([x - phi_bwc["width"] / 2., x + phi_bwc["width"] / 2.], phi_stat["ydn"][i], phi_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_phi.fill_between([x - phi_bwc["width"] / 2., x + phi_bwc["width"] / 2.],  phi_stat["ydn"][i], phi_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(phi_bwc["center"], bins=phi_bwc["binning"], weights=phi_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$\phi$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_phi, xhigh_phi)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_phi.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_phi.png" %(directory,key))
    #plt.close(f_phi)
    plt.close()

    # pt
    pt = results["pt"][{"dataset_pt":key}].values()
    pt_bwc = binning_width_center_1d(nbins=nbins_pt, xlow=xlow_pt, xhigh=xhigh_pt)
    pt_stat = weighted_kinematics(var=pt, width=pt_bwc["width"])

    f_pt, ax_pt = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(pt_bwc["center"]):
        if i == 0:
            ax_pt.fill_between([x - pt_bwc["width"] / 2., x + pt_bwc["width"] / 2.], pt_stat["ydn"][i], pt_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_pt.fill_between([x - pt_bwc["width"] / 2., x + pt_bwc["width"] / 2.],  pt_stat["ydn"][i], pt_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(pt_bwc["center"], bins=pt_bwc["binning"], weights=pt_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$p_{T}$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_pt, xhigh_pt)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_pt.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_pt.png" %(directory,key))
    #plt.close(f_pt)
    plt.close()

    # r9
    r9 = results["r9"][{"dataset_r9":key}].values()
    r9_bwc = binning_width_center_1d(nbins=nbins_r9, xlow=xlow_r9, xhigh=xhigh_r9)
    r9_stat = weighted_kinematics(var=r9, width=r9_bwc["width"])

    f_r9, ax_r9 = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(r9_bwc["center"]):
        if i == 0:
            ax_r9.fill_between([x - r9_bwc["width"] / 2., x + r9_bwc["width"] / 2.], r9_stat["ydn"][i], r9_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_r9.fill_between([x - r9_bwc["width"] / 2., x + r9_bwc["width"] / 2.],  r9_stat["ydn"][i], r9_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(r9_bwc["center"], bins=r9_bwc["binning"], weights=r9_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$R_{9}$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_r9, xhigh_r9)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_r9.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_r9.png" %(directory,key))
    #plt.close(f_r9)
    plt.close()

    # sigma ieta ieta
    sieie = results["sieie"][{"dataset_sieie":key}].values()
    sieie_bwc = binning_width_center_1d(nbins=nbins_sieie, xlow=xlow_sieie, xhigh=xhigh_sieie)
    sieie_stat = weighted_kinematics(var=sieie, width=sieie_bwc["width"])

    f_sieie, ax_sieie = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(sieie_bwc["center"]):
        if i == 0:
            ax_sieie.fill_between([x - sieie_bwc["width"] / 2., x + sieie_bwc["width"] / 2.], sieie_stat["ydn"][i], sieie_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_sieie.fill_between([x - sieie_bwc["width"] / 2., x + sieie_bwc["width"] / 2.],  sieie_stat["ydn"][i], sieie_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(sieie_bwc["center"], bins=sieie_bwc["binning"], weights=sieie_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('$\sigma_{i\eta i\eta}$')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_sieie, xhigh_sieie)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_sieie.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_sieie.png" %(directory,key))
    #plt.close(f_sieie)
    plt.close()

    # mvaID_Fall17V1p1
    mvaID_Fall17V1p1 = results["mvaID_Fall17V1p1"][{"dataset_mvaID_Fall17V1p1":key}].values()
    mvaID_Fall17V1p1_bwc = binning_width_center_1d(nbins=nbins_mvaID_Fall17V1p1, xlow=xlow_mvaID_Fall17V1p1, xhigh=xhigh_mvaID_Fall17V1p1)
    mvaID_Fall17V1p1_stat = weighted_kinematics(var=mvaID_Fall17V1p1, width=mvaID_Fall17V1p1_bwc["width"])

    f_mvaID_Fall17V1p1, ax_mvaID_Fall17V1p1 = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(mvaID_Fall17V1p1_bwc["center"]):
        if i == 0:
            ax_mvaID_Fall17V1p1.fill_between([x - mvaID_Fall17V1p1_bwc["width"] / 2., x + mvaID_Fall17V1p1_bwc["width"] / 2.], mvaID_Fall17V1p1_stat["ydn"][i], mvaID_Fall17V1p1_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_mvaID_Fall17V1p1.fill_between([x - mvaID_Fall17V1p1_bwc["width"] / 2., x + mvaID_Fall17V1p1_bwc["width"] / 2.],  mvaID_Fall17V1p1_stat["ydn"][i], mvaID_Fall17V1p1_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(mvaID_Fall17V1p1_bwc["center"], bins=mvaID_Fall17V1p1_bwc["binning"], weights=mvaID_Fall17V1p1_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('mvaID_Fall17V1p1')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_mvaID_Fall17V1p1, xhigh_mvaID_Fall17V1p1)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_mvaID_Fall17V1p1.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_mvaID_Fall17V1p1.png" %(directory,key))
    #plt.close(f_mvaID_Fall17V1p1)
    plt.close()

    # mvaID_WP80
    mvaID_WP80 = results["mvaID_WP80"][{"dataset_mvaID_WP80":key}].values()
    mvaID_WP80_bwc = binning_width_center_1d(nbins=nbins_mvaID_WP80, xlow=xlow_mvaID_WP80, xhigh=xhigh_mvaID_WP80)
    mvaID_WP80_stat = weighted_kinematics(var=mvaID_WP80, width=mvaID_WP80_bwc["width"])

    f_mvaID_WP80, ax_mvaID_WP80 = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(mvaID_WP80_bwc["center"]):
        if i == 0:
            ax_mvaID_WP80.fill_between([x - mvaID_WP80_bwc["width"] / 2., x + mvaID_WP80_bwc["width"] / 2.], mvaID_WP80_stat["ydn"][i], mvaID_WP80_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_mvaID_WP80.fill_between([x - mvaID_WP80_bwc["width"] / 2., x + mvaID_WP80_bwc["width"] / 2.],  mvaID_WP80_stat["ydn"][i], mvaID_WP80_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''

    #hep.cms.label(loc=2,com=13)
    plt.hist(mvaID_WP80_bwc["center"], bins=mvaID_WP80_bwc["binning"], weights=mvaID_WP80_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('mvaID_WP80')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_mvaID_WP80, xhigh_mvaID_WP80)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_mvaID_WP80.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_mvaID_WP80.png" %(directory,key))
    #plt.close(f_mvaID_WP80)
    plt.close()

    # mvaID_WP90
    mvaID_WP90 = results["mvaID_WP90"][{"dataset_mvaID_WP90":key}].values()
    mvaID_WP90_bwc = binning_width_center_1d(nbins=nbins_mvaID_WP90, xlow=xlow_mvaID_WP90, xhigh=xhigh_mvaID_WP90)
    mvaID_WP90_stat = weighted_kinematics(var=mvaID_WP90, width=mvaID_WP90_bwc["width"])

    f_mvaID_WP90, ax_mvaID_WP90 = plt.subplots(figsize=(xfig,yfig))
    
    # draw error bar for MC
    '''for i, x in enumerate(mvaID_WP90_bwc["center"]):
        if i == 0:
            ax_mvaID_WP90.fill_between([x - mvaID_WP90_bwc["width"] / 2., x + mvaID_WP90_bwc["width"] / 2.], mvaID_WP90_stat["ydn"][i], mvaID_WP90_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label=f"MC stat unc.") # we want just one entry in the legend
        else:
            ax_mvaID_WP90.fill_between([x - mvaID_WP90_bwc["width"] / 2., x + mvaID_WP90_bwc["width"] / 2.],  mvaID_WP90_stat["ydn"][i], mvaID_WP90_stat["yup"][i], facecolor="royalblue", alpha=0.5, edgecolor="royalblue", label="")
    '''
    
    #hep.cms.label(loc=2,com=13)
    plt.hist(mvaID_WP90_bwc["center"], bins=mvaID_WP90_bwc["binning"], weights=mvaID_WP90_stat["weighted"], histtype='step', label='MC', linewidth=2)
    plt.xlabel('mvaID_WP90')
    plt.ylabel('A.U.C.')
    plt.yscale('log')
    plt.legend()
    plt.xlim(xlow_mvaID_WP90, xhigh_mvaID_WP90)
    plt.title(key)#, fontsize = 16, loc = 'left')
    if args.category == "EE":
        plt.savefig("%s/%s_EE_mvaID_WP90.png" %(directory,key))
    if args.category == "EB":
        plt.savefig("%s/%s_EB_mvaID_WP90.png" %(directory,key))
    #plt.close(f_mvaID_WP90)
    plt.close()

