import numpy as np
import awkward as ak
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.nanoevents.methods import vector
import hist
from hist import Hist
import argparse
import os
import json
import pandas as pd
import datetime as dt
import math
#import vector



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
            samplefile = open('../configs/private_nanoaod/GGH_M130_M400_M800_MC_UL2017.json')
            #samplefile = open('../configs/GGH_M130_M400_M800_MC_UL2017.json')
            #samplefile = open('../configs/private_nanoaod/M130_GluGlu_MC_UL2017.json')
        if model == "RSG":
            samplefile = open("../configs/private_nanoaod/RSG_M130_M400_M800_MC_UL2017.json")
            #samplefile = open('../configs/RSG_M130_M400_M800_MC_UL2017.json')
            #samplefile = open('../configs/M130_M400_M800_MC_UL2017.json')
            #samplefile = open('../configs/private_nanoaod/M130_RSG_MC_UL2017.json')
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
current_time = dt.datetime.now()
strtime = current_time.strftime("%Y%m%d_%H%M%S")
#folder="photon_kinematics"
directory="../parquet/%s/%s/diphoton/%s" %(yearinput, model, strtime)
print("Directory: ", directory)
if os.path.exists(directory):
    print("Directory '%s' exists" %(directory))
else:
    os.makedirs(directory)
    print("Directory '%s' just created" %(directory))


bin_dict = {
    "diphotonmass": {"nbins":50, "xlow":0, "xhigh":1000},
    "diphotonpt": {"nbins":50, "xlow":0, "xhigh":1000},
    "lead_pt": {"nbins":50, "xlow":0, "xhigh":1000},
    "sublead_pt": {"nbins":50, "xlow":0, "xhigh":1000},
    "lead_eta": {"nbins":35, "xlow":-3.5, "xhigh":3.5},
    "sublead_eta": {"nbins":35, "xlow":-3.5, "xhigh":3.5},
    "deltaPhi": {"nbins":30, "xlow":-3., "xhigh":3.},
    "deltaEta": {"nbins":30, "xlow":-3., "xhigh":3.},
    "deltaR": {"nbins":60, "xlow":0, "xhigh":6.},
    "lead_SCeta": {"nbins":35, "xlow":-3.5, "xhigh":3.5},
    "sublead_SCeta": {"nbins":35, "xlow":-3.5, "xhigh":3.5},
    "cosThetaStar": {"nbins":50, "xlow":-1.0, "xhigh":1.0},
    "genDiphoton_cosThetaStar": {"nbins":50, "xlow":-1.0, "xhigh":1.0},
    "diphotonInvariantMass": {"nbins":50, "xlow":0, "xhigh":1000}
}

#Processor#    
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
            
        h_diphotonmass = (
            Hist.new.StrCat([], growth=True, name="dataset_diphotonmass", label="Primary dataset")
            .Reg(bin_dict["diphotonmass"]["nbins"], bin_dict["diphotonmass"]["xlow"], bin_dict["diphotonmass"]["xhigh"], overflow=False, underflow=False, name="x_diphotonmass", label = r"m$_{\gamma1 \gamma2}$ (GeV)")
            .Weight()
        )

        h_diphotonpt = (
            Hist.new.StrCat([], growth=True, name="dataset_diphotonpt", label="Primary dataset")
            .Reg(bin_dict["diphotonpt"]["nbins"], bin_dict["diphotonpt"]["xlow"], bin_dict["diphotonpt"]["xhigh"], overflow=False, underflow=False, name="x_diphotonpt", label = r"p$_{T,\gamma1 \gamma2}$ (GeV)")
            .Weight()
        )

        h_lead_pt = (
            Hist.new.StrCat([], growth=True, name="dataset_lead_pt", label="Primary dataset")
            .Reg(bin_dict["lead_pt"]["nbins"], bin_dict["lead_pt"]["xlow"], bin_dict["lead_pt"]["xhigh"], overflow=False, underflow=False, name="x_lead_pt", label = r"p_$T\gamma1$ (GeV)")
            .Weight()
        )

        h_sublead_pt = (
            Hist.new.StrCat([], growth=True, name="dataset_sublead_pt", label="Primary dataset")
            .Reg(bin_dict["sublead_pt"]["nbins"], bin_dict["sublead_pt"]["xlow"], bin_dict["sublead_pt"]["xhigh"], overflow=False, underflow=False, name="x_sublead_pt", label = r"p_$T\gamma2$ (GeV)")
            .Weight()
        )

        h_lead_eta = (
            Hist.new.StrCat([], growth=True, name="dataset_lead_eta", label="Primary dataset")
            .Reg(bin_dict["lead_eta"]["nbins"], bin_dict["lead_eta"]["xlow"], bin_dict["lead_eta"]["xhigh"], overflow=False, underflow=False, name="x_lead_eta", label = r"$\eta_{\gamma1}$")
            .Weight()
        )

        h_sublead_eta = (
            Hist.new.StrCat([], growth=True, name="dataset_sublead_eta", label="Primary dataset")
            .Reg(bin_dict["sublead_eta"]["nbins"], bin_dict["sublead_eta"]["xlow"], bin_dict["sublead_eta"]["xhigh"], overflow=False, underflow=False, name="x_sublead_eta", label = r"$\eta_{\gamma2}$")
            .Weight()
        )

        h_deltaPhi = (
            Hist.new.StrCat([], growth=True, name="dataset_deltaPhi", label="Primary dataset")
            .Reg(bin_dict["deltaPhi"]["nbins"], bin_dict["deltaPhi"]["xlow"], bin_dict["deltaPhi"]["xhigh"], overflow=False, underflow=False, name="x_deltaPhi", label = r"$\Delta\phi_{\gamma1 \gamma2}$")
            .Weight()
        )

        h_deltaEta = (
            Hist.new.StrCat([], growth=True, name="dataset_deltaEta", label="Primary dataset")
            .Reg(bin_dict["deltaEta"]["nbins"], bin_dict["deltaEta"]["xlow"], bin_dict["deltaEta"]["xhigh"], overflow=False, underflow=False, name="x_deltaEta", label = r"$\Delta\eta_{\gamma1 \gamma2}$")
            .Weight()
        )

        h_deltaR = (
            Hist.new.StrCat([], growth=True, name="dataset_deltaR", label="Primary dataset")
            .Reg(bin_dict["deltaR"]["nbins"], bin_dict["deltaR"]["xlow"], bin_dict["deltaR"]["xhigh"], overflow=False, underflow=False, name="x_deltaR", label = r"$\DeltaR_{\gamma1 \gamma2}$")
            .Weight()
        )

        h_lead_SCeta = (
            Hist.new.StrCat([], growth=True, name="dataset_lead_SCeta", label="Primary dataset")
            .Reg(bin_dict["lead_SCeta"]["nbins"], bin_dict["lead_SCeta"]["xlow"], bin_dict["lead_SCeta"]["xhigh"], overflow=False, underflow=False, name="x_lead_SCeta", label = r"$\eta_{SC \gamma1}$")
            .Weight()
        )

        h_sublead_SCeta = (
            Hist.new.StrCat([], growth=True, name="dataset_sublead_SCeta", label="Primary dataset")
            .Reg(bin_dict["sublead_SCeta"]["nbins"], bin_dict["sublead_SCeta"]["xlow"], bin_dict["sublead_SCeta"]["xhigh"], overflow=False, underflow=False, name="x_sublead_SCeta", label = r"$\eta_{SC \gamma2}$")
            .Weight()
        )

        h_cosThetaStar = (
            Hist.new.StrCat([], growth=True, name="dataset_cosThetaStar", label="Primary dataset")
            .Reg(bin_dict["cosThetaStar"]["nbins"], bin_dict["cosThetaStar"]["xlow"], bin_dict["cosThetaStar"]["xhigh"], overflow=False, underflow=False, name="x_cosThetaStar", label = r"cosThetaStar")
            .Weight()
        )

        h_genDiphoton_cosThetaStar = (
            Hist.new.StrCat([], growth=True, name="dataset_genDiphoton_cosThetaStar", label="Primary dataset")
            .Reg(bin_dict["genDiphoton_cosThetaStar"]["nbins"], bin_dict["genDiphoton_cosThetaStar"]["xlow"], bin_dict["genDiphoton_cosThetaStar"]["xhigh"], overflow=False, underflow=False, name="x_genDiphoton_cosThetaStar", label = r"cosThetaStar")
            .Weight()
        )

        h_diphotonInvariantMass = (
            Hist.new.StrCat([], growth=True, name="dataset_diphotonInvariantMass", label="Primary dataset")
            .Reg(bin_dict["diphotonInvariantMass"]["nbins"], bin_dict["diphotonInvariantMass"]["xlow"], bin_dict["diphotonInvariantMass"]["xhigh"], overflow=False, underflow=False, name="x_diphotonInvariantMass", label = r"Diphoton Invariant Mass (GeV)")
            .Weight()
        )


        #get HLT
        # if yearinput == "2016" or yearinput == "2016PreVFP":
        #     hlt_doublephoton60 = events.HLT.DoublePhoton60 #for 2016
        #     events_preselected = events[hlt_doublephoton60]
        # if yearinput == "2017" or yearinput == "2018":
        #     hlt_doublephoton70 = events.HLT.DoublePhoton70 #for 2017 or 2018
        #     events_preselected = events[hlt_doublephoton70]

        # gen level
        genPart = events.GenPart
        #genPhotons = genPart[(genPart.status == 1) & (genPart.pdgId == 22)]
        if ak.any(genPart.pdgId == 22):
            genPhoMomIdx = genPart.genPartIdxMother
        
        if model == "GGH":
            genPhotons = genPart[(genPart.pdgId == 22) & (genPart.pdgId[genPhoMomIdx] == 25)]
        if model == "RSG":
            genPhotons = genPart[(genPart.pdgId == 22) & (genPart.pdgId[genPhoMomIdx] == 5100039)]

        count_number_of_genPhotons = ak.num(genPhotons, axis=1)

        genPhotons_masked = genPhotons.mask[count_number_of_genPhotons > 1]
        # genPhoton_interselection = ak.fill_none(
        #     count_number_of_genPhotons > 1,
        #     False
        # )
        # genPhotons_interselected = genPhotons[genPhoton_interselection]

        # genPhotons_flattened = ak.flatten(genPhotons_interselected, axis=0)

        # genPhotons_masked = genPhotons_interselected.mask[
        # genPhotons_masked = genPhotons_flattened.mask[
        ## genPhotons_masked = genPhotons_masked.mask[
        #     (genPhotons.pt[:,0] > 50) &
        #     (genPhotons.pt[:,1] > 50)
        # ]

        # save only the events pass the selections
        total_genPhoton_selection = ak.fill_none(
            ak.num(genPhotons_masked,axis=1) > 1,
            False
        )
        genPhotons_selected = genPhotons_masked[total_genPhoton_selection]
        genPhotons_selected["charge"] = ak.zeros_like(genPhotons_selected.pt)
        genPhotons_selected_zip = ak.zip(
            {
                "pt": genPhotons_selected.pt,
                "eta": genPhotons_selected.eta,
                "phi": genPhotons_selected.phi,
                "mass": genPhotons_selected.mass,
                "charge": genPhotons_selected.charge,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )

        #print("charge: ", genPhotons_selected["charge"])
        genDiphoton_pairs = ak.combinations(genPhotons_selected, 2, fields=["lead", "sublead"])
        genDiphotons = genDiphoton_pairs.lead+genDiphoton_pairs.sublead

        #print(genDiphotons["charge"])

        genDiphotons_selection = ak.fill_none(
            (genDiphotons.mass[:,0] > 0),
            False
        )

        genDiphotons_selected = genDiphotons[genDiphotons_selection]
        genDiphotons_selected["charge"] = ak.zeros_like(genDiphotons_selected.pt)

        # genDiphoton cosThetaStar is cos(theta*) in the rest frame of the diphoton in gen-level
        genDiphotons_selected_zip = ak.zip(
            {
                "pt": genDiphotons_selected.pt,
                "eta": genDiphotons_selected.eta,
                "phi": genDiphotons_selected.phi,
                "mass": genDiphotons_selected.mass,
                "charge": genDiphotons_selected.charge,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )

        #print("genDiphotons_selected_zip: ", genDiphotons_selected_zip)
        genPhoton1_rest = genPhotons_selected_zip[:,0].boost(-genDiphotons_selected_zip.boostvec)
        #print("genPhoton1_rest: ", genPhoton1_rest)
        genDiphoton_ThetaStar = genPhoton1_rest.theta
        #print("genDiphoton_ThetaStar: ", genDiphoton_ThetaStar)
        genDiphoton_cosThetaStar = np.cos(genDiphoton_ThetaStar)
        #print("genDiphoton_cosThetaStar: ", genDiphoton_cosThetaStar)


        # get photons
        # photons = events_preselected.Photon
        photons = events.Photon

        # add selections

        # high pT photon ID v2
        # pass_pho_id_photons = photons[
        #     ((photons.isScEtaEB) & (photons.hoe < 0.05) & (photons.r9 > 0.8) & (photons.sieie < 0.0105) & (photons.electronVeto)) |
        #     ((photons.isScEtaEE) & (photons.hoe < 0.05) & (photons.r9 > 0.8) & (photons.sieie < 0.0280) & (photons.electronVeto))
        # ]

        pass_pho_id_photons = photons[photons.electronVeto]
        
        count_number_of_photons = ak.num(pass_pho_id_photons, axis=1)

        photons_masked = pass_pho_id_photons.mask[count_number_of_photons > 1]

        # photons_masked = photons_masked.mask[
        #     (photons_masked.pt[:,0] > 125) &
        #     (photons_masked.pt[:,1] > 125)
        # ]

        photons_masked = photons_masked.mask[
            (photons_masked.pt[:,0] > 50) &
            (photons_masked.pt[:,1] > 50)
        ]

        # get both leading and subleading photons in EBEB
        # if args.category == "EBEB":
        #     photons_masked = photons_masked.mask[
        #         (photons_masked.isScEtaEB[:,0]) &
        #         (photons_masked.isScEtaEB[:,1])
        #     ]

        # get leading and subleading photons either in EBEE or EEEB
        # if args.category == "EBEE":
        #     photons_masked = photons_masked.mask[
        #         ((photons_masked.isScEtaEB[:,0]) & (photons_masked.isScEtaEE[:,1])) |
        #         ((photons_masked.isScEtaEE[:,0]) & (photons_masked.isScEtaEB[:,1]))
        #     ]

        # save only the events pass the selections
        total_selection = ak.fill_none(
            ak.num(photons_masked,axis=1) > 1,#1
            False
        )

        photons_selected = photons_masked[total_selection]

        deltaEta = photons_selected.eta[:,0] - photons_selected.eta[:,1]
        deltaPhi = photons_selected.phi[:,0] - photons_selected.phi[:,1]
        # print("deltaphi: ", deltaPhi)
        deltaPhi = ak.where(deltaPhi > np.pi, deltaPhi - (2 * np.pi), deltaPhi)
        deltaPhi = ak.where(deltaPhi <= - np.pi, deltaPhi + (2 * np.pi), deltaPhi)
        deltaR = np.sqrt( (deltaEta)**2 + (deltaPhi)**2 )

        # add charge varibale to the photons: zero for all photons. It will be used in the diphoton 4-momentum calculation
        photons_selected["charge"] = ak.zeros_like(photons_selected.pt)

        photons_selected_zip = ak.zip(
            {
                "pt": photons_selected.pt,
                "eta": photons_selected.eta,
                "phi": photons_selected.phi,
                "mass": photons_selected.mass,
                "charge": photons_selected.charge,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )


        # make the diphoton pair combinations
        diphoton_pairs = ak.combinations(photons_selected, 2, fields=["lead", "sublead"])

        # diphoton four-momentum
        diphotons = diphoton_pairs.lead+diphoton_pairs.sublead

        # apply mass cut to diphoton
        # final_selection = ak.fill_none(
        #     (diphotons.mass[:,0] > 150),
        #     False
        # )
        final_selection = ak.fill_none(
            (diphotons.mass[:,0] > 0),
            False
        )

        diphotons_selected = diphotons[final_selection]

        # diphoton cosThetaStar is cos(theta*) in the rest frame of the diphoton
        diphotons_selected_zip = ak.zip(
            {
                "pt": diphotons_selected.pt,
                "eta": diphotons_selected.eta,
                "phi": diphotons_selected.phi,
                "mass": diphotons_selected.mass,
                "charge": diphotons_selected.charge,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )

        photon1_rest = photons_selected_zip[:,0].boost(-diphotons_selected_zip.boostvec)
        diphoton_ThetaStar = photon1_rest.theta
        #print("diphoton_ThetaStar: ", diphoton_ThetaStar)
        diphoton_cosThetaStar = np.cos(diphoton_ThetaStar)

        h_diphotonmass.fill(dataset_diphotonmass=dataset,x_diphotonmass=diphotons_selected.mass[:,0])
        #h.fill(dataset=dataset,x=photons.mass[:,0])
        results["diphotonmass"] = h_diphotonmass

        h_diphotonpt.fill(dataset_diphotonpt=dataset,x_diphotonpt=diphotons_selected.pt[:,0])
        results["diphotonpt"] = h_diphotonpt

        h_lead_pt.fill(dataset_lead_pt=dataset,x_lead_pt=photons_selected.pt[:,0])
        results["lead_pt"] = h_lead_pt

        h_sublead_pt.fill(dataset_sublead_pt=dataset,x_sublead_pt=photons_selected.pt[:,1])
        results["sublead_pt"] = h_sublead_pt

        h_lead_eta.fill(dataset_lead_eta=dataset,x_lead_eta=photons_selected.eta[:,0])
        results["lead_eta"] = h_lead_eta

        h_sublead_eta.fill(dataset_sublead_eta=dataset,x_sublead_eta=photons_selected.eta[:,1])
        results["sublead_eta"] = h_sublead_eta

        h_deltaEta.fill(dataset_deltaEta=dataset,x_deltaEta=deltaEta)
        results["deltaEta"] = h_deltaEta

        h_deltaPhi.fill(dataset_deltaPhi=dataset,x_deltaPhi=deltaPhi)
        results["deltaPhi"] = h_deltaPhi

        h_deltaR.fill(dataset_deltaR=dataset,x_deltaR=deltaR)
        results["deltaR"] = h_deltaR

        h_lead_SCeta.fill(dataset_lead_SCeta=dataset,x_lead_SCeta=photons_selected.superclusterEta[:,0])
        results["lead_SCeta"] = h_lead_SCeta

        h_sublead_SCeta.fill(dataset_sublead_SCeta=dataset,x_sublead_SCeta=photons_selected.superclusterEta[:,1])
        results["sublead_SCeta"] = h_sublead_SCeta

        h_cosThetaStar.fill(dataset_cosThetaStar=dataset,x_cosThetaStar=diphoton_cosThetaStar[:,0])
        results["cosThetaStar"] = h_cosThetaStar

        h_genDiphoton_cosThetaStar.fill(dataset_genDiphoton_cosThetaStar=dataset,x_genDiphoton_cosThetaStar=genDiphoton_cosThetaStar[:,0])
        results["genDiphoton_cosThetaStar"] = h_genDiphoton_cosThetaStar

        # h_diphotonInvariantMass.fill(dataset_diphotonInvariantMass=dataset,x_diphotonInvariantMass=diphotons_selected_zip.absolute[:,0])
        # results["diphotonInvariantMass"] = h_diphotonInvariantMass

        return results

    def postprocess(self, accumulant):
        pass

# cell 20

sample_dict = json.load(samplefile)

#sample_dict = {
#    "glugluSpin0toGG":[
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
    #executor=processor.IterativeExecutor(workers=4),
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

variables = ["diphotonmass","diphotonpt","lead_pt","sublead_pt","lead_eta","sublead_eta","deltaPhi","deltaEta","deltaR","lead_SCeta","sublead_SCeta","cosThetaStar", "genDiphoton_cosThetaStar"]#, "diphotonInvariantMass"]

for key in sample_dict.keys():
    print("Sample: " + key)
    key_ntuple={}

    for variable in variables:
        dataset_var = "dataset_%s" %(variable) #(f'dataset_{variable}')

        key_var_hist = results[variable][{dataset_var:key}]
        key_var_value = key_var_hist.values()
        key_var_edges = key_var_hist.axes.edges
        
        key_ntuple[variable]=[{
            'bin_center':[(key_var_edges[0][i] + key_var_edges[0][i+1]) / 2 for i in range(len(key_var_edges[0])-1)],
            'value':key_var_value,
            'nbins': [bin_dict[variable]["nbins"]],
            'xlow': [bin_dict[variable]["xlow"]],
            'xhigh': [bin_dict[variable]["xhigh"]]

        #key_ntuple[var_string]['bin_center']=
        }] #if include '[]' it will create bin_center, nbins, xlow, and xhigh as subfields. But, if not, those will being as index at 0 level

    df_key = pd.DataFrame(data=key_ntuple)
    df_key.to_parquet("%s/%s.parquet" %(directory, key))

# import mplhep as hep
# import matplotlib.pyplot as plt

# hep.style.use("CMS")

# for key in sample_dict.keys():
#     print("Sample: " + key)

#     f_diphotonmass, ax_diphotonmass = plt.subplots(figsize=(10,10))
#     ax_diphotonmass.set_ylabel("Counts", fontsize=14)
#     #ax_diphotonmass.set_xlabel(r'm$_{\gamma \gamma}$ [GeV]', fontsize=14)
#     #ax_diphotonmass.set_title(f'{key}', fontsize=14)
#     #ax_diphotonmass.tick_params(axis='x', labelsize=10)
#     #ax_diphotonmass.tick_params(axis='y', labelsize=10)
#     #ax_diphotonmass.legend(prop={'size': 14})
#     #results["mass"][{"dataset":"glugluSpin0toGG"}].plot(ax=ax,label="Spin0toGG")
#     #results["mass"][{"dataset":"DY"}].plot(ax=ax,label="DY")
#     results["mass"][{"dataset_diphotonmass":key}].plot(ax=ax_diphotonmass, label=key, density=True)
#     hep.cms.label("Preliminary",loc=0,com=13)
#     ax_diphotonmass.set_yscale("log")
#     plt.legend()

#     #print(results["mass"].view())

#     if args.category == "EBEB":
#         ax_diphotonmass.figure.savefig("%s/%s_EBEB_mass.png" %(directory,key), bbox_inches="tight")
#         plt.close(f_diphotonmass)
#     if args.category == "EBEE":
#         ax_diphotonmass.figure.savefig("%s/%s_EBEE_mass.png" %(directory,key), bbox_inches="tight")
#         plt.close(f_diphotonmass)
