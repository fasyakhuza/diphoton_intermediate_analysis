##############################################################################################
## Workflow prepared for the high-mass diphoton analysis (Run2,3)                           ##
## Processor: XggHighMassProcessor                                                          ##
## Based on : HggBaseProcessor                                                              ##
## Version1 : Samadhan Kamble, IIT Madras (samadhan.kamble@cern.ch)                         ##
## Run with : run_analysis.py --json-analysis runner.json --executor futures --dump NTuples ##
##############################################################################################


from higgs_dna.workflows.base import HggBaseProcessor

from higgs_dna.tools.chained_quantile import ChainedQuantileRegression
from higgs_dna.tools.SC_eta import add_photon_SC_eta
from higgs_dna.tools.EELeak_region import veto_EEleak_flag
from higgs_dna.tools.EcalBadCalibCrystal_events import remove_EcalBadCalibCrystal_events
from higgs_dna.tools.gen_helpers import get_fiducial_flag, get_genJets, get_higgs_gen_attributes
from higgs_dna.tools.xgb_loader import load_bdt
from higgs_dna.tools.photonid_mva import calculate_photonid_mva, load_photonid_mva
from higgs_dna.tools.photonid_mva import calculate_photonid_mva_run3, load_photonid_mva_run3
from higgs_dna.tools.flow_corrections import calculate_flow_corrections
from higgs_dna.tools.mass_decorrelator import decorrelate_mass_resolution

from higgs_dna.selections.photon_selections import photon_preselection
from higgs_dna.selections.lepton_selections import select_electrons, select_muons
from higgs_dna.selections.jet_selections import select_jets, jetvetomap
from higgs_dna.selections.lumi_selections import select_lumis
from higgs_dna.selections.sv_selections import match_sv

from higgs_dna.metaconditions import photon_id_mva_weights
from higgs_dna.metaconditions import diphoton as diphoton_mva_dir

from higgs_dna.utils.dumping_utils import (
    diphoton_ak_array,
    dump_ak_array,
    diphoton_list_to_pandas,
    dump_pandas,
    get_obj_syst_dict,
    dress_branches
)
from higgs_dna.utils.misc_utils import choose_jet

from higgs_dna.systematics import object_systematics as available_object_systematics
from higgs_dna.systematics import object_corrections as available_object_corrections
from higgs_dna.systematics import weight_systematics as available_weight_systematics
from higgs_dna.systematics import weight_corrections as available_weight_corrections

import functools
import operator
import os
import warnings
import awkward as ak
import numpy
import sys
import vector
from typing import Any, Dict, List, Optional
from coffea import processor
from coffea.analysis_tools import Weights
from copy import deepcopy
from coffea.nanoevents.methods import candidate

import logging
logger = logging.getLogger(__name__)

vector.register_awkward()

class XggHighMassProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        systematics: Dict[str, List[Any]] = None,
        corrections: Dict[str, List[Any]] = None,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
        trigger_group: str = ".*DoubleEG.*",
        analysis: str = "mainAnalysis",   # "highMassAnalysis"
        skipCQR: bool = False,
        skipJetVetoMap: bool = False,
        year: Dict[str, List[str]] = None,
        fiducialCuts: str = "classical",
        doDeco: bool = False,
        Smear_sigma_m: bool = False,
        doFlow_corrections: bool = False,
        output_format: str = "parquet"
    ) -> None:
        super().__init__(   # The super() function is used to call the constructor of
                            # the parent class HggBaseProcessor, passing all
                            # the parameters to it. This ensures that the XggHighMassProcessor
                            # is properly initialized with all the necessary settings inherited
                            # from the base processor, while also allowing for any
                            # additional customization specific to high-mass analyses.
            metaconditions,
            systematics=systematics,
            corrections=corrections,
            apply_trigger=apply_trigger,
            output_location=output_location,
            taggers=taggers,
            trigger_group=trigger_group,
            analysis=analysis,
            skipCQR=skipCQR,
            skipJetVetoMap=skipJetVetoMap,
            year=year,
            fiducialCuts=fiducialCuts,
            doDeco=doDeco,
            Smear_sigma_m=Smear_sigma_m,
            doFlow_corrections=doFlow_corrections,
            output_format=output_format
        )

        # ------- Re-define some cuts specific for Xgg High-Mass analysis -------

        # why the muon and electron selections are not included?
        #maybe I can adjust these values for GGH and RSG analysis, especially for diphoton preselection
        # jet selection cuts
        self.jet_jetId = "tightLepVeto"  # can be "tightLepVeto" or "tight": https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#nanoAOD_Flags
        self.jet_dipho_min_dr = 0.4
        self.jet_pho_min_dr = 0.4
        self.jet_ele_min_dr = 0.4
        self.jet_muo_min_dr = 0.4
        self.jet_pt_threshold = 20
        self.jet_max_eta = 4.7

        self.clean_jet_dipho = False
        self.clean_jet_pho = True
        self.clean_jet_ele = True
        self.clean_jet_muo = True

        # vbf
        self.vbf_lead_jet_min_pt = 40.0

        # from where you get these cut values? do you have any reference?
        # diphoton preselection cuts
        self.min_pt_photon = 25.0 #
        self.min_pt_lead_photon = 35.0
        self.min_pt_sublead_photon = 30.0
        self.min_mvaid = -0.9 #
        self.max_sc_eta = 2.5
        self.gap_barrel_eta = 1.4442
        self.gap_endcap_eta = 1.566
        self.max_hovere = 0.08 #
        self.min_full5x5_r9 = 0.8 #
        self.max_chad_iso = 20.0 #
        self.max_chad_rel_iso = 0.3 #

        self.min_full5x5_r9_EB_high_r9 = 0.85 #
        self.min_full5x5_r9_EE_high_r9 = 0.9 #
        self.min_full5x5_r9_EB_low_r9 = 0.5 #
        self.min_full5x5_r9_EE_low_r9 = 0.8 #
        self.max_trkSumPtHollowConeDR03_EB_low_r9 = 6.0 # # for v11, we cut on Photon_pfChargedIsoPFPV
        self.max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0 # # Leaving the names of the preselection cut variables
                                                         # the same to change as little as possible
        self.max_sieie_EB_low_r9 = 0.015 #
        self.max_sieie_EE_low_r9 = 0.035 #
        self.max_pho_iso_EB_low_r9 = 4.0 #
        self.max_pho_iso_EE_low_r9 = 4.0 #

        self.eta_rho_corr = 1.5 #
        self.low_eta_rho_corr = 0.16544 #
        self.high_eta_rho_corr = 0.13212 #

        # EA values for Run3 from Egamma
        self.EA1_EB1 = 0.102056
        self.EA2_EB1 = -0.000398112
        self.EA1_EB2 = 0.0820317
        self.EA2_EB2 = -0.000286224
        self.EA1_EE1 = 0.0564915
        self.EA2_EE1 = -0.000248591
        self.EA1_EE2 = 0.0428606
        self.EA2_EE2 = -0.000171541
        self.EA1_EE3 = 0.0395282
        self.EA2_EE3 = -0.000121398
        self.EA1_EE4 = 0.0369761
        self.EA2_EE4 = -8.10369e-05
        self.EA1_EE5 = 0.0369417
        self.EA2_EE5 = -2.76885e-05
        self.e_veto = 0.5

        self.num_leptons_to_store = 4 #what is this for?

        logger.debug(f"Setting up processor with metaconditions: {self.meta}")
        logger.debug(f"self.year = {self.year}")
        logger.debug(f"self.apply_trigger = {self.apply_trigger}")
        logger.debug(f"self.trigger_group = {self.trigger_group}")
        logger.debug(f"self.corrections = {self.corrections}")
        logger.debug(f"self.systematics = {self.systematics}")
        logger.debug(f"self.taggers = {self.taggers}")
        logger.debug(f"self.output_location = {self.output_location}")

        # why do you comment this?
        # # ------- build the chained quantile regressions -------
        # if not self.skipCQR:
        #     logger.info(f"\n \tBuilding CQR as required\n")
        #     try:
        #         self.chained_quantile: Optional[ChainedQuantileRegression] = ChainedQuantileRegression(**self.meta["PhoIdInputCorrections"])
        #     except Exception as e:
        #         warnings.warn(f"Could not instantiate ChainedQuantileRegression: {e}")
        #         self.chained_quantile = None
        # else:
        #     logger.info("\n \tSkipping CQR as required\n")
        #     self.chained_quantile = None

        # # ------- initialize photonid_mva -------
        # photon_id_mva_dir = os.path.dirname(photon_id_mva_weights.__file__)
        # try:
        #     logger.info(f"\n \tInitializing photonid_mva\n")
        #     logger.debug(
        #         f"\t \tLooking for {self.meta['flashggPhotons']['photonIdMVAweightfile_EB']} in {photon_id_mva_dir}"
        #     )
        #     self.photonid_mva_EB = load_photonid_mva(
        #         os.path.join(
        #             photon_id_mva_dir,
        #             self.meta["flashggPhotons"]["photonIdMVAweightfile_EB"],
        #         )
        #     )
        #     self.photonid_mva_EE = load_photonid_mva(
        #         os.path.join(
        #             photon_id_mva_dir,
        #             self.meta["flashggPhotons"]["photonIdMVAweightfile_EE"],
        #         )
        #     )
        # except Exception as e:
        #     warnings.warn(f"\t \tCould not instantiate PhotonID MVA on the fly: {e}")
        #     self.photonid_mva_EB = None
        #     self.photonid_mva_EE = None

        # # ------- initialize diphoton mva -------
        # diphoton_weights_dir = os.path.dirname(diphoton_mva_dir.__file__)
        # logger.info(f"\n \tInitializing diphoton mva\n")
        # logger.debug(f"\t \tBase path to look for diphoton IDMVA weight files: {diphoton_weights_dir}")

        # try:
        #     self.diphoton_mva = load_bdt(
        #         os.path.join(
        #             diphoton_weights_dir, self.meta["customDiPhotonMVA"]["weightFile"]
        #         )
        #     )
        # except Exception as e:
        #     warnings.warn(f"\t \tCould not instantiate diphoton MVA: {e}")
        #     self.diphoton_mva = None

        logger.info(
            f"\n"
            f"\tself.photonid_mva_EB: {self.photonid_mva_EB}\n"
            f"\tself.photonid_mva_EE: {self.photonid_mva_EE}\n"
            f"\tself.diphoton_mva: {self.diphoton_mva}\n"
        )

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events, {} #what does this mean?

    def process(self, events: ak.Array) -> Dict[Any, Any]:

        logger.info("\n \t \t#------- Running the Xgg High-Mass Processor -------#\n")

        dataset_name = events.metadata["dataset"]
        logger.info(f"\n \tdataset_name: {dataset_name} \n")

        filename = events.metadata["filename"]
        logger.info(f"\n \tfilename: {filename} \n")

        # data or monte carlo?
        self.data_kind = "mc" if hasattr(events, "GenPart") else "data"

        # here we start recording possible coffea accumulators
        # most likely histograms, could be counters, arrays, ...
        histos_etc = {}
        histos_etc[dataset_name] = {}

        if self.data_kind == "mc":
            histos_etc[dataset_name]["nTot"] = int(ak.num(events.genWeight, axis=0))
            histos_etc[dataset_name]["nPos"] = int(ak.sum(events.genWeight > 0))
            histos_etc[dataset_name]["nNeg"] = int(ak.sum(events.genWeight < 0))
            histos_etc[dataset_name]["nEff"] = int(
                histos_etc[dataset_name]["nPos"] - histos_etc[dataset_name]["nNeg"]
            )
            histos_etc[dataset_name]["genWeightSum"] = float(ak.sum(events.genWeight))
        else:
            histos_etc[dataset_name]["nTot"] = int(len(events))
            histos_etc[dataset_name]["nPos"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["nNeg"] = int(0)
            histos_etc[dataset_name]["nEff"] = int(histos_etc[dataset_name]["nTot"])
            histos_etc[dataset_name]["genWeightSum"] = float(len(events))

        # ------ lumi mask -------
        if self.data_kind == "data":
            try:
                lumimask = select_lumis(self.year[dataset_name][0], events, logger)
                logger.info("\n \tApplying lumi-mask\n")
                events = events[lumimask]
            except:
                logger.info(f"[ lumimask ] Skip now! Unable to find year info of {dataset_name}")


        # ------ jetvetomap: only retain events that has no jets in the EE leakage region --------
        if not self.skipJetVetoMap:
            logger.info("\n \tApplying jetveto-map\n")
            events = jetvetomap(events, logger, dataset_name, year=self.year[dataset_name][0])

        # metadata array to append to HiggsDNA output
        metadata = {}

        if self.data_kind == "mc":
            logger.debug(f'genWeight: {events["genWeight"]}')
            # Add sum of gen weights before selection for normalisation in postprocessing
            metadata["sum_genw_presel"] = str(ak.sum(events.genWeight))
        else:
            metadata["sum_genw_presel"] = "Data"

        # ------ apply filters and triggers -------
        #maybe I can adjust this trigger for GGH and RSG analysis
        logger.debug(f"\n \tBefore filters and triggers: Number of events (for {dataset_name}): {len(events.Photon)}")
        logger.info("\n \tApplying filters and triggers\n")
        events = self.apply_filters_and_triggers(events)
        logger.debug(f"\n \tAfters filters and triggers: Number of events (for {dataset_name}): {len(events.Photon)}")

        # ------ remove events affected by EcalBadCalibCrystal -------
        #maybe I can adjust this crystal callibration for data for GGH and RSG analysis
        if self.data_kind == "data":
            logger.info("\n \tRemoving EcalBadCalibCrystal\n")
            events = remove_EcalBadCalibCrystal_events(events)

        # ------ add veto EE leak branch for photons, could also be used for electrons -------
        if (self.year[dataset_name][0] == "2022EE" or self.year[dataset_name][0] == "2022postEE"):
            logger.info("\n \tAdding veto EE leak branch for photons\n")
            events.Photon = veto_EEleak_flag(self, events.Photon)

        # ------ add SC_eta (we need it for corrections and systematics, it is present in NanoAODv13+
        # and can be calculated using PV for older versions) -------
        logger.info("\n \tAdding photon SC_eta\n")
        events.Photon = add_photon_SC_eta(events.Photon, events.PV)

        # read which systematics and corrections to process
        try:
            correction_names = self.corrections[dataset_name]
        except KeyError:
            correction_names = []
        try:
            systematic_names = self.systematics[dataset_name]
        except KeyError:
            systematic_names = []

        # ------ Compatibility check for applying 'Smearing' correction and Performing the mass resolution smearing (Smear_sigma_m_m) ------
        # for MC: expliciltly needs 'Smearing' correction.
        # If --Smear_sigma_m == True and no 'Smearing' correction in .json -->
        # throw an error, since the pt spectrum needs to be smeared in order to properly calculate the smeared sigma_m_m)
        if (
            self.data_kind == "mc"
            and self.Smear_sigma_m
            and "Smearing" not in correction_names
        ):
            warnings.warn(
                "Smearing should be specified in the corrections field in .json in order to smear the mass!"
            )
            sys.exit(0)

        # ------ for Data: apply 'Smearing' correction if --Smear_sigma_m == True ------
        # since now we are applying Smearing term to the sigma_m_over_m, I added this portion of code
        # specially for the estimation of smearing terms for the data events [data pt/energy] are not smeared!
        #maybe I can add "smearing" in the correction inside runner.json
        if self.data_kind == "data" and self.Smear_sigma_m:
            correction_name = "Smearing"   # doubt: wouldn't 'Smearing' correction for data repeat in the corrections loop ahead ?

            logger.info(f"\n \tApplying correction: {correction_name}, to dataset: {dataset_name}\n")
            varying_function = available_object_corrections[correction_name] ## for getting run Smearing function from photon_systematics.py 
            events = varying_function(events=events, year=self.year[dataset_name][0]) ## for run Smearing function from photon_systematics.py 

        ## ------ Applying corrections -------

        logger.info("\n \tApplying corrections...")
        logger.info(f"\n[Correction_names]: {correction_names}\n")

        logger.debug(f"\nPhoton pt before corrections: {events.Photon.pt}\n")
        for correction_name in correction_names:
            if correction_name in available_object_corrections.keys():
                logger.info(
                    f"\n \tApplying object correction: {correction_name}, to dataset: {dataset_name} \n"
                )
                varying_function = available_object_corrections[correction_name]
                events = varying_function(
                    events=events, year=self.year[dataset_name][0]
                )
            elif correction_name in available_weight_corrections:
                # event weight corrections will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(f"Could not process the corrections: {correction_name}, continuing..")
                continue

        original_photons = events.Photon   # these are in-fact after the corrections, not sure why is it called 'original',
                                           # perhaps as these are before systematics
        logger.debug(f"Photon pt after correction: {original_photons.pt}")
        # NOTE: jet jerc systematics are added in the correction functions and handled later
        original_jets = events.Jet

        # ------- computing the normalizing flow correction -------
        if self.data_kind == "mc" and self.doFlow_corrections:
            logger.info("\n \tApplying normalizing flow corrections\n")

            # applying the flow corrections to all photons before pre-selection
            counts = ak.num(original_photons)
            corrected_inputs,var_list = calculate_flow_corrections(
                original_photons, events, self.meta["flashggPhotons"]["flow_inputs"],
                self.meta["flashggPhotons"]["Isolation_transform_order"],
                year=self.year[dataset_name][0])

            # store the raw nanoAOD value and update the photon ID MVA value for preselection
            original_photons["mvaID_nano"] = original_photons["mvaID"]

            # store the raw values of the inputs and update the input values with the corrections
            # since some variables will be used in the preselection
            for i in range(len(var_list)):
                original_photons["raw_" + str(var_list[i])] = original_photons[str(var_list[i])] ## the raw values of the inputs
                original_photons[str(var_list[i])] = ak.unflatten(corrected_inputs[:,i] , counts) ## the input values with the corrections

            original_photons["mvaID"] = ak.unflatten(self.add_photonid_mva_run3(original_photons, events), counts) ## only for run3 ? should be calculate_photonid_ or load_photonid_ ?

        ## ------- Systematic object variations -------

        # ------- 'adding' systematic variations to original_photons
        # (original_photons-->systematics-->fields(scale, smearing, etc.)-->variations (up, down, etc.))-------
        logger.info("\nAdding systematic variations...")
        logger.info(f"\n[systematic_names]: {systematic_names}\n") ## systematic_names will contain dictionary/list of several systematics
        for systematic_name in systematic_names:
            if systematic_name in available_object_systematics.keys():
                systematic_dct = available_object_systematics[systematic_name]
                if systematic_dct["object"] == "Photon":
                    logger.info(
                        f"\n \tAdding systematic: {systematic_name}, to photons collection of dataset: {dataset_name}"
                    )
                    original_photons.add_systematic(
                        # passing the arguments here explicitly since I want to pass the events to the varying function.
                        # If there is a more elegant / flexible way, just change it!
                        name=systematic_name,
                        kind=systematic_dct["args"]["kind"],
                        what=systematic_dct["args"]["what"],
                        varying_function=functools.partial(
                            systematic_dct["args"]["varying_function"],
                            events=events,
                            year=self.year[dataset_name][0],
                        )
                        # name=systematic_name, **systematic_dct["args"]
                    ) ## add systematics from smearing, scale, energyerrshift, etc of poton to original_photons
                # to be implemented for other objects here
            elif systematic_name in available_weight_systematics:
                # event weight systematics will be applied after photon preselection / application of further taggers
                continue
            else:
                # may want to throw an error instead, needs to be discussed
                warnings.warn(
                    f"Could not process systematic variation: {systematic_name}."
                )
                continue
        logger.debug(f"After addition of systematic variations, original_photons.systematics: {original_photons.systematics}")

        # -------  'applying' systematic variations: creating 'systematic variations' deepcopies (scale_up, smearing_down, etc.)
        # in photons_dct from original_photons.systematics[systematic][variation] -------
        photons_dct = {}
        photons_dct["nominal"] = original_photons

        logger.debug(f"Systematics: original_photon.systematics.fields: {original_photons.systematics.fields}")
        for systematic in original_photons.systematics.fields:

            logger.debug(f"Variations: original_photon.systematic[systematic].fields: {original_photons.systematics[systematic].fields}")
            for variation in original_photons.systematics[systematic].fields:

                logger.debug("Variation fields: original_photon.systematic[systematic][variation].fields:"
                            f"{original_photons.systematics[systematic][variation].fields}")
                photons_dct[f"{systematic}_{variation}"] = deepcopy(original_photons.systematics[systematic][variation])   # deepcopy to allow for independent calculations on photon variables with CQR

        logger.debug(f"\n \t[ photons_dct ]: {photons_dct}\n")

        # NOTE: jet jerc systematics are added in the corrections, now extract those variations and create the dictionary
        logger.debug(f"\n \t[ original_jets.fields ]: {original_jets.fields}\n")

        jerc_syst_list, jets_dct = get_obj_syst_dict(original_jets, ["pt", "mass"])

        # object systematics dictionary
        logger.debug(f"\n \t[ jerc systematics ]: {jerc_syst_list}\n")
        logger.debug(f"\n \t[ jets_dct ]: {jets_dct}\n")

        ## ------- Build the flattened array of all possible variations -------

        variations_combined = []

        variations_combined.append(original_photons.systematics.fields)
        logger.debug(f"Photons variations_combined: {variations_combined}")

        # NOTE: jet jerc systematics are not added with add_systematics
        variations_combined.append(jerc_syst_list)
        logger.debug(f"Photon + Jet variations_combined: {variations_combined}")

        # ------- flatten -------
        variations_flattened = sum(variations_combined, [])  # begin with empty list and keep concatenating

        # ------- attach _down and _up -------
        variations = [item + suffix for item in variations_flattened for suffix in ['_down', '_up']]

        # ------- add nominal to the list -------
        variations.append('nominal')
        logger.info(f"\nTotal Systematics Variations: {variations}\n")

        ## ------- Processing block for each variation -------

        logger.info("\n \t#------- Processing block for each variation starts here -------#\n")

        for variation in variations:
            logger.info(f"\n \t#------- Variation: {variation} -------#\n")

            photons, jets = photons_dct["nominal"], events.Jet  # doubt: why don't we also take jets from jet_dct ?

            logger.debug(f"[photons fields at the beginning of the variations loop]: {photons.fields}")
            logger.debug(f"[Jets fields at the beginning of the variations loop]: {jets.fields}")

            if variation == "nominal":
                pass  # Do nothing since we already get the unvaried, but nominally corrected objets above
            elif variation in [*photons_dct]:   # [*dict] gets the keys of the dict since Python >= 3.5
                photons = photons_dct[variation]   # assigns to 'photons' the photons collection that corresponds to the specific variation
                logger.debug(f"[photons.fields] after {variation} variation assignment: {photons.fields}")
            elif variation in [*jets_dct]:
                jets = jets_dct[variation]   # assigns to jets the jet collection that corresponds to the specific variation
                logger.debug(f"jets fields after {variation} variation assignment: {jets.fields}")

            do_variation = variation  # We can also simplify this a bit but for now it works

            # ------- apply CQR -------
            if self.chained_quantile is not None:
                logger.info(f"\n \tApplying CQR to photons in {variation}\n")
                photons = self.chained_quantile.apply(photons, events)
            # ------- recompute photonid_mva on the fly and add to photons -------
            if self.photonid_mva_EB and self.photonid_mva_EE:
                logger.info(f"\n \tAdding recomputed photonid_mva on the fly to photons in {variation}\n")
                photons = self.add_photonid_mva(photons, events) ## photonid_mva not yet applied

            # ------- photon preselection -------
            # logger.debug(f"Before preselction: Number of events (for {dataset_name}, {variation}): {len(events)}")
            # logger.debug(f"Before preselection: Number of photons (for {dataset_name}, {variation}): {variation}: {ak.num(photons)}")
            # logger.debug(f"Before preselection: Total number of photons across events (for {dataset_name}, {variation}): {variation}: {ak.sum(ak.num(photons))}")
            logger.info(f"Before preselction: Number of events (for {dataset_name}, {variation}): {len(events)}")
            logger.info(f"Before preselection: Number of photons (for {dataset_name}, {variation}): {variation}: {ak.num(photons)}")
            logger.info(f"Before preselection: Total number of photons across events (for {dataset_name}, {variation}): {variation}: {ak.sum(ak.num(photons))}")

            logger.info(f"\n \tApplying photon preselection to photons in {variation}\n")

            photons = photon_preselection(self, photons, events, year=self.year[dataset_name][0])

            # logger.debug(f"After preselction: Number of events (for {dataset_name} & {variation}): {len(events)}")   # number of events retained even if no photon in a event passes preselection
            # logger.debug(f"After preselection: Number of photons per event (for {dataset_name}, {variation}): {ak.num(photons)}")
            # logger.debug(f"After preselection: Total number of photons acroos events (for {dataset_name}, {variation}): {ak.sum(ak.num(photons))}")
            logger.info(f"After preselction: Number of events (for {dataset_name} & {variation}): {len(events)}")   # number of events retained even if no photon in a event passes preselection
            logger.info(f"After preselection: Number of photons per event (for {dataset_name}, {variation}): {ak.num(photons)}")
            logger.info(f"After preselection: Total number of photons acroos events (for {dataset_name}, {variation}): {ak.sum(ak.num(photons))}")

            ## ------- Making diphoton candidates -------

            # sort photons in each event descending in pt
            photons = photons[ak.argsort(photons.pt, ascending=False)]
            photons["charge"] = ak.zeros_like(photons.pt)  # added this because 'charge'
                                                           # is not a property of photons in nanoAOD v11.
                                                           # We just assume every photon has charge zero...
            
            photons_zip = ak.zip(
            {
                "pt": photons.pt,
                "eta": photons.eta,
                "phi": photons.phi,
                "mass": photons.mass,
                "charge": photons.charge,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
            )

            logger.info(f"\n \tMaking diphoton candidates using photons in {variation}\n")
            diphotons = ak.combinations(
                photons, 2, fields=["pho_lead", "pho_sublead"]
            )

            # apply diphoton specific cuts
            logger.debug(f"Number of diphotons per event before diphoton specific cuts: {ak.num(diphotons)}")
            logger.debug(f"Total number of diphotons across events before diphoton specific cuts: {ak.sum(ak.num(diphotons))}")

            diphotons = diphotons[
                (diphotons["pho_lead"].pt > self.min_pt_lead_photon)
                & (diphotons["pho_sublead"].pt > self.min_pt_sublead_photon)
                # & ((diphotons.pho_lead.pt / diphotons.mass) >= 0.4)
                # & ((diphotons.pho_sublead.pt / diphotons.mass) >= 0.)
                # & (diphotons["mass"] >= self.min_mass_range)
                # & (diphotons["mass"] <= self.max_mass_range)
            ]
            logger.info(f"Number of diphotons per event after diphoton specific cuts: {ak.num(diphotons)}")
            logger.info(f"Total number of diphotons across events after diphoton specific cuts: {ak.sum(ak.num(diphotons))}")

            # Adding leading and subleading photons' pt in each event, and then, sort these values from big to small.
            # This will be used the moment we want to pick out the diphoton pair having the largest pt sum in each event.
            lead_and_sublead_pt_sum = diphotons.pho_lead.pt + diphotons.pho_sublead.pt
            diphotons = diphotons[ak.argsort(lead_and_sublead_pt_sum, axis = 1, ascending = False)]

            # now turn the diphotons into candidates with four momenta and such
            diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
            diphotons["pt"] = diphoton_4mom.pt
            diphotons["eta"] = diphoton_4mom.eta
            diphotons["phi"] = diphoton_4mom.phi
            diphotons["mass"] = diphoton_4mom.mass
            diphotons["charge"] = diphoton_4mom.charge

            diphoton_pz = diphoton_4mom.z
            diphoton_e = diphoton_4mom.energy
            diphotons["rapidity"] = 0.5 * numpy.log((diphoton_e + diphoton_pz) / (diphoton_e - diphoton_pz))

            diphotons = ak.with_name(diphotons, "PtEtaPhiMCandidate")

            ## Additional for run-2 studies
            # diphoton_4mom_zip = ak.zip(
            # {
            #     "pt": diphoton_4mom.pt,
            #     "eta": diphoton_4mom.eta,
            #     "phi": diphoton_4mom.phi,
            #     "mass": diphoton_4mom.mass,
            #     "charge": diphoton_4mom.charge,
            # },
            # with_name="PtEtaPhiMLorentzVector",
            # behavior=vector.behavior,
            # )
            
            # Leading photon is transformed using Lorentz transformation to the C.O.M. of the resonance
            # Then, the theta and cos(theta) between the leading photon in the resonance's C.O.M and the z-direction is calculated
            photon1_rest = photons_zip[:,0].boost(-diphoton_4mom.boostvec)
            diphoton_ThetaStar = photon1_rest.theta
            #print("diphoton_ThetaStar: ", diphoton_ThetaStar)
            diphotons["cosThetaStar"] = np.cos(diphoton_ThetaStar)

            # sort diphotons by pT
            diphotons = diphotons[
                ak.argsort(diphotons.pt, ascending=False)
            ] #ask Dmitry if we still have to aplly this diphoton sort again after the storing?

            logger.info(f"\n[diphotons.fields] at the beginning, for {variation}: {diphotons.fields}")

            # ------- determine if event passes fiducial Hgg cuts at the detector-level -------
            logger.info(f"\n \tApplying fiducial Hgg cuts and adding flags to diphotons in {variation}\n")
            if self.fiducialCuts == 'classical':
                # fid_det_passed = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & \
                #                  (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & \
                #                  (diphotons.pho_lead.pfRelIso03_all_quadratic * diphotons.pho_lead.pt < 10) & \
                #                  (diphotons.pho_sublead.pfRelIso03_all_quadratic * diphotons.pho_sublead.pt < 10) & \
                #                  (numpy.abs(diphotons.pho_lead.eta) < 2.5) & \
                #                  (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
                fid_det_passed = (diphotons.pho_lead.pt / diphotons.mass > 1 / 3) & \
                                 (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & \
                                 (diphotons.pho_lead.pfRelIso03_all * diphotons.pho_lead.pt < 10) & \
                                 (diphotons.pho_sublead.pfRelIso03_all * diphotons.pho_sublead.pt < 10) & \
                                 (numpy.abs(diphotons.pho_lead.eta) < 2.5) & \
                                 (numpy.abs(diphotons.pho_sublead.eta) < 2.5) ## the error coming from pho_lead.pfRelIso03_all_quadratic, there's no pho_lead.pfRelIso03_all_quadratic in events.photons.fields
            elif self.fiducialCuts == 'geometric':
                # fid_det_passed = (numpy.sqrt(diphotons.pho_lead.pt * diphotons.pho_sublead.pt) / diphotons.mass > 1 / 3) & \
                #                  (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & \
                #                  (diphotons.pho_lead.pfRelIso03_all_quadratic * diphotons.pho_lead.pt < 10) & \ 
                #                  (diphotons.pho_sublead.pfRelIso03_all_quadratic * diphotons.pho_sublead.pt < 10) & \
                #                  (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5) ## the error coming from pho_lead.pfRelIso03_all_quadratic, there's no pho_lead.pfRelIso03_all_quadratic in events.photons.fields
                fid_det_passed = (numpy.sqrt(diphotons.pho_lead.pt * diphotons.pho_sublead.pt) / diphotons.mass > 1 / 3) & \
                                 (diphotons.pho_sublead.pt / diphotons.mass > 1 / 4) & \
                                 (diphotons.pho_lead.pfRelIso03_all * diphotons.pho_lead.pt < 10) & \
                                 (diphotons.pho_sublead.pfRelIso03_all * diphotons.pho_sublead.pt < 10) & \
                                 (numpy.abs(diphotons.pho_lead.eta) < 2.5) & (numpy.abs(diphotons.pho_sublead.eta) < 2.5)
            elif self.fiducialCuts == 'none':
                fid_det_passed = diphotons.pho_lead.pt > -10  # This is a very dummy way but I do not know how to make a true array of outer shape of diphotons
            else:
                warnings.warn(f"You chose {self.fiducialCuts} the fiducialCuts mode, but this is currently not supported."
                              "You should check your settings. For this run, no fiducial selection at"
                              "detector level is applied.")
                fid_det_passed = diphotons.pho_lead.pt > -10

            diphotons = diphotons[fid_det_passed] #another selection for diphoton
            logger.debug(f"Total number of diphotons across events after fiducial cuts: {ak.sum(ak.num(diphotons))}")   # could use '--fiducialCuts none'

            if self.data_kind == "mc":
                # ------- add the fiducial flags for particle level -------
                diphotons['fiducialClassicalFlag'] = get_fiducial_flag(events, flavour='Classical')
                diphotons['fiducialGeometricFlag'] = get_fiducial_flag(events, flavour='Geometric')

                diphotons['GenPTH'], GenYH, GenPhiH = get_higgs_gen_attributes(events) ## check get_higgs_gen_attributes function in tools/gen_helpers.py

                genJets = get_genJets(self, events, pt_cut=30., eta_cut=2.5)
                diphotons['GenNJ'] = ak.num(genJets)
                GenPTJ0 = choose_jet(genJets.pt, 0, -999.0)  # Choose zero (leading) jet and pad with -999 if none
                diphotons['GenPTJ0'] = GenPTJ0

                gen_first_jet_eta = choose_jet(genJets.eta, 0, -999.0)
                gen_first_jet_mass = choose_jet(genJets.mass, 0, -999.0)
                gen_first_jet_phi = choose_jet(genJets.phi, 0, -999.0)

                gen_first_jet_pz = GenPTJ0 * numpy.sinh(gen_first_jet_eta)
                gen_first_jet_energy = numpy.sqrt((GenPTJ0**2 * numpy.cosh(gen_first_jet_eta)**2) + gen_first_jet_mass**2)

                with numpy.errstate(divide='ignore', invalid='ignore'):
                    GenYJ0 = 0.5 * numpy.log((gen_first_jet_energy + gen_first_jet_pz) / (gen_first_jet_energy - gen_first_jet_pz))

                GenYJ0 = ak.fill_none(GenYJ0, -999)
                GenYJ0 = ak.where(numpy.isnan(GenYJ0), -999, GenYJ0)
                diphotons['GenYJ0'] = GenYJ0

                GenYH = ak.fill_none(GenYH, -999)
                GenYH = ak.where(numpy.isnan(GenYH), -999, GenYH)
                diphotons['GenYH'] = GenYH

                GenAbsPhiHJ0 = numpy.abs(gen_first_jet_phi - GenPhiH)

                # Set all entries above 2*pi to -999
                GenAbsPhiHJ0 = ak.where(
                    GenAbsPhiHJ0 > 2 * numpy.pi,
                    -999,
                    GenAbsPhiHJ0
                )
                GenAbsPhiHJ0_pi_array = ak.full_like(GenAbsPhiHJ0, 2 * numpy.pi)

                # Select the smallest angle
                GenAbsPhiHJ0 = ak.where(
                    GenAbsPhiHJ0 > numpy.pi,
                    GenAbsPhiHJ0_pi_array - GenAbsPhiHJ0,
                    GenAbsPhiHJ0
                )

                diphotons["GenDPhiHJ0"] = GenAbsPhiHJ0

                GenAbsYHJ0 = numpy.abs(GenYJ0 - GenYH)

                # Set all entries above 500 to -999
                GenAbsYHJ0 = ak.where(
                    GenAbsYHJ0 > 500,
                    -999,
                    GenAbsYHJ0
                )

                diphotons["GenDYHJ0"] = GenAbsYHJ0

                logger.debug(f"\n \t[diphotons.fields] after gen block for {variation}: {diphotons.fields}")

            # ------- Baseline modifications to diphotons -------
            # ------- adding diphoton mva -------
            if self.diphoton_mva is not None:
                logger.info(f"\n \tAdding diphoton mva to diphotons of {variation}")
                diphotons = self.add_diphoton_mva(diphotons, events)
            #     dipho_bdt_cut = diphotons["bdt_score"] > 0.7
            #     dipho_bdt_cut_low = (diphotons["bdt_score"] > 0.7) & (
            #         diphotons["bdt_score"] <= 0.8
            #     )
            #     dipho_bdt_cut_med = (diphotons["bdt_score"] > 0.8) & (
            #         diphotons["bdt_score"] <= 0.9
            #     )
            #     dipho_bdt_cut_hig = diphotons["bdt_score"] > 0.9
            #     diphotons["bdt_low"] = dipho_bdt_cut_low
            #     diphotons["bdt_med"] = dipho_bdt_cut_med
            #     diphotons["bdt_hig"] = dipho_bdt_cut_hig
            # else:
            #     logger.info(
            #         "weight_file for diphoton_mva based cut not found, skipping..."
            #     )
            #     dipho_bdt_cut = dipho_events.event >= 0
            #     diphotons["bdt_low"] = dipho_events.event >= 0
            #     diphotons["bdt_med"] = dipho_events.event >= 0
            #     diphotons["bdt_hig"] = dipho_events.event >= 0

            # # select events falling with diphotonID > 0.8, probably there's a smarter way but anyhow...
            # (https://gitlab.cern.ch/tbevilac/higgs-dna-tiziano-bevilacqua/-/blob/HpC_dev/higgs_dna/workflows/HplusCharm_systematics.py?ref_type=heads)
            # diphotons = diphotons[dipho_bdt_cut]
            # sel_jets = sel_jets[dipho_bdt_cut]
            # n_jets = n_jets[dipho_bdt_cut]
            # dipho_events = dipho_events[dipho_bdt_cut]

            # ------- Workflow specific processing -------
            events, process_extra = self.process_extra(events)
            histos_etc.update(process_extra)

            # ------- jets_variables -------

            logger.debug(f"[jets.fields] before zipping: {jets.fields}")
            jets = ak.zip(        # Doubt: but we already have created this earlier: jets = jets_dct[variation], why are we redefining it again?
                {
                    "pt": jets.pt,
                    "eta": jets.eta,
                    "phi": jets.phi,
                    "mass": jets.mass,
                    "charge": ak.zeros_like(jets.pt),  # added this because jet charge is not a property of photons in nanoAOD v11. 
                                                       # We just need the charge to build jet collection.
                    "hFlav": jets.hadronFlavour if self.data_kind == "mc" else ak.ones_like(jets.pt) * -1.,
                    "btagDeepFlav_B": jets.btagDeepFlavB if hasattr(jets, "btagDeepFlavB") else ak.ones_like(jets.pt) * -1.,
                    "btagDeepFlav_CvB": jets.btagDeepFlavCvB if hasattr(jets, "btagDeepFlavCvB") else ak.ones_like(jets.pt) * -1.,
                    "btagDeepFlav_CvL": jets.btagDeepFlavCvL if hasattr(jets, "btagDeepFlavCvL") else ak.ones_like(jets.pt) * -1.,
                    "btagDeepFlav_QG": jets.btagDeepFlavQG if hasattr(jets, "btagDeepFlavQG") else ak.ones_like(jets.pt) * -1.,
                    "btagPNetB": jets.btagPNetB,
                    "btagPNetQvG": jets.btagPNetQvG,
                    "PNetRegPtRawCorr": jets.PNetRegPtRawCorr,
                    "PNetRegPtRawCorrNeutrino": jets.PNetRegPtRawCorrNeutrino,
                    "PNetRegPtRawRes": jets.PNetRegPtRawRes,
                    "btagRobustParTAK4B": jets.btagRobustParTAK4B,
                    "jetId": jets.jetId,
                    "n_sv": jets.nSVs if hasattr(jets, "nSVs") else ak.ones_like(jets.pt) * -1.,
                    "n_muons": jets.nMuons if hasattr(jets, "nMuons") else ak.ones_like(jets.pt) * -1.,
                    "n_electrons": jets.nElectrons if hasattr(jets, "nElectrons") else ak.ones_like(jets.pt) * -1.,
                    "n_const": jets.nConstituents if hasattr(jets, "nConstituents") else ak.ones_like(jets.pt) * -1.,
                }
            )
            jets = ak.with_name(jets, "PtEtaPhiMCandidate")
            logger.debug(f"[jets.fields] after zipping: {jets.fields}")

            # ------- electrons_variables -------
            electrons = ak.zip(
                {
                    "pt": events.Electron.pt,
                    "eta": events.Electron.eta,
                    "phi": events.Electron.phi,
                    "mass": events.Electron.mass,
                    "charge": events.Electron.charge,
                    "cutBased": events.Electron.cutBased,
                    "mvaIso_WP90": events.Electron.mvaIso_WP90,
                    "mvaIso_WP80": events.Electron.mvaIso_WP80,
                    # "mvaIso_Fall17V2_WP90": events.Electron.mvaIso_Fall17V2_WP90 if not self.nAODv10 else events.Electron.mvaIso_WP90,
                    # "mvaIso_Fall17V2_WP80": events.Electron.mvaIso_Fall17V2_WP80 if not self.nAODv10 else events.Electron.mvaIso_WP80,
                    # "mvaIso_Fall17V2_WPL": events.Electron.mvaIso_Fall17V2_WPL if not self.nAODv10 else events.Electron.mvaIso_WPL
                }
            )
            electrons = ak.with_name(electrons, "PtEtaPhiMCandidate")

            # ------- special cut for the base workflow to replicate electrons iso cut in case of muons too -------

            logger.debug(f"len(events) before special muon cut: {len(events)}")
            logger.debug(f"Total number of muons per event before special muon cut: {ak.num(events.Muon)}")
            logger.debug(f"Total  number of muons before special muon cut: {ak.sum(ak.num(events.Muon))}")
            events['Muon'] = events.Muon[events.Muon.pfRelIso03_all < 0.2]
            logger.debug(f"len(events) after special muon cut: {len(events)}")
            logger.debug(f"Total  number of muons per event after special muon cut: {ak.num(events.Muon)}")
            logger.debug(f"Total  number of muons after special muon cut: {ak.sum(ak.num(events.Muon))}")

            # ------- muons_variables -------
            muons = ak.zip(
                {
                    "pt": events.Muon.pt,
                    "eta": events.Muon.eta,
                    "phi": events.Muon.phi,
                    "mass": events.Muon.mass,
                    "charge": events.Muon.charge,
                    "tightId": events.Muon.tightId,
                    "mediumId": events.Muon.mediumId,
                    "looseId": events.Muon.looseId,
                    "isGlobal": events.Muon.isGlobal,
                }
            )
            muons = ak.with_name(muons, "PtEtaPhiMCandidate")

            ## ------- lepton cleaning -------

            logger.debug(f"Number of electrons per event before ele cuts: {ak.num(electrons)}")
            logger.debug(f"Total number of electrons across events before ele cuts: {ak.sum(ak.num(electrons))}")

            logger.info("\n \tApplying lepton cleaning selections\n")

            sel_electrons = electrons[select_electrons(self, electrons, diphotons)]
            events["sel_electrons"] = sel_electrons

            logger.debug(f"Number of electrons per event after ele cuts: {ak.num(sel_electrons)}")
            logger.debug(f"Total number of electrons across events after ele cuts: {ak.sum(ak.num(sel_electrons))}")

            sel_muons = muons[select_muons(self, muons, diphotons)]
            events["sel_muons"] = sel_muons

            ## ------- Combine electrons and muons into a single leptons collection -------
            logger.info("\n \tCombining electrons and muons into a single 'leptons' collection\n")
            leptons = ak.concatenate([sel_electrons, sel_muons], axis=1)
            leptons = ak.with_name(leptons, "PtEtaPhiMCandidate")

            # ------- sort leptons by pt in descending order -------
            leptons = leptons[ak.argsort(leptons.pt, ascending=False)]

            n_leptons = ak.num(leptons)
            diphotons["n_leptons"] = n_leptons

            # ------- annotate diphotons with selected leptons properties -------
            lepton_properties = ["pt", "eta", "phi", "mass", "charge"]
            for i in range(self.num_leptons_to_store):  # Number of leptons to select
                for prop in lepton_properties:
                    key = f"lepton{i+1}_{prop}"
                    # Retrieve the value using the choose_jet function (which can be used for leptons as well)
                    value = choose_jet(getattr(leptons, prop), i, -999.0)
                    # Store the value in the diphotons dictionary
                    diphotons[key] = value

            ## ------- Jet selections -------
            logger.info("\n \tApplying jet selections\n")

            jets = jets[select_jets(self, jets, diphotons, sel_muons, sel_electrons)]
            jets = jets[ak.argsort(jets.pt, ascending=False)]

            events["sel_jets"] = jets

            # ------- add first and second jet props to diphotons -------

            logger.info(f"\n \tAdding first and second jet properties to diphotons of {variation}\n")
            n_jets = ak.num(jets)
            Njets2p5 = ak.num(jets[(jets.pt > 30) & (numpy.abs(jets.eta) < 2.5)])

            first_jet_pt = choose_jet(jets.pt, 0, -999.0)
            first_jet_eta = choose_jet(jets.eta, 0, -999.0)
            first_jet_phi = choose_jet(jets.phi, 0, -999.0)
            first_jet_mass = choose_jet(jets.mass, 0, -999.0)
            first_jet_charge = choose_jet(jets.charge, 0, -999.0)

            second_jet_pt = choose_jet(jets.pt, 1, -999.0)
            second_jet_eta = choose_jet(jets.eta, 1, -999.0)
            second_jet_phi = choose_jet(jets.phi, 1, -999.0)
            second_jet_mass = choose_jet(jets.mass, 1, -999.0)
            second_jet_charge = choose_jet(jets.charge, 1, -999.0)

            diphotons["first_jet_pt"] = first_jet_pt
            diphotons["PTJ0"] = first_jet_pt
            diphotons["first_jet_eta"] = first_jet_eta
            diphotons["first_jet_phi"] = first_jet_phi
            diphotons["first_jet_mass"] = first_jet_mass
            diphotons["first_jet_charge"] = first_jet_charge

            diphotons["second_jet_pt"] = second_jet_pt
            diphotons["PTJ1"] = second_jet_pt
            diphotons["second_jet_eta"] = second_jet_eta
            diphotons["second_jet_phi"] = second_jet_phi
            diphotons["second_jet_mass"] = second_jet_mass
            diphotons["second_jet_charge"] = second_jet_charge

            diphotons["n_jets"] = n_jets
            diphotons["Njets2p5"] = Njets2p5

            diphotons["NJ"] = Njets2p5

            first_jet_pz = first_jet_pt * numpy.sinh(first_jet_eta)
            first_jet_energy = numpy.sqrt((first_jet_pt**2 * numpy.cosh(first_jet_eta)**2) + first_jet_mass**2)

            first_jet_y = 0.5 * numpy.log((first_jet_energy + first_jet_pz) / (first_jet_energy - first_jet_pz))
            first_jet_y = ak.fill_none(first_jet_y, -999)
            first_jet_y = ak.where(numpy.isnan(first_jet_y), -999, first_jet_y)
            diphotons["YJ0"] = first_jet_y

            AbsPhiHJ0 = numpy.abs(first_jet_phi - diphotons["phi"])

            AbsPhiHJ0_pi_array = ak.full_like(AbsPhiHJ0, 2 * numpy.pi)

            # Select the smallest angle
            AbsPhiHJ0 = ak.where(
                AbsPhiHJ0 > numpy.pi,
                AbsPhiHJ0_pi_array - AbsPhiHJ0,
                AbsPhiHJ0
            )
            AbsPhiHJ0 = ak.where(
                AbsPhiHJ0 > 2 * numpy.pi,
                -999,
                ak.where(
                    AbsPhiHJ0 < 0,
                    -999,
                    AbsPhiHJ0
                )
            )
            diphotons["DPhiHJ0"] = AbsPhiHJ0

            AbsYHJ0 = numpy.abs(first_jet_y - diphotons["rapidity"])

            # Set all entries above 500 to -999
            AbsYHJ0 = ak.where(
                AbsYHJ0 > 500,
                -999,
                AbsYHJ0
            )

            diphotons["DYHJ0"] = AbsYHJ0

            # ------- Run taggers on the events list with added diphotons -------
            logger.info(f"\n \tRunning taggers on diphotons in the event for variation, {variation} (creating fields in diphotons for each tagger with its prio)\n")

            # Note: the shape here is ensured to be broadcastable
            for tagger in self.taggers:
                (
                    diphotons["_".join([tagger.name, str(tagger.priority)])],
                    tagger_extra,
                ) = tagger(
                    events, diphotons
                )  # creates new column in diphotons - tagger priority, or 0, also return list of histrograms here?
                histos_etc.update(tagger_extra)

            # ------- decide the best tag -------
            # if there are taggers to run, arbitrate by them first
            # Deal with order of tagger priorities
            # Turn from diphoton jagged array to whether or not an event was selected
            if self.taggers:
                counts = ak.num(diphotons.pt, axis=1)
                # flatten all the tags in diphotons
                flat_tags = numpy.stack(
                    (
                        ak.flatten(
                            diphotons[
                                "_".join([tagger.name, str(tagger.priority)])
                            ]
                        )
                        for tagger in self.taggers
                    ),
                    axis=1,
                )
                tags = ak.from_regular(
                    ak.unflatten(flat_tags, counts), axis=2
                )
                winner = ak.min(tags[tags != 0], axis=2)
                diphotons["best_tag"] = winner
                logger.info(f"\n \tThe best tag is: {diphotons.best_tag}\n")

                # ------- choose diphoton with the highest pT in case > 1 diphoton survives
                # lowest priority is most important (ascending sort)
                # leave in order of diphoton pT in case of ties (stable sort)
                sorted_tagwise = ak.argsort(diphotons.best_tag, stable=True)
                diphotons = diphotons[sorted_tagwise]

            logger.debug(f"\n \tNumber of diphotons per event after all selections: {ak.num(diphotons)}\n")
            logger.debug(f"\n \tTotal number of diphotons across events after all selections: {ak.sum(ak.num(diphotons))}\n")

            # Calculate the number of events that have at least one diphoton
            events_with_diphotons_before = ak.count(diphotons, axis=0)
            logger.debug(f"\n \tNumber of events with at least one diphoton: {events_with_diphotons_before}\n")

            diphotons = ak.firsts(diphotons)

            logger.info(f"\n[diphotons.fields] after the tagger loop: {diphotons.fields}\n")

            # ------- set diphotons as part of the event record -------
            events[f"diphotons_{do_variation}"] = diphotons

            # ------- annotate diphotons with event information -------
            diphotons["event"] = events.event
            diphotons["lumi"] = events.luminosityBlock
            diphotons["run"] = events.run
            diphotons["nPV"] = events.PV.npvs       # nPV just for validation of pileup reweighting
            diphotons["fixedGridRhoAll"] = events.Rho.fixedGridRhoAll
            diphotons = dress_branches(diphotons, events.PV, "PV")
            diphotons = dress_branches(diphotons, events.Rho, "Rho")

            logger.info(f"\n \t[diphotons.fields] after annotating with event information: {diphotons.fields}\n")

            # ------- annotate diphotons with dZ information (difference between z position of GenVtx and PV) as required by flashggfinalfits -------
            if self.data_kind == "mc":
                diphotons["dZ"] = events.GenVtx.z - events.PV.z
                diphotons["genWeight"] = events.genWeight
            # fill zeros for data because there is no GenVtx for data, obviously
            else:
                diphotons["dZ"] = ak.zeros_like(events.PV.z)

            # ------- drop events without a preselected diphoton candidate -------
            # ------- drop events without a tag, if there are tags -------

            logger.info("\n \tApplying selection_mask to drop events without a preselected"
                        "diphoton candidate and without a tag, if there are tags\n")

            if self.taggers:
                selection_mask = ~(
                    ak.is_none(diphotons)
                    | ak.is_none(diphotons.best_tag)
                )
                diphotons = diphotons[selection_mask]
            else:
                selection_mask = ~ak.is_none(diphotons)
                diphotons = diphotons[selection_mask]
            logger.debug(f"[len(diphotons)] after selection_mask: {len(diphotons)})")

            # return if there is no surviving events
            if len(diphotons) == 0:
                logger.debug("No surviving events in this run, return now!")
                return histos_etc

            ## ------- Event weight corrections and systematics -------

            logger.info("\n \tApplying event weight corrections and systematics\n")

            if self.data_kind == "mc":
                # initiate Weight container here, after selection, since event selection cannot easily be applied to weight container afterwards
                event_weights = Weights(size=len(events[selection_mask]))

                logger.debug(f"event_weights.weight before: {event_weights.weight()}, for variation: {variation}")
                logger.debug(f"event_weights.variations before: {event_weights.variations}, for variation: {variation}")

                # ------- applying corrections to event weights -------
                for correction_name in correction_names:
                    if correction_name in available_weight_corrections:
                        logger.info(
                            f"\n \tAdding correction: {correction_name}, to the weight collection of dataset: {dataset_name}, for variation {variation}\n"
                        )
                        varying_function = available_weight_corrections[correction_name]
                        event_weights = varying_function(
                            events=events[selection_mask],
                            photons=events[f"diphotons_{do_variation}"][selection_mask],
                            weights=event_weights,
                            dataset_name=dataset_name,
                            year=self.year[dataset_name][0],
                        )

                logger.debug(f"[event_weights.weight] after corrections: {event_weights.weight()}, for variation: {variation}")
                logger.debug(f"[event_weights.variations] after corrections: {event_weights.variations}, for variation: {variation}")

                # ------- adding systematic variations to event weights (these go to nominal output dataframe) -------
                if do_variation == "nominal":
                    for systematic_name in systematic_names:
                        if systematic_name in available_weight_systematics:
                            logger.info(
                                f"Adding systematic: {systematic_name}, to weight collection of dataset: {dataset_name}, for variation {variation}\n"
                            )
                            if systematic_name == "LHEScale":
                                if hasattr(events, "LHEScaleWeight"):
                                    diphotons["nweight_LHEScale"] = ak.num(
                                        events.LHEScaleWeight[selection_mask],
                                        axis=1,
                                    )
                                    diphotons["weight_LHEScale"] = events.LHEScaleWeight[selection_mask]
                                else:
                                    logger.info(
                                        f"No {systematic_name} Weights in dataset {dataset_name}"
                                    )
                            elif systematic_name == "LHEPdf":
                                if hasattr(events, "LHEPdfWeight"):
                                    # two AlphaS weights are removed
                                    diphotons["nweight_LHEPdf"] = (
                                        ak.num(
                                            events.LHEPdfWeight[selection_mask],
                                            axis=1,
                                        )
                                        - 2
                                    )
                                    diphotons[
                                        "weight_LHEPdf"
                                    ] = events.LHEPdfWeight[selection_mask][
                                        :, :-2
                                    ]
                                else:
                                    logger.info(
                                        f"No {systematic_name} Weights in dataset {dataset_name}"
                                    )
                            else:
                                varying_function = available_weight_systematics[
                                    systematic_name
                                ]
                                event_weights = varying_function(
                                    events=events[selection_mask],
                                    photons=events[f"diphotons_{do_variation}"][selection_mask],
                                    weights=event_weights,
                                    dataset_name=dataset_name,
                                    year=self.year[dataset_name][0],
                                )

                logger.debug(f"[event_weights.weight] after corrections + systs: {event_weights.weight()}, for variation: {variation}")
                logger.debug(f"[event_weights.variations] after corrections + systs: {event_weights.variations}, for variation: {variation}")

                logger.info("\n \tAdding weight variables to diphotons")

                diphotons["weight_central"] = event_weights.weight()

                # ------- store variations with respect to central weight -------
                if do_variation == "nominal":
                    if event_weights.variations:
                        logger.info("\n \tAdding systematic weight variations to nominal output file\n")
                    for modifier in event_weights.variations:
                        diphotons["weight_" + modifier] = event_weights.weight(
                            modifier=modifier
                        )

                # ------- multiply weight by genWeight for normalisation in post-processing chain -------
                event_weights._weight = (
                    events["genWeight"][selection_mask]
                    * diphotons["weight_central"]
                )
                diphotons["weight"] = event_weights.weight()

            # ------- add weight variables (=1) for data for consistent datasets -------
            else:
                diphotons["weight_central"] = ak.ones_like(diphotons["event"])
                diphotons["weight"] = ak.ones_like(diphotons["event"])

            logger.info(f"\n[diphotons.fields] after weights block: {diphotons.fields}")

            ## ------- Add mass resolution uncertainty (sigma_m_over_m) -------
            logger.info("\n \t#------- Adding mass resolution uncertainty (sigma_m_over_m) -------#\n")

            if (self.data_kind == "mc" and self.doFlow_corrections):
                logger.info("\n \tAs self.data_kind == mc and self.doFlow_corrections;"
                            "adding mass resolution uncertainties: sigma_m_over_m, sigma_m_over_m_corr\n")

                diphotons["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                    (
                        diphotons["pho_lead"].raw_energyErr
                        / (
                            diphotons["pho_lead"].pt                  # Note: 'pt*cosh(eta)' is equal to the energy of a four vector
                            * numpy.cosh(diphotons["pho_lead"].eta)   # Note that one needs to call it slightly differently than in
                                                                      # the output of HiggsDNA; as 'pho_lead' -> 'lead' is only done in dumping utils
                        )
                    )
                    ** 2
                    + (
                        diphotons["pho_sublead"].raw_energyErr
                        / (
                            diphotons["pho_sublead"].pt
                            * numpy.cosh(diphotons["pho_sublead"].eta)
                        )
                    )
                    ** 2
                )

                diphotons["sigma_m_over_m_corr"] = 0.5 * numpy.sqrt(
                    (
                        diphotons["pho_lead"].energyErr
                        / (
                            diphotons["pho_lead"].pt
                            * numpy.cosh(diphotons["pho_lead"].eta)
                        )
                    )
                    ** 2
                    + (
                        diphotons["pho_sublead"].energyErr
                        / (
                            diphotons["pho_sublead"].pt
                            * numpy.cosh(diphotons["pho_sublead"].eta)
                        )
                    )
                    ** 2
                )

            else:
                logger.info("\n \tAdding mass resolution uncertainty: sigma_m_over_m\n")
                diphotons["sigma_m_over_m"] = 0.5 * numpy.sqrt(
                    (
                        diphotons["pho_lead"].energyErr
                        / (
                            diphotons["pho_lead"].pt
                            * numpy.cosh(diphotons["pho_lead"].eta)
                        )
                    )
                    ** 2
                    + (
                        diphotons["pho_sublead"].energyErr
                        / (
                            diphotons["pho_sublead"].pt
                            * numpy.cosh(diphotons["pho_sublead"].eta)
                        )
                    )
                    ** 2
                )

            ## ------- Add smeared mass resolution uncertainty (sigma_m_over_m_Smeared) -------
            logger.info("\n \tAdding smeared mass resolution uncertainty (sigma_m_over_m_Smeared)\n")

            # The mass SigmaM/M value including the smearing term from the Scale and smearing
            # The implementation follows the flashGG implementation
            # -> https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/DataFormats/src/Photon.cc#L293
            # adittional flashGG link when the smearing of the SigmaE/E is called
            # -> https://github.com/cms-analysis/flashgg/blob/4edea8897e2a4b0518dca76ba6c9909c20c40ae7/Systematics/plugins/PhotonSigEoverESmearingEGMTool.cc#L83C40-L83C45
            # Just a reminder, the pt/energy of the data is not smearing, but the smearing term is added to the data sigma_m_over_m
            if self.Smear_sigma_m:

                if (self.doFlow_corrections and self.data_kind == "mc"):
                    # adding the smeared BDT error to the ntuples!
                    diphotons["pho_lead","energyErr_Smeared"] = numpy.sqrt(
                    (diphotons["pho_lead"].raw_energyErr)**2 +
                    (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2
                    )

                    diphotons["pho_sublead","energyErr_Smeared"] = numpy.sqrt(
                    (diphotons["pho_sublead"].raw_energyErr) ** 2 +
                    (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2
                    )


                    logger.info("\n \tAs self.data_kind == mc and self.doFlow_corrections;"
                                "adding smeared mass resolution uncertainties: sigma_m_over_m_Smeared, sigma_m_over_m_Smeared_corr\n")
                    diphotons["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                        (
                            numpy.sqrt((diphotons["pho_lead"].raw_energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                            / (
                                diphotons["pho_lead"].pt
                                * numpy.cosh(diphotons["pho_lead"].eta)
                            )
                        )
                        ** 2
                        + (
                            numpy.sqrt((diphotons["pho_sublead"].raw_energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                            / (
                                diphotons["pho_sublead"].pt
                                * numpy.cosh(diphotons["pho_sublead"].eta)
                            )
                        )
                        ** 2
                    )

                    diphotons["sigma_m_over_m_Smeared_corr"] = 0.5 * numpy.sqrt(
                        (
                            numpy.sqrt((diphotons["pho_lead"].energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                            / (
                                diphotons["pho_lead"].pt
                                * numpy.cosh(diphotons["pho_lead"].eta)
                            )
                        )
                        ** 2
                        + (
                            numpy.sqrt((diphotons["pho_sublead"].energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                            / (
                                diphotons["pho_sublead"].pt
                                * numpy.cosh(diphotons["pho_sublead"].eta)
                            )
                        )
                        ** 2
                    )

                else:
                    # adding the smeared BDT error to the ntuples!
                    diphotons["pho_lead", "energyErr_Smeared"] = numpy.sqrt(
                    (diphotons["pho_lead"].energyErr)**2 +
                    (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2
                    )

                    diphotons["pho_sublead","energyErr_Smeared"] = numpy.sqrt(
                    (diphotons["pho_sublead"].energyErr) ** 2 +
                    (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)

                    logger.info("\n \tAdding smeared mass resolution uncertainty: sigma_m_over_m_Smeared")
                    diphotons["sigma_m_over_m_Smeared"] = 0.5 * numpy.sqrt(
                        (
                            numpy.sqrt((diphotons["pho_lead"].energyErr)**2 + (diphotons["pho_lead"].rho_smear * ((diphotons["pho_lead"].pt * numpy.cosh(diphotons["pho_lead"].eta)))) ** 2)
                            / (
                                diphotons["pho_lead"].pt
                                * numpy.cosh(diphotons["pho_lead"].eta)
                            )
                        )
                        ** 2
                        + (
                            numpy.sqrt((diphotons["pho_sublead"].energyErr) ** 2 + (diphotons["pho_sublead"].rho_smear * ((diphotons["pho_sublead"].pt * numpy.cosh(diphotons["pho_sublead"].eta)))) ** 2)
                            / (
                                diphotons["pho_sublead"].pt
                                * numpy.cosh(diphotons["pho_sublead"].eta)
                            )
                        )
                        ** 2
                    )

                logger.info(f"\n[diphotons.fields] after Smear_sigma_m block: {diphotons.fields}")

            ## ------- Decorrelate the mass resolution (still need to supress the decorrelator noises) -------

            if self.doDeco:
                logger.info("\n \tDecorrelating the mass resolution\n")
                logger.info("\n \tAdding various (nominal, nominal_smeared, corr, corr_smeared) decorrelated mass resolutions\n")

                # ------- 1) decorrelate nominal sigma_m_over_m -------
                diphotons["sigma_m_over_m_nominal_decorr"] = decorrelate_mass_resolution(diphotons, type="nominal", year=self.year[dataset_name][0])

                # ------- 2) decorrelate smeared nominal sigma_m_overm_m -------
                if self.Smear_sigma_m:
                    diphotons["sigma_m_over_m_smeared_decorr"] = decorrelate_mass_resolution(diphotons, type="smeared", year=self.year[dataset_name][0])

                # ------- 3) decorrelate flow corrected sigma_m_over_m -------
                if self.doFlow_corrections:
                    diphotons["sigma_m_over_m_corr_decorr"] = decorrelate_mass_resolution(diphotons, type="corr", year=self.year[dataset_name][0])   # doubt/bug: for 'Data', we dont calculate/have 'Sigma_m_over_m_corr' in events/diphotons; so, for Data, 'sigma_m_over_m_corr_decorr' throws error..

                # ------- 4) decorrelate flow smeared corrected sigma_m_over_m -------
                if (self.doFlow_corrections and self.Smear_sigma_m):
                    diphotons["sigma_m_over_m_corr_smeared_decorr"] = decorrelate_mass_resolution(diphotons, type="corr_smeared", year=self.year[dataset_name][0])

                # Instead of the nominal sigma_m_over_m, we will use the smeared version of it -> (https://indico.cern.ch/event/1319585/#169-update-on-the-run-3-mass-r)
                # else:
                #    warnings.warn("Smearing need to be applied in order to decorrelate the (Smeared) mass resolution. -- Exiting!")
                #    sys.exit(0)

                logger.info(f"\n[diphotons.fields] after doDeco block: {diphotons.fields}")

            logger.info(f"After all (pre)selction: Number of events (for {dataset_name} & {variation}): {len(events)}")   # number of events retained even if no photon in a event passes all (pre)selection
            logger.info(f"\n \t#----- diphotons at the end: {ak.type(diphotons)}")

            ## ------- Output processing -------

            if self.output_location is not None:
                if self.output_format == "root":
                    df = diphoton_list_to_pandas(self, diphotons)
                else:
                    akarr = diphoton_ak_array(self, diphotons)

                    # Remove fixedGridRhoAll from photons to avoid having event-level info per photon
                    akarr = akarr[
                        [
                            field
                            for field in akarr.fields
                            if "lead_fixedGridRhoAll" not in field
                        ]
                    ]

                fname = (
                    events.behavior[
                        "__events_factory__"
                    ]._partition_key.replace("/", "_")
                    + ".%s" % self.output_format
                )
                subdirs = []
                if "dataset" in events.metadata:
                    subdirs.append(events.metadata["dataset"])
                subdirs.append(do_variation)
                if self.output_format == "root":
                    dump_pandas(self, df, fname, self.output_location, subdirs)
                else:
                    dump_ak_array(
                        self, akarr, fname, self.output_location, metadata, subdirs,
                    )

        return histos_etc

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
