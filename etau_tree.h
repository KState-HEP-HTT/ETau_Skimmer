#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <cmath>

class etau_tree {
private:
  TTree *tree;
  TTree *original;
  std::vector<Int_t> good_events;

public:
  // Member variables

  // variables copied directly tree to tree
  ULong64_t evt;
  Float_t GenWeight, genpX, genpY, vispX, vispY, genpT, genM, type1_pfMetEt, type1_pfMetPhi, metSig, metcov00, metcov01, metcov10, metcov11, met_EESDown, met_EESUp, met_JESUp, met_JESDown;
  Float_t met_MESDown, met_MESUp, met_PESUp, met_PESDown, met_TESUp, met_TESDown, met_UESUp, met_UESDown, met_JERDown, met_JERUp, metphi_EESDown, metphi_EESUp, metphi_JESUp;
  Float_t metphi_JESDown, metphi_MESDown, metphi_MESUp, metphi_PESUp, metphi_PESDown, metphi_TESUp, metphi_TESDown, metphi_UESUp, metphi_UESDown, metphi_JERDown, metphi_JERUp;
  Float_t mvaMet, mvaMetcov00, mvaMetcov11, mvaMetcov10, mvaMetcov01, mvaMetphi, dphi_12, dphi_emet, dphi_taumet, singleE25eta2p1TightPass, passEle27, eMatchesEle25TightFilter, eMatchesEle25eta2p1TightPath;
  Float_t topQuarkPt1, topQuarkPt2, NUP, jetVeto30, bjetCISVVeto20Medium, jetVeto20, eZTTGenMatching, tZTTGenMatching, tZTTGenDR;
  Float_t tDecayMode, ePVDZ, ePVDXY, eIsoDB03, eCharge, tPVDZ, tPVDXY, tByIsolationMVArun2v1DBoldDMwLTraw, tCharge;
  Float_t tAgainstMuonTight3, tAgainstMuonLoose3, tAgainstElectronVLooseMVA6, tAgainstElectronLooseMVA6, tAgainstElectronMediumMVA6, tAgainstElectronTightMVA6;
  Float_t tAgainstElectronVTightMVA6, tByLooseCombinedIsolationDeltaBetaCorr3Hits, tByMediumCombinedIsolationDeltaBetaCorr3Hits, tByTightCombinedIsolationDeltaBetaCorr3Hits;
  Float_t tByCombinedIsolationDeltaBetaCorrRaw3Hits, tByIsolationMVA3oldDMwLTraw, tByIsolationMVA3newDMwLTraw, tByVLooseIsolationMVArun2v1DBoldDMwLT, tByLooseIsolationMVArun2v1DBoldDMwLT;
  Float_t tByMediumIsolationMVArun2v1DBoldDMwLT, tByTightIsolationMVArun2v1DBoldDMwLT, tByVTightIsolationMVArun2v1DBoldDMwLT, tByVVTightIsolationMVArun2v1DBoldDMwLT;
  Float_t tNeutralIsoPtSum, tChargedIsoPtSum, tPuCorrPtSum, tDecayModeFinding, tDecayModeFindingNewDMs, j1pt, j2pt, j1eta, j2eta, j1phi, j2phi, j1csv, j2csv, jb1pt, jb2pt;
  Float_t jb1eta, jb2eta, jb1phi, jb2phi, jb1csv, jb2csv, nTruePU, nvtx, rho, extratau_veto, isZtt, idisoweight_2, tByIsolationMVArun2v1DBnewDMwLTraw;

  Float_t vbfMass_JetAbsoluteFlavMapUp, vbfMass_JetAbsoluteMPFBiasUp, vbfMass_JetAbsoluteScaleUp, vbfMass_JetAbsoluteStatUp;
  Float_t vbfMass_JetEnUp, vbfMass_JetFlavorQCDUp, vbfMass_JetFragmentationUp, vbfMass_JetPileUpDataMCUp, vbfMass_JetPileUpPtBBUp;
  Float_t vbfMass_JetPileUpPtEC1Up, vbfMass_JetPileUpPtEC2Up, vbfMass_JetPileUpPtHFUp, vbfMass_JetPileUpPtRefUp, vbfMass_JetRelativeBalUp;
  Float_t vbfMass_JetRelativeFSRUp, vbfMass_JetRelativeJEREC1Up, vbfMass_JetRelativeJEREC2Up, vbfMass_JetRelativeJERHFUp;
  Float_t vbfMass_JetRelativePtBBUp, vbfMass_JetRelativePtEC1Up, vbfMass_JetRelativePtEC2Up, vbfMass_JetRelativePtHFUp;
  Float_t vbfMass_JetRelativeStatECUp, vbfMass_JetRelativeStatFSRUp, vbfMass_JetRelativeStatHFUp, vbfMass_JetSinglePionECALUp;
  Float_t vbfMass_JetSinglePionHCALUp, vbfMass_JetTimePtEtaUp, vbfMass_JetAbsoluteFlavMapDown, vbfMass_JetAbsoluteMPFBiasDown;
  Float_t vbfMass_JetAbsoluteScaleDown, vbfMass_JetAbsoluteStatDown, vbfMass_JetEnDown, vbfMass_JetFlavorQCDDown;
  Float_t vbfMass_JetFragmentationDown, vbfMass_JetPileUpDataMCDown, vbfMass_JetPileUpPtBBDown, vbfMass_JetPileUpPtEC1Down;
  Float_t vbfMass_JetPileUpPtEC2Down, vbfMass_JetPileUpPtHFDown, vbfMass_JetPileUpPtRefDown, vbfMass_JetRelativeBalDown;
  Float_t vbfMass_JetRelativeFSRDown, vbfMass_JetRelativeJEREC1Down, vbfMass_JetRelativeJEREC2Down, vbfMass_JetRelativeJERHFDown;
  Float_t vbfMass_JetRelativePtBBDown, vbfMass_JetRelativePtEC1Down, vbfMass_JetRelativePtEC2Down, vbfMass_JetRelativePtHFDown;
  Float_t vbfMass_JetRelativeStatECDown, vbfMass_JetRelativeStatFSRDown, vbfMass_JetRelativeStatHFDown, vbfMass_JetSinglePionECALDown;
  Float_t vbfMass_JetSinglePionHCALDown, vbfMass_JetTimePtEtaDown;

  Float_t jetVeto30_JetAbsoluteFlavMapUp, jetVeto30_JetAbsoluteMPFBiasUp, jetVeto30_JetAbsoluteScaleUp, jetVeto30_JetAbsoluteStatUp;
  Float_t jetVeto30_JetFlavorQCDUp, jetVeto30_JetFragmentationUp, jetVeto30_JetPileUpDataMCUp, jetVeto30_JetPileUpPtBBUp;
  Float_t jetVeto30_JetPileUpPtEC1Up, jetVeto30_JetPileUpPtEC2Up, jetVeto30_JetPileUpPtHFUp, jetVeto30_JetPileUpPtRefUp;
  Float_t jetVeto30_JetRelativeBalUp, jetVeto30_JetRelativeFSRUp, jetVeto30_JetRelativeJEREC1Up, jetVeto30_JetRelativeJEREC2Up;
  Float_t jetVeto30_JetRelativeJERHFUp, jetVeto30_JetRelativePtBBUp, jetVeto30_JetRelativePtEC1Up, jetVeto30_JetRelativePtEC2Up;
  Float_t jetVeto30_JetRelativePtHFUp, jetVeto30_JetRelativeStatECUp, jetVeto30_JetRelativeStatFSRUp, jetVeto30_JetRelativeStatHFUp;
  Float_t jetVeto30_JetSinglePionECALUp, jetVeto30_JetSinglePionHCALUp, jetVeto30_JetTimePtEtaUp, jetVeto30_JetAbsoluteFlavMapDown;
  Float_t jetVeto30_JetAbsoluteMPFBiasDown, jetVeto30_JetAbsoluteScaleDown, jetVeto30_JetAbsoluteStatDown;
  Float_t jetVeto30_JetFlavorQCDDown, jetVeto30_JetFragmentationDown, jetVeto30_JetPileUpDataMCDown, jetVeto30_JetPileUpPtBBDown;
  Float_t jetVeto30_JetPileUpPtEC1Down, jetVeto30_JetPileUpPtEC2Down, jetVeto30_JetPileUpPtHFDown, jetVeto30_JetPileUpPtRefDown;
  Float_t jetVeto30_JetRelativeBalDown, jetVeto30_JetRelativeFSRDown, jetVeto30_JetRelativeJEREC1Down, jetVeto30_JetRelativeJEREC2Down;
  Float_t jetVeto30_JetRelativeJERHFDown, jetVeto30_JetRelativePtBBDown, jetVeto30_JetRelativePtEC1Down, jetVeto30_JetRelativePtEC2Down;
  Float_t jetVeto30_JetRelativePtHFDown, jetVeto30_JetRelativeStatECDown, jetVeto30_JetRelativeStatFSRDown, jetVeto30_JetRelativeStatHFDown;
  Float_t jetVeto30_JetSinglePionECALDown, jetVeto30_JetSinglePionHCALDown, jetVeto30_JetTimePtEtaDown;

  Float_t jetVeto30_JetEnUp, jetVeto30_JetEnDown, jetVeto20_JetEnUp, jetVeto20_JetEnDown;

  Float_t type1_pfMet_shiftedPhi_ElectronEnDown, type1_pfMet_shiftedPhi_ElectronEnUp, type1_pfMet_shiftedPhi_JetEnDown;
  Float_t type1_pfMet_shiftedPhi_JetEnUp, type1_pfMet_shiftedPhi_JetResDown, type1_pfMet_shiftedPhi_JetResUp;
  Float_t type1_pfMet_shiftedPhi_MuonEnDown, type1_pfMet_shiftedPhi_MuonEnUp, type1_pfMet_shiftedPhi_PhotonEnDown;
  Float_t type1_pfMet_shiftedPhi_PhotonEnUp, type1_pfMet_shiftedPhi_TauEnDown, type1_pfMet_shiftedPhi_TauEnUp;
  Float_t type1_pfMet_shiftedPhi_UnclusteredEnDown, type1_pfMet_shiftedPhi_UnclusteredEnUp, type1_pfMet_shiftedPt_ElectronEnDown;
  Float_t type1_pfMet_shiftedPt_ElectronEnUp, type1_pfMet_shiftedPt_JetEnDown, type1_pfMet_shiftedPt_JetEnUp;
  Float_t type1_pfMet_shiftedPt_JetResDown, type1_pfMet_shiftedPt_JetResUp, type1_pfMet_shiftedPt_MuonEnDown;
  Float_t type1_pfMet_shiftedPt_MuonEnUp, type1_pfMet_shiftedPt_PhotonEnDown, type1_pfMet_shiftedPt_PhotonEnUp;
  Float_t type1_pfMet_shiftedPt_TauEnUp, type1_pfMet_shiftedPt_UnclusteredEnDown, type1_pfMet_shiftedPt_UnclusteredEnUp;
  Float_t type1_pfMet_shiftedPt_TauEnDown;

  Float_t m_coll, m_coll_uesU, m_coll_uesD, m_coll_jesU, m_coll_jesD, m_coll_tesU, m_coll_tesD;

  Int_t run, lumi, njetingap, njetingap20;
  Int_t njetingap_JESUp, njetingap20_JESUp, njetingap_JESDown, njetingap20_JESDown, njets_JESUp, njetspt20_JESUp, njets_JESDown, njetspt20_JESDown;

  // temporary storage variables
  Float_t eVetoZTTp001dxyzR0, muVetoZTTp001dxyzR0, dielectronVeto;

  // new variables to store in new tree
  Float_t met_px, met_py, extraelec_veto, dilepton_veto, m_1, pt_1, eta_1, phi_1, e_1, px_1, py_1, pz_1, m_2, pt_2, eta_2, phi_2, e_2, px_2, py_2, pz_2, dijetphi, hdijetphi, visjeteta, isZet;
  Float_t jdeta, mjj, dijetpt, mjj_JESUp, jdeta_JESUp, mjj_JESDown, jdeta_JESDown;
  Float_t weight;

  // forgotten
  Float_t byIsolationMVA3oldDMwoLTraw_2, trigweight_2, byIsolationMVA3newDMwoLTraw_2, filterEle27, ePt, eMass, ePhi, eEta, tPhi, tEta, tMass, tPt, numGenJets, vbfDeta_JetEnDown, vbfDeta_JetEnUp;
  Float_t vbfMass, vbfJetVeto20, vbfJetVeto20_JetEnDown, vbfJetVeto20_JetEnUp, vbfJetVeto30, vbfJetVeto30_JetEnDown, vbfJetVeto30_JetEnUp, eGenPdgId, vbfDeta, extramuon_veto;
  Float_t eMVANonTrigWP80, ePassesConversionVeto, eMissingHits, e_t_DR;

  Float_t Flag_BadChargedCandidateFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_HBHENoiseFilter;
  Float_t Flag_badCloneMuonFilter, Flag_badGlobalMuonFilter;
  Float_t Flag_eeBadScFilter, Flag_globalTightHalo2016Filter, Flag_goodVertices, Flag_BadPFMuonFilter, Flag_HBHENoiseIsoFilter;
  Float_t e_t_MvaMetCovMatrix00, e_t_MvaMetCovMatrix10, e_t_MvaMetCovMatrix01, e_t_MvaMetCovMatrix11;
  
  // Member functions
  etau_tree (TTree* orig, TTree* itree);
  virtual ~etau_tree () {};
  void do_skimming();
  void set_branches();
  TTree* fill_tree();
};

//////////////////////////////////////////////////////////////////
// Purpose: Initialize tree and original, then read branches    //
//          needed for skimming/sorting from Original           //
//////////////////////////////////////////////////////////////////
// Parameters                                                   //
//   - Original: A tree read from an input root file.           //
//               At this point, only the branches necessary     //
//               for skim selection and sorting are read        //
//   - itree: A newly constructed tree. This tree will be       //
//            filled for all events passing the skim selection  //
//////////////////////////////////////////////////////////////////
etau_tree::etau_tree(TTree* Original, TTree* itree) :
tree(itree),
original(Original)
{
  // read only what is needed for skimming and sorting
  original->SetBranchAddress("evt", &evt);
  original->SetBranchAddress("ePt", &ePt);
  original->SetBranchAddress("eEta", &eEta);
  original->SetBranchAddress("ePhi", &ePhi);
  original->SetBranchAddress("ePassesConversionVeto", &ePassesConversionVeto);
  original->SetBranchAddress("ePVDZ", &ePVDZ);
  original->SetBranchAddress("ePVDXY", &ePVDXY);
  original->SetBranchAddress("tPt", &tPt);
  original->SetBranchAddress("tEta", &tEta);
  original->SetBranchAddress("tPhi", &tPhi);
  original->SetBranchAddress("tDecayModeFinding", &tDecayModeFinding);
  original->SetBranchAddress("tPVDXY", &tPVDXY);
  original->SetBranchAddress("tPVDZ", &tPVDZ);
  original->SetBranchAddress("tByVLooseIsolationMVArun2v1DBoldDMwLT", &tByVLooseIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByTightIsolationMVArun2v1DBoldDMwLT", &tByTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("dielectronVeto", &dielectronVeto);
  original->SetBranchAddress("tAgainstMuonLoose3", &tAgainstMuonLoose3);
  original->SetBranchAddress("tAgainstElectronTightMVA6", &tAgainstElectronTightMVA6);
  original->SetBranchAddress("eCharge", &eCharge);
  original->SetBranchAddress("tCharge", &tCharge);
  original->SetBranchAddress("eIsoDB03", &eIsoDB03);
  original->SetBranchAddress("eMVANonTrigWP80", &eMVANonTrigWP80);
  original->SetBranchAddress("singleE25eta2p1TightPass", &singleE25eta2p1TightPass);
  original->SetBranchAddress("eMatchesEle25TightFilter", &eMatchesEle25TightFilter);
  original->SetBranchAddress("eMatchesEle25eta2p1TightPath", &eMatchesEle25eta2p1TightPath);
  original->SetBranchAddress("eMissingHits", &eMissingHits);
  original->SetBranchAddress("eVetoZTTp001dxyzR0", &eVetoZTTp001dxyzR0);
  original->SetBranchAddress("tByMediumIsolationMVArun2v1DBoldDMwLT", &tByMediumIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("e_t_DR", &e_t_DR);
  original->SetBranchAddress("muVetoZTTp001dxyzR0", &muVetoZTTp001dxyzR0);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void etau_tree::do_skimming() {

  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> eleCandidate, tauCandidate;

  Int_t nevt = (Int_t)original->GetEntries();
  for (auto ievt = 0; ievt < nevt; ievt++) {
    original->GetEntry(ievt);
    evt_now = evt;

    // apply event selection 
    if (!singleE25eta2p1TightPass || !eMatchesEle25TightFilter || !eMatchesEle25eta2p1TightPath) // apply trigger HLT Ele25 eta2p1 WPTight Gsf with matching
      continue;

    if (ePt < 24 || fabs(eEta) > 2.1 || fabs(ePVDZ) > 0.2 || fabs(ePVDXY) > 0.045 || !eMVANonTrigWP80 || !ePassesConversionVeto || eMissingHits > 1 || eVetoZTTp001dxyzR0 > 1) // electron selection
      continue;

    if (tPt < 27 || fabs(tEta) > 2.3 || fabs(tPVDZ) > 0.2 || !tByVLooseIsolationMVArun2v1DBoldDMwLT || !tDecayModeFinding || fabs(tCharge) > 1) // tau selection
      continue;

    if (muVetoZTTp001dxyzR0 > 0 || dielectronVeto > 0)
      continue;

    if (e_t_DR < 0.5)
      continue;

    // implement new sorting per 
	  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
	  if (evt_now != evt_before) { // new event, save the tau candidates
      //   since it is new event, do we have the best entry to save? If yes, save it!
      if ( best_evt > -1  )
        good_events.push_back(best_evt);
      
      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      eleCandidate = std::make_pair(ePt, eIsoDB03);
      tauCandidate  = std::make_pair(tPt,  tByVLooseIsolationMVArun2v1DBoldDMwLT);
    } 
    else { // not a new event

      std::pair<float, float> currEleCandidate(ePt, eIsoDB03);
      std::pair<float, float> currTauCandidate(tPt, tByVLooseIsolationMVArun2v1DBoldDMwLT);

      // clause 1, select the pair that has most isolated tau lepton 1
      if (currEleCandidate.second - eleCandidate.second  < 0.0001 ) best_evt = ievt;

      // check if the first tau is the same, and if so - move to clause 2
      if ( fabs(currEleCandidate.second - eleCandidate.second)  <  0.0001 ) {
        // pick up  the pair with the highest pT of the first candidate
        if (currEleCandidate.first - eleCandidate.first > 0.0001 ) best_evt = ievt;
        if ( fabs(currEleCandidate.first -eleCandidate.first) < 0.0001 ) { 
          // same pT, same iso, move to clause 3
          if (currTauCandidate.second - tauCandidate.second > 0.0001 ) best_evt = ievt;
          if ( fabs(currTauCandidate.second - tauCandidate.second) < 0.0001 ) {
            // same iso - pick the pair with the highest pT
            if ( currTauCandidate.first - tauCandidate.first  > 0.0001 ) best_evt = ievt;
          } // tau2 has the same isolation
        } // tau1 has the same pT
      } // tau1 has the same isolation
    } // not a new event
    evt_before = evt_now;
  }
}

//////////////////////////////////////////////////////////////////
// Purpose: Fill tree with variables from original and new      //
//          variables. Only events that have been skimmed,      //
//          sorted, and stored in the good_events vector will   //
//          be stored in the tree.                              //
//////////////////////////////////////////////////////////////////
// Return: The same TTree passed to the constructor and stored  //
//         in original, but now it is filled with good events   //
//////////////////////////////////////////////////////////////////
TTree* etau_tree::fill_tree() {

  set_branches();  // get all the branches set up

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    met_px = type1_pfMetEt * cos(type1_pfMetPhi);
    met_py = type1_pfMetEt * sin(type1_pfMetPhi);

    extraelec_veto = eVetoZTTp001dxyzR0 > 1;
    extramuon_veto = muVetoZTTp001dxyzR0 > 0;
    dilepton_veto = dielectronVeto > 0;

    TLorentzVector ele;
    TLorentzVector tau;
    ele.SetPtEtaPhiM(ePt, eEta, ePhi, eMass);
    tau.SetPtEtaPhiM(tPt, tEta, tPhi, tMass);

    m_1 = ele.M();
    px_1 = ele.Px();
    py_1 = ele.Py();
    pz_1 = ele.Pz();
    e_1 = ele.E();
    pt_1 = ele.Pt();
    phi_1 = ele.Phi();
    eta_1 = ele.Eta();
    m_2 = tau.M();
    px_2 = tau.Px();
    py_2 = tau.Py();
    pz_2 = tau.Pz();
    e_2 = tau.E();
    pt_2 = tau.Pt();
    phi_2 = tau.Phi();
    eta_2 = tau.Eta();

    TLorentzVector h = ele + tau;
    TLorentzVector jet1;
    if (jetVeto20 > 0 && j1pt > 0)
      jet1.SetPtEtaPhiM(j1pt, j1eta, j1phi, 0);
    TLorentzVector jet2;
    if (jetVeto20 > 1 && j2pt > 0)
      jet2.SetPtEtaPhiM(j2pt, j2eta, j2phi, 0);
    TLorentzVector dijet=jet1+jet2;

    if (jetVeto20 > 1)
    {
      jdeta = vbfDeta;
      dijetphi = dijet.Phi();
      hdijetphi = h.DeltaPhi(dijet);
      visjeteta = h.Eta() - dijet.Eta();
      mjj = vbfMass;
      njetingap20 = vbfJetVeto20;
      njetingap = vbfJetVeto30;
    }
    else{
       jdeta = -10000;
       dijetpt = -10000;
       dijetphi = -10000;
       hdijetphi = -10000;
       visjeteta = -10000;
       mjj = -10000;
       njetingap20 = -10000;
       njetingap = -100000;
    }
    if (njetspt20_JESUp > 1){
       njetingap20_JESUp = vbfJetVeto20_JetEnUp;
       njetingap_JESUp = vbfJetVeto30_JetEnUp;
       mjj_JESUp = vbfMass_JetEnUp;
       jdeta_JESUp = vbfDeta_JetEnUp;
    }
    else{
       jdeta_JESUp = -10000;
       mjj_JESUp = -10000;
       njetingap20_JESUp = -10000;
       njetingap_JESUp = -100000;
    }
    if (njetspt20_JESDown > 1){
       njetingap20_JESDown = vbfJetVeto20_JetEnDown;
       njetingap_JESDown = vbfJetVeto30_JetEnDown;
       mjj_JESDown = vbfMass_JetEnDown;
       jdeta_JESDown = vbfDeta_JetEnDown;
    }
    else{
       jdeta_JESDown  = -10000;
       mjj_JESDown  = -10000;
       njetingap20_JESDown  = -10000;
       njetingap_JESDown  = -100000;
    }

    isZet = isZtt && fabs(eGenPdgId) == 11;

    tree->Fill();
  }
  return tree;
}

//////////////////////////////////////////////////////////////////
// Purpose:                                                     //
//   - Read all branches we need from original. These branches  //
//          branches may be directly stored, be copied with a   //
//          new name, or used to construct new variables.       //
//   - Create branches in tree to store any variable we want    //
//////////////////////////////////////////////////////////////////
void etau_tree::set_branches() {

  // straight from input tree
  tree->Branch("run", &run, "run/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("evt", &evt);
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("metcov00", &metcov00, "metcov00/F");
  tree->Branch("metcov10", &metcov10, "metcov10/F");
  tree->Branch("metcov11", &metcov11, "metcov11/F");
  tree->Branch("metcov01", &metcov01, "metcov01/F");
  tree->Branch("NUP", &NUP, "NUP/I");
  tree->Branch("genpX", &genpX, "genpX/F");
  tree->Branch("genpY", &genpY, "genpY/F");
  tree->Branch("genM", &genM, "genM/F");
  tree->Branch("genpT", &genpT, "genpT/F");
  tree->Branch("numGenJets", &numGenJets, "numGenJets/F");
  tree->Branch("tZTTGenDR", &tZTTGenDR, "tZTTGenDR/F");

  // from input tree, but change name
  tree->Branch("q_1", &eCharge, "q_1/F");
  tree->Branch("d0_1", &ePVDXY, "d0_1/F");
  tree->Branch("dZ_1", &ePVDZ, "dZ_1/F");
  tree->Branch("iso_1", &eIsoDB03, "iso_1/F");
  tree->Branch("q_2", &tCharge, "q_2/F");
  tree->Branch("d0_2", &tPVDXY, "d0_2/F");
  tree->Branch("dZ_2", &tPVDZ, "dZ_2/F");
  tree->Branch("iso_2", &tByIsolationMVArun2v1DBoldDMwLTraw, "iso_2/F");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/F");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec/F");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon/F");
  tree->Branch("gen_match_1", &eZTTGenMatching, "gen_match_1/I");
  tree->Branch("gen_match_2", &tZTTGenMatching, "gen_match_2/I");
  tree->Branch("againstElectronLooseMVA6_2", &tAgainstElectronLooseMVA6, "againstElectronLooseMVA6_2/F");
  tree->Branch("againstElectronMediumMVA6_2", &tAgainstElectronMediumMVA6, "againstElectronMediumMVA6_2/F");
  tree->Branch("againstElectronTightMVA6_2", &tAgainstElectronTightMVA6, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &tAgainstElectronVLooseMVA6, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstElectronVTightMVA6_2", &tAgainstElectronVTightMVA6, "againstElectronVTightMVA6_2/F");
  tree->Branch("againstMuonLoose3_2", &tAgainstMuonLoose3, "againstMuonLoose3_2/F");
  tree->Branch("againstMuonTight3_2", &tAgainstMuonTight3, "againstMuonTight3_2/F");
  tree->Branch("npv", &nTruePU, "npv/F");
  tree->Branch("npu", &nvtx, "npu/F");
  tree->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &tByLooseCombinedIsolationDeltaBetaCorr3Hits, "byLooseCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &tByMediumCombinedIsolationDeltaBetaCorr3Hits, "byMediumCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &tByTightCombinedIsolationDeltaBetaCorr3Hits, "byTightCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &tByCombinedIsolationDeltaBetaCorrRaw3Hits, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
  tree->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, "byIsolationMVA3oldDMwoLTraw_2/F");
  tree->Branch("byIsolationMVA3oldDMwLTraw_2", &tByIsolationMVA3oldDMwLTraw, "byIsolationMVA3oldDMwLTraw_2/F");
  tree->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, "byIsolationMVA3newDMwoLTraw_2/F");
  tree->Branch("byIsolationMVA3newDMwLTraw_2", &tByIsolationMVA3newDMwLTraw, "byIsolationMVA3newDMwLTraw_2/F");
  tree->Branch("chargedIsoPtSum_2", &tChargedIsoPtSum, "chargedIsoPtSum_2/F");
  tree->Branch("decayModeFinding_2", &tDecayModeFinding, "decayModeFinding_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &tDecayModeFindingNewDMs, "decayModeFindingNewDMs_2/F");
  tree->Branch("neutralIsoPtSum_2", &tNeutralIsoPtSum, "neutralIsoPtSum_2/F");
  tree->Branch("puCorrPtSum_2", &tPuCorrPtSum, "puCorrPtSum_2/F");
  tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  tree->Branch("idisoweight_2", &idisoweight_2, "idisoweight_2/F");
  tree->Branch("met", &type1_pfMetEt, "met/F");
  tree->Branch("metphi", &type1_pfMetPhi, "metphi/F");
  tree->Branch("nbtag", &bjetCISVVeto20Medium, "nbtag/F");
  tree->Branch("njets", &jetVeto30, "njets/F");
  tree->Branch("njetspt20", &jetVeto20, "njetspt20/F");
  tree->Branch("jpt_1", &j1pt, "jpt_1/F");
  tree->Branch("jeta_1", &j1eta, "jeta_1/F");
  tree->Branch("jphi_1", &j1phi, "jphi_1/F");
  tree->Branch("jcsv_1", &j1csv, "jcsv_1/F");
  tree->Branch("jpt_2", &j2pt, "jpt_2/F");
  tree->Branch("jeta_2", &j2eta, "jeta_2/F");
  tree->Branch("jphi_2", &j2phi, "jphi_2/F");
  tree->Branch("jcsv_2", &j2csv, "jcsv_2/F");
  tree->Branch("bpt_1", &jb1pt, "bpt_1/F");
  tree->Branch("beta_1", &jb1eta, "beta_1/F");
  tree->Branch("bphi_1", &jb1phi, "bphi_1/F");
  tree->Branch("bcsv_1", &jb1csv, "bcsv_1/F");
  tree->Branch("bpt_2", &jb2pt, "bpt_2/F");
  tree->Branch("beta_2", &jb2eta, "beta_2/F");
  tree->Branch("bphi_2", &jb2phi, "bphi_2/F");
  tree->Branch("bcsv_2", &jb2csv, "bcsv_2/F");
  tree->Branch("MVANonTrigWP80_1", &eMVANonTrigWP80, "MVANonTrigWP80_1/F");
  tree->Branch("passEle25", &singleE25eta2p1TightPass, "passEle25/F");
  tree->Branch("filterEle25", &eMatchesEle25TightFilter, "filterEle25/F");
  tree->Branch("matchEle25", &eMatchesEle25eta2p1TightPath, "matchEle25/F");
  tree->Branch("l2_decayMode", &tDecayMode, "l2_decayMode/F");
  tree->Branch("genweight", &GenWeight, "genweight/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &tByTightIsolationMVArun2v1DBoldDMwLT, "byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2", &tByMediumIsolationMVArun2v1DBoldDMwLT, "byMediumIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("pt_top1", &topQuarkPt1, "pt_top1/F");
  tree->Branch("pt_top2", &topQuarkPt2, "pt_top2/F");
  tree->Branch("metSig", &metSig, "metSig/F");
  tree->Branch("flag_BadChargedCandidate", &Flag_BadChargedCandidateFilter);
  tree->Branch("flag_BadPFMuon", &Flag_BadPFMuonFilter);
  tree->Branch("flag_EcalDeadCellTriggerPrimitive", &Flag_EcalDeadCellTriggerPrimitiveFilter);
  tree->Branch("flag_HBHENoise", &Flag_HBHENoiseFilter);
  tree->Branch("flag_HBHENoiseIso", &Flag_HBHENoiseIsoFilter);
  tree->Branch("flag_badCloneMuon", &Flag_badCloneMuonFilter);
  tree->Branch("flag_badGlobalMuon", &Flag_badGlobalMuonFilter);
  tree->Branch("flag_eeBadSc", &Flag_eeBadScFilter);
  tree->Branch("flag_globalTightHalo2016", &Flag_globalTightHalo2016Filter);
  tree->Branch("flag_goodVertices", &Flag_goodVertices);
  tree->Branch("metcov00_v2", &e_t_MvaMetCovMatrix00, "metcov00_v2/F");
  tree->Branch("metcov10_v2", &e_t_MvaMetCovMatrix10, "metcov10_v2/F");
  tree->Branch("metcov11_v2", &e_t_MvaMetCovMatrix01, "metcov11_v2/F");
  tree->Branch("metcov01_v2", &e_t_MvaMetCovMatrix11, "metcov01_v2/F");


  tree->Branch("njets_JetAbsoluteFlavMapUp", &jetVeto30_JetAbsoluteFlavMapUp);
  tree->Branch("njets_JetAbsoluteMPFBiasUp", &jetVeto30_JetAbsoluteMPFBiasUp);
  tree->Branch("njets_JetAbsoluteScaleUp", &jetVeto30_JetAbsoluteScaleUp);
  tree->Branch("njets_JetAbsoluteStatUp", &jetVeto30_JetAbsoluteStatUp);
  tree->Branch("njets_JetEnUp", &jetVeto30_JetEnUp);
  tree->Branch("njets_JetFlavorQCDUp", &jetVeto30_JetFlavorQCDUp);
  tree->Branch("njets_JetFragmentationUp", &jetVeto30_JetFragmentationUp);
  tree->Branch("njets_JetPileUpDataMCUp", &jetVeto30_JetPileUpDataMCUp);
  tree->Branch("njets_JetPileUpPtBBUp", &jetVeto30_JetPileUpPtBBUp);
  tree->Branch("njets_JetPileUpPtEC1Up", &jetVeto30_JetPileUpPtEC1Up);
  tree->Branch("njets_JetPileUpPtEC2Up", &jetVeto30_JetPileUpPtEC2Up);
  tree->Branch("njets_JetPileUpPtHFUp", &jetVeto30_JetPileUpPtHFUp);
  tree->Branch("njets_JetPileUpPtRefUp", &jetVeto30_JetPileUpPtRefUp);
  tree->Branch("njets_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  tree->Branch("njets_JetRelativeFSRUp", &jetVeto30_JetRelativeFSRUp);
  tree->Branch("njets_JetRelativeJEREC1Up", &jetVeto30_JetRelativeJEREC1Up);
  tree->Branch("njets_JetRelativeJEREC2Up", &jetVeto30_JetRelativeJEREC2Up);
  tree->Branch("njets_JetRelativeJERHFUp", &jetVeto30_JetRelativeJERHFUp);
  tree->Branch("njets_JetRelativePtBBUp", &jetVeto30_JetRelativePtBBUp);
  tree->Branch("njets_JetRelativePtEC1Up", &jetVeto30_JetRelativePtEC1Up);
  tree->Branch("njets_JetRelativePtEC2Up", &jetVeto30_JetRelativePtEC2Up);
  tree->Branch("njets_JetRelativePtHFUp", &jetVeto30_JetRelativePtHFUp);
  tree->Branch("njets_JetRelativeStatECUp", &jetVeto30_JetRelativeStatECUp);
  tree->Branch("njets_JetRelativeStatFSRUp", &jetVeto30_JetRelativeStatFSRUp);
  tree->Branch("njets_JetRelativeStatHFUp", &jetVeto30_JetRelativeStatHFUp);
  tree->Branch("njets_JetSinglePionECALUp", &jetVeto30_JetSinglePionECALUp);
  tree->Branch("njets_JetSinglePionHCALUp", &jetVeto30_JetSinglePionHCALUp);
  tree->Branch("njets_JetTimePtEtaUp", &jetVeto30_JetTimePtEtaUp);

  tree->Branch("njets_JetAbsoluteFlavMapDown", &jetVeto30_JetAbsoluteFlavMapDown);
  tree->Branch("njets_JetAbsoluteMPFBiasDown", &jetVeto30_JetAbsoluteMPFBiasDown);
  tree->Branch("njets_JetAbsoluteScaleDown", &jetVeto30_JetAbsoluteScaleDown);
  tree->Branch("njets_JetAbsoluteStatDown", &jetVeto30_JetAbsoluteStatDown);
  tree->Branch("njets_JetEnDown", &jetVeto30_JetEnDown);
  tree->Branch("njets_JetFlavorQCDDown", &jetVeto30_JetFlavorQCDDown);
  tree->Branch("njets_JetFragmentationDown", &jetVeto30_JetFragmentationDown);
  tree->Branch("njets_JetPileUpDataMCDown", &jetVeto30_JetPileUpDataMCDown);
  tree->Branch("njets_JetPileUpPtBBDown", &jetVeto30_JetPileUpPtBBDown);
  tree->Branch("njets_JetPileUpPtEC1Down", &jetVeto30_JetPileUpPtEC1Down);
  tree->Branch("njets_JetPileUpPtEC2Down", &jetVeto30_JetPileUpPtEC2Down);
  tree->Branch("njets_JetPileUpPtHFDown", &jetVeto30_JetPileUpPtHFDown);
  tree->Branch("njets_JetPileUpPtRefDown", &jetVeto30_JetPileUpPtRefDown);
  tree->Branch("njets_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  tree->Branch("njets_JetRelativeFSRDown", &jetVeto30_JetRelativeFSRDown);
  tree->Branch("njets_JetRelativeJEREC1Down", &jetVeto30_JetRelativeJEREC1Down);
  tree->Branch("njets_JetRelativeJEREC2Down", &jetVeto30_JetRelativeJEREC2Down);
  tree->Branch("njets_JetRelativeJERHFDown", &jetVeto30_JetRelativeJERHFDown);
  tree->Branch("njets_JetRelativePtBBDown", &jetVeto30_JetRelativePtBBDown);
  tree->Branch("njets_JetRelativePtEC1Down", &jetVeto30_JetRelativePtEC1Down);
  tree->Branch("njets_JetRelativePtEC2Down", &jetVeto30_JetRelativePtEC2Down);
  tree->Branch("njets_JetRelativePtHFDown", &jetVeto30_JetRelativePtHFDown);
  tree->Branch("njets_JetRelativeStatECDown", &jetVeto30_JetRelativeStatECDown);
  tree->Branch("njets_JetRelativeStatFSRDown", &jetVeto30_JetRelativeStatFSRDown);
  tree->Branch("njets_JetRelativeStatHFDown", &jetVeto30_JetRelativeStatHFDown);
  tree->Branch("njets_JetSinglePionECALDown", &jetVeto30_JetSinglePionECALDown);
  tree->Branch("njets_JetSinglePionHCALDown", &jetVeto30_JetSinglePionHCALDown);
  tree->Branch("njets_JetTimePtEtaDown", &jetVeto30_JetTimePtEtaDown);

  tree->Branch("mjj_JetAbsoluteFlavMapUp", &vbfMass_JetAbsoluteFlavMapUp);
  tree->Branch("mjj_JetAbsoluteMPFBiasUp", &vbfMass_JetAbsoluteMPFBiasUp);
  tree->Branch("mjj_JetAbsoluteScaleUp", &vbfMass_JetAbsoluteScaleUp);
  tree->Branch("mjj_JetAbsoluteStatUp", &vbfMass_JetAbsoluteStatUp);
  tree->Branch("mjj_JetEnUp", &vbfMass_JetEnUp);
  tree->Branch("mjj_JetFlavorQCDUp", &vbfMass_JetFlavorQCDUp);
  tree->Branch("mjj_JetFragmentationUp", &vbfMass_JetFragmentationUp);
  tree->Branch("mjj_JetPileUpDataMCUp", &vbfMass_JetPileUpDataMCUp);
  tree->Branch("mjj_JetPileUpPtBBUp", &vbfMass_JetPileUpPtBBUp);
  tree->Branch("mjj_JetPileUpPtEC1Up", &vbfMass_JetPileUpPtEC1Up);
  tree->Branch("mjj_JetPileUpPtEC2Up", &vbfMass_JetPileUpPtEC2Up);
  tree->Branch("mjj_JetPileUpPtHFUp", &vbfMass_JetPileUpPtHFUp);
  tree->Branch("mjj_JetPileUpPtRefUp", &vbfMass_JetPileUpPtRefUp);
  tree->Branch("mjj_JetRelativeBalUp", &vbfMass_JetRelativeBalUp);
  tree->Branch("mjj_JetRelativeFSRUp", &vbfMass_JetRelativeFSRUp);
  tree->Branch("mjj_JetRelativeJEREC1Up", &vbfMass_JetRelativeJEREC1Up);
  tree->Branch("mjj_JetRelativeJEREC2Up", &vbfMass_JetRelativeJEREC2Up);
  tree->Branch("mjj_JetRelativeJERHFUp", &vbfMass_JetRelativeJERHFUp);
  tree->Branch("mjj_JetRelativePtBBUp", &vbfMass_JetRelativePtBBUp);
  tree->Branch("mjj_JetRelativePtEC1Up", &vbfMass_JetRelativePtEC1Up);
  tree->Branch("mjj_JetRelativePtEC2Up", &vbfMass_JetRelativePtEC2Up);
  tree->Branch("mjj_JetRelativePtHFUp", &vbfMass_JetRelativePtHFUp);
  tree->Branch("mjj_JetRelativeStatECUp", &vbfMass_JetRelativeStatECUp);
  tree->Branch("mjj_JetRelativeStatFSRUp", &vbfMass_JetRelativeStatFSRUp);
  tree->Branch("mjj_JetRelativeStatHFUp", &vbfMass_JetRelativeStatHFUp);
  tree->Branch("mjj_JetSinglePionECALUp", &vbfMass_JetSinglePionECALUp);
  tree->Branch("mjj_JetSinglePionHCALUp", &vbfMass_JetSinglePionHCALUp);
  tree->Branch("mjj_JetTimePtEtaUp", &vbfMass_JetTimePtEtaUp);

  tree->Branch("mjj_JetAbsoluteFlavMapDown", &vbfMass_JetAbsoluteFlavMapDown);
  tree->Branch("mjj_JetAbsoluteMPFBiasDown", &vbfMass_JetAbsoluteMPFBiasDown);
  tree->Branch("mjj_JetAbsoluteScaleDown", &vbfMass_JetAbsoluteScaleDown);
  tree->Branch("mjj_JetAbsoluteStatDown", &vbfMass_JetAbsoluteStatDown);
  tree->Branch("mjj_JetEnDown", &vbfMass_JetEnDown);
  tree->Branch("mjj_JetFlavorQCDDown", &vbfMass_JetFlavorQCDDown);
  tree->Branch("mjj_JetFragmentationDown", &vbfMass_JetFragmentationDown);
  tree->Branch("mjj_JetPileUpDataMCDown", &vbfMass_JetPileUpDataMCDown);
  tree->Branch("mjj_JetPileUpPtBBDown", &vbfMass_JetPileUpPtBBDown);
  tree->Branch("mjj_JetPileUpPtEC1Down", &vbfMass_JetPileUpPtEC1Down);
  tree->Branch("mjj_JetPileUpPtEC2Down", &vbfMass_JetPileUpPtEC2Down);
  tree->Branch("mjj_JetPileUpPtHFDown", &vbfMass_JetPileUpPtHFDown);
  tree->Branch("mjj_JetPileUpPtRefDown", &vbfMass_JetPileUpPtRefDown);
  tree->Branch("mjj_JetRelativeBalDown", &vbfMass_JetRelativeBalDown);
  tree->Branch("mjj_JetRelativeFSRDown", &vbfMass_JetRelativeFSRDown);
  tree->Branch("mjj_JetRelativeJEREC1Down", &vbfMass_JetRelativeJEREC1Down);
  tree->Branch("mjj_JetRelativeJEREC2Down", &vbfMass_JetRelativeJEREC2Down);
  tree->Branch("mjj_JetRelativeJERHFDown", &vbfMass_JetRelativeJERHFDown);
  tree->Branch("mjj_JetRelativePtBBDown", &vbfMass_JetRelativePtBBDown);
  tree->Branch("mjj_JetRelativePtEC1Down", &vbfMass_JetRelativePtEC1Down);
  tree->Branch("mjj_JetRelativePtEC2Down", &vbfMass_JetRelativePtEC2Down);
  tree->Branch("mjj_JetRelativePtHFDown", &vbfMass_JetRelativePtHFDown);
  tree->Branch("mjj_JetRelativeStatECDown", &vbfMass_JetRelativeStatECDown);
  tree->Branch("mjj_JetRelativeStatFSRDown", &vbfMass_JetRelativeStatFSRDown);
  tree->Branch("mjj_JetRelativeStatHFDown", &vbfMass_JetRelativeStatHFDown);
  tree->Branch("mjj_JetSinglePionECALDown", &vbfMass_JetSinglePionECALDown);
  tree->Branch("mjj_JetSinglePionHCALDown", &vbfMass_JetSinglePionHCALDown);
  tree->Branch("mjj_JetTimePtEtaDown", &vbfMass_JetTimePtEtaDown);
  tree->Branch("e_t_collinearmass", &m_coll);
  tree->Branch("e_t_collinearmass_UnclusteredEnUp", &m_coll_uesU);
  tree->Branch("e_t_collinearmass_UnclusteredEnDown", &m_coll_uesD);
  tree->Branch("e_t_collinearmass_JetEnUp", &m_coll_jesU);
  tree->Branch("e_t_collinearmass_JetEnDown", &m_coll_jesD);
  tree->Branch("e_t_collinearmass_TauEnUp", &m_coll_tesU);
  tree->Branch("e_t_collinearmass_TauEnDown", &m_coll_tesD);
  tree->Branch("njets_JESDown", &jetVeto20_JetEnDown);
  tree->Branch("njetspt20_JESDown", &jetVeto20_JetEnUp);
  tree->Branch("njets_JESUp", &jetVeto30_JetEnDown);
  tree->Branch("njetspt20_JESUp", &jetVeto20_JetEnUp);

  tree->Branch("metphi_EESDown", &type1_pfMet_shiftedPhi_ElectronEnDown);
  tree->Branch("metphi_EESUp", &type1_pfMet_shiftedPhi_ElectronEnUp);
  tree->Branch("metphi_JESDown", &type1_pfMet_shiftedPhi_JetEnDown);
  tree->Branch("metphi_JESUp", &type1_pfMet_shiftedPhi_JetEnUp);
  tree->Branch("metphi_JERDown", &type1_pfMet_shiftedPhi_JetResDown);
  tree->Branch("metphi_JERUp", &type1_pfMet_shiftedPhi_JetResUp);
  tree->Branch("metphi_MESDown", &type1_pfMet_shiftedPhi_MuonEnDown);
  tree->Branch("metphi_MESUp", &type1_pfMet_shiftedPhi_MuonEnUp);
  tree->Branch("metphi_PESDown", &type1_pfMet_shiftedPhi_PhotonEnDown);
  tree->Branch("metphi_PESUp", &type1_pfMet_shiftedPhi_PhotonEnUp);
  tree->Branch("metphi_TESDown", &type1_pfMet_shiftedPhi_TauEnDown);
  tree->Branch("metphi_TESUp", &type1_pfMet_shiftedPhi_TauEnUp);
  tree->Branch("metphi_UESDown", &type1_pfMet_shiftedPhi_UnclusteredEnDown);
  tree->Branch("metphi_UESUp", &type1_pfMet_shiftedPhi_UnclusteredEnUp);
  tree->Branch("met_EESDown", &type1_pfMet_shiftedPt_ElectronEnDown);
  tree->Branch("met_EESUp", &type1_pfMet_shiftedPt_ElectronEnUp);
  tree->Branch("met_JESDown", &type1_pfMet_shiftedPt_JetEnDown);
  tree->Branch("met_JESUp", &type1_pfMet_shiftedPt_JetEnUp);
  tree->Branch("met_JERDown", &type1_pfMet_shiftedPt_JetResDown);
  tree->Branch("met_JERUp", &type1_pfMet_shiftedPt_JetResUp);
  tree->Branch("met_MESDown", &type1_pfMet_shiftedPt_MuonEnDown);
  tree->Branch("met_MESUp", &type1_pfMet_shiftedPt_MuonEnUp);
  tree->Branch("met_PESDown", &type1_pfMet_shiftedPt_PhotonEnDown);
  tree->Branch("met_PESUp", &type1_pfMet_shiftedPt_PhotonEnUp);
  tree->Branch("met_TESDown", &type1_pfMet_shiftedPt_TauEnDown);
  tree->Branch("met_TESUp", &type1_pfMet_shiftedPt_TauEnUp);
  tree->Branch("met_UESDown", &type1_pfMet_shiftedPt_UnclusteredEnDown);
  tree->Branch("met_UESUp", &type1_pfMet_shiftedPt_UnclusteredEnUp);

  // created during skimming
  tree->Branch("pt_1", &pt_1, "pt_1/F");
  tree->Branch("phi_1", &phi_1, "phi_1/F");
  tree->Branch("eta_1", &eta_1, "eta_1/F");
  tree->Branch("m_1", &m_1, "m_1/F");
  tree->Branch("px_1", &px_1, "px_1/F");
  tree->Branch("py_1", &py_1, "py_1/F");
  tree->Branch("pz_1", &pz_1, "pz_1/F");
  tree->Branch("e_1", &e_1, "e_1/F");
  tree->Branch("pt_2", &pt_2, "pt_2/F");
  tree->Branch("phi_2", &phi_2, "phi_2/F");
  tree->Branch("eta_2", &eta_2, "eta_2/F");
  tree->Branch("px_2", &px_2, "px_2/F");
  tree->Branch("py_2", &py_2, "py_2/F");
  tree->Branch("pz_2", &pz_2, "pz_2/F");
  tree->Branch("m_2", &m_2, "m_2/F");
  tree->Branch("e_2", &e_2, "e_2/F");
  tree->Branch("mjj", &mjj, "mjj/F");
  tree->Branch("jdeta", &jdeta, "jdeta/F");
  tree->Branch("njetingap", &njetingap, "njetingap/I");
  tree->Branch("njetingap20", &njetingap20, "njetingap20/I");
  tree->Branch("dijetpt", &dijetpt, "dijetpt/F");
  tree->Branch("dijetphi", &dijetphi, "dijetphi/F");
  tree->Branch("hdijetphi", &hdijetphi, "hdijetphi/F");
  tree->Branch("visjeteta", &visjeteta, "visjeteta/F");
  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");

  // read input tree

  // straight from input tree
  original->SetBranchAddress("run", &run); 
  original->SetBranchAddress("lumi", &lumi);
  original->SetBranchAddress("rho", &rho);
  original->SetBranchAddress("metcov00", &metcov00);
  original->SetBranchAddress("metcov01", &metcov01);
  original->SetBranchAddress("metcov10", &metcov10);
  original->SetBranchAddress("metcov11", &metcov11);
  original->SetBranchAddress("NUP", &NUP);
  original->SetBranchAddress("tZTTGenDR", &tZTTGenDR);

  // read from tree and change name
  original->SetBranchAddress("GenWeight", &GenWeight);
  original->SetBranchAddress("genpX", &genpX);
  original->SetBranchAddress("genpY", &genpY);
  original->SetBranchAddress("vispX", &vispX);
  original->SetBranchAddress("vispY", &vispY);
  original->SetBranchAddress("genpT", &genpT);
  original->SetBranchAddress("genM", &genM);
  original->SetBranchAddress("type1_pfMetEt", &type1_pfMetEt);
  original->SetBranchAddress("type1_pfMetPhi", &type1_pfMetPhi);
  original->SetBranchAddress("metSig", &metSig);
  original->SetBranchAddress("jetVeto30", &jetVeto30);
  original->SetBranchAddress("bjetCISVVeto20Medium", &bjetCISVVeto20Medium);
  original->SetBranchAddress("jetVeto20", &jetVeto20);
  original->SetBranchAddress("eZTTGenMatching", &eZTTGenMatching);
  original->SetBranchAddress("tZTTGenMatching", &tZTTGenMatching);
  original->SetBranchAddress("tDecayMode", &tDecayMode);
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &tByIsolationMVArun2v1DBoldDMwLTraw);
  original->SetBranchAddress("tAgainstMuonTight3", &tAgainstMuonTight3);
  original->SetBranchAddress("tAgainstElectronVLooseMVA6", &tAgainstElectronVLooseMVA6);
  original->SetBranchAddress("tAgainstElectronLooseMVA6", &tAgainstElectronLooseMVA6);
  original->SetBranchAddress("tAgainstElectronMediumMVA6", &tAgainstElectronMediumMVA6);
  original->SetBranchAddress("tAgainstElectronVTightMVA6", &tAgainstElectronVTightMVA6);
  original->SetBranchAddress("tByLooseCombinedIsolationDeltaBetaCorr3Hits", &tByLooseCombinedIsolationDeltaBetaCorr3Hits);
  original->SetBranchAddress("tByMediumCombinedIsolationDeltaBetaCorr3Hits", &tByMediumCombinedIsolationDeltaBetaCorr3Hits);
  original->SetBranchAddress("tByTightCombinedIsolationDeltaBetaCorr3Hits", &tByTightCombinedIsolationDeltaBetaCorr3Hits);
  original->SetBranchAddress("tByCombinedIsolationDeltaBetaCorrRaw3Hits", &tByCombinedIsolationDeltaBetaCorrRaw3Hits);
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &tByIsolationMVArun2v1DBoldDMwLTraw);
  original->SetBranchAddress("tByIsolationMVArun2v1DBnewDMwLTraw", &tByIsolationMVArun2v1DBnewDMwLTraw);
  original->SetBranchAddress("tByLooseIsolationMVArun2v1DBoldDMwLT", &tByLooseIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByVTightIsolationMVArun2v1DBoldDMwLT", &tByVTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByVVTightIsolationMVArun2v1DBoldDMwLT", &tByVVTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tNeutralIsoPtSum", &tNeutralIsoPtSum);
  original->SetBranchAddress("tChargedIsoPtSum", &tChargedIsoPtSum);
  original->SetBranchAddress("tPuCorrPtSum", &tPuCorrPtSum);
  original->SetBranchAddress("tDecayModeFindingNewDMs", &tDecayModeFindingNewDMs);
  original->SetBranchAddress("j1eta", &j1eta);
  original->SetBranchAddress("j1csv", &j1csv);
  original->SetBranchAddress("j1phi", &j1phi);
  original->SetBranchAddress("j1pt", &j1pt);
  original->SetBranchAddress("j2csv", &j1csv);
  original->SetBranchAddress("j2eta", &j1eta);
  original->SetBranchAddress("j2phi", &j1phi);
  original->SetBranchAddress("j2pt", &j1pt);
  original->SetBranchAddress("jb1csv", &jb1csv);
  original->SetBranchAddress("jb1eta", &jb1eta);
  original->SetBranchAddress("jb1phi", &jb1phi);
  original->SetBranchAddress("jb1pt", &jb1pt);
  original->SetBranchAddress("jb2csv", &jb2csv);
  original->SetBranchAddress("jb2eta", &jb2eta);
  original->SetBranchAddress("jb2phi", &jb2phi);
  original->SetBranchAddress("jb2pt", &jb2pt);
  original->SetBranchAddress("nTruePU", &nTruePU);
  original->SetBranchAddress("numGenJets", &numGenJets);
  original->SetBranchAddress("nvtx", &nvtx);
  original->SetBranchAddress("topQuarkPt1", &topQuarkPt1);
  original->SetBranchAddress("topQuarkPt2", &topQuarkPt2);
  original->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
  original->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
  original->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
  original->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
  original->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
  original->SetBranchAddress("Flag_badCloneMuonFilter", &Flag_badCloneMuonFilter);
  original->SetBranchAddress("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter);
  original->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
  original->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
  original->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
  original->SetBranchAddress("e_t_MvaMetCovMatrix00", &e_t_MvaMetCovMatrix00);
  original->SetBranchAddress("e_t_MvaMetCovMatrix10", &e_t_MvaMetCovMatrix10);
  original->SetBranchAddress("e_t_MvaMetCovMatrix01", &e_t_MvaMetCovMatrix01);
  original->SetBranchAddress("e_t_MvaMetCovMatrix11", &e_t_MvaMetCovMatrix11);

  // used to construct something
  original->SetBranchAddress("eMass", &eMass);
  original->SetBranchAddress("tMass", &tMass);
  original->SetBranchAddress("vbfDeta", &vbfDeta);
  original->SetBranchAddress("vbfMass", &vbfMass);
  original->SetBranchAddress("vbfJetVeto20", &vbfJetVeto20);
  original->SetBranchAddress("vbfJetVeto30", &vbfJetVeto30);
  original->SetBranchAddress("eGenPdgId", &eGenPdgId);

  // Systematics
  original->SetBranchAddress("jetVeto30_JetAbsoluteFlavMapUp", &jetVeto30_JetAbsoluteFlavMapUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteMPFBiasUp", &jetVeto30_JetAbsoluteMPFBiasUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteScaleUp", &jetVeto30_JetAbsoluteScaleUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteStatUp", &jetVeto30_JetAbsoluteStatUp);
  original->SetBranchAddress("jetVeto30_JetEnUp", &jetVeto30_JetEnUp);
  original->SetBranchAddress("jetVeto30_JetFlavorQCDUp", &jetVeto30_JetFlavorQCDUp);
  original->SetBranchAddress("jetVeto30_JetFragmentationUp", &jetVeto30_JetFragmentationUp);
  original->SetBranchAddress("jetVeto30_JetPileUpDataMCUp", &jetVeto30_JetPileUpDataMCUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtBBUp", &jetVeto30_JetPileUpPtBBUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC1Up", &jetVeto30_JetPileUpPtEC1Up);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC2Up", &jetVeto30_JetPileUpPtEC2Up);
  original->SetBranchAddress("jetVeto30_JetPileUpPtHFUp", &jetVeto30_JetPileUpPtHFUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtRefUp", &jetVeto30_JetPileUpPtRefUp);
  original->SetBranchAddress("jetVeto30_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  original->SetBranchAddress("jetVeto30_JetRelativeFSRUp", &jetVeto30_JetRelativeFSRUp);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC1Up", &jetVeto30_JetRelativeJEREC1Up);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC2Up", &jetVeto30_JetRelativeJEREC2Up);
  original->SetBranchAddress("jetVeto30_JetRelativeJERHFUp", &jetVeto30_JetRelativeJERHFUp);
  original->SetBranchAddress("jetVeto30_JetRelativePtBBUp", &jetVeto30_JetRelativePtBBUp);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC1Up", &jetVeto30_JetRelativePtEC1Up);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC2Up", &jetVeto30_JetRelativePtEC2Up);
  original->SetBranchAddress("jetVeto30_JetRelativePtHFUp", &jetVeto30_JetRelativePtHFUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatECUp", &jetVeto30_JetRelativeStatECUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatFSRUp", &jetVeto30_JetRelativeStatFSRUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatHFUp", &jetVeto30_JetRelativeStatHFUp);
  original->SetBranchAddress("jetVeto30_JetSinglePionECALUp", &jetVeto30_JetSinglePionECALUp);
  original->SetBranchAddress("jetVeto30_JetSinglePionHCALUp", &jetVeto30_JetSinglePionHCALUp);
  original->SetBranchAddress("jetVeto30_JetTimePtEtaUp", &jetVeto30_JetTimePtEtaUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteFlavMapDown", &jetVeto30_JetAbsoluteFlavMapDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteMPFBiasDown", &jetVeto30_JetAbsoluteMPFBiasDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteScaleDown", &jetVeto30_JetAbsoluteScaleDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteStatDown", &jetVeto30_JetAbsoluteStatDown);
  original->SetBranchAddress("jetVeto30_JetEnDown", &jetVeto30_JetEnDown);
  original->SetBranchAddress("jetVeto30_JetFlavorQCDDown", &jetVeto30_JetFlavorQCDDown);
  original->SetBranchAddress("jetVeto30_JetFragmentationDown", &jetVeto30_JetFragmentationDown);
  original->SetBranchAddress("jetVeto30_JetPileUpDataMCDown", &jetVeto30_JetPileUpDataMCDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtBBDown", &jetVeto30_JetPileUpPtBBDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC1Down", &jetVeto30_JetPileUpPtEC1Down);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC2Down", &jetVeto30_JetPileUpPtEC2Down);
  original->SetBranchAddress("jetVeto30_JetPileUpPtHFDown", &jetVeto30_JetPileUpPtHFDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtRefDown", &jetVeto30_JetPileUpPtRefDown);
  original->SetBranchAddress("jetVeto30_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  original->SetBranchAddress("jetVeto30_JetRelativeFSRDown", &jetVeto30_JetRelativeFSRDown);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC1Down", &jetVeto30_JetRelativeJEREC1Down);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC2Down", &jetVeto30_JetRelativeJEREC2Down);
  original->SetBranchAddress("jetVeto30_JetRelativeJERHFDown", &jetVeto30_JetRelativeJERHFDown);
  original->SetBranchAddress("jetVeto30_JetRelativePtBBDown", &jetVeto30_JetRelativePtBBDown);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC1Down", &jetVeto30_JetRelativePtEC1Down);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC2Down", &jetVeto30_JetRelativePtEC2Down);
  original->SetBranchAddress("jetVeto30_JetRelativePtHFDown", &jetVeto30_JetRelativePtHFDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatECDown", &jetVeto30_JetRelativeStatECDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatFSRDown", &jetVeto30_JetRelativeStatFSRDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatHFDown", &jetVeto30_JetRelativeStatHFDown);
  original->SetBranchAddress("jetVeto30_JetSinglePionECALDown", &jetVeto30_JetSinglePionECALDown);
  original->SetBranchAddress("jetVeto30_JetSinglePionHCALDown", &jetVeto30_JetSinglePionHCALDown);
  original->SetBranchAddress("jetVeto30_JetTimePtEtaDown", &jetVeto30_JetTimePtEtaDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteFlavMapUp", &vbfMass_JetAbsoluteFlavMapUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteMPFBiasUp", &vbfMass_JetAbsoluteMPFBiasUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteScaleUp", &vbfMass_JetAbsoluteScaleUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteStatUp", &vbfMass_JetAbsoluteStatUp);
  original->SetBranchAddress("vbfMass_JetEnUp", &vbfMass_JetEnUp);
  original->SetBranchAddress("vbfMass_JetFlavorQCDUp", &vbfMass_JetFlavorQCDUp);
  original->SetBranchAddress("vbfMass_JetFragmentationUp", &vbfMass_JetFragmentationUp);
  original->SetBranchAddress("vbfMass_JetPileUpDataMCUp", &vbfMass_JetPileUpDataMCUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtBBUp", &vbfMass_JetPileUpPtBBUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC1Up", &vbfMass_JetPileUpPtEC1Up);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC2Up", &vbfMass_JetPileUpPtEC2Up);
  original->SetBranchAddress("vbfMass_JetPileUpPtHFUp", &vbfMass_JetPileUpPtHFUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtRefUp", &vbfMass_JetPileUpPtRefUp);
  original->SetBranchAddress("vbfMass_JetRelativeBalUp", &vbfMass_JetRelativeBalUp);
  original->SetBranchAddress("vbfMass_JetRelativeFSRUp", &vbfMass_JetRelativeFSRUp);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC1Up", &vbfMass_JetRelativeJEREC1Up);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC2Up", &vbfMass_JetRelativeJEREC2Up);
  original->SetBranchAddress("vbfMass_JetRelativeJERHFUp", &vbfMass_JetRelativeJERHFUp);
  original->SetBranchAddress("vbfMass_JetRelativePtBBUp", &vbfMass_JetRelativePtBBUp);
  original->SetBranchAddress("vbfMass_JetRelativePtEC1Up", &vbfMass_JetRelativePtEC1Up);
  original->SetBranchAddress("vbfMass_JetRelativePtEC2Up", &vbfMass_JetRelativePtEC2Up);
  original->SetBranchAddress("vbfMass_JetRelativePtHFUp", &vbfMass_JetRelativePtHFUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatECUp", &vbfMass_JetRelativeStatECUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatFSRUp", &vbfMass_JetRelativeStatFSRUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatHFUp", &vbfMass_JetRelativeStatHFUp);
  original->SetBranchAddress("vbfMass_JetSinglePionECALUp", &vbfMass_JetSinglePionECALUp);
  original->SetBranchAddress("vbfMass_JetSinglePionHCALUp", &vbfMass_JetSinglePionHCALUp);
  original->SetBranchAddress("vbfMass_JetTimePtEtaUp", &vbfMass_JetTimePtEtaUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteFlavMapDown", &vbfMass_JetAbsoluteFlavMapDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteMPFBiasDown", &vbfMass_JetAbsoluteMPFBiasDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteScaleDown", &vbfMass_JetAbsoluteScaleDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteStatDown", &vbfMass_JetAbsoluteStatDown);
  original->SetBranchAddress("vbfMass_JetEnDown", &vbfMass_JetEnDown);
  original->SetBranchAddress("vbfMass_JetFlavorQCDDown", &vbfMass_JetFlavorQCDDown);
  original->SetBranchAddress("vbfMass_JetFragmentationDown", &vbfMass_JetFragmentationDown);
  original->SetBranchAddress("vbfMass_JetPileUpDataMCDown", &vbfMass_JetPileUpDataMCDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtBBDown", &vbfMass_JetPileUpPtBBDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC1Down", &vbfMass_JetPileUpPtEC1Down);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC2Down", &vbfMass_JetPileUpPtEC2Down);
  original->SetBranchAddress("vbfMass_JetPileUpPtHFDown", &vbfMass_JetPileUpPtHFDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtRefDown", &vbfMass_JetPileUpPtRefDown);
  original->SetBranchAddress("vbfMass_JetRelativeBalDown", &vbfMass_JetRelativeBalDown);
  original->SetBranchAddress("vbfMass_JetRelativeFSRDown", &vbfMass_JetRelativeFSRDown);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC1Down", &vbfMass_JetRelativeJEREC1Down);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC2Down", &vbfMass_JetRelativeJEREC2Down);
  original->SetBranchAddress("vbfMass_JetRelativeJERHFDown", &vbfMass_JetRelativeJERHFDown);
  original->SetBranchAddress("vbfMass_JetRelativePtBBDown", &vbfMass_JetRelativePtBBDown);
  original->SetBranchAddress("vbfMass_JetRelativePtEC1Down", &vbfMass_JetRelativePtEC1Down);
  original->SetBranchAddress("vbfMass_JetRelativePtEC2Down", &vbfMass_JetRelativePtEC2Down);
  original->SetBranchAddress("vbfMass_JetRelativePtHFDown", &vbfMass_JetRelativePtHFDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatECDown", &vbfMass_JetRelativeStatECDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatFSRDown", &vbfMass_JetRelativeStatFSRDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatHFDown", &vbfMass_JetRelativeStatHFDown);
  original->SetBranchAddress("vbfMass_JetSinglePionECALDown", &vbfMass_JetSinglePionECALDown);
  original->SetBranchAddress("vbfMass_JetSinglePionHCALDown", &vbfMass_JetSinglePionHCALDown);
  original->SetBranchAddress("vbfMass_JetTimePtEtaDown", &vbfMass_JetTimePtEtaDown);
  original->SetBranchAddress("e_t_collinearmass", &m_coll);
  original->SetBranchAddress("e_t_collinearmass_UnclusteredEnUp", &m_coll_uesU);
  original->SetBranchAddress("e_t_collinearmass_UnclusteredEnDown", &m_coll_uesD);
  original->SetBranchAddress("e_t_collinearmass_JetEnUp", &m_coll_jesU);
  original->SetBranchAddress("e_t_collinearmass_JetEnDown", &m_coll_jesD);
  original->SetBranchAddress("e_t_collinearmass_TauEnUp", &m_coll_tesU);
  original->SetBranchAddress("e_t_collinearmass_TauEnDown", &m_coll_tesD);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_ElectronEnDown", &type1_pfMet_shiftedPhi_ElectronEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_ElectronEnUp", &type1_pfMet_shiftedPhi_ElectronEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown", &type1_pfMet_shiftedPhi_JetEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp", &type1_pfMet_shiftedPhi_JetEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetResDown", &type1_pfMet_shiftedPhi_JetResDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetResUp", &type1_pfMet_shiftedPhi_JetResUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_MuonEnDown", &type1_pfMet_shiftedPhi_MuonEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_MuonEnUp", &type1_pfMet_shiftedPhi_MuonEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_PhotonEnDown", &type1_pfMet_shiftedPhi_PhotonEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_PhotonEnUp", &type1_pfMet_shiftedPhi_PhotonEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_TauEnDown", &type1_pfMet_shiftedPhi_TauEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_TauEnUp", &type1_pfMet_shiftedPhi_TauEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown", &type1_pfMet_shiftedPhi_UnclusteredEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp", &type1_pfMet_shiftedPhi_UnclusteredEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_ElectronEnDown", &type1_pfMet_shiftedPt_ElectronEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_ElectronEnUp", &type1_pfMet_shiftedPt_ElectronEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown", &type1_pfMet_shiftedPt_JetEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp", &type1_pfMet_shiftedPt_JetEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetResDown", &type1_pfMet_shiftedPt_JetResDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetResUp", &type1_pfMet_shiftedPt_JetResUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_MuonEnDown", &type1_pfMet_shiftedPt_MuonEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_MuonEnUp", &type1_pfMet_shiftedPt_MuonEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_PhotonEnDown", &type1_pfMet_shiftedPt_PhotonEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_PhotonEnUp", &type1_pfMet_shiftedPt_PhotonEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_TauEnDown", &type1_pfMet_shiftedPt_TauEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_TauEnUp", &type1_pfMet_shiftedPt_TauEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown", &type1_pfMet_shiftedPt_UnclusteredEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp", &type1_pfMet_shiftedPt_UnclusteredEnUp);

  //    // not needed for sync
  //    original->SetBranchAddress("vbfMass_JetEnUp", &vbfMass_JetEnUp);
  //    original->SetBranchAddress("vbfMass_JetEnDown", &vbfMass_JetEnDown);
  //    original->SetBranchAddress("vbfJetVeto20_JetEnDown", &vbfJetVeto20_JetEnDown);
  //    original->SetBranchAddress("vbfJetVeto20_JetEnUp", &vbfJetVeto20_JetEnUp);
  //    original->SetBranchAddress("vbfJetVeto30_JetEnDown", &vbfJetVeto30_JetEnDown);
  //    original->SetBranchAddress("vbfJetVeto30_JetEnUp", &vbfJetVeto30_JetEnUp);
  //    original->SetBranchAddress("vbfDeta_JetEnDown", &vbfDeta_JetEnDown);
  //    original->SetBranchAddress("vbfDeta_JetEnUp", &vbfDeta_JetEnUp);
  //    original->SetBranchAddress("tauVetoPt20Loose3HitsVtx", &extratau_veto);
  //    original->SetBranchAddress("isZtautau", &isZtt);
  //    original->SetBranchAddress("e_t_MvaMet", &mvaMet);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix00", &mvaMetcov00);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix01", &mvaMetcov01);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix10", &mvaMetcov10);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix11", &mvaMetcov11);
  //    original->SetBranchAddress("e_t_MvaMetPhi", &mvaMetphi);
  //    original->SetBranchAddress("e_t_DPhi", &dphi_12);
  //    original->SetBranchAddress("eDPhiToPfMet_type1", &dphi_emet);
  //    original->SetBranchAddress("tDPhiToPfMet_type1", &dphi_taumet);
}
