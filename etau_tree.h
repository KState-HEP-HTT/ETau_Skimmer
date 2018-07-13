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
  Float_t genweight, genpX, genpY, vispX, vispY, genpT, genM, met, metphi, metSig, metcov00, metcov01, metcov10, metcov11, met_EESDown, met_EESUp, met_JESUp, met_JESDown;
  Float_t met_MESDown, met_MESUp, met_PESUp, met_PESDown, met_TESUp, met_TESDown, met_UESUp, met_UESDown, met_JERDown, met_JERUp, metphi_EESDown, metphi_EESUp, metphi_JESUp;
  Float_t metphi_JESDown, metphi_MESDown, metphi_MESUp, metphi_PESUp, metphi_PESDown, metphi_TESUp, metphi_TESDown, metphi_UESUp, metphi_UESDown, metphi_JERDown, metphi_JERUp;
  Float_t mvaMet, mvaMetcov00, mvaMetcov11, mvaMetcov10, mvaMetcov01, mvaMetphi, dphi_12, dphi_emet, dphi_taumet, passEle25, passEle27, filterEle25;
  Float_t pt_top1, pt_top2, NUP, njets, nbtag, njetspt20, gen_match_1, gen_match_2;
  Float_t m_vis, l2_decayMode, dZ_1, d0_1, iso_1, q_1, dZ_2, d0_2, iso_2, q_2, m_coll, m_coll_uesU, m_coll_uesD, m_coll_jesU, m_coll_jesD, m_coll_tesU, m_coll_tesD;
  Float_t againstMuonTight3_2, againstMuonLoose3_2, againstElectronVLooseMVA6_2, againstElectronLooseMVA6_2, againstElectronMediumMVA6_2, againstElectronTightMVA6_2;
  Float_t againstElectronVTightMVA6_2, byLooseCombinedIsolationDeltaBetaCorr3Hits_2, byMediumCombinedIsolationDeltaBetaCorr3Hits_2, byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  Float_t byCombinedIsolationDeltaBetaCorrRaw3Hits_2, byIsolationMVA3oldDMwLTraw_2, byIsolationMVA3newDMwLTraw_2, byVLooseIsolationMVArun2v1DBoldDMwLT_2, byLooseIsolationMVArun2v1DBoldDMwLT_2;
  Float_t byMediumIsolationMVArun2v1DBoldDMwLT_2, byTightIsolationMVArun2v1DBoldDMwLT_2, byVTightIsolationMVArun2v1DBoldDMwLT_2, byVVTightIsolationMVArun2v1DBoldDMwLT_2;
  Float_t neutralIsoPtSum_2, chargedIsoPtSum_2, puCorrPtSum_2, decayModeFinding_2, decayModeFindingNewDMs_2, jpt_1, jpt_2, jeta_1, jeta_2, jphi_1, jphi_2, jcsv_1, jcsv_2, bpt_1, bpt_2;
  Float_t beta_1, beta_2, bphi_1, bphi_2, bcsv_1, bcsv_2, npu, npv, rho, extratau_veto, isZtt, idisoweight_2;

  Int_t run, lumi, njetingap, njetingap20;
  Int_t njetingap_JESUp, njetingap20_JESUp, njetingap_JESDown, njetingap20_JESDown, njets_JESUp, njetspt20_JESUp, njets_JESDown, njetspt20_JESDown;

  // temporary storage variables
  Float_t eVetoZTTp001dxyzR0, muVetoZTTp001dxyzR0, dielectronVeto, vbfMass_JetEnUp, vbfMass_JetEnDown;

  // new variables to store in new tree
  Float_t met_px, met_py, extraelec_veto, dilepton_veto, m_1, pt_1, eta_1, phi_1, e_1, px_1, py_1, pz_1, m_2, pt_2, eta_2, phi_2, e_2, px_2, py_2, pz_2, dijetphi, hdijetphi, visjeteta, isZet;
  Float_t jdeta, mjj, dijetpt, mjj_JESUp, jdeta_JESUp, mjj_JESDown, jdeta_JESDown;
  Float_t gen_Higgs_pt, gen_Higgs_mass, weight;

  // forgotten
  Float_t byIsolationMVA3oldDMwoLTraw_2, trigweight_2, byIsolationMVA3newDMwoLTraw_2, filterEle27, ePt, eMass, ePhi, eEta, tPhi, tEta, tMass, tPt, numGenJets, vbfDeta_JetEnDown, vbfDeta_JetEnUp;
  Float_t vbfMass, vbfJetVeto20, vbfJetVeto20_JetEnDown, vbfJetVeto20_JetEnUp, vbfJetVeto30, vbfJetVeto30_JetEnDown, vbfJetVeto30_JetEnUp, eGenPdgId, vbfDeta, extramuon_veto;
  Float_t eMVANonTrigWP80, ePassesConversionVeto, eMissingHits, e_t_DR;

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
  original->SetBranchAddress("ePVDZ", &dZ_1);
  original->SetBranchAddress("ePVDXY", &d0_1);
  original->SetBranchAddress("tPt", &tPt);
  original->SetBranchAddress("tEta", &tEta);
  original->SetBranchAddress("tPhi", &tPhi);
  original->SetBranchAddress("tDecayModeFinding", &decayModeFinding_2);
  original->SetBranchAddress("tPVDXY", &d0_2);
  original->SetBranchAddress("tPVDZ", &dZ_2);
  original->SetBranchAddress("tByVLooseIsolationMVArun2v1DBoldDMwLT", &byVLooseIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("tByTightIsolationMVArun2v1DBoldDMwLT", &byTightIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("dielectronVeto", &dielectronVeto);
  original->SetBranchAddress("tAgainstMuonLoose3", &againstMuonLoose3_2);
  original->SetBranchAddress("tAgainstElectronTightMVA6", &againstElectronTightMVA6_2);
  original->SetBranchAddress("eCharge", &q_1);
  original->SetBranchAddress("tCharge", &q_2);
  original->SetBranchAddress("eIsoDB03", &iso_1);
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
    if (!passEle25 || !filterEle25) // apply trigger HLT Ele25 eta2p1 WPTight Gsf with matching
      continue;

    if (ePt < 25 || fabs(eEta) > 2.5 || fabs(dZ_1) > 0.2 || fabs(d0_1) > 0.045 || !eMVANonTrigWP80 || !ePassesConversionVeto || eMissingHits > 1 || eVetoZTTp001dxyzR0 > 1) // electron selection
      continue;

    if (tPt < 19 || fabs(tEta) > 2.3 || fabs(dZ_2) > 0.2 || !byMediumIsolationMVArun2v1DBoldDMwLT_2 || !decayModeFinding_2 || fabs(q_2) > 1) // tau selection
      continue;

    if (muVetoZTTp001dxyzR0 > 0 || dielectronVeto > 0)
      continue;

    double dR = sqrt( pow(eEta - tEta, 2) + pow(ePhi - tPhi, 2) ); // pair selection
    if (dR < 0.5)
      continue;

    // implement new sorting per 
	  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
	  if (evt_now != evt_before) { // new event, save the tau candidates
      //   since it is new event, do we have the best entry to save? If yes, save it!
      if ( best_evt > -1  )
        good_events.push_back(best_evt);
      
      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      eleCandidate  = std::make_pair(ePt,  iso_1);
      tauCandidate  = std::make_pair(tPt,  byVLooseIsolationMVArun2v1DBoldDMwLT_2);
    } 
    else { // not a new event

      std::pair<float, float> currEleCandidate(ePt, iso_1);
      std::pair<float, float> currTauCandidate(tPt, byVLooseIsolationMVArun2v1DBoldDMwLT_2);

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

    met_px = met*cos(metphi);
    met_py = met*sin(metphi);

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
    if (njetspt20 > 0 && jpt_1>0)
       jet1.SetPtEtaPhiM(jpt_1,jeta_1,jphi_1,0);
    TLorentzVector jet2;
    if (njetspt20 > 1 && jpt_2>0)
       jet2.SetPtEtaPhiM(jpt_2,jeta_2,jphi_2,0);
    TLorentzVector dijet=jet1+jet2;

    if (njetspt20 > 1){
       jdeta = vbfDeta;
       dijetphi = dijet.Phi();
       hdijetphi = h.DeltaPhi(dijet);
       visjeteta = h.Eta()-dijet.Eta();
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

  // from input tree, but change name
  tree->Branch("q_1", &q_1, "q_1/F");
  tree->Branch("d0_1", &d0_1, "d0_1/F");
  tree->Branch("dZ_1", &dZ_1, "dZ_1/F");
  tree->Branch("iso_1", &iso_1, "iso_1/F");
  tree->Branch("q_2", &q_2, "q_2/F");
  tree->Branch("d0_2", &d0_2, "d0_2/F");
  tree->Branch("dZ_2", &dZ_2, "dZ_2/F");
  tree->Branch("iso_2", &iso_2, "iso_2/F");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/F");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec/F");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon/F");
  tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
  tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
  tree->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, "againstElectronLooseMVA6_2/F");
  tree->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, "againstElectronMediumMVA6_2/F");
  tree->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, "againstElectronVTightMVA6_2/F");
  tree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/F");
  tree->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/F");
  tree->Branch("npv", &npv, "npv/F");
  tree->Branch("npu", &npu, "npu/F");
  tree->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2, "byLooseCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2, "byMediumCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2, "byTightCombinedIsolationDeltaBetaCorr3Hits_2/F");
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
  tree->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, "byIsolationMVA3oldDMwoLTraw_2/F");
  tree->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, "byIsolationMVA3oldDMwLTraw_2/F");
  tree->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, "byIsolationMVA3newDMwoLTraw_2/F");
  tree->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, "byIsolationMVA3newDMwLTraw_2/F");
  tree->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2, "chargedIsoPtSum_2/F");
  tree->Branch("decayModeFinding_2", &decayModeFinding_2, "decayModeFinding_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2, "decayModeFindingNewDMs_2/F");
  tree->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2, "neutralIsoPtSum_2/F");
  tree->Branch("puCorrPtSum_2", &puCorrPtSum_2, "puCorrPtSum_2/F");
  tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
  tree->Branch("idisoweight_2", &idisoweight_2, "idisoweight_2/F");
  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
  tree->Branch("jpt_1", &jpt_1, "jpt_1/F");
  tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
  tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
  tree->Branch("jcsv_1", &jcsv_1, "jcsv_1/F");
  tree->Branch("jpt_2", &jpt_2, "jpt_2/F");
  tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
  tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
  tree->Branch("jcsv_2", &jcsv_2, "jcsv_2/F");
  tree->Branch("bpt_1", &bpt_1, "bpt_1/F");
  tree->Branch("beta_1", &beta_1, "beta_1/F");
  tree->Branch("bphi_1", &bphi_1, "bphi_1/F");
  tree->Branch("bcsv_1", &bcsv_1, "bcsv_1/F");
  tree->Branch("bpt_2", &bpt_2, "bpt_2/F");
  tree->Branch("beta_2", &beta_2, "beta_2/F");
  tree->Branch("bphi_2", &bphi_2, "bphi_2/F");
  tree->Branch("bcsv_2", &bcsv_2, "bcsv_2/F");
  tree->Branch("MVANonTrigWP80_1", &eMVANonTrigWP80, "MVANonTrigWP80_1/F");
  tree->Branch("passEle25", &passEle25, "passEle25/F");
  tree->Branch("filterEle25", &filterEle25, "filterEle25/F");
  tree->Branch("l2_decayMode", &l2_decayMode, "l2_decayMode/F");

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




  //    // others not for sync
  //    tree->Branch("isZtt", &isZtt, "isZtt/O");
  //    tree->Branch("isZet", &isZet, "isZet/O");
  //    tree->Branch("vispX", &vispX, "vispX/F");
  //    tree->Branch("vispY", &vispY, "vispY/F");
  //    tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2", &byMediumIsolationMVArun2v1DBoldDMwLT_2, "byMediumIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &byTightIsolationMVArun2v1DBoldDMwLT_2, "byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  //    tree->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2", &byVTightIsolationMVArun2v1DBoldDMwLT_2, "byVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  //    tree->Branch("byVVTightIsolationMVArun2v1DBoldDMwLT_2", &byVVTightIsolationMVArun2v1DBoldDMwLT_2, "byVVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  //    tree->Branch("metSig", &metSig, "metSig/F");
  //    tree->Branch("mvaMet", &mvaMet, "mvaMet/F");
  //    tree->Branch("mvaMetcov00", &mvaMetcov00, "mvaMetcov00/F");
  //    tree->Branch("mvaMetcov10", &mvaMetcov10, "mvaMetcov10/F");
  //    tree->Branch("mvaMetcov11", &mvaMetcov11, "mvaMetcov11/F");
  //    tree->Branch("mvaMetcov01", &mvaMetcov01, "mvaMetcov01/F");
  //    tree->Branch("mvaMetphi", &mvaMetphi, "mvaMetphi/F");
  //    tree->Branch("mjj_JESUp", &mjj_JESUp, "mjj_JESUp/F");
  //    tree->Branch("jdeta_JESUp", &jdeta_JESUp, "jdeta_JESUp/F");
  //    tree->Branch("njetingap_JESUp", &njetingap_JESUp, "njetingap_JESUp/I");
  //    tree->Branch("njetingap20_JESUp", &njetingap20_JESUp, "njetingap20_JESUp/I");
  //    tree->Branch("mjj_JESDown", &mjj_JESDown, "mjj_JESDown/F");
  //    tree->Branch("jdeta_JESDown", &jdeta_JESDown, "jdeta_JESDown/F");
  //    tree->Branch("njetingap_JESDown", &njetingap_JESDown, "njetingap_JESDown/I");
  //    tree->Branch("njetingap20_JESDown", &njetingap20_JESDown, "njetingap20_JESDown/I");
  //    tree->Branch("njets_JESUp", &njets_JESUp, "njets_JESUp/I");
  //    tree->Branch("njetspt20_JESUp", &njetspt20_JESUp, "njetspt20_JESUp/I");
  //    tree->Branch("njets_JESDown", &njets_JESDown, "njets_JESDown/I");
  //    tree->Branch("njetspt20_JESDown", &njetspt20_JESDown, "njetspt20_JESDown/I");
  //    tree->Branch("pt_top1", &pt_top1, "pt_top1/F");
  //    tree->Branch("pt_top2", &pt_top2, "pt_top2/F");
  //    tree->Branch("genweight", &genweight, "genweight/F");
  //    tree->Branch("gen_Higgs_pt", &gen_Higgs_pt, "gen_Higgs_pt/F");
  //    tree->Branch("gen_Higgs_mass", &gen_Higgs_mass, "gen_Higgs_mass/F");


  // read input tree

  // already read for skimming
  // original->SetBranchAddress("evt", &evt);
  // original->SetBranchAddress("ePt", &ePt);
  // original->SetBranchAddress("eEta", &eEta);
  // original->SetBranchAddress("ePassesConversionVeto", &ePassesConversionVeto);
  // original->SetBranchAddress("ePVDZ", &dZ_1);
  // original->SetBranchAddress("ePVDXY", &d0_1);
  // original->SetBranchAddress("tPt", &tPt);
  // original->SetBranchAddress("tEta", &tEta);
  // original->SetBranchAddress("tDecayModeFinding", &decayModeFinding_2);
  // original->SetBranchAddress("tPVDXY", &dZ_2);
  // original->SetBranchAddress("tPVDZ", &d0_2);
  // original->SetBranchAddress("tByTightIsolationMVArun2v1DBoldDMwLT", &byTightIsolationMVArun2v1DBoldDMwLT_2);
  // original->SetBranchAddress("dielectronVeto", &dielectronVeto);
  // original->SetBranchAddress("tAgainstMuonLoose3", &againstMuonLoose3_2);
  // original->SetBranchAddress("tAgainstElectronTightMVA6", &againstElectronTightMVA6_2);
  // original->SetBranchAddress("eCharge", &q_1);
  // original->SetBranchAddress("tCharge", &q_2);
  // original->SetBranchAddress("eIsoDB03", &iso_1);
  // original->SetBranchAddress("ePhi", &ePhi);
  // original->SetBranchAddress("tPhi", &tPhi);

  // straight from input tree
  original->SetBranchAddress("run", &run); 
  original->SetBranchAddress("lumi", &lumi);
  original->SetBranchAddress("rho", &rho);
  original->SetBranchAddress("metcov00", &metcov00);
  original->SetBranchAddress("metcov01", &metcov01);
  original->SetBranchAddress("metcov10", &metcov10);
  original->SetBranchAddress("metcov11", &metcov11);
  original->SetBranchAddress("NUP", &NUP);

  // read from tree and change name
  original->SetBranchAddress("GenWeight", &genweight);
  original->SetBranchAddress("genpX", &genpX);
  original->SetBranchAddress("genpY", &genpY);
  original->SetBranchAddress("vispX", &vispX);
  original->SetBranchAddress("vispY", &vispY);
  original->SetBranchAddress("genpT", &genpT);
  original->SetBranchAddress("genM", &genM);
  original->SetBranchAddress("type1_pfMetEt", &met);
  original->SetBranchAddress("type1_pfMetPhi", &metphi);
  original->SetBranchAddress("metSig", &metSig);
  original->SetBranchAddress("GenWeight", &weight);
  original->SetBranchAddress("jetVeto30", &njets);
  original->SetBranchAddress("bjetCISVVeto20Medium", &nbtag);
  original->SetBranchAddress("jetVeto20", &njetspt20);
  original->SetBranchAddress("eZTTGenMatching", &gen_match_1);
  original->SetBranchAddress("tZTTGenMatching", &gen_match_2);
  original->SetBranchAddress("tDecayMode", &l2_decayMode);
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &iso_2);
  original->SetBranchAddress("tAgainstMuonTight3", &againstMuonTight3_2);
  original->SetBranchAddress("tAgainstElectronVLooseMVA6", &againstElectronVLooseMVA6_2);
  original->SetBranchAddress("tAgainstElectronLooseMVA6", &againstElectronLooseMVA6_2);
  original->SetBranchAddress("tAgainstElectronMediumMVA6", &againstElectronMediumMVA6_2);
  original->SetBranchAddress("tAgainstElectronVTightMVA6", &againstElectronVTightMVA6_2);
  original->SetBranchAddress("tByLooseCombinedIsolationDeltaBetaCorr3Hits", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
  original->SetBranchAddress("tByMediumCombinedIsolationDeltaBetaCorr3Hits", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
  original->SetBranchAddress("tByTightCombinedIsolationDeltaBetaCorr3Hits", &byTightCombinedIsolationDeltaBetaCorr3Hits_2);
  original->SetBranchAddress("tByCombinedIsolationDeltaBetaCorrRaw3Hits", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &byIsolationMVA3oldDMwLTraw_2);
  original->SetBranchAddress("tByIsolationMVArun2v1DBnewDMwLTraw", &byIsolationMVA3newDMwLTraw_2);
  original->SetBranchAddress("tByLooseIsolationMVArun2v1DBoldDMwLT", &byLooseIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("tByMediumIsolationMVArun2v1DBoldDMwLT", &byMediumIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("tByVTightIsolationMVArun2v1DBoldDMwLT", &byVTightIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("tByVVTightIsolationMVArun2v1DBoldDMwLT", &byVVTightIsolationMVArun2v1DBoldDMwLT_2);
  original->SetBranchAddress("tNeutralIsoPtSum", &neutralIsoPtSum_2);
  original->SetBranchAddress("tChargedIsoPtSum", &chargedIsoPtSum_2);
  original->SetBranchAddress("tPuCorrPtSum", &puCorrPtSum_2);
  original->SetBranchAddress("tDecayModeFindingNewDMs", &decayModeFindingNewDMs_2);
  original->SetBranchAddress("j1eta", &jeta_1);
  original->SetBranchAddress("j1csv", &jcsv_1);
  original->SetBranchAddress("j1phi", &jphi_1);
  original->SetBranchAddress("j1pt", &jpt_1);
  original->SetBranchAddress("j2csv", &jcsv_2);
  original->SetBranchAddress("j2eta", &jeta_2);
  original->SetBranchAddress("j2phi", &jphi_2);
  original->SetBranchAddress("j2pt", &jpt_2);
  original->SetBranchAddress("jb1csv", &bcsv_1);
  original->SetBranchAddress("jb1eta", &beta_1);
  original->SetBranchAddress("jb1phi", &bphi_1);
  original->SetBranchAddress("jb1pt", &bpt_1);
  original->SetBranchAddress("jb2csv", &bcsv_2);
  original->SetBranchAddress("jb2eta", &beta_2);
  original->SetBranchAddress("jb2phi", &bphi_2);
  original->SetBranchAddress("jb2pt", &bpt_2);
  original->SetBranchAddress("nTruePU", &npu);
  original->SetBranchAddress("numGenJets", &numGenJets);
  original->SetBranchAddress("nvtx", &npv);
  original->SetBranchAddress("eMVANonTrigWP80", &eMVANonTrigWP80);
  original->SetBranchAddress("singleE25eta2p1TightPass", &passEle25);
  original->SetBranchAddress("eMatchesEle25TightFilter", &filterEle25);
  

  // used to construct something
  original->SetBranchAddress("eVetoZTTp001dxyzR0", &eVetoZTTp001dxyzR0);
  original->SetBranchAddress("muVetoZTTp001dxyzR0", &muVetoZTTp001dxyzR0);
  original->SetBranchAddress("eMass", &eMass);
  original->SetBranchAddress("tMass", &tMass);
  original->SetBranchAddress("vbfDeta", &vbfDeta);
  original->SetBranchAddress("vbfMass", &vbfMass);
  original->SetBranchAddress("vbfJetVeto20", &vbfJetVeto20);
  original->SetBranchAddress("vbfJetVeto30", &vbfJetVeto30);
  original->SetBranchAddress("eGenPdgId", &eGenPdgId);
  original->SetBranchAddress("eMissingHits", &eMissingHits);
  // original->SetBranchAddress("e_t_DR", &e_t_DR);

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
  //    original->SetBranchAddress("e_t_collinearmass", &m_coll);
  //    original->SetBranchAddress("e_t_collinearmass_UnclusteredEnUp", &m_coll_uesU);
  //    original->SetBranchAddress("e_t_collinearmass_UnclusteredEnDown", &m_coll_uesD);
  //    original->SetBranchAddress("e_t_collinearmass_JetEnUp", &m_coll_jesU);
  //    original->SetBranchAddress("e_t_collinearmass_JetEnDown", &m_coll_jesD);
  //    original->SetBranchAddress("e_t_collinearmass_TauEnUp", &m_coll_tesU);
  //    original->SetBranchAddress("e_t_collinearmass_TauEnDown", &m_coll_tesD);
  //    original->SetBranchAddress("jetVeto20_JetEnDown", &njets_JESDown);
  //    original->SetBranchAddress("jetVeto20_JetEnUp", &njetspt20_JESDown);
  //    original->SetBranchAddress("jetVeto30_JetEnDown", &njets_JESUp);
  //    original->SetBranchAddress("jetVeto30_JetEnUp", &njetspt20_JESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_ElectronEnDown", &metphi_EESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_ElectronEnUp", &metphi_EESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown", &metphi_JESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp", &metphi_JESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_JetResDown", &metphi_JERDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_JetResUp", &metphi_JERUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_MuonEnDown", &metphi_MESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_MuonEnUp", &metphi_MESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_PhotonEnDown", &metphi_PESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_PhotonEnUp", &metphi_PESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_TauEnDown", &metphi_TESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_TauEnUp", &metphi_TESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown", &metphi_UESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp", &metphi_UESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_ElectronEnDown", &met_EESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_ElectronEnUp", &met_EESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown", &met_JESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp", &met_JESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_JetResDown", &met_JERDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_JetResUp", &met_JERUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_MuonEnDown", &met_MESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_MuonEnUp", &met_MESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_PhotonEnDown", &met_PESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_PhotonEnUp", &met_PESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_TauEnDown", &met_TESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_TauEnUp", &met_TESUp);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown", &met_UESDown);
  //    original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp", &met_UESUp);
  //    original->SetBranchAddress("e_t_MvaMet", &mvaMet);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix00", &mvaMetcov00);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix01", &mvaMetcov01);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix10", &mvaMetcov10);
  //    original->SetBranchAddress("e_t_MvaMetCovMatrix11", &mvaMetcov11);
  //    original->SetBranchAddress("e_t_MvaMetPhi", &mvaMetphi);
  //    original->SetBranchAddress("genpT", &gen_Higgs_pt);
  //    original->SetBranchAddress("genM", &gen_Higgs_mass);
  //    original->SetBranchAddress("e_t_DPhi", &dphi_12);
  //    original->SetBranchAddress("eDPhiToPfMet_type1", &dphi_emet);
  //    original->SetBranchAddress("tDPhiToPfMet_type1", &dphi_taumet);
  //    original->SetBranchAddress("topQuarkPt1", &pt_top1);
  //    original->SetBranchAddress("topQuarkPt2", &pt_top2);
}
