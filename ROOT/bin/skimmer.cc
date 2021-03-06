// general includes
#include <dirent.h>
#include <sys/types.h>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

// user includes
#include "ETau_Skimmer/ROOT/src/etau_tree.h"
#include "ETau_Skimmer/ROOT/src/CLParser.h"

static unsigned events(0);
int main(int argc, char *argv[]) {

  CLParser parser(argc, argv);
  std::string dir_name = parser.Option("-d");
  std::string job_type = parser.Option("-j");
  std::string inrecoil = parser.Option("-r");
  std::string ifile = parser.Option("-i");
  std::string ofile = parser.Option("-o");

  // recoil corrections
  int recoil(0);
  if (inrecoil.find("W") != std::string::npos) {
    recoil = 1;
  } else if (inrecoil.find("Z") != std::string::npos) {
    recoil = 2;
  }

  // set flag for MC or data
  bool isMC(true);
  if (job_type == "data")
    isMC = false;

  bool isEmbed(false);
  if (job_type == "embed")
    isEmbed = true;

  std::cout << job_type << " " << isMC << " " << isEmbed << std::endl;

  RecoilCorrector recoilPFMetCorrector("SMH_ettau/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root");

  TH1F *nevents = new TH1F("nevents", "N(events)", 2, 0.5, 2.5);
  TH1F* cutflow = new TH1F("cutflow", "cutflow", 10, 0.5, 10.5);

    auto open_file = new TFile(ifile.c_str(), "READ");
    auto ntuple = (TTree *)open_file->Get("et/final/Ntuple");
    auto evt_count = (TH1F *)open_file->Get("et/eventCount")->Clone();
    auto wt_count = (TH1F *)open_file->Get("et/summedWeights")->Clone();

    nevents->SetBinContent(1, evt_count->Integral());
    nevents->SetBinContent(2, wt_count->Integral());

    auto fout = new TFile(ofile.c_str(), "RECREATE");

    TTree *newtree = new TTree("etau_tree", "etau_tree");
    etau_tree *skimmer = new etau_tree(ntuple, newtree, isMC, isEmbed, recoil);
    skimmer->do_skimming(cutflow);
    auto skimmed_tree = skimmer->fill_tree(recoilPFMetCorrector);
    events += skimmed_tree->GetEntries();

    open_file->Close();
    fout->cd();
    nevents->Write();
    cutflow->Write();
    skimmed_tree->Write();
    fout->Close();

  std::cout << std::endl;
  std::cout << "\n" << events << " events saved in the tree.\n" << std::endl;
  return 0;
}


