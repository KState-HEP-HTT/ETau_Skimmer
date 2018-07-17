#include <sys/types.h>
#include <dirent.h>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"

#include "etau_tree.h"
#include "util.h"

void read_directory(const std::string &name, std::vector<std::string> &v) {
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != 0) {
    if (static_cast<std::string>(dp->d_name).find("root") != std::string::npos)
      v.push_back(dp->d_name);
  }
  closedir(dirp);
}
static unsigned events(0);
int main(int argc, char* argv[]) {


  std::vector<std::string> all_files;
  std::string dir_name = "test";
  if (argc > 1)
    dir_name = argv[1];

  std::string in = samples[dir_name];
  const int dir_err = system(("mkdir -p /nfs_scratch/tmitchel/mc2016_ntuples_July13_skim/"+dir_name).c_str());
  read_directory(in, all_files);

  TH1F* nevents = new TH1F("nevents", "N(events)", 2, 0.5, 2.5);

  // I/O files and the ntuple
  std::cout << "Begin loading files..." << std::endl;
  unsigned i = 0;
  for (auto& ifile : all_files) {
    std::cout << "Loading file: " << i+1 << " out of " << all_files.size() << " files.\r" << std::flush;
    i += 1;
    auto ntuple = new TChain("et/final/Ntuple");
    auto open_file = new TFile((in+"/"+ifile).c_str(), "READ");
    auto evt_count = (TH1F*)open_file->Get("et/eventCount")->Clone();
    auto wt_count = (TH1F*)open_file->Get("et/summedWeights")->Clone();

    nevents->SetBinContent(1, evt_count->Integral());
    nevents->SetBinContent(2, wt_count->Integral());

    open_file->Close();
    ntuple->Add((in+"/"+ifile).c_str());
    std::string suffix = "/nfs_scratch/tmitchel/mc2016_ntuples_July13_skim/"+dir_name+"/Skim_";
    auto fout = new TFile((suffix+ifile).c_str(), "RECREATE");

    TTree* newtree = new TTree("etau_tree","etau_tree");
    etau_tree* skimmer = new etau_tree(ntuple, newtree);
    skimmer->do_skimming();
    auto skimmed_tree = skimmer->fill_tree();
    events += skimmed_tree->GetEntries();

    fout->cd();
    nevents->Write();
    skimmed_tree->Write();
    fout->Close();
  }
  std::cout << std::endl;
  std::cout << events << "events saved in the tree." << std::endl;

  return 0;
}
