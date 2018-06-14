#include <sys/types.h>
// #include <sys/stat.h>
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

int main(int argc, char* argv[]) {


  std::vector<std::string> all_files;
  std::string dir_name = "test";
  if (argc > 1)
    dir_name = argv[1];

  std::string in = samples[dir_name];
  auto outprefix = "/hdfs/store/user/tmitchel/skims_2016/"+dir_name;
  const int dir_err = system(("gsido mkdir "+outprefix).c_str());
  // const int dir_err = mkdir(outprefix, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  read_directory(in, all_files);

  TH1F* nevents = new TH1F("nevents", "N(events)", 1, 0.5, 1.5);
  TH1F* nweights = new TH1F("nweights", "N(weights)", 1, 0.5, 1.5);

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

    nevents->SetBinContent(1, nevents->GetBinContent(1)+evt_count->Integral());
    nweights->SetBinContent(1, nweights->GetBinContent(1)+wt_count->Integral());

    open_file->Close();
    ntuple->Add((in+"/"+ifile).c_str());
    std::string suffix = "Skim_";
    auto fout = new TFile((outprefix+suffix+ifile).c_str(), "RECREATE");

    TTree* newtree = new TTree("skim","skim");
    etau_tree* skimmer = new etau_tree(ntuple, newtree);
    auto skimmed_tree = skimmer->do_skimming();

    fout->cd();
    nevents->Write();
    nweights->Write();
    skimmed_tree->Write();
    fout->Close();
  }

  return 0;
}
