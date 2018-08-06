#include <sys/types.h>
#include <dirent.h>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"

#include "etau_tree.h"
#include "util.h"

void read_directory(const std::string &name, std::vector<std::string> &v, std::vector<std::string> &dir) {
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != 0) {
    if (static_cast<std::string>(dp->d_name).find("root") != std::string::npos) {
      v.push_back(dp->d_name);
      dir.push_back(name);
    }
  }
  closedir(dirp);
}

void add_pref(std::vector<std::string> &files, std::string prefix) {
  for (auto& ifile : files) {
    ifile = prefix + ifile;
  }
}

static unsigned events(0);
int main(int argc, char* argv[]) {

  // get the short-cut name for files to process
  std::vector<std::string> all_files, all_dirs;
  std::string dir_name = "test";
  if (argc > 1)
    dir_name = argv[1];

  // get the jobType for the files
  std::string jobType = "bkg";
  if (argc > 2)
    jobType = argv[2];

  // decide which map to look-up files from and which prefix to prepend
  std::vector<std::string> files;
  if (jobType == "bkg") {
    files = bkg_samples[dir_name];
    add_pref(files, bkg_pref);
  } else if (jobType == "sig") {
    files = sig_samples[dir_name];
    add_pref(files, sig_pref);
  } else if (jobType == "data") {
    files = data_samples[dir_name];
    add_pref(files, data_pref);
  } else {
    std::cerr << "jobType must be bkg, sig, or data. You gave " << jobType << std::endl;
    return -1;
  }

  // read all root files in the given directory
  const int dir_err = system(("mkdir "+dir_name).c_str());
  for (auto& ifile : files) {
    read_directory(ifile, all_files, all_dirs);
  }

  // set flag for MC or data
  bool isMC(true);
  if (jobType == "data")
    isMC = false;

  TH1F* nevents = new TH1F("nevents", "N(events)", 2, 0.5, 2.5);

  // I/O files and the ntuple
  std::cout << "Begin loading files " << all_files.size() << " ..." << std::endl;
  unsigned i = 0;
  for (auto& ifile : all_files) {
    std::cout << "Loading file: " << i+1 << " out of " << all_files.size() << " files.\r" << std::flush;
    
    auto ntuple = new TChain("et/final/Ntuple");
    auto idir = all_dirs.at(i);
    auto open_file = new TFile((idir+"/"+ifile).c_str(), "READ");
    auto evt_count = (TH1F*)open_file->Get("et/eventCount")->Clone();
    auto wt_count = (TH1F*)open_file->Get("et/summedWeights")->Clone();

    nevents->SetBinContent(1, evt_count->Integral());
    nevents->SetBinContent(2, wt_count->Integral());

    open_file->Close();
    ntuple->Add((idir+"/"+ifile).c_str());
    std::string suffix = dir_name+"/Skim_";
    auto fout = new TFile((suffix+ifile).c_str(), "RECREATE");

    TTree* newtree = new TTree("etau_tree","etau_tree");
    etau_tree* skimmer = new etau_tree(ntuple, newtree, isMC); // hard-code true for now
    skimmer->do_skimming();
    auto skimmed_tree = skimmer->fill_tree();
    events += skimmed_tree->GetEntries();

    fout->cd();
    nevents->Write();
    skimmed_tree->Write();
    fout->Close();

    i += 1;
  }
  std::cout << std::endl;
  std::cout << "\n"  << events << " events saved in the tree.\n" << std::endl;

  return 0;
}
