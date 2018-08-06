#include <map>

std::string bkg_pref("/hdfs/store/user/caillol/SMHTT_mc_feb13/");
std::string sig_pref("/hdfs/store/user/truggles/SMHTT_signals_may30/");
std::string data_pref("/hdfs/store/user/caillol/SMHTT_reminiaod_feb14/");

static std::map<std::string, std::vector<std::string>> bkg_samples = {
  {"test", std::vector<std::string>{"root_files/csync"}},
  
  {"DYJets"    , std::vector<std::string>{"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v2/", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1/"}},
  {"DYJets1"   , std::vector<std::string>{"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"DYJets2"   , std::vector<std::string>{"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"DYJets3"   , std::vector<std::string>{"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"DYJets4"   , std::vector<std::string>{"DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"WJets"     , std::vector<std::string>{"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"WJets1"    , std::vector<std::string>{"W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/"}},
  {"WJets2"    , std::vector<std::string>{"W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/", "W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1/"}},
  {"WJets3"    , std::vector<std::string>{"W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/", "W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1/"}},
  {"WJets4"    , std::vector<std::string>{"W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1/", "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1/", "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1/"}},
  {"TT"        , std::vector<std::string>{"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1/", "VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1/"}},
  {"VV2l2nu"   , std::vector<std::string>{"VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6_ext1-v1/"}},
  {"Tbar-tchan", std::vector<std::string>{"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1/"}},
  {"T-tchan"   , std::vector<std::string>{"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1/"}},
  {"Tbar-tW"   , std::vector<std::string>{"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1/"}},
  {"T-tW"      , std::vector<std::string>{"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1/"}},
  {"WW1l1nu2q" , std::vector<std::string>{"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1/"}},
  {"ZZ2l2q"    , std::vector<std::string>{"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1/"}},
  {"ZZ4l"      , std::vector<std::string>{"ZZTo4L_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1/"}},
  {"EWKWMinus" , std::vector<std::string>{"EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6-v1", "EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext1-v1/", "EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext2-v1/"}},
  {"EWKWPlus"  , std::vector<std::string>{"EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6-v1/", "EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext1-v1/", "EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext2-v1/"}},
  {"EWKZ2l"    , std::vector<std::string>{"EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6-v1/", "EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6_ext1-v1/", "EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6_ext2-v1/"}},
  {"EWKZ2nu"   , std::vector<std::string>{"EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6-v1/", "EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext1-v1/", "EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext2-v1/"}},
  {"ggH_WW125" , std::vector<std::string>{"GluGluHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"VBF_WW125" , std::vector<std::string>{"VBFHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"WGLNu"     , std::vector<std::string>{"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1/", "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v1/"}},
  {"WGstarEE"  , std::vector<std::string>{"WGstarToLNuEE_012Jets_13TeV-madgraph_v6-v1/"}},
  {"WGstarMuMu", std::vector<std::string>{"WGstarToLNuMuMu_012Jets_13TeV-madgraph_v6-v1/"}},
  {"WWW"       , std::vector<std::string>{"WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_v6-v1/"}},
  {"WZ3l1nu"   , std::vector<std::string>{"WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8_v6-v1/"}},
  {"WZ1l1nu2q" , std::vector<std::string>{"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v3/"}},
  {"WZ1l3nu"   , std::vector<std::string>{"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1/"}},
  {"WZ2l2q"    , std::vector<std::string>{"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1/"}},
};

static std::map<std::string, std::vector<std::string>> sig_samples = {
  {"ggHtoTauTau125"  , std::vector<std::string>{"GluGluHToTauTau_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"VBFHtoTauTau125" , std::vector<std::string>{"VBFHToTauTau_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"WPlusHTauTau125" , std::vector<std::string>{"WplusHToTauTau_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"WMinusHTauTau125", std::vector<std::string>{"WminusHToTauTau_M125_13TeV_powheg_pythia8_v6-v1/"}},
  {"ZHTauTau125"     , std::vector<std::string>{"ZHToTauTau_M125_13TeV_powheg_pythia8_v6-v1/"}},
};

static std::map<std::string, std::vector<std::string>> data_samples = {
    {"dataET-B", std::vector<std::string>{"data_Tau_Run2016B_v1/", "data_Tau_Run2016B_v2/"}},
    {"dataET-C", std::vector<std::string>{"data_Tau_Run2016C/"}},
    {"dataET-D", std::vector<std::string>{"data_Tau_Run2016D/"}},
    {"dataET-E", std::vector<std::string>{"data_Tau_Run2016E/"}},
    {"dataET-F", std::vector<std::string>{"data_Tau_Run2016F/"}},
    {"dataET-G", std::vector<std::string>{"data_Tau_Run2016G/"}},
    {"dataET-H", std::vector<std::string>{"data_Tau_Run2016H_v2/", "data_Tau_Run2016H_v3/"}},
};

