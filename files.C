// Run final selection on data and all MC sets
{
  gSystem->Load("monoPhotonAnalysis_C.so");
  // Load up files used
  TFile * fdata = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/data.root");
  TFile * fqcd = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/qcdMc.root");
  TFile * fznunu = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/zNuNuGammaMc.root");
  TFile * fwJets = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/wJetsMc.root");
  TFile * fsig = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/addMc.root");
  TTree * data = fdata->Get("EventTree");
  TTree * qcd = fqcd->Get("EventTree");
  TTree * znunu = fznunu->Get("EventTree");
  TTree * wJets = fwJets->Get("EventTree");
  TTree * signal = fsig->Get("EventTree");

  monoPhotonAnalysis adata(data, "data_selection.root");
  monoPhotonAnalysis aqcd(qcd, "qcd_selection.root");
  monoPhotonAnalysis aznunu(znunu, "znunu_selection.root");
  monoPhotonAnalysis awJets(wJets, "wJets_selection.root");
  monoPhotonAnalysis asig(signal, "signal_selection.root");
}
