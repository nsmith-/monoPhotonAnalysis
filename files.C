{
  // Load up files used
  TFile * fdata = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/data.root");
  TFile * fqcd = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/qcdMc.root");
  TFile * fznunu = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/zNuNuGammaMc.root");
  TFile * fwJets = TFile::Open("/eos/uscms/store/user/weinberg/cmsdas2015/wJetsMc.root");
  TTree * data = fdata->Get("EventTree");
  TTree * qcd = fqcd->Get("EventTree");
  TTree * znunu = fznunu->Get("EventTree");
  TTree * wJets = fwJets->Get("EventTree");

  // Now, use
  // .L monoPhotonAnalysis.C+
  // .x files.C
  // monoPhotonAnalysis adata(data);
}
