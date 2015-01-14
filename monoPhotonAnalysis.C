#define monoPhotonAnalysis_cxx
#include "monoPhotonAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void monoPhotonAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

bool monoPhotonAnalysis::HasMediumPhoton(int& photonNo)
{
  int nMediumPhotons = 0;

  for (int i = 0; i < nPho; i++) {

    bool etCut            = phoEt->at           (i)  > 100.0;
    bool scEtaCut         = fabs(phoSCEta->at   (i)) <   1.4442;
    bool hOverECut        = phoHoverE->at       (i)  <   0.05;
    bool sigmaIEtaIEtaCut = phoSigmaIEtaIEta->at(i)  >   0.001 &&
                            phoSigmaIEtaIEta->at(i)  <   0.011;
    bool sigmaIPhiIPhiCut = phoSigmaIPhiIPhi->at(i)  >   0.001;
    bool pixelSeedCut     = phohasPixelSeed->at (i)  ==  0;
    bool r9Cut            = phoR9->at           (i)  <   1.0;
    auto rhoCorrection = [&](float var, float eta)
    {
      // Just barrel, only need the two
      if ( fabs(eta) < 1. )
        return max(var-rho*0.148, 0.);
      else if ( fabs(eta) < 1.479 )
        return max(var-rho*0.130, 0.);
      return 0.;
    };
    bool rhoCorrPFchi     = rhoCorrection(phoPFChIso->at(i), phoEta->at(i)) < 1.2;
    bool rhoCorrPFnhi     = rhoCorrection(phoPFNeuIso->at(i), phoEta->at(i)) < 1.+0.04*phoEt->at(i);
    bool rhoCorrPFphoi    = rhoCorrection(phoPFPhoIso->at(i), phoEta->at(i)) < 0.7+0.005*phoEt->at(i);

    bool isMediumPhoton = etCut            && scEtaCut         &&  hOverECut   &&
                          sigmaIEtaIEtaCut && sigmaIPhiIPhiCut && pixelSeedCut && 
                          r9Cut && rhoCorrPFchi && rhoCorrPFnhi && rhoCorrPFphoi;

    if (isMediumPhoton) {
      nMediumPhotons++;
      photonNo = i;
    }
  }

  return nMediumPhotons == 1;
}
