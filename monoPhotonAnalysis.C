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

      // first medium photon cut
      int selectedPhoton = 0;
      if ( !HasMediumPhoton(selectedPhoton) ) continue;

      // Make sure MET is significant
      if ( pfMET < 90. ) continue;

      // photon and MET should be back to back
      auto deltaPhi = [](float phi1, float phi2)
      {
        // Thanks to reco::deltaPhi in CMSSW
        double result = phi1 - phi2;
        while (result > M_PI) result -= 2*M_PI;
        while (result <= -M_PI) result += 2*M_PI;
        return result;
      };
      if ( deltaPhi(pfMETPhi, phoPhi->at(selectedPhoton)) < 2. ) continue;

      // Veto event if electron or muon
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

    // Only barrel photons, so only two bins in effective area
    bool rhoCorrPFchi     = ( (phoEta->at(i) < 1.) ? max(phoPFChIso->at(i)-rho*0.012, 0.) : max(phoPFChIso->at(i)-rho*0.010, 0.) ) < 1.2;
    bool rhoCorrPFnhi     = ( (phoEta->at(i) < 1.) ? max(phoPFNeuIso->at(i)-rho*0.030, 0.) : max(phoPFNeuIso->at(i)-rho*0.057, 0.) ) < 1.+0.04*phoEt->at(i);
    bool rhoCorrPFphoi    = ( (phoEta->at(i) < 1.) ? max(phoPFPhoIso->at(i)-rho*0.148, 0.) : max(phoPFPhoIso->at(i)-rho*0.130, 0.) ) < 0.7+0.005*phoEt->at(i);

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
