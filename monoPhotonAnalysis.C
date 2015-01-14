#define monoPhotonAnalysis_cxx
#include "monoPhotonAnalysis.h"
#include <cstdio>
#include <iostream>
#include <map>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void monoPhotonAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t npassed = 0;
   std::map<std::string, long> cutFlow;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if ( ientry % 10000 == 0 ) printf("Processed %7d / %7d events (% 2.1f\%)\n", jentry, nentries, jentry*100./nentries);
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // first medium photon cut
      int selectedPhoton = 0;
      if ( !HasMediumPhoton(selectedPhoton) ) continue;
      cutFlow["medium photon"]++;

      // Make sure MET is significant
      if ( pfMET < 90. ) continue;
      cutFlow["pfMET < 90"]++;

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
      cutFlow["deltaPhi"]++;

      // Veto event if electron or muon
      // if ( electronVeto(selectedPhoton) ) continue;
      cutFlow["electron veto"]++;
      // if ( muonVeto(selectedPhoton) ) continue;
      cutFlow["muon veto"]++;

      // -----  End Cut Selection -----
      npassed++;
   }
   std::cout << "Passed " << npassed << " events out of " << nentries << std::endl;
   std::cout << "Cut flow summary --------" << std::endl;
   for(const auto cut : cutFlow)
   {
     std::cout << setw(30) << cut.first << " : " << cut.second << " events passed." << std::endl;
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
    float chi = ( kUseWorstChIso ) ? phoPFChWorstIso->at(i) : phoPFChIso->at(i);
    bool rhoCorrPFchi     = ( (phoEta->at(i) < 1.) ? max(chi-rho*0.012, 0.) : max(chi-rho*0.010, 0.) ) < 1.2;
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

bool monoPhotonAnalysis::electronVeto(const int photonNo)
{
  return true;
}

bool monoPhotonAnalysis::muonVeto(const int photonNo)
{
  return true;
}
