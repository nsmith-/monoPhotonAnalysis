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
      if ( ientry % 10000 == 0 ) printf("Processed %7lld / %7lld events (% 2.1f%%)\n", jentry, nentries, jentry*100./nentries);
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // first medium photon cut
      int selectedPhoton = 0;
      if ( !HasMediumPhoton(selectedPhoton) ) continue;
      cutFlow["medium photon"]++;

      // Make sure MET is significant
      if ( pfMET < 90. ) continue;
      cutFlow["pfMET > 90"]++;

      // photon and MET should be back to back
      if ( fabs(deltaPhi(pfMETPhi, phoPhi->at(selectedPhoton))) < 2. ) continue;
      cutFlow["deltaPhi"]++;

      // Veto event if electron or muon
      if ( electronVeto(selectedPhoton) ) continue;
      cutFlow["electron veto"]++;
      if ( muonVeto(selectedPhoton) ) continue;
      cutFlow["muon veto"]++;

      // -----  End Cut Selection -----
      npassed++;
   }
   std::cout << "Passed " << npassed << " events out of " << nentries << std::endl;
   std::cout << "Cut flow summary --------" << std::endl;
   for ( const auto& cut : cutFlow ) std::cout << setw(30) << cut.first << " : " << cut.second << " events passed." << std::endl;
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
    float effectiveAreaLowEta = ( kUseWorstChIso ) ? 0.075 : 0.012;
    float effectiveAreaHighEta = ( kUseWorstChIso ) ? 0.0617 : 0.010;
    bool rhoCorrPFchi     = ( (phoSCEta->at(i) < 1.) ? max(chi-rho*effectiveAreaLowEta, 0.f) : max(chi-rho*effectiveAreaHighEta, 0.f) ) < 1.2;
    bool rhoCorrPFnhi     = ( (phoSCEta->at(i) < 1.) ? max(phoPFNeuIso->at(i)-rho*0.030, 0.) : max(phoPFNeuIso->at(i)-rho*0.057, 0.) ) < 1.+0.04*phoEt->at(i);
    bool rhoCorrPFphoi    = ( (phoSCEta->at(i) < 1.) ? max(phoPFPhoIso->at(i)-rho*0.148, 0.) : max(phoPFPhoIso->at(i)-rho*0.130, 0.) ) < 0.7+0.005*phoEt->at(i);

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
  for (int i=0; i < nEle; i++)
  {
    // at least 10GeV electron
    if ( elePt->at(i) < 10. ) continue;

    // electron PF Charged Hadron Isolation + max (0.0, electron PF Photon Isolation + electron PF Neutral Hadron Isolation - 0.5*electron PF PU Isolation)
    float eleAbsIso = elePFChIso->at(i) + max(0., elePFPhoIso->at(i) + elePFNeuIso->at(i) - 0.5*elePFPUIso->at(i));

    // Barrel and endcap quality criteria
    if ( eleEta->at(i) < 1.479 ) 
    {
      // i.   check  absolute electron dEtaAtVtx < 0.007
      // ii.  check absolute electron dPhiAtVtx < 0.15 
      // iii. check electron SigmaIEtaIEta_2012 < 0.01
      // iv.  electron H/E < 0.12 
      // v.   absolute value of electron D0 < 0.02
      // vi.  absolute value of electron Dz < 0.2
      // vii. electron EoverPInv < 0.05
      // viii.electron MissHits <= 1   
      // ix.  absolute isolation/electron Pt < 0.15 
      if ( eledEtaAtVtx->at(i) > 0.007 ) continue;
      if ( eledPhiAtVtx->at(i) > 0.15 ) continue;
      if ( eleSigmaIEtaIEta_2012->at(i) > 0.01 ) continue;
      if ( eleHoverE->at(i) > 0.12 ) continue;
      if ( fabs(eleD0->at(i)) > 0.02 ) continue;
      if ( fabs(eleDz->at(i)) > 0.2 ) continue;
      if ( eleEoverPInv->at(i) > 0.05 ) continue;
      if ( eleMissHits->at(i) > 1 ) continue;
      if ( eleAbsIso/elePt->at(i) > 0.15 ) continue;
    }
    else
    {
      // i.   absolute value of electron dEtaAtVtx < 0.009
      // ii.  absolute value of electron dPhiAtVtx < 0.10
      // iii. electron SigmaIEtaIEta_2012 < 0.03
      // iv.  electron H/E < 0.10
      // v.   absolute value of electron D0 < 0.02
      // vi.  absolute value of electron Dz < 0.02
      // vii. electron EoverPInv < 0.05
      // viii.electron MissHits <=1 
      // ix.  absolute isolation/electron Pt < 0.10
      if ( eledEtaAtVtx->at(i) > 0.009 ) continue;
      if ( eledPhiAtVtx->at(i) > 0.10 ) continue;
      if ( eleSigmaIEtaIEta_2012->at(i) > 0.03 ) continue;
      if ( eleHoverE->at(i) > 0.10 ) continue;
      if ( fabs(eleD0->at(i)) > 0.02 ) continue;
      if ( fabs(eleDz->at(i)) > 0.02 ) continue;
      if ( eleEoverPInv->at(i) > 0.05 ) continue;
      if ( eleMissHits->at(i) > 1 ) continue;
      if ( eleAbsIso/elePt->at(i) > 0.10 ) continue;
    }

    // If electron overlaps good photon, skip
    float deltaR = sqrt( pow(deltaPhi(elePhi->at(i), phoPhi->at(photonNo)),2) + pow(eleEta->at(i)-phoSCEta->at(photonNo),2) );
    if ( deltaR < 0.5 ) continue;

    // We found a good electron that doesn't overlap the selected photon, veto
    return true;
  }

  // We found no good electrons, don't veto event
  return false;
}

bool monoPhotonAnalysis::muonVeto(const int photonNo)
{
  for ( int i=0; i < nMu; i++ )
  {
    // a) muon Pt > 10
    // b) muon Isolation Track / muon Pt < 0.10
    // c) muon Chi2/NDF < 10
    // d) muon hits > 0
    // e) muon pixel hits > 0
    // f) muon Stations > 1
    // g) muon D0 < 0.2
    // h) muon Dz < 0.5
    // i) muon Track Layers > 5
    if ( muPt->at(i) < 10. ) continue;
    if ( muIsoTrk->at(i) / muPt->at(i) > 0.1 ) continue;
    if ( muChi2NDF->at(i) >= 10 ) continue;
    if ( muMuonHits->at(i) == 0 ) continue;
    if ( muPixelHits->at(i) == 0 ) continue;
    if ( muStations->at(i) <= 1 ) continue;
    if ( fabs(muD0->at(i)) > 0.2 ) continue;
    if ( fabs(muDz->at(i)) > 0.5 ) continue;
    if ( muTrkLayers->at(i) <= 5 ) continue;

    // If muon overlaps good photon, skip
    float deltaR = sqrt( pow(deltaPhi(muPhi->at(i), phoPhi->at(photonNo)),2) + pow(muEta->at(i)-phoSCEta->at(photonNo),2) );
    if ( deltaR < 0.5 ) continue;

    // We found a good muon that doesn't overlap the selected photon, veto
    return true;
  }

  // No good muons
  return false;
}

float monoPhotonAnalysis::deltaPhi(float phi1, float phi2)
{
  // Thanks to reco::deltaPhi in CMSSW
  double result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

