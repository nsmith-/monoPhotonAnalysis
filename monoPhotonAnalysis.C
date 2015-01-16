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

   TFile * output = new TFile(_outFilename.c_str(), "recreate");
   output->cd();
   TTree * outTree = fChain->GetTree()->CloneTree(0); // 0 means do not copy entries
   outTree->CopyAddresses(fChain->GetTree());
   int selectedPhoton = 0;
   outTree->Branch("selectedPhoton", &selectedPhoton);
   float deltaPhi_PhoMET = 0;
   outTree->Branch("deltaPhi_phoMET", &deltaPhi_PhoMET);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t npassed = 0;
   double nQCD = 0.;
   std::map<std::string, long> cutFlow;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if ( ientry % 10000 == 0 )
      {
        printf("Processed %7lld / %7lld events (% 2.1f%%)\r", jentry, nentries, jentry*100./nentries);
        cout.flush();
      }
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // first medium photon cut
      if ( kDoQCDBackground )
      {
        if ( !HasQCD(selectedPhoton) ) continue;
        cutFlow["qcd selection"]++;
      }
      else
      {
        if ( !HasMediumPhoton(selectedPhoton) ) continue;
        cutFlow["medium photon"]++;
      }

      // Make sure MET is significant
      if ( pfMET < 90. ) continue;
      cutFlow["pfMET > 90"]++;

      // photon and MET should be back to back
      deltaPhi_PhoMET = deltaPhi(pfMETPhi, phoPhi->at(selectedPhoton));
      if ( fabs(deltaPhi_PhoMET) < 2. ) continue;
      cutFlow["deltaPhi"]++;

      // Veto event if electron or muon
      if ( electronVeto(selectedPhoton) ) continue;
      cutFlow["electron veto"]++;
      if ( muonVeto(selectedPhoton) ) continue;
      cutFlow["muon veto"]++;

      // -----  End Cut Selection -----
      npassed++;
      nQCD += 4.3865e-02+7.0223e-04*phoEt->at(selectedPhoton);

      outTree->Fill();
   }
   std::cout << "\nPassed " << npassed << " events out of " << nentries << std::endl;
   if ( kDoQCDBackground )
     std::cout << "nQCD = " << nQCD << std::endl;
   std::cout << "Cut flow summary --------" << std::endl;
   for ( const auto& cut : cutFlow ) std::cout << setw(30) << cut.first << " : " << cut.second << " events passed." << std::endl;

   outTree->Write();
   output->Close();
   delete output;
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
    bool isMediumPhoton = etCut            && scEtaCut         &&  hOverECut   &&
                          sigmaIEtaIEtaCut && sigmaIPhiIPhiCut && pixelSeedCut && 
                          r9Cut && isIsolatedPhoton(i);
                            
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
    if ( fabs(eleEta->at(i)) < 1.479 )
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
      if ( ! eledEtaAtVtx->at(i) < 0.007 ) continue;
      if ( ! eledPhiAtVtx->at(i) < 0.15 ) continue;
      if ( ! eleSigmaIEtaIEta_2012->at(i) < 0.01 ) continue;
      if ( ! eleHoverE->at(i) < 0.12 ) continue;
      if ( ! fabs(eleD0->at(i)) < 0.02 ) continue;
      if ( ! fabs(eleDz->at(i)) < 0.2 ) continue;
      if ( ! eleEoverPInv->at(i) < 0.05 ) continue;
      if ( ! eleMissHits->at(i) <= 1 ) continue;
      if ( ! eleAbsIso/elePt->at(i) < 0.15 ) continue;
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
      if ( ! eledEtaAtVtx->at(i) < 0.009 ) continue;
      if ( ! eledPhiAtVtx->at(i) < 0.10 ) continue;
      if ( ! eleSigmaIEtaIEta_2012->at(i) < 0.03 ) continue;
      if ( ! eleHoverE->at(i) < 0.10 ) continue;
      if ( ! fabs(eleD0->at(i)) < 0.02 ) continue;
      if ( ! fabs(eleDz->at(i)) < 0.02 ) continue;
      if ( ! eleEoverPInv->at(i) < 0.05 ) continue;
      if ( ! eleMissHits->at(i) <= 1 ) continue;
      if ( ! eleAbsIso/elePt->at(i) < 0.10 ) continue;
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

bool monoPhotonAnalysis::isIsolatedPhoton(int i)
{
  const float MED_CH_HADRON_ISO = 1.5;
  const float MED_NEU_HADRON_ISO = 1.0 + 0.04*phoEt->at(i);
  const float MED_PHOTON_ISO = 0.7 + 0.005*phoEt->at(i);

  auto correctedIso = [this, i](double naiveIso, std::string isoParticle)
  {
    return max(naiveIso - this->rho*getPhotonEffectiveArea(isoParticle, i), 0.);
  };

  bool isolatedFromCH = correctedIso(phoPFChWorstIso->at(i), ( kUseWorstChIso ) ? "worst charged hadron" : "charged hadron" ) < MED_CH_HADRON_ISO;
  bool isolatedFromNH = correctedIso(phoPFNeuIso->at(i), "neutral hadron") < MED_NEU_HADRON_ISO;
  bool isolatedFromPho = correctedIso(phoPFPhoIso->at(i), "photon") < MED_PHOTON_ISO;

  return (isolatedFromCH && isolatedFromNH && isolatedFromPho);
}

double monoPhotonAnalysis::getPhotonEffectiveArea(std::string isoParticle, int i)
{
  std::map<std::string, std::vector<float>> effAreaLookup;
  effAreaLookup["charged hadron"] = { 0.012, 0.010 };
  effAreaLookup["worst charged hadron"] = { 0.075, 0.0617 };
  effAreaLookup["neutral hadron"] = { 0.030, 0.057 };
  effAreaLookup["photon"] = { 0.148, 0.130 };

  auto etaRegion = [](float photonEta)
  {
    if (fabs(photonEta) < 1.0)
      return 0;
    else if (fabs(photonEta) < 1.479)
      return 1;
    else
      return 1;
  };
  return effAreaLookup[isoParticle][etaRegion(phoSCEta->at(i))];
}

// ///////// QCD selection ///////////
//     1.  maximum PF Charged Hadron Isolation = Min(5.0*2.6 , 0.20*photonPt) 
//     2.  maximum PF Photon Isolation         = Min(5.0*(1.3+0.005*photonPt) , 0.20*photonPt)
//     3.  maximum PF Neutral Hadron Isolation = Min(5.0*(3.5+0.04*photonPt) , 0.20*photonPt)
//     4.  Define an upperBound using ALL these cuts (look at the bottom for area values): 
// 
//       a)  (Max(photon PF Worst Charged Hadron Isolation - rho* WorstChargedHadronsAreas ,0.0) < maximum PF Charged Hadron Isolation ) 
//       b)  (photon SigmaIEtaIEta < 0.013 ) 
//       c)  (Max(photon PF Photon Isolation - rho *PhotonAreas) ,0.0) < maximum PF Photon Isolation ) 
//       d)  (Max(photon PF Neutral Hadron Isolation - rho*NeutralHadronAreas) ,0.0) < maximum PF Neutral Hadron Isolation)
// 
//     5.  Define a lowerBound where ANY ONE of these cuts are satisfied:
//     
//       a) (Max(photon PF Photon Isolation - rho*PhotonAreas) ,0.0) > 1.3+0.005*photonPt
//       b) (Max(photon PF Neutral Isolation - rho*NeutralHadronAreas) ,0.0) > 3.5+0.04*photonPt
//       c) (Max(photon PF Worst Charged Hadron Isolation - rho*WorstChargedHadronAreas) ,0.0) > 2.6
bool monoPhotonAnalysis::isQCDLike(int i)
{
  const float CH_HADRON_ISO = min(5.*2.6, 0.2*phoEt->at(i));
  const float PHOTON_ISO = min(5.*(1.3+0.005*phoEt->at(i)), 0.2*phoEt->at(i));
  const float NEU_HADRON_ISO = min(5.*(3.5+0.04*phoEt->at(i)), 0.2*phoEt->at(i));

  auto correctedIso = [&](double naiveIso, std::string isoParticle)
  {
    return max(naiveIso - this->rho*getPhotonEffectiveArea(isoParticle, i), 0.);
  };

  bool tooNotIsolatedFromCH = correctedIso(phoPFChWorstIso->at(i), "worst charged hadron") < CH_HADRON_ISO;
  bool tooNotIsolatedFromNH = correctedIso(phoPFNeuIso->at(i), "neutral hadron") < NEU_HADRON_ISO;
  bool tooNotIsolatedFromPho = correctedIso(phoPFPhoIso->at(i), "photon") < PHOTON_ISO;

  bool upper = tooNotIsolatedFromCH && tooNotIsolatedFromNH && tooNotIsolatedFromPho;

  bool notIsolatedFromCH =  correctedIso(phoPFChWorstIso->at(i), "worst charged hadron") > 2.6;
  bool notIsolatedFromNH =  correctedIso(phoPFNeuIso->at(i), "neutral hadron") > 3.5+0.04*phoEt->at(i);
  bool notIsolatedFromPho = correctedIso(phoPFPhoIso->at(i), "photon") > 1.3+0.005*phoEt->at(i);

  bool lower = notIsolatedFromCH || notIsolatedFromNH || notIsolatedFromPho;

  return upper && lower;
}

bool monoPhotonAnalysis::HasQCD(int& qcdNo)
{
  int nQCD = 0;

  for (int i = 0; i < nPho; i++) {

    bool etCut            = phoEt->at           (i)  > 100.0;
    bool scEtaCut         = fabs(phoSCEta->at   (i)) <   1.4442;
    bool hOverECut        = phoHoverE->at       (i)  <   0.05;
    bool sigmaIEtaIEtaCut = phoSigmaIEtaIEta->at(i)  >   0.001 &&
                            phoSigmaIEtaIEta->at(i)  <   0.013;
    bool sigmaIPhiIPhiCut = phoSigmaIPhiIPhi->at(i)  >   0.001;
    bool pixelSeedCut     = phohasPixelSeed->at (i)  ==  0;
    bool r9Cut            = phoR9->at           (i)  <   1.0;

    // Only barrel photons, so only two bins in effective area
    bool isMediumPhoton = etCut            && scEtaCut         &&  hOverECut   &&
                          sigmaIEtaIEtaCut && sigmaIPhiIPhiCut && pixelSeedCut && 
                          r9Cut && isQCDLike(i);
                            
    if (isMediumPhoton) {
      nQCD++;
      qcdNo = i;
    }
  }

  return nQCD == 1;
}
