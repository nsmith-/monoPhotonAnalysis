//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 14 08:59:44 2015 by ROOT version 5.34/18
// from TTree EventTree/Event data
// found on file: /eos/uscms/store/user/weinberg/cmsdas2015/data.root
//////////////////////////////////////////////////////////

#ifndef monoPhotonAnalyzer_h
#define monoPhotonAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iomanip>

// Fixed size dimensions of array or collections stored in the TTree if any.

class monoPhotonAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nTrks;
   Float_t         rho;
   std::vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   std::vector<int>     *mcPID;
   std::vector<float>   *mcVtx_x;
   std::vector<float>   *mcVtx_y;
   std::vector<float>   *mcVtx_z;
   std::vector<float>   *mcPt;
   std::vector<float>   *mcMass;
   std::vector<float>   *mcEta;
   std::vector<float>   *mcPhi;
   std::vector<float>   *mcE;
   std::vector<float>   *mcEt;
   std::vector<int>     *mcGMomPID;
   std::vector<int>     *mcMomPID;
   std::vector<float>   *mcMomPt;
   std::vector<float>   *mcMomMass;
   std::vector<float>   *mcMomEta;
   std::vector<float>   *mcMomPhi;
   std::vector<int>     *mcIndex;
   std::vector<int>     *mcDecayType;
   std::vector<int>     *mcParentage;
   std::vector<int>     *mcStatus;
   std::vector<float>   *mcCalIsoDR03;
   std::vector<float>   *mcTrkIsoDR03;
   std::vector<float>   *mcCalIsoDR04;
   std::vector<float>   *mcTrkIsoDR04;
   Int_t           nPUInfo;
   std::vector<int>     *nPU;
   std::vector<int>     *puBX;
   std::vector<float>   *puTrue;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Int_t           nEle;
   std::vector<int>     *eleCharge;
   std::vector<int>     *eleChargeConsistent;
   std::vector<float>   *eleEn;
   std::vector<float>   *eleSCEn;
   std::vector<float>   *eleESEn;
   std::vector<float>   *eleD0;
   std::vector<float>   *eleDz;
   std::vector<float>   *elePt;
   std::vector<float>   *eleEta;
   std::vector<float>   *elePhi;
   std::vector<float>   *eleSCEta;
   std::vector<float>   *eleSCPhi;
   std::vector<float>   *eleSCRawEn;
   std::vector<float>   *eleSCEtaWidth;
   std::vector<float>   *eleSCPhiWidth;
   std::vector<float>   *eleHoverE;
   std::vector<float>   *eleEoverP;
   std::vector<float>   *eleEoverPInv;
   std::vector<float>   *eleBrem;
   std::vector<float>   *eledEtaAtVtx;
   std::vector<float>   *eledPhiAtVtx;
   std::vector<float>   *eleSigmaIEtaIEta;
   std::vector<float>   *eleSigmaIEtaIPhi;
   std::vector<float>   *eleSigmaIPhiIPhi;
   std::vector<float>   *eleSigmaIEtaIEta_2012;
   std::vector<int>     *eleConvVeto;
   std::vector<int>     *eleMissHits;
   std::vector<float>   *eleESEffSigmaRR;
   std::vector<float>   *elePFChIso;
   std::vector<float>   *elePFPhoIso;
   std::vector<float>   *elePFNeuIso;
   std::vector<float>   *elePFPUIso;
   std::vector<float>   *eleBC1E;
   std::vector<float>   *eleBC1Eta;
   std::vector<float>   *eleBC2E;
   std::vector<float>   *eleBC2Eta;
   Int_t           nPho;
   std::vector<float>   *phoE;
   std::vector<float>   *phoEt;
   std::vector<float>   *phoEta;
   std::vector<float>   *phoPhi;
   std::vector<float>   *phoSCE;
   std::vector<float>   *phoSCRawE;
   std::vector<float>   *phoESEn;
   std::vector<float>   *phoSCEta;
   std::vector<float>   *phoSCPhi;
   std::vector<float>   *phoSCEtaWidth;
   std::vector<float>   *phoSCPhiWidth;
   std::vector<float>   *phoSCBrem;
   std::vector<int>     *phohasPixelSeed;
   std::vector<int>     *phoEleVeto;
   std::vector<float>   *phoR9;
   std::vector<float>   *phoHoverE;
   std::vector<float>   *phoSigmaIEtaIEta;
   std::vector<float>   *phoSigmaIEtaIPhi;
   std::vector<float>   *phoSigmaIPhiIPhi;
   std::vector<float>   *phoE1x3;
   std::vector<float>   *phoE2x2;
   std::vector<float>   *phoE2x5Max;
   std::vector<float>   *phoE5x5;
   std::vector<float>   *phoESEffSigmaRR;
   std::vector<float>   *phoSigmaIEtaIEta_2012;
   std::vector<float>   *phoSigmaIEtaIPhi_2012;
   std::vector<float>   *phoSigmaIPhiIPhi_2012;
   std::vector<float>   *phoE1x3_2012;
   std::vector<float>   *phoE2x2_2012;
   std::vector<float>   *phoE2x5Max_2012;
   std::vector<float>   *phoE5x5_2012;
   std::vector<float>   *phoPFChIso;
   std::vector<float>   *phoPFPhoIso;
   std::vector<float>   *phoPFNeuIso;
   std::vector<float>   *phoPFChWorstIso;
   std::vector<float>   *phoPFChIsoFrix1;
   std::vector<float>   *phoPFChIsoFrix2;
   std::vector<float>   *phoPFChIsoFrix3;
   std::vector<float>   *phoPFChIsoFrix4;
   std::vector<float>   *phoPFChIsoFrix5;
   std::vector<float>   *phoPFChIsoFrix6;
   std::vector<float>   *phoPFChIsoFrix7;
   std::vector<float>   *phoPFChIsoFrix8;
   std::vector<float>   *phoPFPhoIsoFrix1;
   std::vector<float>   *phoPFPhoIsoFrix2;
   std::vector<float>   *phoPFPhoIsoFrix3;
   std::vector<float>   *phoPFPhoIsoFrix4;
   std::vector<float>   *phoPFPhoIsoFrix5;
   std::vector<float>   *phoPFPhoIsoFrix6;
   std::vector<float>   *phoPFPhoIsoFrix7;
   std::vector<float>   *phoPFPhoIsoFrix8;
   std::vector<float>   *phoPFNeuIsoFrix1;
   std::vector<float>   *phoPFNeuIsoFrix2;
   std::vector<float>   *phoPFNeuIsoFrix3;
   std::vector<float>   *phoPFNeuIsoFrix4;
   std::vector<float>   *phoPFNeuIsoFrix5;
   std::vector<float>   *phoPFNeuIsoFrix6;
   std::vector<float>   *phoPFNeuIsoFrix7;
   std::vector<float>   *phoPFNeuIsoFrix8;
   std::vector<float>   *phoBC1E;
   std::vector<float>   *phoBC1Eta;
   std::vector<float>   *phoBC2E;
   std::vector<float>   *phoBC2Eta;
   Int_t           nMu;
   std::vector<float>   *muPt;
   std::vector<float>   *muEta;
   std::vector<float>   *muPhi;
   std::vector<int>     *muCharge;
   std::vector<int>     *muType;
   std::vector<int>     *muIsGood;
   std::vector<float>   *muD0;
   std::vector<float>   *muDz;
   std::vector<float>   *muChi2NDF;
   std::vector<float>   *muInnerD0;
   std::vector<float>   *muInnerDz;
   std::vector<int>     *muTrkLayers;
   std::vector<int>     *muPixelLayers;
   std::vector<int>     *muPixelHits;
   std::vector<int>     *muMuonHits;
   std::vector<int>     *muStations;
   std::vector<int>     *muTrkQuality;
   std::vector<float>   *muIsoTrk;
   std::vector<float>   *muPFChIso;
   std::vector<float>   *muPFPhoIso;
   std::vector<float>   *muPFNeuIso;
   std::vector<float>   *muPFPUIso;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nTrks;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleSigmaIEtaIEta_2012;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_eleBC1E;   //!
   TBranch        *b_eleBC1Eta;   //!
   TBranch        *b_eleBC2E;   //!
   TBranch        *b_eleBC2Eta;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEta_2012;   //!
   TBranch        *b_phoSigmaIEtaIPhi_2012;   //!
   TBranch        *b_phoSigmaIPhiIPhi_2012;   //!
   TBranch        *b_phoE1x3_2012;   //!
   TBranch        *b_phoE2x2_2012;   //!
   TBranch        *b_phoE2x5Max_2012;   //!
   TBranch        *b_phoE5x5_2012;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFChIsoFrix1;   //!
   TBranch        *b_phoPFChIsoFrix2;   //!
   TBranch        *b_phoPFChIsoFrix3;   //!
   TBranch        *b_phoPFChIsoFrix4;   //!
   TBranch        *b_phoPFChIsoFrix5;   //!
   TBranch        *b_phoPFChIsoFrix6;   //!
   TBranch        *b_phoPFChIsoFrix7;   //!
   TBranch        *b_phoPFChIsoFrix8;   //!
   TBranch        *b_phoPFPhoIsoFrix1;   //!
   TBranch        *b_phoPFPhoIsoFrix2;   //!
   TBranch        *b_phoPFPhoIsoFrix3;   //!
   TBranch        *b_phoPFPhoIsoFrix4;   //!
   TBranch        *b_phoPFPhoIsoFrix5;   //!
   TBranch        *b_phoPFPhoIsoFrix6;   //!
   TBranch        *b_phoPFPhoIsoFrix7;   //!
   TBranch        *b_phoPFPhoIsoFrix8;   //!
   TBranch        *b_phoPFNeuIsoFrix1;   //!
   TBranch        *b_phoPFNeuIsoFrix2;   //!
   TBranch        *b_phoPFNeuIsoFrix3;   //!
   TBranch        *b_phoPFNeuIsoFrix4;   //!
   TBranch        *b_phoPFNeuIsoFrix5;   //!
   TBranch        *b_phoPFNeuIsoFrix6;   //!
   TBranch        *b_phoPFNeuIsoFrix7;   //!
   TBranch        *b_phoPFNeuIsoFrix8;   //!
   TBranch        *b_phoBC1E;   //!
   TBranch        *b_phoBC1Eta;   //!
   TBranch        *b_phoBC2E;   //!
   TBranch        *b_phoBC2Eta;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIsGood;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!

   monoPhotonAnalyzer(TTree *tree=0);
   virtual ~monoPhotonAnalyzer();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   bool             HasMediumPhoton(int& photonNo);
   bool             electronVeto(const int photonNo);
   bool             muonVeto(const int photonNo);
   float            deltaPhi(float phi1, float phi2);

   // options
   bool kUseWorstChIso = true;
};

#endif
