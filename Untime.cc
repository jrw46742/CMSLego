// -*- C++ -*-
//
// Package:    Askew/Untime
// Class:      Untime
// 
/**\class Untime Untime.cc Askew/Untime/plugins/Untime.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Askew
//         Created:  Thu, 17 Sep 2015 15:25:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrack.h"
//​DataFormats/​ParticleFlowCandidate/​interface/​PFCandidate.h 
#include "TH1.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//
// class declaration
//
#include "TTree.h"
#include "TFile.h"
#include <map>
#include <math.h>
#include <TMath.h>
float dRCalc(float etaLead, float phiLead, float etaTrail, float phiTrail){
  
  float dphi = fabs(phiLead - phiTrail);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  float deta = fabs(etaLead - etaTrail);
  float dR = sqrt(deta*deta + dphi*dphi);
  return dR;
  
}

class Untime : public edm::EDAnalyzer {
public:
  explicit Untime(const edm::ParameterSet&);
  ~Untime();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  edm::EDGetTokenT<reco::PFMETCollection> PFMET_token;
  edm::EDGetTokenT<reco::PhotonCollection> photons_old;
  edm::EDGetTokenT<reco::VertexCollection> _vtxTag;
  edm::EDGetTokenT<double> _rhoTag;
  edm::EDGetTokenT<reco::TrackCollection> _trackTag;
  edm::EDGetTokenT<EcalRecHitCollection> _ecalHits;
  edm::EDGetTokenT<HBHERecHitCollection> _hcalHits;
  edm::EDGetTokenT<reco::PFJetCollection> pfjettoken;
  edm::EDGetTokenT<std::vector<reco::Muon> > muToken;
  
  edm::EDGetTokenT<reco::PFCandidateCollection> pfcandtoken;
  edm::EDGetTokenT<std::vector<reco::PFBlock> > pfblock;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<GenEventInfoProduct> genEvToken;
  edm::EDGetTokenT<LHEEventProduct> LHEToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > mcpartToken;
  
  TFile*  rootFile_;
  TTree *tinyTree;
  Int_t Run;
  Long64_t Event;
  Int_t LumiSec;
  Float_t rho;
  

  Float_t GENpdf[150];
  Float_t GENpthat;
  Float_t GENprocessID;
  Float_t GENweight;
  Float_t GENht;
  Float_t GENpdfweight;
  Float_t GENpdfSysweight[150];
  
  Int_t nPUInfo;
  Int_t nPU[200];
  Int_t puBX[200];
  Float_t puTrue[200];
  
  Int_t nMC;
  Int_t MCPid[10000];
  Int_t MCCharge[10000];
  Float_t MCVtx_x[10000];
  Float_t MCVtx_y[10000];
  Float_t MCVtx_z[10000];
  Float_t MCPt[10000];
  Float_t MCMass[10000];
  Float_t MCPhi[10000];
  Float_t MCEta[10000];
  Float_t MCE[10000];
  Int_t MCGMomPID[10000];
  Int_t MCMomPID[10000];
  Float_t MCMomPt[10000];
  Float_t MCMomPhi[10000];
  Float_t MCMomEta[10000];
  Int_t MCStatFlag[10000];
  Int_t MCParentage[10000];
  Int_t MCStatus[10000];


  Float_t PFMET;
  Float_t PFMETPhi;
  Float_t PFMET_old;
  Float_t PFMETPhi_old;
  


  // Int_t nAllCellsEE;
  // Float_t AllCellsEtaEE[30000];
  // Float_t AllCellsPhiEE[30000];
  // Int_t AllCellsIXEE[30000];
  // Int_t AllCellsIYEE[30000];
  // Float_t AllCellsEE_X[30000];
  // Float_t AllCellsEE_Y[30000];
  // Float_t AllCellsEE_Z[30000];
  // Float_t AllCellsE_EE[30000];
  // Float_t AllCellsEt_EE[30000];
  // Int_t AllClusteredEE[30000];
  // Int_t AllClusteredEEpf[30000];
  // Int_t AllFlagEE[30000];
  // Float_t AllTimeEE[30000];
  // Float_t AllTimeErrEE[30000]; 
  // Int_t AllSideEE[30000];
  

  Int_t nAllCellsEB;
  Float_t AllCellsEtaEB[30000];
  Float_t AllCellsPhiEB[30000];
  Int_t AllCellsIEtaEB[30000];
  Int_t AllCellsIPhiEB[30000];
  Float_t AllCellsEB_X[30000];
  Float_t AllCellsEB_Y[30000];
  Float_t AllCellsEB_Z[30000];
  Float_t AllCellsE_EB[30000];
  Float_t AllCellsEt_EB[30000];
  Int_t AllClusteredEB[30000];
  Int_t AllClusteredEBUt[30000];
  Int_t AllFlagEB[30000];
  Float_t AllTimeEB[30000];
  Float_t AllTimeErrEB[30000]; 
  
  Int_t nAllCellsHE;
  Int_t AllCellsIEtaHE[30000];
  Int_t AllCellsIPhiHE[30000];
  Float_t AllCellsE_HE[30000];
  Int_t AllCellsIETAHI[30000];
  Int_t AllCellsIETALO[30000];
  Int_t AllCellsIPHIHI[30000];
  Int_t AllCellsIPHILO[30000];
  // // Int_t nPFPhoton;
  // // Float_t PFPhotonE[100];
  // // Float_t PFPhotonET[100];
  // // Float_t PFPhotonHadEM[100];
  // // Float_t PFPhotonEta[100];
  // // Float_t PFPhotonPhi[100];
  // // Float_t PFPhotonSigmaIetaIeta[100];
  // // Float_t PFPhotonR9[100];
  // // Float_t PFPhotonZ[100];
  // // Float_t PFChargedHadronIso[100];
  // // Float_t PFNeutralHadronIso[100];
  // // Float_t PFPhotonIso[100];
  // // Float_t PFMIPTotE[100];
  // // Int_t PFMIPNhit[100];
  // // Float_t PFEcalIso03[100];
  // // Float_t PFHcalIso03[100];
  // // Float_t PFTrkIso03[100];
  // // Int_t PFPixMatch[100];

  Int_t nucPFPhoton;
  Float_t ucPFPhotonE[100];
  Float_t ucPFPhotonET[100];
  Float_t ucPFPhotonHadEM[100];
  Float_t ucPFPhotonEta[100];
  Float_t ucPFPhotonPhi[100];
  Float_t ucPFPhotonSigmaIetaIeta[100];
  Float_t ucPFPhotonSigmaIphiIphi[100];
  Float_t ucPFPhotonR9[100];
  Float_t ucPFPhotonZ[100];
  Float_t ucPFChargedHadronIso[100];
  Float_t ucPFNeutralHadronIso[100];
  Float_t ucPFPhotonIso[100];
  Float_t ucPFMIPTotE[100];
  Int_t ucPFMIPNhit[100];
  Float_t ucPFEcalIso03[100];
  Float_t ucPFHcalIso03[100];
  Float_t ucPFTrkIso03[100];
  Int_t ucPFPixMatch[100];
  
  Int_t nPFJet;
  Float_t PFJetE[500];
  Float_t PFJetEta[500];
  Float_t PFJetPhi[500];
  Float_t PFJetZ[500];
  Float_t PFNeutralFraction[500];
  Float_t PFChargedFraction[500];
  Float_t PFPhotonFraction[500];
  Float_t PFElectronFraction[500];
  Float_t PFMuonFraction[500];
  Float_t PFTrackDxy[500];
  Float_t PFTrackDsz[500];
  Int_t PFJetNCon[500];

  Int_t Vertex_n;
  Float_t Vertex_x[200];
  Float_t Vertex_y[200];
  Float_t Vertex_z[200];
  Int_t Vertex_tracksize[200];
  Int_t Vertex_ndof[200];
  Float_t Vertex_chi2[200];
  Float_t Vertex_d0[200];
  Bool_t Vertex_isFake[200];

  Int_t nAllTracks;
  Float_t AllTracksPt[10000];
  Float_t AllTracksEta[10000];
  Float_t AllTracksPhi[10000];
  Int_t AllTracksNhit[10000];
  Float_t AllTracksChi2[10000];
  Float_t AllTracksPtErr[10000];
  Float_t AllTracksDxy[10000];
  Float_t AllTracksZ0[10000];
  Int_t AllTracksQ[10000];


  Int_t nPFBlock;
  Int_t nBlockEle[1000];
  Int_t nBlockEleN[1000][100];
  Float_t BlockElePt[1000][100];
  Float_t BlockElePhi[1000][100];
  Float_t BlockEleEta[1000][100];
  Int_t BlockEleType[1000][100];

  Int_t nPFClusECALEB;
  Float_t PFClusECAL_Eta[1000];
  Float_t PFClusECAL_Phi[1000];
  Float_t PFClusECAL_E[1000];
  Int_t PFClusECAL_Nhit[1000];
  Float_t PFClusECAL_Frac[1000][100];
  Int_t PFClusECAL_iEta[1000][100];
  Int_t PFClusECAL_iPhi[1000][100];
  Float_t PFClusECAL_hitE[1000][100];

  Int_t nPFClusHCALHB;
  Float_t PFClusHCAL_Eta[1000];
  Float_t PFClusHCAL_Phi[1000];
  Float_t PFClusHCAL_E[1000];
  Int_t PFClusHCAL_Nhit[1000];
  Float_t PFClusHCAL_Frac[1000][100];
  Int_t PFClusHCAL_iEta[1000][100];
  Int_t PFClusHCAL_iPhi[1000][100];
  Float_t PFClusHCAL_hitE[1000][100];


  Int_t nPFCand;
  Int_t PFid[15000];
  Float_t PFCandE[15000];
  Float_t PFCandEta[15000];
  Float_t PFCandPhi[15000];
  Float_t PFCandPt[15000];
  Float_t PFCandPositionZ[15000];
  Float_t PFCandCharge[15000];
  Float_t PFCandHcalE[15000];
  Float_t PFCandEcalE[15000];

  Int_t     nItems[15000];

  Int_t   PFBlock1[15000];
  Int_t PFElement1[15000];
  Int_t PFEleType1[15000];
  Float_t PFItem1Eta[15000];
  Float_t PFItem1Phi[15000];
  Float_t PFItem1Pt[15000];
  Float_t PFItem1EtaEcal[15000];
  Float_t PFItem1PhiEcal[15000];

  Int_t   PFBlock2[15000];
  Int_t PFElement2[15000];
  Int_t PFEleType2[15000];
  Float_t PFItem2Eta[15000];
  Float_t PFItem2Phi[15000];
  Float_t PFItem2Pt[15000];
  Float_t PFItem2EtaEcal[15000];
  Float_t PFItem2PhiEcal[15000];

  Int_t   PFBlock3[15000];
  Int_t PFElement3[15000];
  Int_t PFEleType3[15000];
  Float_t PFItem3Eta[15000];
  Float_t PFItem3Phi[15000];
  Float_t PFItem3Pt[15000];
  Float_t PFItem3EtaEcal[15000];
  Float_t PFItem3PhiEcal[15000];


  Int_t   PFBlock4[15000];  
  Int_t PFElement4[15000];
  Int_t PFEleType4[15000];
  Float_t PFItem4Eta[15000];
  Float_t PFItem4Phi[15000];
  Float_t PFItem4Pt[15000];
  Float_t PFItem4EtaEcal[15000];
  Float_t PFItem4PhiEcal[15000];


  Int_t   PFBlock5[15000];
  Int_t PFElement5[15000];
  Int_t PFEleType5[15000];
  Float_t PFItem5Eta[15000];
  Float_t PFItem5Phi[15000];
  Float_t PFItem5Pt[15000];
  Float_t PFItem5EtaEcal[15000];
  Float_t PFItem5PhiEcal[15000];
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Untime::Untime(const edm::ParameterSet& iConfig)
{

  pfjettoken = consumes<reco::PFJetCollection>(edm::InputTag("ak4PFJets"));
  pfcandtoken = consumes<reco::PFCandidateCollection>(edm::InputTag("particleFlow"));
  pfblock = consumes<std::vector<reco::PFBlock> > (edm::InputTag("particleFlowBlock"));
  photons_old = consumes<reco::PhotonCollection>(edm::InputTag("gedPhotons"));
  _ecalHits = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit:EcalRecHitsEB"));
  _hcalHits = consumes<HBHERecHitCollection>( edm::InputTag("hbhereco::RECO"));
  _rhoTag =consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  _vtxTag =consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  _trackTag = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
  PFMET_token =consumes<reco::PFMETCollection>(edm::InputTag("pfMet"));
  muToken =consumes<std::vector<reco::Muon> > (edm::InputTag("muons"));
  
  
  vtxToken = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  
  genEvToken = consumes<GenEventInfoProduct> (edm::InputTag("generator"));
  LHEToken = consumes<LHEEventProduct> (edm::InputTag("externalLHEProducer"));
  puToken = consumes<std::vector<PileupSummaryInfo> > (edm::InputTag("addPileupInfo"));
  mcpartToken = consumes<std::vector<reco::GenParticle> > (edm::InputTag("genParticles"));
  
  
  rootFile_ = new TFile("EvtDump.root","RECREATE");
  //rootFile_ = new TFile("check.root","RECREATE");
  rootFile_->cd();
  tinyTree = new TTree("tinyTree","Tiny Photon Tree");
  tinyTree->Branch("Run", &Run, "Run/I");
  tinyTree->Branch("Event", &Event, "Event/L");
  tinyTree->Branch("LumiSec", &LumiSec, "LumiSec/I");
  tinyTree->Branch("rho",&rho,"rho/F");
  
  tinyTree->Branch("GENpdf", GENpdf, "GENpdf[150]/F");
  tinyTree->Branch("GENpthat", &GENpthat,"GENpthat/F");
  tinyTree->Branch("GENprocessID", &GENprocessID, "GENprocessID/F");
  tinyTree->Branch("GENweight", &GENweight, "GENweight/F");
  tinyTree->Branch("GENht", &GENht, "GENht/F");
  tinyTree->Branch("GENpdfweight", &GENpdfweight, "GENpdfweight/F");
  tinyTree->Branch("GENpdfSysweght", GENpdfSysweight, "GENpdfSysweight[150]/F");

  tinyTree->Branch("nPUInfo", &nPUInfo, "nPUInfo/I");
  tinyTree->Branch("nPU", nPU, "nPU[nPUInfo]/I");
  tinyTree->Branch("puBX", puBX, "puBX[nPUInfo]/I");
  tinyTree->Branch("puTrue", puTrue, "puTrue[nPUInfo]/F");

  tinyTree->Branch("nMC", &nMC, "nMC/I");
  tinyTree->Branch("MCPid", MCPid, "MCPid[nMC]/I");
  tinyTree->Branch("MCCharge", MCCharge, "MCCharge[nMC]/I");
  tinyTree->Branch("MCVtx_x", MCVtx_x, "MCVtx_x[nMC]/F");
  tinyTree->Branch("MCVtx_y", MCVtx_y, "MCVtx_y[nMC]/F");
  tinyTree->Branch("MCVtx_z", MCVtx_z, "MCVtx_z[nMC]/F");
  tinyTree->Branch("MCpt", MCPt, "MCPt[nMC]/F");
  tinyTree->Branch("MCMass", MCMass, "MCMass[nMC]/F");
  tinyTree->Branch("MCPhi", MCPhi, "MCPhi[nMC]/F");
  tinyTree->Branch("MCEta", MCEta, "MCEta[nMC]/F");
  tinyTree->Branch("MCE", MCE, "MCE[nMC]/F");
  tinyTree->Branch("MCGMomPID", MCGMomPID, "MCGMomPID[nMC]/I");
  tinyTree->Branch("MCMomPID", MCMomPID, "MCMomPID[nMC]/I");
  tinyTree->Branch("MCMomPt", MCMomPt, "MCMomPt[nMC]/F");
  tinyTree->Branch("MCMomPhi", MCMomPhi, "MCMomPhi[nMC]/F");
  tinyTree->Branch("MCMomEta", MCMomEta, "MCMomEta[nMC]/F");
  tinyTree->Branch("MCStatFlag", MCStatFlag, "MCStatFlag[nMC]/I");
  tinyTree->Branch("MCParentage", MCParentage, "MCParentage[nMC]/I");
  tinyTree->Branch("MCStatus", MCStatus, "MCStatus[nMC]/I");

  tinyTree->Branch("Vertex_n", &Vertex_n, "Vertex_n/I");
  tinyTree->Branch("Vertex_x", Vertex_x, "Vertex_x[Vertex_n]/F");
  tinyTree->Branch("Vertex_y", Vertex_y, "Vertex_y[Vertex_n]/F");
  tinyTree->Branch("Vertex_z", Vertex_z, "Vertex_z[Vertex_n]/F");
  tinyTree->Branch("Vertex_tracksize", Vertex_tracksize, "Vertex_tracksize[Vertex_n]/I");
  tinyTree->Branch("Vertex_ndof", Vertex_ndof, "Vertex_ndof[Vertex_n]/I");
  tinyTree->Branch("Vertex_chi2", Vertex_chi2, "Vertex_chi2[Vertex_n]/F");
  tinyTree->Branch("Vertex_d0", Vertex_d0, "Vertex_d0[Vertex_n]/F");
  tinyTree->Branch("Vertex_isFake", Vertex_isFake, "Vertex_isFake[Vertex_n]/O");

  tinyTree->Branch("nAllTracks", &nAllTracks,"nAllTracks/I");
  tinyTree->Branch("AllTracksPt", AllTracksPt,"AllTracksPt[nAllTracks]/F");
  tinyTree->Branch("AllTracksEta", AllTracksEta,"AllTracksEta[nAllTracks]/F");
  tinyTree->Branch("AllTracksPhi", AllTracksPhi,"AllTracksPhi[nAllTracks]/F");
  tinyTree->Branch("AllTracksNhit", AllTracksNhit,"AllTracksNhit[nAllTracks]/I");
  tinyTree->Branch("AllTracksChi2", AllTracksChi2, "AllTracksChi2[nAllTracks]/F");
  tinyTree->Branch("AllTracksPtErr", AllTracksPtErr,"AllTracksPtErr[nAllTracks]/F");
  tinyTree->Branch("AllTracksDxy", AllTracksDxy, "AllTracksDxy[nAllTracks]/F");
  tinyTree->Branch("AllTracksZ0", AllTracksZ0, "AllTracksZ0[nAllTracks]/F");
  tinyTree->Branch("AllTracksQ", AllTracksQ, "AllTracksQ[nAllTracks]/I");

  tinyTree->Branch("nPFCand", &nPFCand, "nPFCand/I");
  tinyTree->Branch("PFid", PFid, "PFid[nPFCand]/I");
  tinyTree->Branch("PFCandE", PFCandE, "PFCandE[nPFCand]/F");
  tinyTree->Branch("PFCandEta",PFCandEta, "PFCandEta[nPFCand]/F");
  tinyTree->Branch("PFCandPhi", PFCandPhi, "PFCandPhi[nPFCand]/F");
  tinyTree->Branch("PFCandPt", PFCandPt, "PFCandPt[nPFCand]/F");
  tinyTree->Branch("PFCandPositionZ", PFCandPositionZ, "PFCandPositionZ[nPFCand]/F");
  tinyTree->Branch("PFCandCharge", PFCandCharge, "PFCandCharge[nPFCand]/F");
  tinyTree->Branch("PFCandHcalE", PFCandHcalE, "PFCandHcalE[nPFCand]/F");
  tinyTree->Branch("PFCandEcalE", PFCandEcalE, "PFCandEcalE[nPFCand]/F");
  
  tinyTree->Branch("nItems", nItems, "nItems[nPFCand]/I");
  tinyTree->Branch("PFBlock1", PFBlock1, "PFBlock1[nPFCand]/I");
  tinyTree->Branch("PFElement1", PFElement1, "PFElement1[nPFCand]/I");
  tinyTree->Branch("PFEleType1", PFEleType1, "PFEleType1[nPFCand]/I");
  tinyTree->Branch("PFItem1Eta", PFItem1Eta, "PFItem1Eta[nPFCand]/F");;
  tinyTree->Branch("PFItem1Phi", PFItem1Phi, "PFItem1Phi[nPFCand]/F");;
  tinyTree->Branch("PFItem1Pt", PFItem1Pt, "PFItem1Pt[nPFCand]/F");;
  tinyTree->Branch("PFItem1EtaEcal", PFItem1EtaEcal, "PFItem1EtaEcal[nPFCand]/F");
  tinyTree->Branch("PFItem1PhiEcal", PFItem1PhiEcal, "PFItem1PhiEcal[nPFCand]/F");

  tinyTree->Branch("PFBlock2", PFBlock2, "PFBlock2[nPFCand]/I");
  tinyTree->Branch("PFElement2", PFElement2, "PFElement2[nPFCand]/I");
  tinyTree->Branch("PFEleType2", PFEleType2, "PFEleType2[nPFCand]/I");
  tinyTree->Branch("PFItem2Eta", PFItem2Eta, "PFItem2Eta[nPFCand]/F");;
  tinyTree->Branch("PFItem2Phi", PFItem2Phi, "PFItem2Phi[nPFCand]/F");;
  tinyTree->Branch("PFItem2Pt", PFItem2Pt, "PFItem2Pt[nPFCand]/F");;
  tinyTree->Branch("PFItem2EtaEcal", PFItem2EtaEcal, "PFItem2EtaEcal[nPFCand]/F");
  tinyTree->Branch("PFItem2PhiEcal", PFItem2PhiEcal, "PFItem2PhiEcal[nPFCand]/F");

  tinyTree->Branch("PFBlock3", PFBlock3, "PFBlock3[nPFCand]/I");
  tinyTree->Branch("PFElement3", PFElement3,"PFElement3[nPFCand]/I");
  tinyTree->Branch("PFEleType3", PFEleType3, "PFEleType3[nPFCand]/I");
  tinyTree->Branch("PFItem3Eta", PFItem3Eta, "PFItem3Eta[nPFCand]/F");;
  tinyTree->Branch("PFItem3Phi", PFItem3Phi, "PFItem3Phi[nPFCand]/F");;
  tinyTree->Branch("PFItem3Pt", PFItem3Pt, "PFItem3Pt[nPFCand]/F");;
  tinyTree->Branch("PFItem3EtaEcal", PFItem3EtaEcal, "PFItem3EtaEcal[nPFCand]/F");
  tinyTree->Branch("PFItem3PhiEcal", PFItem3PhiEcal, "PFItem3PhiEcal[nPFCand]/F");

  tinyTree->Branch("PFBlock4", PFBlock4,"PFBlock4[nPFCand]/I");
  tinyTree->Branch("PFElement4", PFElement4,"PFElement4[nPFCand]/I");
  tinyTree->Branch("PFEleType4", PFEleType4, "PFEleType4[nPFCand]/I");
  tinyTree->Branch("PFItem4Eta", PFItem4Eta, "PFItem4Eta[nPFCand]/F");;
  tinyTree->Branch("PFItem4Phi", PFItem4Phi, "PFItem4Phi[nPFCand]/F");;
  tinyTree->Branch("PFItem4Pt", PFItem4Pt, "PFItem4Pt[nPFCand]/F");;
  tinyTree->Branch("PFItem4EtaEcal", PFItem4EtaEcal, "PFItem4EtaEcal[nPFCand]/F");
  tinyTree->Branch("PFItem4PhiEcal", PFItem4PhiEcal, "PFItem4PhiEcal[nPFCand]/F");

  tinyTree->Branch("PFBlock5", PFBlock5, "PFBlock5[nPFCand]/I");
  tinyTree->Branch("PFElement5", PFElement5,"PFElement5[nPFCand]/I");
  tinyTree->Branch("PFEleType5", PFEleType5, "PFEleType5[nPFCand]/I");
  tinyTree->Branch("PFItem5Eta", PFItem5Eta, "PFItem5Eta[nPFCand]/F");;
  tinyTree->Branch("PFItem5Phi", PFItem5Phi, "PFItem5Phi[nPFCand]/F");;
  tinyTree->Branch("PFItem5Pt", PFItem5Pt, "PFItem5Pt[nPFCand]/F");;
  tinyTree->Branch("PFItem5EtaEcal", PFItem5EtaEcal, "PFItem5EtaEcal[nPFCand]/F");
  tinyTree->Branch("PFItem5PhiEcal", PFItem5PhiEcal, "PFItem5PhiEcal[nPFCand]/F");

  tinyTree->Branch("nPFBlock", &nPFBlock, "nPFBlock/I");
  tinyTree->Branch("nBlockEle", nBlockEle, "nBlockEle[nPFBlock]/I");
  tinyTree->Branch("nBlockEleN", nBlockEleN,"nBlockEleN[nPFBlock][100]/I");
  tinyTree->Branch("BlockElePt", BlockElePt, "BlockElePt[nPFBlock][100]/F");
  tinyTree->Branch("BlockElePhi", BlockElePhi,"BlockElePhi[nPFBlock][100]/F");
  tinyTree->Branch("BlockEleEta", BlockEleEta,"BlockEleEta[nPFBlock][100]/F");
  tinyTree->Branch("BlockEleType", BlockEleType, "BlockEleType[nPFBlock][100]/I");


  tinyTree->Branch("nPFClusECALEB", &nPFClusECALEB, "nPFClusECALEB/I");
  tinyTree->Branch("PFClusECAL_Eta", PFClusECAL_Eta, "PFClusECAL_Eta[nPFClusECALEB]/F");
  tinyTree->Branch("PFClusECAL_Phi", PFClusECAL_Phi, "PFClusECAL_Phi[nPFClusECALEB]/F");;
  tinyTree->Branch("PFClusECAL_E", PFClusECAL_E, "PFClusECAL_E[nPFClusECALEB]/F");
  tinyTree->Branch("PFClusECAL_Nhit", PFClusECAL_Nhit, "PFClusECAL_Nhit[nPFClusECALEB]/I");
  tinyTree->Branch("PFClusECAL_Frac", PFClusECAL_Frac, "PFClusECAL_Frac[nPFClusECALEB][100]/F");
  tinyTree->Branch("PFClusECAL_iEta", PFClusECAL_iEta, "PFClusECAL_iEta[nPFClusECALEB][100]/I");
  tinyTree->Branch("PFClusECAL_iPhi", PFClusECAL_iPhi, "PFClusECAL_iPhi[nPFClusECALEB][100]/I");
  tinyTree->Branch("PFClusECAL_hitE", PFClusECAL_hitE, "PFClusECAL_hitE[nPFClusECALEB][100]/F");

  tinyTree->Branch("nPFClusHCALHB", &nPFClusHCALHB, "nPFClusHCALHB/I");
  tinyTree->Branch("PFClusHCAL_Eta", PFClusHCAL_Eta, "PFClusHCAL_Eta[nPFClusHCALHB]/F");
  tinyTree->Branch("PFClusHCAL_Phi", PFClusHCAL_Phi, "PFClusHCAL_Phi[nPFClusHCALHB]/F");;
  tinyTree->Branch("PFClusHCAL_E", PFClusHCAL_E, "PFClusHCAL_E[nPFClusHCALHB]/F");
  tinyTree->Branch("PFClusHCAL_Nhit", PFClusHCAL_Nhit, "PFClusHCAL_Nhit[nPFClusHCALHB]/I");
  tinyTree->Branch("PFClusHCAL_Frac", PFClusHCAL_Frac, "PFClusHCAL_Frac[nPFClusHCALHB][100]/F");
  tinyTree->Branch("PFClusHCAL_iEta", PFClusHCAL_iEta, "PFClusHCAL_iEta[nPFClusHCALHB][100]/I");
  tinyTree->Branch("PFClusHCAL_iPhi", PFClusHCAL_iPhi, "PFClusHCAL_iPhi[nPFClusHCALHB][100]/I");
  tinyTree->Branch("PFClusHCAL_hitE", PFClusHCAL_hitE, "PFClusHCAL_hitE[nPFClusHCALHB][100]/F");



  tinyTree->Branch("PFMET", &PFMET,"PFMET/F");
  tinyTree->Branch("PFMETPhi", &PFMETPhi, "PFMETPhi/F");


  tinyTree->Branch("nAllCellsEB", &nAllCellsEB, "nAllCellsEB/I");
  tinyTree->Branch("AllCellsEtaEB", AllCellsEtaEB, "AllCellsEtaEB[nAllCellsEB]/F");
  tinyTree->Branch("AllCellsPhiEB", AllCellsPhiEB, "AllCellsPhiEB[nAllCellsEB]/F");
  tinyTree->Branch("AllCellsIEtaEB", AllCellsIEtaEB, "AllCellsIEtaEB[nAllCellsEB]/I");
  tinyTree->Branch("AllCellsIPhiEB", AllCellsIPhiEB, "AllCellsIPhiEB[nAllCellsEB]/I");  
  tinyTree->Branch("AllCellsE_EB", AllCellsE_EB, "AllCellsE_EB[nAllCellsEB]/F");
  tinyTree->Branch("AllCellsEt_EB", AllCellsEt_EB, "AllCellsEt_EB[nAllCellsEB]/F");
  tinyTree->Branch("AllClusteredEB", AllClusteredEB, "AllClusteredEB[nAllCellsEB]/I");
  tinyTree->Branch("AllClusteredEBUt", AllClusteredEBUt, "AllClusteredEBUt[nAllCellsEB]/I");
  tinyTree->Branch("AllCellsEB_X", AllCellsEB_X, "AllCellsEB_X[nAllCellsEB]/F");
  tinyTree->Branch("AllCellsEB_Y", AllCellsEB_Y, "AllCellsEB_Y[nAllCellsEB]/F");
  tinyTree->Branch("AllCellsEB_Z", AllCellsEB_Z, "AllCellsEB_Z[nAllCellsEB]/F");

  tinyTree->Branch("AllFlagEB", AllFlagEB,"AllFlagEB[nAllCellsEB]/I");
  tinyTree->Branch("AllTimeEB", AllTimeEB,"AllTimeEB[nAllCellsEB]/F");
  tinyTree->Branch("AllTimeErrEB", AllTimeErrEB, "AllTimeErrEB[nAllCellsEB]/F");

  tinyTree->Branch("nAllCellsHE", &nAllCellsHE, "nAllCellsHE/I");
  tinyTree->Branch("AllCellsIEtaHE",AllCellsIEtaHE, "AllCellsIEtaHE[nAllCellsHE]/I");
  tinyTree->Branch("AllCellsIPhiHE", AllCellsIPhiHE, "AllCellsIPhiHE[nAllCellsHE]/I");
  tinyTree->Branch("AllCellsE_HE", AllCellsE_HE,"AllCellsE_HE[nAllCellsHE]/F");
  tinyTree->Branch("AllCellsIETAHI", AllCellsIETAHI, "AllCellsIETAHI[nAllCellsHE]/I");
  tinyTree->Branch("AllCellsIETALO", AllCellsIETALO, "AllCellsIETALO[nAllCellsHE]/I");
  tinyTree->Branch("AllCellsIPHIHI", AllCellsIPHIHI, "AllCellsIPHIHI[nAllCellsHE]/I");
  tinyTree->Branch("AllCellsIPHILO", AllCellsIPHILO, "AllCellsIPHILO[nAllCellsHE]/I");
 

  tinyTree->Branch("nucPFPhoton", &nucPFPhoton, "nucPFPhoton/I");
  tinyTree->Branch("ucPFPhotonE", ucPFPhotonE, "ucPFPhotonE[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonHadEM", ucPFPhotonHadEM, "ucPFPhotonHadEM[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonEta", ucPFPhotonEta, "ucPFPhotonEta[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonPhi",ucPFPhotonPhi, "ucPFPhotonPhi[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonSigmaIetaIeta", ucPFPhotonSigmaIetaIeta,"ucPFPhotonSigmaIetaIeta[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonSigmaIphiIphi", ucPFPhotonSigmaIphiIphi,"ucPFPhotonSigmaIphiIphi[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonR9", ucPFPhotonR9, "ucPFPhotonR9[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonZ", ucPFPhotonZ, "ucPFPhotonZ[nucPFPhoton]/F");
  tinyTree->Branch("ucPFChargedHadronIso", ucPFChargedHadronIso, "ucPFChargedHadronIso[nucPFPhoton]/F");
  tinyTree->Branch("ucPFNeutralHadronIso", ucPFNeutralHadronIso, "ucPFNeutralHadronIso[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonIso", ucPFPhotonIso, "ucPFPhotonIso[nucPFPhoton]/F");

  tinyTree->Branch("ucPFMIPTotE", ucPFMIPTotE, "ucPFMIPTotE[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPhotonET", ucPFPhotonET, "ucPFPhotonET[nucPFPhoton]/F");
  tinyTree->Branch("ucPFMIPNhit", ucPFMIPNhit, "ucPFMIPNhit[nucPFPhoton]/I");
  tinyTree->Branch("ucPFEcalIso03", ucPFEcalIso03, "ucPFEcalIso03[nucPFPhoton]/F");
  tinyTree->Branch("ucPFHcalIso03", ucPFHcalIso03, "ucPFHcalIso03[nucPFPhoton]/F");
  tinyTree->Branch("ucPFTrkIso03", ucPFTrkIso03, "ucPFTrkIso03[nucPFPhoton]/F");
  tinyTree->Branch("ucPFPixMatch", ucPFPixMatch, "ucPFPixMatch[nucPFPhoton]/I");


  tinyTree->Branch("nPFJet", &nPFJet, "nPFJet/I");
  tinyTree->Branch("PFJetE", PFJetE, "PFJetE[nPFJet]/F");
  tinyTree->Branch("PFJetEta", PFJetEta, "PFJetEta[nPFJet]/F");
  tinyTree->Branch("PFJetPhi", PFJetPhi, "PFJetPhi[nPFJet]/F");
  tinyTree->Branch("PFJetZ", PFJetZ, "PFJetZ[nPFJet]/F");
  tinyTree->Branch("PFNeutralFraction", PFNeutralFraction, "PFNeutralFraction[nPFJet]/F");
  tinyTree->Branch("PFChargedFraction", PFChargedFraction, "PFChargedFraction[nPFJet]/F");
  tinyTree->Branch("PFPhotonFraction", PFPhotonFraction, "PFPhotonFraction[nPFJet]/F");
  tinyTree->Branch("PFElectronFraction", PFElectronFraction, "PFElectronFraction[nPFJet]/F");
  tinyTree->Branch("PFMuonFraction", PFMuonFraction, "PFMuonFraction[nPFJet]/F");
  tinyTree->Branch("PFTrackDxy", PFTrackDxy, "PFTrackDxy[nPFJet]/F");
  tinyTree->Branch("PFTrackDsz", PFTrackDsz, "PFTrackDsz[nPFJet]/F");
  tinyTree->Branch("PFJetNCon", &PFJetNCon, "PFJetNCon[nPFJet]/I");



}


Untime::~Untime()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Untime::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  Run      = iEvent.id().run();
  Event    = iEvent.id().event();
  LumiSec  = iEvent.id().luminosityBlock();
  
  if(Event == 72){
    
    cout << "============== EVENT " << Event << " ====================" << endl;
    
    edm::Handle<double> rH;
    iEvent.getByToken(_rhoTag,rH);
    rho = *rH;
    
    
    nPUInfo=0;
    
    edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
    iEvent.getByToken(puToken, genPileupHandle);
    
    
    for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end()&&nPUInfo<200; ++pu) {
      //  if (pu->getBunchCrossing() == 0) {
      // hPU_->Fill(pu->getPU_NumInteractions());
      //  hPUTrue_->Fill(pu->getTrueNumInteractions());
      // }
      
      nPU[nPUInfo]= pu->getPU_NumInteractions();
      puTrue[nPUInfo]=(pu->getTrueNumInteractions());
      puBX[nPUInfo] = (pu->getBunchCrossing());
      
      nPUInfo++;
    }
    /////////// gen info //////////////
    
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByToken(genEvToken, genEventInfoHandle);
    if (genEventInfoHandle.isValid()){
      if (genEventInfoHandle->pdf()){
	GENpdf[0] = genEventInfoHandle->pdf()->id.first;
	GENpdf[1] = genEventInfoHandle->pdf()->id.second;
	GENpdf[2] = genEventInfoHandle->pdf()->x.first;
	GENpdf[3] = genEventInfoHandle->pdf()->x.second;
	GENpdf[4] = genEventInfoHandle->pdf()->xPDF.first;
	GENpdf[5] = genEventInfoHandle->pdf()->xPDF.second;
	GENpdf[6] = genEventInfoHandle->pdf()->scalePDF;
      }
      
      
      if (genEventInfoHandle->hasBinningValues()) GENpthat = genEventInfoHandle->binningValues()[0];
      
      GENprocessID= genEventInfoHandle->signalProcessID();
      GENweight = genEventInfoHandle->weight();
      
    }
    
    
    edm::Handle<vector <reco::GenParticle> > genParticlesHandle;
    iEvent.getByToken(mcpartToken, genParticlesHandle);
    cout << "-----------------------------------------" << endl;
    cout << "----------- MC INFORMATION --------------" << endl;
    nMC=0;
    for (vector< reco::GenParticle>::const_iterator ip = genParticlesHandle->begin();  ip!=genParticlesHandle->end() && nMC<10000; ++ip){
      
      // int status = ip->status();
      bool quarks = abs(ip->pdgId()) <7;
      bool heavyParticle =
				      ((    ip->pdgId()  == 23 && ip->isHardProcess()) || 
				       (abs(ip->pdgId()) == 24 && ip->isHardProcess()) || 
				       (    ip->pdgId()  == 25 && ip->isHardProcess()) ||
				       (abs(ip->pdgId()) ==  6 && ip->isHardProcess()) || 
				       (abs(ip->pdgId()) ==  5 && ip->isHardProcess()));
      bool IsLastPT = ip->pt()>1. && ip->isPromptFinalState();
      
      if(abs(ip->pdgId()) <= 5){ 
	cout << "*** Quark with (eta,phi,pt) = " << ip->eta() << ",\t" << ip->phi() << ",\t" << ip->pt() << endl;
	if(ip->eta() < 0){ 
	  cout << "    HCAL (ieta,iphi) = " << (ip->eta()/0.087)-0.5 << ",\t" << (ip->phi()/0.087)+0.5 << endl;
	  cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)-0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
	}
	else{
	  cout << "    HCAL (ieta,iphi) = " << (ip->eta()/0.087)+0.5 << ",\t" << (ip->phi()/0.087)+0.5 << endl;
	  cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)+0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
	}
      }
      if(abs(ip->pdgId()) == 11){ 
	cout << "*** ELECTRON with (eta,phi,pt) = " << ip->eta() << ",\t" << ip->phi() << ",\t" << ip->pt() << endl;
	if(ip->eta() < 0) cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)-0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
	else cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)+0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
      }
      if(abs(ip->pdgId()) == 13){ 
	cout << "*** MUON with (eta,phi,pt) = " << ip->eta() << ",\t" << ip->phi() << ",\t" << ip->pt() << endl;
	if(ip->eta() < 0) cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)-0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
	else cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)+0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
      }
      
      if(abs(ip->pdgId()) == 18){
	cout << "*** DARK MATTER with (eta,phi,pt) = " << ip->eta() << ",\t" << ip->phi() << ",\t" << ip->pt() << endl;
	if(ip->eta() < 0) cout << "    ECAL (ieta,iphi) = " << (ip->eta()/0.0175)-0.5 << ",\t" << (ip->phi()/0.0175)+10.5 << endl;
      }	
      
      if (IsLastPT || heavyParticle || quarks){
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	
	
	MCPid[nMC]=p->pdgId();
	MCCharge[nMC]=p->charge();
	MCVtx_x[nMC]=p->vx();
	MCVtx_y[nMC]=p->vy();
	MCVtx_z[nMC]=p->vz();
	MCPt[nMC]=p->pt();
	MCMass[nMC]=p->mass();
	MCPhi[nMC]=p->phi();
	MCEta[nMC]=p->eta();
	MCE[nMC]= p->energy();
	MCStatus[nMC]=ip->status();
	//     if ( abs(ip->pdgId())==6 && ip->isHardProcess()==1){
	//      if ( abs(ip->pdgId())==11){
	//	cout << "top: " << p->mass() << endl;
	//	cout << "eta: " << p->eta() << endl;
	//	cout << "phi: " << p->phi() << endl;
	//	cout << "status: " << ip->status() << endl;
	//	cout << "mom: " << p->mother()->pdgId() << endl;;
	//	cout << "isPromptDecayed: " << ip->isPromptDecayed() << endl;
	//	cout << "isHardProcess: " << ip->isHardProcess() << endl;
	//}
	Int_t tmpStatusFlag = 0;
	if (ip->fromHardProcessFinalState()) tmpStatusFlag+=1;
	if (ip->isPromptFinalState())        tmpStatusFlag+=10;
	if (ip->isHardProcess())             tmpStatusFlag+=100;
	
	// if genParticle is W or Z, check its decay type
	if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
	  for (size_t k=0; k < p->numberOfDaughters(); ++k) {
	    const reco::Candidate *dp = p->daughter(k);
	    if (abs(dp->pdgId())<=6)                             tmpStatusFlag+=1000;
	    else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12) tmpStatusFlag+=10000;
	    else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14) tmpStatusFlag+=100000;
	    else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16) tmpStatusFlag+=1000000;
	  }
	}
	MCStatFlag[nMC] = tmpStatusFlag;
	int mcGMomPID_ = -999;
	int mcMomPID_  = -999;
	float mcMomPt_    = -999.;
	//      float mcMomMass_  = -999.;
	float mcMomEta_   = -999.;
	float mcMomPhi_   = -999.;
	
	MCGMomPID[nMC] = mcGMomPID_;
	MCMomPID[nMC]=mcMomPID_;
	MCMomPt[nMC] = mcMomPt_;
	MCMomEta[nMC] = mcMomEta_;
	MCMomPhi[nMC] = mcMomPhi_;
	
	
	nMC++;
      }
    }
    
    cout << "-----------------------------------------\n\n" << endl;    

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByToken(_vtxTag, vtxH);
    reco::VertexRef primVtxRef(vtxH, 0);
    
    //cout << "Got vertex coll." <<endl;
    Vertex_n = 0;
    for(reco::VertexCollection::const_iterator v=vtxH->begin();v!=vtxH->end() && Vertex_n<200; ++v){
      
      Vertex_x[Vertex_n]     = v->x();
      Vertex_y[Vertex_n]     = v->y();
      Vertex_z[Vertex_n]     = v->z();
      Vertex_chi2[Vertex_n]   = v->chi2();
      Vertex_tracksize[Vertex_n] = v->tracksSize();
      Vertex_ndof[Vertex_n] = v->ndof();
      Vertex_isFake[Vertex_n] = v->isFake();
      Vertex_d0[Vertex_n] = v->position().rho();
      Vertex_n++;
      
    }

    

    edm::Handle<reco::TrackCollection> trackCollection;
    iEvent.getByToken(_trackTag, trackCollection);
    nAllTracks=0;
    //Fill track bank.
    for(reco::TrackCollection::const_iterator trItr = trackCollection->begin(); trItr != trackCollection->end() && nAllTracks<10000; ++trItr){
      //      if (trItr->pt()>1.){	
      AllTracksPt[nAllTracks]=trItr->pt();
      // if (trItr->pt()>25) n20Trk++;
      AllTracksPhi[nAllTracks]=trItr->phi();
      AllTracksEta[nAllTracks]=trItr->eta();
      AllTracksNhit[nAllTracks]=trItr->numberOfValidHits();
      AllTracksChi2[nAllTracks]=trItr->normalizedChi2();
      AllTracksPtErr[nAllTracks]=trItr->ptError();
      AllTracksQ[nAllTracks]=trItr->charge();
      //AllTracksDxy[nAllTracks]=trItr->dxy(math::XYZPoint(vertexBeamSpot.x0(
      //),vertexBeamSpot.y0(),vertexBeamSpot.z0()));
      //AllTracksZ0[nAllTracks]=trItr->dz(math::XYZPoint(vertexBeamSpot.x0(),
      //                                               vertexBeamSpot.y0(),vertexBeamSpot.z0()));
      AllTracksZ0[nAllTracks]=trItr->vz();
      nAllTracks++;
      
      //	if(trItr->eta() < 0 && trItr->phi() > 0 && trItr->phi() < TMath::Pi()*3/4){
      //	  cout << "Track: (eta,phi,pt) = " << trItr->eta() << ",\t" << trItr->phi() << ",\t" << trItr->pt() << endl;
      //	}
    }//Loop over tracks
    
	
     // Muons
   edm::Handle<std::vector<reco::Muon> > muonHandle;
   iEvent.getByToken(muToken, muonHandle);
   int count = 0;
   for (std::vector<reco::Muon>::const_iterator mu = muonHandle->begin(); mu != muonHandle->end(); ++mu) {
     float itracketa;
     try {
	itracketa = mu->innerTrack()->eta();
     }
     catch (...) {
	cout << "itracketa failed... breaking" << endl;
	continue;
	}
     cout << "1" << endl;
     float itrackphi = mu->innerTrack()->phi();
     cout << "2" << endl;
     float itrackpt = mu->innerTrack()->pt(); //trackRef()->pt()
     cout << "3" << endl;
     //float otracketa = mu->outerTrack()->eta();
     cout << "4" << endl;
     float otracketa;
     count++;
     try {
       otracketa = mu->outerTrack()->eta();
     }
     catch (...) {
       cout << "otracketa failed... breaking" << endl;
       continue;
     }
	cout << "5" << endl;
     float otrackphi;
     try {
	otrackphi = mu->outerTrack()->phi();
	}
     catch (...) {
	cout << "otrackphi failed... breaking" << endl;
	continue;
	}
	cout << "6" << endl;
     float otrackpt;
     try {
	otrackpt = mu->outerTrack()->pt(); //trackRef()->pt()
     }
     catch (...) { 
	cout << "otrackpt failed... breaking" << endl;
	continue;
     }
	cout << "7" << endl;
     //std::vector<MuonChamberMatch> matches = mu->matches();
     //int nmatches = matches.size();
     cout << "-----------------------------------------" << endl;
     cout << "----- MUON -----" << endl;
     cout << "Inner track (eta,phi,pt) = " << itracketa << ",\t" << itrackphi << ",\t" << itrackpt << endl;
     cout << "Onner track (eta,phi,pt) = " << otracketa << ",\t" << otrackphi << ",\t" << otrackpt << endl; 
   //  cout << "N matched stations = " << mu->numberOfMatchedStations() << ",\t N matches = " << nmatches << endl;
     // for(unsigned int imatch = 0; imatch < matches.size(); imatch++){
      //	//	DetId mudet = matches.at(imatch).id;
      //	int detector = matches.at(imatch).detector();
      //	cout << "Match " << imatch << ": detector = " << detector << endl;
      //}
   }
    
    cout << "-----------------------------------------\n\n" << endl;
    
    cout << "-----------------------------------------" << endl;
    cout << "--------------- MET ---------------------" << endl;
    Handle<reco::PFMETCollection> pfMET;
    
    iEvent.getByToken(PFMET_token, pfMET);
    PFMET = pfMET->begin()->et();
    PFMETPhi = pfMET->begin()->phi();  
    
    cout << "MET pt = " << pfMET->begin()->et() << ",\t phi = " << pfMET->begin()->phi() << endl;
    
    cout << "-----------------------------------------\n\n" << endl;
    
    Handle<reco::PFCandidateCollection> pfcandcoll;
    iEvent.getByToken(pfcandtoken, pfcandcoll);
    
    Handle<reco::PFBlockCollection> pfblockcoll;
    iEvent.getByToken(pfblock, pfblockcoll);
    nPFBlock=0;
    nPFClusECALEB=0;
    nPFClusHCALHB=0;
    
    Handle<HBHERecHitCollection> hcalHits;
    iEvent.getByToken(_hcalHits, hcalHits);
    const HBHERecHitCollection* rechitsHcal = hcalHits.product();
    Handle<EcalRecHitCollection> ecalhitsCollHEB;
    iEvent.getByToken(_ecalHits, ecalhitsCollHEB);
    const EcalRecHitCollection* rechitsCollectionEB_ = ecalhitsCollHEB.product();
    
    for (int ik=0;ik< int((*pfblockcoll).size()); ++ik){
      //if(ik != 22 && ik != 23 && ik != 82 && ik != 26 && ik != 113 && ik != 109 && ik != 55 && ik != 25){
      const edm::OwnVector< reco::PFBlockElement > elee = (*pfblockcoll)[ik].elements();  
      cout << "-----------------------------------------" << endl;      
      cout << "----------- BLOCK " << ik << " with " << elee.size() << " Elements -------------" << endl;
      
      nBlockEle[nPFBlock]= elee.size();
      for (int ij=0;ij<int(elee.size()) && ij<100;++ij){
	cout << "*** ELEMENT " << ij << " of type " << elee[ij].type() << endl;
	nBlockEleN[nPFBlock][ij]=ij;
	
	BlockEleType[nPFBlock][ij]=elee[ij].type();
	BlockElePt[nPFBlock][ij]=-999;
	BlockElePhi[nPFBlock][ij]=-999;
	BlockEleEta[nPFBlock][ij]=-999;
	
	// Check for a muon:
	//      if((elee[ij]).muonRef()->isMuon()){	
	//	cout << "*** Check for MUON ref with (pt,eta,phi) = " << (elee[ij]).muonRef()->pt() << ",\t" << (elee[ij]).muonRef()->eta() << ",\t" << (elee[ij]).muonRef()->phi() << endl;
	//      }
	
	if (elee[ij].type()==4 || elee[ij].type()==5){
	  //ECAL or HCAL, which means this is a cluster element
	  //PFClusterRef reffer = (elee[ij]).clusterRef();
	  BlockElePt[nPFBlock][ij]=(elee[ij]).clusterRef()->pt();
	  BlockEleEta[nPFBlock][ij] = (elee[ij]).clusterRef()->eta();
	  BlockElePhi[nPFBlock][ij] = (elee[ij]).clusterRef()->phi();
	  if (elee[ij].type()==4 && fabs((elee[ij]).clusterRef()->eta())<1.2 && nPFClusECALEB<1000){
	    PFClusECAL_Eta[nPFClusECALEB] = (elee[ij]).clusterRef()->eta();
	    PFClusECAL_Phi[nPFClusECALEB] = (elee[ij]).clusterRef()->phi();
	    PFClusECAL_E[nPFClusECALEB] = (elee[ij]).clusterRef()->energy();
	    //	  const vector< reco::PFRecHitFraction> PFVec = (elee[ij]).clusterRef()->recHitFractions();
	    const vector< std::pair<DetId, float> > PFVec = (elee[ij]).clusterRef()->hitsAndFractions();
	    PFClusECAL_Nhit[nPFClusECALEB]= PFVec.size();
	    cout << "****** ECAL Element with " << PFVec.size() << " hits" << ",\t eta = " << (elee[ij]).clusterRef()->eta() << ",\t phi = " << (elee[ij]).clusterRef()->phi() << ",\t E = " << (elee[ij]).clusterRef()->energy() << endl;
	    for (int kk=0;kk<PFClusECAL_Nhit[nPFClusECALEB] && kk<100;++kk){
	      PFClusECAL_Frac[nPFClusECALEB][kk]= PFVec[kk].second;
	      EBDetId det = PFVec[kk].first;
	      PFClusECAL_iEta[nPFClusECALEB][kk]= det.ieta();
	      PFClusECAL_iPhi[nPFClusECALEB][kk] = det.iphi();
	       float hitenergy = -1;
	       for (EcalRecHitCollection::const_iterator it = rechitsCollectionEB_->begin();it!=rechitsCollectionEB_->end();++it){
	       	EBDetId dit = it->detid();
	       	if(dit.ieta() != det.ieta()) continue;
	       	if(dit.iphi() != det.iphi()) continue;
	       	hitenergy = it->energy();
	       }
	       bool significant = false;
	       if(hitenergy > 10) significant = true;
	       if(!significant) cout << "********* ECAL hit with (ieta,iphi,hitE) = " << det.ieta() << ",\t" << det.iphi() << ",\t" << hitenergy << endl;
	       else cout << "********* ECAL hit with (ieta,iphi,hitE) = " << det.ieta() << ",\t" << det.iphi() << ",\t" << hitenergy << " ******** SIGNIFICANT *******" << endl;
	   //   	    PFClusECAL_hitE[nPFClusECALEB][kk]= PFVec[kk].recHitRef()->energy();
	    }
	    nPFClusECALEB++;
	  }
	  
	  
	  if (elee[ij].type()==5 && fabs((elee[ij]).clusterRef()->eta())<1.44
	      && nPFClusHCALHB<1000){
	    PFClusHCAL_Eta[nPFClusHCALHB] = (elee[ij]).clusterRef()->eta();
	    PFClusHCAL_Phi[nPFClusHCALHB] = (elee[ij]).clusterRef()->phi();
	    PFClusHCAL_E[nPFClusHCALHB] = (elee[ij]).clusterRef()->energy();
	    //	  const vector< reco::PFRecHitFraction> PFVec = (elee[ij]).clusterRef()->recHitFractions();
	    const vector< std::pair<DetId, float> > PFVec = (elee[ij]).clusterRef()->hitsAndFractions();
	    PFClusHCAL_Nhit[nPFClusHCALHB]= PFVec.size();
	    cout << "****** HCAL Element with " << PFVec.size() << " hits" << ",\t eta = " << (elee[ij]).clusterRef()->eta() << ",\t phi = " << (elee[ij]).clusterRef()->phi() << ",\t E = " << (elee[ij]).clusterRef()->energy() << endl;
	    for (int kk=0;kk<PFClusHCAL_Nhit[nPFClusHCALHB] && kk<100;++kk){
	      PFClusHCAL_Frac[nPFClusHCALHB][kk]= PFVec[kk].second;
	      HcalDetId det = PFVec[kk].first;
	      PFClusHCAL_iEta[nPFClusHCALHB][kk]= det.ieta();
	      PFClusHCAL_iPhi[nPFClusHCALHB][kk] = det.iphi();
	       float hitenergy = -1;
	       for (HBHERecHitCollection::const_iterator it=rechitsHcal->begin();it!=rechitsHcal->end();++it){
	       	HcalDetId dit = it->id();
	       	if(dit.ieta() != det.ieta()) continue;
	       	if(dit.iphi() != det.iphi()) continue;
	       	hitenergy = it->energy();
	       }
	       bool significant = false;
	       if(hitenergy > 10) significant = true;
	       if(!significant) cout << "********* HCAL hit with (ieta,iphi,hitE) = " << det.ieta() << ",\t" << det.iphi() << ",\t" << hitenergy << endl;
	       else cout << "********* HCAL hit with (ieta,iphi,hitE) = " << det.ieta() << ",\t" << det.iphi() << ",\t" << hitenergy << " ******** SIGNIFICANT *******" << endl;
	      // //	    PFClusECAL_hitE[nPFClusECALEB][kk]=PFVec[kk].recHitRef()->energy();
	    }
	    nPFClusHCALHB++;
	  }	         	  
	}
	if (elee[ij].type()==1 && fabs((elee[ij]).trackRef()->eta()) < 0.73){
	  //Track
	  //TrackRef reffer = (elee[ij]).trackRef();
	  BlockElePt[nPFBlock][ij]= (elee[ij]).trackRef()->pt();
	  BlockEleEta[nPFBlock][ij]= (elee[ij]).trackRef()->eta();
	  BlockElePhi[nPFBlock][ij]= (elee[ij]).trackRef()->phi();
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ij]);
	  //	PFBlockElementTrack trk = PFBlockElementTrack((elee[ij]).trackRefPF());
	  float etaecal = trk->positionAtECALEntrance().eta();
	  float phiecal = trk->positionAtECALEntrance().phi();
	  cout << "****** TRACK Element with (eta,phi,pt,etaEcal,phiEcal) = " << (elee[ij]).trackRef()->eta() << ",\t" << (elee[ij]).trackRef()->phi() << ",\t" 
	       << (elee[ij]).trackRef()->pt() << ",\t" << etaecal << ",\t" << phiecal << endl; 
	  if(etaecal < 0) {
	    if(phiecal > 0) cout << "    ECAL (ieta,iphi) = " << (etaecal/0.0175)-0.5 << ",\t" << (phiecal/0.0175)+10.5 << endl;
	    else cout << "    ECAL (ieta,iphi) = " << (etaecal/0.0175)-0.5 << ",\t" << ((phiecal+2*TMath::Pi())/0.0175)+10.5 << endl;
	  }else{
	    if(phiecal > 0) cout << "    ECAL (ieta,iphi) = " << (etaecal/0.0175)+0.5 << ",\t" << (phiecal/0.0175)+10.5 << endl;
	    else cout << "    ECAL (ieta,iphi) = " << (etaecal/0.0175)+0.5 << ",\t" << ((phiecal+2*TMath::Pi())/0.0175)+10.5 << endl;
	  }
	}
	// if (elee[ij].type()==7){
	// 	BlockElePt[nPFBlock][ij] = (elee[ij]).trackRefPF()->trackRef()->pt();
	// 	BlockEleEta[nPFBlock][ij] = (elee[ij]).trackRefPF()->trackRef()->eta();
	// 	BlockElePhi[nPFBlock][ij] = (elee[ij]).trackRefPF()->trackRef()->phi();
	
	
	// }
	
	
      }
      
      nPFBlock++;
    }
    //}
    cout << "-----------------------------------------\n\n" << endl;    
    
    const reco::PFCandidateCollection *pfcands = pfcandcoll.product();
    //reco::PFCandidateCollection::const_iterator pfcand = pfcands->begin();
    nPFCand = 0;
    
    for (int iii =0;iii<15000;iii++){
      nItems[iii]=-1;
      PFBlock1[iii]=-999;
      PFBlock2[iii]=-999;
      PFBlock3[iii]=-999;
      PFBlock4[iii]=-999;
      PFBlock5[iii]=-999;
      PFElement1[iii]=-999;
      PFElement2[iii]=-999;
      PFElement3[iii]=-999;
      PFElement4[iii]=-999;
      PFElement5[iii]=-999;
      PFEleType1[iii]=-999;
      PFEleType2[iii]=-999;
      PFEleType3[iii]=-999;
      PFEleType4[iii]=-999;
      PFEleType5[iii]=-999;
    }
    
    for(reco::PFCandidateCollection::const_iterator pfcand = pfcands->begin(); pfcand!=pfcands->end() && nPFCand < 15000;++pfcand){
      PFid[nPFCand] = pfcand->pdgId();
      cout << "--------------------------------------" << endl;
      cout << "------------- PF CAND ----------------" << endl;
      cout << "CANDIDATE: type = " << pfcand->pdgId() << ",\t pt = " << pfcand->pt() << endl;
      PFCandE[nPFCand] = pfcand->energy();
      PFCandEta[nPFCand] = pfcand->eta();
      PFCandPhi[nPFCand] = pfcand->phi();
      PFCandPt[nPFCand] = pfcand->pt();
      PFCandPositionZ[nPFCand] = pfcand->vertex().Z();//dz(primVtxRef->position());
      PFCandCharge[nPFCand] = pfcand->charge();
      PFCandHcalE[nPFCand] = pfcand->rawHcalEnergy();
      PFCandEcalE[nPFCand] = pfcand->rawEcalEnergy();
      //cout << "Get the elements in blocks" << endl;
      const PFCandidate::ElementsInBlocks &thisPart = pfcand->elementsInBlocks();
      //cout << "How many elements? " << endl;
      nItems[nPFCand] = thisPart.size();
      //cout << "nItems: " << nItems[nPFCand] << endl;
      if (nItems[nPFCand]==0) cout << "NO ELEMENTS !!!!!!!!!!" << endl;
      if (nItems[nPFCand]>0){
	//   cout << "Get the first element in the block" << endl;
	PFCandidate::ElementInBlock ele1 = thisPart.at(0);
	//PFBlock1[nPFCand] = ele1.first;
	PFElement1[nPFCand] = ele1.second;
	//cout << "What element in this block? " << PFElement1[nPFCand] << endl;
	PFBlockRef thisBlockVec = ele1.first;
	// cout << "What is this size?  Is it ever non 1?: " << (*thisBlockVec).size() << endl;
	const PFBlock thisBlock = (*thisBlockVec);
	const edm::OwnVector< reco::PFBlockElement > elee = (thisBlock).elements();
	PFEleType1[nPFCand] =elee[ele1.second].type();
	cout << "*** Element 0, type = " << elee[ele1.second].type() << endl;
	if (PFEleType1[nPFCand]==1){
	  
	  PFItem1Eta[nPFCand] = elee[ele1.second].trackRef()->eta();
	  PFItem1Phi[nPFCand] = elee[ele1.second].trackRef()->phi();
	  PFItem1Pt[nPFCand] = elee[ele1.second].trackRef()->pt();
	  cout << "*** Element 0, pt = " << elee[ele1.second].trackRef()->pt() << endl;
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ele1.second]);
	  PFItem1EtaEcal[nPFCand]=trk->positionAtECALEntrance().eta();
	  PFItem1PhiEcal[nPFCand] = trk->positionAtECALEntrance().phi();
	  
	}
	if (PFEleType1[nPFCand]==4){
	  PFItem1Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem1Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem1Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 0, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	  
	}
	if (PFEleType1[nPFCand]==5){
	  PFItem1Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem1Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem1Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 0, pt = " << elee[ele1.second].clusterRef()->pt() << endl;

	}
	
      }
      if (nItems[nPFCand]>1){
	//   cout << "Get the first element in the block" << endl;
	PFCandidate::ElementInBlock ele1 = thisPart.at(1);
	//PFBlock2[nPFCand] = ele1.first;
	PFElement2[nPFCand] = ele1.second;
	//cout << "What element in this block? " << PFElement1[nPFCand] << endl;
	PFBlockRef thisBlockVec = ele1.first;
	// cout << "What is this size?  Is it ever non 1?: " << (*thisBlockVec).size() << endl;
	const PFBlock thisBlock = (*thisBlockVec);
	const edm::OwnVector< reco::PFBlockElement > elee = (thisBlock).elements();
	PFEleType2[nPFCand] =elee[ele1.second].type();
	cout << "*** Element 1, type = " << elee[ele1.second].type() << endl;
	if (PFEleType2[nPFCand]==1){
	  
	  PFItem2Eta[nPFCand] = elee[ele1.second].trackRef()->eta();
	  PFItem2Phi[nPFCand] = elee[ele1.second].trackRef()->phi();
	  PFItem2Pt[nPFCand] = elee[ele1.second].trackRef()->pt();
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ele1.second]);
	  PFItem2EtaEcal[nPFCand]=trk->positionAtECALEntrance().eta();
	  PFItem2PhiEcal[nPFCand] = trk->positionAtECALEntrance().phi();
	  cout << "*** Element 1, pt = " << elee[ele1.second].trackRef()->pt() << endl;
	}
	if (PFEleType2[nPFCand]==4){
	  PFItem2Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem2Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem2Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 1, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	}
	if (PFEleType2[nPFCand]==5){
	  PFItem2Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem2Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem2Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 1, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	}
	
      }
      if (nItems[nPFCand]>2){
      //   cout << "Get the first element in the block" << endl;
	PFCandidate::ElementInBlock ele1 = thisPart.at(2);
	//PFBlock3[nPFCand] = ele1.first;
	PFElement3[nPFCand] = ele1.second;
      //cout << "What element in this block? " << PFElement1[nPFCand] << endl;
	PFBlockRef thisBlockVec = ele1.first;
	// cout << "What is this size?  Is it ever non 1?: " << (*thisBlockVec).size() << endl;
	const PFBlock thisBlock = (*thisBlockVec);
	const edm::OwnVector< reco::PFBlockElement > elee = (thisBlock).elements();
	PFEleType3[nPFCand] =elee[ele1.second].type();
	cout << "*** Element 2, type = " << elee[ele1.second].type() << endl;
	if (PFEleType3[nPFCand]==1){
	  
	  PFItem3Eta[nPFCand] = elee[ele1.second].trackRef()->eta();
	  PFItem3Phi[nPFCand] = elee[ele1.second].trackRef()->phi();
	  PFItem3Pt[nPFCand] = elee[ele1.second].trackRef()->pt();
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ele1.second]);
	  PFItem3EtaEcal[nPFCand]=trk->positionAtECALEntrance().eta();
	  PFItem3PhiEcal[nPFCand] = trk->positionAtECALEntrance().phi();
	  cout << "*** Element 2, pt = " << elee[ele1.second].trackRef()->pt() << endl;
	}
	if (PFEleType3[nPFCand]==4){
	  PFItem3Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem3Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem3Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 2, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	}
	if (PFEleType3[nPFCand]==5){
	  PFItem3Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem3Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem3Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 2, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	}
      }
      if (nItems[nPFCand]>3){
	//   cout << "Get the first element in the block" << endl;
	PFCandidate::ElementInBlock ele1 = thisPart.at(3);
	//PFBlock4[nPFCand] = ele1.first;
	PFElement4[nPFCand] = ele1.second;
	//cout << "What element in this block? " << PFElement1[nPFCand] << endl;
	PFBlockRef thisBlockVec = ele1.first;
      // cout << "What is this size?  Is it ever non 1?: " << (*thisBlockVec).size() << endl;
	const PFBlock thisBlock = (*thisBlockVec);
	const edm::OwnVector< reco::PFBlockElement > elee = (thisBlock).elements();
	PFEleType4[nPFCand] =elee[ele1.second].type();
	cout << "*** Element 3, type = " << elee[ele1.second].type() << endl;
	if (PFEleType4[nPFCand]==1){
	  
	  PFItem4Eta[nPFCand] = elee[ele1.second].trackRef()->eta();
	  PFItem4Phi[nPFCand] = elee[ele1.second].trackRef()->phi();
	  PFItem4Pt[nPFCand] = elee[ele1.second].trackRef()->pt();
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ele1.second]);
	  PFItem4EtaEcal[nPFCand]=trk->positionAtECALEntrance().eta();
	  PFItem4PhiEcal[nPFCand] = trk->positionAtECALEntrance().phi();
	  cout << "*** Element 3, pt = " << elee[ele1.second].trackRef()->pt() << endl;
	}
	if (PFEleType4[nPFCand]==4){
	  PFItem4Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem4Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem4Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 3, pt = " << elee[ele1.second].clusterRef()->pt() << endl;	  
	}
	if (PFEleType4[nPFCand]==5){
	  PFItem4Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem4Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem4Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 3, pt = " << elee[ele1.second].clusterRef()->pt() << endl;
	}
      }
      if (nItems[nPFCand]>4){
	//   cout << "Get the first element in the block" << endl;
	PFCandidate::ElementInBlock ele1 = thisPart.at(4);
	//PFBlock5[nPFCand] = ele1.first;
	PFElement5[nPFCand] = ele1.second;
	//cout << "What element in this block? " << PFElement1[nPFCand] << endl;
	PFBlockRef thisBlockVec = ele1.first;
	// cout << "What is this size?  Is it ever non 1?: " << (*thisBlockVec).size() << endl;
	const PFBlock thisBlock = (*thisBlockVec);
	const edm::OwnVector< reco::PFBlockElement > elee = (thisBlock).elements();
	PFEleType5[nPFCand] =elee[ele1.second].type();
	cout << "*** Element 4, type = " << elee[ele1.second].type() << endl;
	if (PFEleType5[nPFCand]==1){
	  
	  PFItem5Eta[nPFCand] = elee[ele1.second].trackRef()->eta();
	  PFItem5Phi[nPFCand] = elee[ele1.second].trackRef()->phi();
	  PFItem5Pt[nPFCand] = elee[ele1.second].trackRef()->pt();
	  const PFBlockElementTrack *trk = dynamic_cast<const PFBlockElementTrack*>(&elee[ele1.second]);
	  PFItem5EtaEcal[nPFCand]=trk->positionAtECALEntrance().eta();
	  PFItem5PhiEcal[nPFCand] = trk->positionAtECALEntrance().phi();
	  cout << "*** Element 4, pt = " << elee[ele1.second].trackRef()->pt() << endl;
	}
	if (PFEleType5[nPFCand]==4){
	  PFItem5Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem5Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem5Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 4, pt = " << elee[ele1.second].clusterRef()->pt() << endl;	  
	}
	if (PFEleType5[nPFCand]==5){
	  PFItem5Eta[nPFCand] = elee[ele1.second].clusterRef()->eta();
	  PFItem5Phi[nPFCand] = elee[ele1.second].clusterRef()->phi();
	  PFItem5Pt[nPFCand] = elee[ele1.second].clusterRef()->pt();
	  cout << "*** Element 4, pt = " << elee[ele1.second].clusterRef()->pt() << endl;	  
	}
      }
      
      
      
      nPFCand++;
      
    }
    
    cout << "-----------------------------------\n\n" << endl;
    
    Handle<reco::PhotonCollection> ucphotoncoll;
    iEvent.getByToken(photons_old, ucphotoncoll);
    const reco::PhotonCollection *ucphotons = ucphotoncoll.product();
    reco::PhotonCollection::const_iterator ucphotonit = ucphotons->begin();
    std::map<DetId, int> crysclusEB;
    std::map<DetId, int> uccrysclusEB;
    
    
    
    nucPFPhoton=0;
    
    for(ucphotonit = ucphotons->begin(); ucphotonit!= ucphotons->end()&&nucPFPhoton<100; ++ucphotonit){
      
      std::vector<std::pair<DetId, float> > clusdet = ucphotonit->superCluster()->hitsAndFractions();
      for (int iii=0;iii<int(clusdet.size());++iii){
	uccrysclusEB.insert(make_pair(clusdet[iii].first, nucPFPhoton));
      }
      ucPFPhotonHadEM[nucPFPhoton] = ucphotonit->hadTowOverEm();
      ucPFPhotonE[nucPFPhoton]=ucphotonit->energy();
      ucPFPhotonET[nucPFPhoton]=ucphotonit->pt();
      //if (ucphotonit->pt()>40. && fabs(ucphotonit->eta())<1.6) writehipho=true;
      ucPFMIPTotE[nucPFPhoton]=ucphotonit->mipTotEnergy();
      ucPFMIPNhit[nucPFPhoton]=ucphotonit->mipNhitCone();
      //if(nPFPhoton<20)cout << "PhotonE["<< nPFPhoton << "]: " << PFPhotonE[nPFPhoton] << endl; 
      ucPFPhotonEta[nucPFPhoton]=ucphotonit->eta();
      //if(nPFPhoton<20)cout << "PhotonEta["<< nPFPhoton << "]: " << PFPhotonEta[nPFPhoton] << endl; 
      ucPFPhotonPhi[nucPFPhoton]=ucphotonit->phi();
      //if(nPFPhoton<20)cout << "PhotonPhi["<< nPFPhoton << "]: " << PFPhotonPhi[nPFPhoton] << endl; 
      ucPFPhotonSigmaIetaIeta[nucPFPhoton]=ucphotonit->full5x5_sigmaIetaIeta();
      //ucPFPhotonSigmaIphiIphi[nucPFPhoton]=ucphotonit->full5x5_sigmaIphiIphi();
      //if(nPFPhoton<20)cout << "PhotonSigmaietaIeta["<< nPFPhoton << "]: " << PFPhotonSigmaIetaIeta[nPFPhoton] << endl; 
      ucPFPhotonR9[nucPFPhoton]=ucphotonit->r9();
      //if(nPFPhoton<20)cout << "PhotonR9["<< nPFPhoton << "]: " << PFPhotonR9[nPFPhoton] << endl;
      ucPFPhotonZ[nucPFPhoton]=(*primVtxRef).z();
      //cout << "PFPhoZ["<< nPFPhoton << "]: " << PFPhotonZ[nPFPhoton] << endl;
      
      // isolator03_.fGetIsolation(&*photonit,pfH.product(),primVtxRef,vtxH);
      ucPFChargedHadronIso[nucPFPhoton]= ucphotonit->chargedHadronIso();
      //if(nPFPhoton<20)cout << "ChargedIso["<< nPFPhoton << "]: " << PFChargedHadronIso[nPFPhoton] << endl;
      ucPFNeutralHadronIso[nucPFPhoton]= ucphotonit->neutralHadronIso();
      //if(nPFPhoton<20)cout << "NeutralIso["<< nPFPhoton << "]: " << PFNeutralHadronIso[nPFPhoton] << endl;
      ucPFPhotonIso[nucPFPhoton]= ucphotonit->photonIso();
      //if(nPFPhoton<20)cout << "PhotonIso["<< nPFPhoton << "]: " << PFPhotonIso[nPFPhoton] << endl; 
      ucPFEcalIso03[nucPFPhoton] = ucphotonit->ecalRecHitSumEtConeDR03();
      ucPFHcalIso03[nucPFPhoton] = ucphotonit->hcalTowerSumEtConeDR03();
      ucPFTrkIso03[nucPFPhoton] = ucphotonit->trkSumPtHollowConeDR03();
      ucPFPixMatch[nucPFPhoton] = ucphotonit->hasPixelSeed();
      
      //PFChargedHadronIso[nPFPhoton]=photonit->chargedHadronIso();
      
      //PFNeutralHadronIso[nPFPhoton]=photonit->neutralHadronIso();
      
      //PFPhotonIso[nPFPhoton]=photonit->photonIso();
      
      //PFPUChargedHadronIso[nPFPhoton]=photonit->puChargedHadronIso();
      //if(nPFPhoton<20)cout << "PUChargedIso["<< nPFPhoton << "]: " << PFPUChargedHadronIso[nPFPhoton] << endl; 
      
      //if(nPFPhoton<20)cout << "R9["<< nPFPhoton << "]: " << PFPhotonR9[nPFPhoton] << endl; 
      nucPFPhoton++;
      
    }
    nAllCellsHE=0;
    //  Handle<HBHERecHitCollection> hcalHits;
    //  iEvent.getByToken(_hcalHits, hcalHits);
    //  const HBHERecHitCollection* rechitsHcal = hcalHits.product();
    for (HBHERecHitCollection::const_iterator it=rechitsHcal->begin();it!=rechitsHcal->end() && nAllCellsHE<30000;++it){
      HcalDetId det = it->id();
      AllCellsIEtaHE[nAllCellsHE]=det.ieta();
      AllCellsIPhiHE[nAllCellsHE]=det.iphi();
      AllCellsE_HE[nAllCellsHE]=it->energy();
      AllCellsIETAHI[nAllCellsHE]=det.crystal_ieta_high();
      AllCellsIETALO[nAllCellsHE]=det.crystal_ieta_low();
      AllCellsIPHIHI[nAllCellsHE]=det.crystal_iphi_high();
      AllCellsIPHILO[nAllCellsHE]=det.crystal_iphi_low();
      nAllCellsHE++;

       if(abs(det.ieta()) <= 8 && it->energy() >= 0.8){
       	cout << "********* HCAL hit with (ieta,iphi,hitE) = " << det.ieta() << ",\t" << det.iphi() << ",\t" << it->energy() << endl;
       }
      
    }
    //cout << "Got EE hits" <<endl;
    // ECAL clusters
    //std::cout << nGenPartjkdg << std::endl;
    //std::map<DetId, int> crysclusEB;
    
    
    nAllCellsEB=0;
    for (EcalRecHitCollection::const_iterator it = rechitsCollectionEB_->begin();it!=rechitsCollectionEB_->end()&&nAllCellsEB<30000;++it){
      DetId blarg = it->detid();
      //Here's where the map is used to put what photon this crystal goes with. 
      std::map<DetId, int>::const_iterator selcrys = crysclusEB.find(blarg);
      if (selcrys==crysclusEB.end()){
	AllClusteredEB[nAllCellsEB]=-1;
      }
      else{
	AllClusteredEB[nAllCellsEB]=selcrys->second;
      }
      std::map<DetId, int>::const_iterator ucselcrys = uccrysclusEB.find(blarg);
      if (ucselcrys==uccrysclusEB.end()){
	AllClusteredEBUt[nAllCellsEB]=-1;
      }
      else{
      AllClusteredEBUt[nAllCellsEB]=ucselcrys->second;
      }
      
      
      EBDetId dit = it->detid();
      AllCellsIEtaEB[nAllCellsEB]=dit.ieta();
      AllCellsIPhiEB[nAllCellsEB]=dit.iphi();
      AllCellsE_EB[nAllCellsEB]=it->energy();
      AllTimeEB[nAllCellsEB]=it->time();
      AllTimeErrEB[nAllCellsEB]=it->timeError();

      if(abs(dit.ieta()) < 40 && it->energy() >= 0.08){
       	cout << "********* ECAL hit with (ieta,iphi,hitE) = " << dit.ieta() << ",\t" << dit.iphi() << ",\t" << it->energy() << endl;
       }

      nAllCellsEB++;
    }
    //cout << "Got EB hits" <<endl;
    
    /////////// ak5pf jets //////////////
    Handle<reco::PFJetCollection> pfjetcoll;
    iEvent.getByToken(pfjettoken,//InputTag("ak4PFJets")
		      pfjetcoll);
    const reco::PFJetCollection *pfjets = pfjetcoll.product();
    reco::PFJetCollection::const_iterator pfjetclus = pfjets->begin();

    cout << "-----------------------------------------" << endl;
    cout << "---------------- JETS -------------------" << endl;
    
    nPFJet=0;
    //int nPFCloseJet = 0 ;
    //int nPFFarJet = 0;
    for(pfjetclus = pfjets->begin(); pfjetclus!= pfjets->end() && nPFJet<500; ++pfjetclus){
      //int flag(0);
      //if((pfjetclus->energy() < 290) || (pfjetclus->energy() > 1050) || (fabs(pfjetclus->eta()) < 1.39) || (fabs(pfjetclus->eta()) > 2.41))continue;
      //for(int l=0 ; l<nGenJet; ++l){
      //float dRGen = dRCalc(GenJetEta[l], GenJetPhi[l], pfjetclus->eta(), pfjetclus->phi());
      //float GenJetE = GenJetEt[l] * cosh(GenJetEta[l]);
      //if((GenJetE < 200) || (GenJetE > 1050) || (fabs(GenJetEta[l]) < 1.39) || (fabs(GenJetEta[l]) > 2.41))continue;
      //if (dRGen < 0.2)nPFCloseJet++;
      //if (dRGen > 0.2)nPFFarJet++;
      //if(dRGen < 0.2){
      //flag=1;
      PFJetE[nPFJet]=pfjetclus->energy();
      //cout << "JetE["<< nPFJet << "]: " << PFJetE[nPFJet] << endl;
      PFJetEta[nPFJet]=pfjetclus->eta();
      //cout << "JetEta["<< nPFJet << "]: " << PFJetEta[nPFJet] << endl;
      PFJetPhi[nPFJet]=pfjetclus->phi();
      //cout << "JetPhi["<< nPFJet << "]: " << PFJetPhi[nPFJet] << endl;
      PFJetZ[nPFJet]=(*primVtxRef).z();
      //cout << pfjetclus->vertex().Z() << endl;
      //cout << pfjetclus->track()->position().z() << endl;
      //cout << "JetZ[" << nPFJet << "]: " << PFJetZ[nPFJet] << endl;
      PFChargedFraction[nPFJet]=pfjetclus->chargedHadronEnergyFraction();
      //cout << "ChargedFrac["<< nPFJet << "]: " << PFChargedFraction[nPFJet] << endl;
      PFNeutralFraction[nPFJet]=pfjetclus->neutralHadronEnergyFraction();
      //cout << "NeutralFrac["<< nPFJet << "]: " << PFNeutralFraction[nPFJet] << endl;
      PFPhotonFraction[nPFJet]=pfjetclus->photonEnergyFraction();
      //cout << "PhotonFrac["<< nPFJet << "]: " << PFPhotonFraction[nPFJet] << endl;
      PFElectronFraction[nPFJet]=pfjetclus->electronEnergyFraction();
      //cout << "ElectronFrac["<< nPFJet << "]: " << PFElectronFraction[nPFJet] << endl;
      PFMuonFraction[nPFJet]=pfjetclus->muonEnergyFraction();
      PFJetNCon[nPFJet]=pfjetclus->chargedHadronMultiplicity()+pfjetclus->neutralHadronMultiplicity()+pfjetclus->photonMultiplicity();
      //PFTrackDxy[nPFJet] = pfjetclus->dxy(thebs->position());
      //PFTrackDsz[nPFJet] = pfjetclus->dsz(thebs->position());
      //if(nPFJet<20)cout << "MuonFrac["<< nPFJet << "]: " << PFMuonFraction[nPFJet] << endl;
      nPFJet++;
      //}
      //}
      cout << "JET: (eta,phi,E) = " << pfjetclus->eta() << ".\t" << pfjetclus->phi() << ",\t" << pfjetclus->energy() << endl;
      
    }
    cout << "-----------------------------------------------\n\n" << endl;

    cout << "nPUInfo: "<< nPUInfo << endl;
    cout << "nMC: " << nMC << endl;
    cout << "Vertex_n: " << Vertex_n << endl;
    cout << "nAllTracks: " << nAllTracks << endl;
    cout << "nPFCand: " << nPFCand << endl;
    cout << "nAllCellsEB: " << nAllCellsEB << endl;
    cout << "nAllCellsHE: " << nAllCellsHE << endl;
    cout << "nPFJet: " << nPFJet << endl;
    cout << "nucPFPhoton: " << nucPFPhoton << endl;
    
    //if (writehipho)
    tinyTree->Fill();
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
Untime::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Untime::endJob() 
{
  // go to *OUR* root file and store histograms
  rootFile_->cd();
  tinyTree->Write();

  rootFile_->Write();
  rootFile_->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
Untime::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Untime::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Untime::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Untime::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Untime::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

}

//define this as a plug-in
DEFINE_FWK_MODULE(Untime);
