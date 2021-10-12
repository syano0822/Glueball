#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TRefArray.h"
#include "TLorentzVector.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolMuon.h"
#include "AliEventPoolManager.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliPIDCombined.h"
#include "AliAODDimuon.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliAODPid.h"
#include "AliTOFHeader.h"

#include "AliMultSelection.h"
#include "AliAnalysisMuonUtility.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskAODTrackPair.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)       

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair() : 
  AliAnalysisTaskSE(),
  fEvent(NULL),
  fPoolKaonTrackMgr(NULL),
  fPoolPionTrackMgr(NULL),
  fPoolK0sMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),  
  fOnKaonAna(false),
  fOnPionAna(false),
  fOnK0sAna(false),
  fOnTrackQA(false),
  fRunNumber(-99999),
  fTrackDepth(1000),
  fPoolSize(1),
  fReadyFraction(0.1),
  fTriggerMaskForSame(AliVEvent::kINT7),
  fTriggerMaskForMixing(AliVEvent::kINT7),
  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  fOutputList(NULL),
  fEventCounter(NULL),

  fHistArmenterosPodolanskiPlot(NULL),
  fHistArmenterosPodolanskiPlot_K0s(NULL),

  fHistULSKKPairMassPt(NULL),
  fHistLSppKKPairMassPt(NULL),
  fHistLSmmKKPairMassPt(NULL),

  fHistMixULSKKPairMassPt(NULL),
  fHistMixLSppKKPairMassPt(NULL),
  fHistMixLSmmKKPairMassPt(NULL),

  fHistULSKKPairMassPt_PionID(NULL),
  fHistLSppKKPairMassPt_PionID(NULL),
  fHistLSmmKKPairMassPt_PionID(NULL),

  fHistMixULSKKPairMassPt_PionID(NULL),
  fHistMixLSppKKPairMassPt_PionID(NULL),
  fHistMixLSmmKKPairMassPt_PionID(NULL),

  fHistULSPiPiPairMassPt(NULL),
  fHistLSppPiPiPairMassPt(NULL),
  fHistLSmmPiPiPairMassPt(NULL),

  fHistMixULSPiPiPairMassPt(NULL),
  fHistMixLSppPiPiPairMassPt(NULL),
  fHistMixLSmmPiPiPairMassPt(NULL),
  
  fHistK0sMassPt(NULL),
  fHistK0sK0sMassPt(NULL),
  fHistMixK0sK0sMassPt(NULL),

  fHistPrimaryTrackNSPDCluster(NULL),
  fHistPrimaryTrackNITSCluster(NULL),
  fHistPrimaryTrackNTPCCluster(NULL),
  fHistPrimaryTrackEta(NULL),
  fHistPrimaryTrackChi2perNDF(NULL),
  fHistPrimaryTrackITSChi2perNDF(NULL),
  fHistPrimaryTrackTPCChi2perNDF(NULL),
  fHistPrimaryTrackTPCCrossRaws(NULL),
  fHistPrimaryTrackTPCFindableCrossRaws(NULL),
  fHistPrimaryTrackTPCCrossRawsOverFindable(NULL),
  fHistPrimaryTrackTOFBunchCrossing(NULL),
  fHistPrimaryTrackKinkIndex(NULL),
  fHistPrimaryTrackDCAxy(NULL),
  fHistPrimaryTrackDCAz(NULL),

  fHistPrimaryTrackNSPDClusterQuality(NULL),
  fHistPrimaryTrackNITSClusterQuality(NULL),
  fHistPrimaryTrackNTPCClusterQuality(NULL),
  fHistPrimaryTrackEtaQuality(NULL),
  fHistPrimaryTrackChi2perNDFQuality(NULL),
  fHistPrimaryTrackITSChi2perNDFQuality(NULL),
  fHistPrimaryTrackTPCChi2perNDFQuality(NULL),
  fHistPrimaryTrackTPCCrossRawsQuality(NULL),
  fHistPrimaryTrackTPCFindableCrossRawsQuality(NULL),
  fHistPrimaryTrackTPCCrossRawsOverFindableQuality(NULL),
  fHistPrimaryTrackTOFBunchCrossingQuality(NULL),
  fHistPrimaryTrackKinkIndexQuality(NULL),
  fHistPrimaryTrackDCAxyQuality(NULL),
  fHistPrimaryTrackDCAzQuality(NULL),


  fHistMomTPCSignalTrackQuality(NULL),
  fHistMomTOFBetaTrackQuality(NULL),
  fHistMomTPCSignalTrackQualityPionPID(NULL),
  fHistMomTOFBetaTrackQualityPionPID(NULL),
  fHistMomTPCSignalTrackQualityKaonPID(NULL),
  fHistMomTOFBetaTrackQualityKaonPID(NULL),
  fHistMomTPCSignalTrackQualityProtonPID(NULL),
  fHistMomTOFBetaTrackQualityProtonPID(NULL),
  
  fHistMomTPCSigmaTrackQuality(NULL),
  fHistMomTOFSigmaTrackQuality(NULL),
  fHistMomTPCSigmaTrackQualityPionPID(NULL),
  fHistMomTOFSigmaTrackQualityPionPID(NULL),
  fHistMomTPCSigmaTrackQualityKaonPID(NULL),
  fHistMomTOFSigmaTrackQualityKaonPID(NULL),
  fHistMomTPCSigmaTrackQualityProtonPID(NULL),
  fHistMomTOFSigmaTrackQualityProtonPID(NULL)

{ 
  
}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char* name) : 
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fPoolKaonTrackMgr(NULL),
  fPoolPionTrackMgr(NULL),
  fPoolK0sMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),
  fOnKaonAna(false),
  fOnPionAna(false),
  fOnK0sAna(false),
  fOnTrackQA(false),
  fRunNumber(-99999),
  fTrackDepth(1000),
  fPoolSize(1),
  fReadyFraction(0.1),
  fTriggerMaskForSame(AliVEvent::kINT7),
  fTriggerMaskForMixing(AliVEvent::kINT7),
  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  fOutputList(NULL),
  fEventCounter(NULL),
  fHistArmenterosPodolanskiPlot(NULL),
  fHistArmenterosPodolanskiPlot_K0s(NULL),

  fHistULSKKPairMassPt(NULL),
  fHistLSppKKPairMassPt(NULL),
  fHistLSmmKKPairMassPt(NULL),

  fHistMixULSKKPairMassPt(NULL),
  fHistMixLSppKKPairMassPt(NULL),
  fHistMixLSmmKKPairMassPt(NULL),

  fHistULSKKPairMassPt_PionID(NULL),
  fHistLSppKKPairMassPt_PionID(NULL),
  fHistLSmmKKPairMassPt_PionID(NULL),

  fHistMixULSKKPairMassPt_PionID(NULL),
  fHistMixLSppKKPairMassPt_PionID(NULL),
  fHistMixLSmmKKPairMassPt_PionID(NULL),
  
  fHistULSPiPiPairMassPt(NULL),
  fHistLSppPiPiPairMassPt(NULL),
  fHistLSmmPiPiPairMassPt(NULL),

  fHistMixULSPiPiPairMassPt(NULL),
  fHistMixLSppPiPiPairMassPt(NULL),
  fHistMixLSmmPiPiPairMassPt(NULL),

  fHistK0sMassPt(NULL),
  fHistK0sK0sMassPt(NULL),
  fHistMixK0sK0sMassPt(NULL),
  
  fHistPrimaryTrackNSPDCluster(NULL),
  fHistPrimaryTrackNITSCluster(NULL),
  fHistPrimaryTrackNTPCCluster(NULL),
  fHistPrimaryTrackEta(NULL),
  fHistPrimaryTrackChi2perNDF(NULL),
  fHistPrimaryTrackITSChi2perNDF(NULL),
  fHistPrimaryTrackTPCChi2perNDF(NULL),
  fHistPrimaryTrackTPCCrossRaws(NULL),
  fHistPrimaryTrackTPCFindableCrossRaws(NULL),
  fHistPrimaryTrackTPCCrossRawsOverFindable(NULL),
  fHistPrimaryTrackTOFBunchCrossing(NULL),
  fHistPrimaryTrackKinkIndex(NULL),
  fHistPrimaryTrackDCAxy(NULL),
  fHistPrimaryTrackDCAz(NULL),

  fHistPrimaryTrackNSPDClusterQuality(NULL),
  fHistPrimaryTrackNITSClusterQuality(NULL),
  fHistPrimaryTrackNTPCClusterQuality(NULL),
  fHistPrimaryTrackEtaQuality(NULL),
  fHistPrimaryTrackChi2perNDFQuality(NULL),
  fHistPrimaryTrackITSChi2perNDFQuality(NULL),
  fHistPrimaryTrackTPCChi2perNDFQuality(NULL),
  fHistPrimaryTrackTPCCrossRawsQuality(NULL),
  fHistPrimaryTrackTPCFindableCrossRawsQuality(NULL),
  fHistPrimaryTrackTPCCrossRawsOverFindableQuality(NULL),
  fHistPrimaryTrackTOFBunchCrossingQuality(NULL),
  fHistPrimaryTrackKinkIndexQuality(NULL),
  fHistPrimaryTrackDCAxyQuality(NULL),
  fHistPrimaryTrackDCAzQuality(NULL),

  fHistMomTPCSignalTrackQuality(NULL),
  fHistMomTOFBetaTrackQuality(NULL),
  fHistMomTPCSignalTrackQualityPionPID(NULL),
  fHistMomTOFBetaTrackQualityPionPID(NULL),
  fHistMomTPCSignalTrackQualityKaonPID(NULL),
  fHistMomTOFBetaTrackQualityKaonPID(NULL),
  fHistMomTPCSignalTrackQualityProtonPID(NULL),
  fHistMomTOFBetaTrackQualityProtonPID(NULL),
  fHistMomTPCSigmaTrackQuality(NULL),
  fHistMomTOFSigmaTrackQuality(NULL),
  fHistMomTPCSigmaTrackQualityPionPID(NULL),
  fHistMomTOFSigmaTrackQualityPionPID(NULL),
  fHistMomTPCSigmaTrackQualityKaonPID(NULL),
  fHistMomTOFSigmaTrackQualityKaonPID(NULL),
  fHistMomTPCSigmaTrackQualityProtonPID(NULL),
  fHistMomTOFSigmaTrackQualityProtonPID(NULL)
{ 
  
  
  double fCentBins[] = {-1,10,20,30,40,50,60,70,80,90,101};
  double fVtxBins[] = {-50,-10.5,-6,-2,0,2,6,10.5,50};
  double fPsiBins[] = {-10,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,10};
  
  int fNCentBins = sizeof(fCentBins)/sizeof(double)-1;
  int fNVtxZBins = sizeof(fVtxBins)/sizeof(double)-1;
  int fNPsiBins = sizeof(fPsiBins)/sizeof(double)-1;

  fPoolKaonTrackMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,
					      fNCentBins,(double*)fCentBins,
					      fNVtxZBins,(double*)fVtxBins,
					      fNPsiBins,(double*)fPsiBins);  
  fPoolKaonTrackMgr->SetTargetValues(fTrackDepth,(double)fReadyFraction,fPoolSize);  

  fPoolPionTrackMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,
					      fNCentBins,(double*)fCentBins,
					      fNVtxZBins,(double*)fVtxBins,
					      fNPsiBins,(double*)fPsiBins);  
  fPoolPionTrackMgr->SetTargetValues(fTrackDepth,(double)fReadyFraction,fPoolSize);  

  fPoolK0sMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,
					      fNCentBins,(double*)fCentBins,
					      fNVtxZBins,(double*)fVtxBins,
					      fNPsiBins,(double*)fPsiBins);  
  fPoolK0sMgr->SetTargetValues(fTrackDepth,(double)fReadyFraction,fPoolSize);  
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPair::~AliAnalysisTaskAODTrackPair() 
{
  cout<<"Destructure"<<endl;
}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPair::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  double bins_cent_hist[]={0.,1.,2.,3.,4.,5.,10.,20.,30.,50.,70.,100.};
  int binnum_cent_hist = sizeof(bins_cent_hist)/sizeof(double) - 1;
  
  double bins_event_hist[]={0,0.1,0.5,1,5,10,20,40,60,90,100};
  int binnum_event_hist = sizeof(bins_event_hist)/sizeof(double) - 1;

  fEventCounter = new TH1F("fEventCounter","",binnum_event_hist,bins_event_hist);
  fOutputList->Add(fEventCounter);
  
  if(fOnK0sAna){
    fHistK0sMassPt = new TH2F("fHistK0sMassPt","fHistK0sMassPt",5000,0.0,5.0,100,0,10);
    fHistK0sK0sMassPt = new TH2F("fHistK0sK0sMassPt","fHistK0sK0sMassPt",600,1.0,4.0,100,0,10);
    fHistMixK0sK0sMassPt = new TH2F("fHistMixK0sK0sMassPt","fHistMixK0sK0sMassPt",600,1.0,4.0,100,0,10);
    fOutputList->Add(fHistK0sMassPt);  
    fOutputList->Add(fHistK0sK0sMassPt);
    fOutputList->Add(fHistMixK0sK0sMassPt);
    
    fHistArmenterosPodolanskiPlot = new TH2F("fHistArmenterosPodolanskiPlot","",220,-1.1,1.1,300,0,0.3);
    fHistArmenterosPodolanskiPlot_K0s = new TH2F("fHistArmenterosPodolanskiPlot_K0s","",220,-1.1,1.1,300,0,0.3);
    fOutputList->Add(fHistArmenterosPodolanskiPlot);
    fOutputList->Add(fHistArmenterosPodolanskiPlot_K0s);
  }

  if(fOnKaonAna){
    fHistULSKKPairMassPt = new TH2F("fHistULSKKPairMassPt","",700,0.5,4.,100,0,10);
    fHistLSppKKPairMassPt = new TH2F("fHistLSppKKPairMassPt","",700,0.5,4.,100,0,10);
    fHistLSmmKKPairMassPt = new TH2F("fHistLSmmKKPairMassPt","",700,0.5,4.,100,0,10);
    fOutputList->Add(fHistULSKKPairMassPt);
    fOutputList->Add(fHistLSppKKPairMassPt);
    fOutputList->Add(fHistLSmmKKPairMassPt);

    fHistMixULSKKPairMassPt = new TH2F("fHistMixULSKKPairMassPt","",700,0.5,4.,100,0,10);
    fHistMixLSppKKPairMassPt = new TH2F("fHistMixLSppKKPairMassPt","",700,0.5,4.,100,0,10);
    fHistMixLSmmKKPairMassPt = new TH2F("fHistMixLSmmKKPairMassPt","",700,0.5,4.,100,0,10);
    fOutputList->Add(fHistMixULSKKPairMassPt);
    fOutputList->Add(fHistMixLSppKKPairMassPt);
    fOutputList->Add(fHistMixLSmmKKPairMassPt);
  }

  if(fOnPionAna){
    fHistULSPiPiPairMassPt = new TH2F("fHistULSPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistLSppPiPiPairMassPt = new TH2F("fHistLSppPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistLSmmPiPiPairMassPt = new TH2F("fHistLSmmPiPiPairMassPt","",800,0.,4.,100,0,10);
    fOutputList->Add(fHistULSPiPiPairMassPt);
    fOutputList->Add(fHistLSppPiPiPairMassPt);
    fOutputList->Add(fHistLSmmPiPiPairMassPt);
    
    fHistMixULSPiPiPairMassPt = new TH2F("fHistMixULSPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistMixLSppPiPiPairMassPt = new TH2F("fHistMixLSppPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistMixLSmmPiPiPairMassPt = new TH2F("fHistMixLSmmPiPiPairMassPt","",800,0.,4.,100,0,10);
    fOutputList->Add(fHistMixULSPiPiPairMassPt);
    fOutputList->Add(fHistMixLSppPiPiPairMassPt);
    fOutputList->Add(fHistMixLSmmPiPiPairMassPt);
  }

  if(fOnTrackQA){
    fHistPrimaryTrackNSPDClusterQuality = new TH2F("fHistPrimaryTrackNSPDClusterQuality","",100,0,10,3,0,3);
    fHistPrimaryTrackNITSClusterQuality = new TH2F("fHistPrimaryTrackNITSClusterQuality","",100,0,10,7,0,7);
    fHistPrimaryTrackNTPCClusterQuality = new TH2F("fHistPrimaryTrackNTPCClusterQuality","",100,0,10,160,0,160);
    fHistPrimaryTrackEtaQuality = new TH2F("fHistPrimaryTrackEtaQuality","",100,0,10,200,-1.,1.);
    fHistPrimaryTrackChi2perNDFQuality = new TH2F("fHistPrimaryTrackChi2perNDFQuality","",100,0,10,100,0,10);
    fHistPrimaryTrackITSChi2perNDFQuality = new TH2F("fHistPrimaryTrackITSChi2perNDFQuality","",100,0,10,50,0,50);
    fHistPrimaryTrackTPCChi2perNDFQuality = new TH2F("fHistPrimaryTrackTPCChi2perNDFQuality","",100,0,10,100,0,10);
    fHistPrimaryTrackTPCCrossRawsQuality = new TH2F("fHistPrimaryTrackTPCCrossRawsQuality","",100,0,10,160,0,160);
    fHistPrimaryTrackTPCFindableCrossRawsQuality = new TH2F("fHistPrimaryTrackTPCFindableCrossRawsQuality","",100,0,10,160,0,160);
    fHistPrimaryTrackTPCCrossRawsOverFindableQuality = new TH2F("fHistPrimaryTrackTPCCrossRawsOverFindableQuality","",100,0,10,100,0,1);
    fHistPrimaryTrackTOFBunchCrossingQuality = new TH2F("fHistPrimaryTrackTOFBunchCrossingQuality","",100,0,10,20,0,20);
    fHistPrimaryTrackKinkIndexQuality = new TH2F("fHistPrimaryTrackKinkIndexQuality","",100,0,10,2,0,2);
    fHistPrimaryTrackDCAxyQuality = new TH2F("fHistPrimaryTrackDCAxyQuality","",100,0,10,200,-1,1);
    fHistPrimaryTrackDCAzQuality = new TH2F("fHistPrimaryTrackDCAzQuality","",100,0,10,50,-2.5,2.5);
    fOutputList->Add(fHistPrimaryTrackNSPDClusterQuality);
    fOutputList->Add(fHistPrimaryTrackNITSClusterQuality);
    fOutputList->Add(fHistPrimaryTrackNTPCClusterQuality);
    fOutputList->Add(fHistPrimaryTrackEtaQuality);
    fOutputList->Add(fHistPrimaryTrackChi2perNDFQuality);
    fOutputList->Add(fHistPrimaryTrackITSChi2perNDFQuality);
    fOutputList->Add(fHistPrimaryTrackTPCChi2perNDFQuality);
    fOutputList->Add(fHistPrimaryTrackTPCCrossRawsQuality);
    fOutputList->Add(fHistPrimaryTrackTPCFindableCrossRawsQuality);
    fOutputList->Add(fHistPrimaryTrackTPCCrossRawsOverFindableQuality);
    fOutputList->Add(fHistPrimaryTrackTOFBunchCrossingQuality);
    fOutputList->Add(fHistPrimaryTrackKinkIndexQuality);
    fOutputList->Add(fHistPrimaryTrackDCAxyQuality);
    fOutputList->Add(fHistPrimaryTrackDCAzQuality);
    
  
    fHistMomTPCSignalTrackQuality = new TH2F("fHistMomTPCSignalTrackQuality","",1000,0,10,1000,0,1000);
    fHistMomTOFBetaTrackQuality = new TH2F("fHistMomTOFBetaTrackQuality","",1000,0,10,200,0,2);
    fHistMomTPCSignalTrackQualityPionPID = new TH2F("fHistMomTPCSignalTrackQualityPionPID","",1000,0,10,1000,0,1000);
    fHistMomTOFBetaTrackQualityPionPID = new TH2F("fHistMomTOFBetaTrackQualityPionPID","",1000,0,10,200,0,2);
    fHistMomTPCSignalTrackQualityKaonPID = new TH2F("fHistMomTPCSignalTrackQualityKaonPID","",1000,0,10,1000,0,1000);
    fHistMomTOFBetaTrackQualityKaonPID = new TH2F("fHistMomTOFBetaTrackQualityKaonPID","",1000,0,10,200,0,2);
    fHistMomTPCSignalTrackQualityProtonPID = new TH2F("fHistMomTPCSignalTrackQualityProtonPID","",1000,0,10,1000,0,1000);
    fHistMomTOFBetaTrackQualityProtonPID = new TH2F("fHistMomTOFBetaTrackQualityProtonPID","",1000,0,10,200,0,2);  
    fOutputList->Add(fHistMomTPCSignalTrackQuality);
    fOutputList->Add(fHistMomTOFBetaTrackQuality);
    fOutputList->Add(fHistMomTPCSignalTrackQualityPionPID);
    fOutputList->Add(fHistMomTOFBetaTrackQualityPionPID);
    fOutputList->Add(fHistMomTPCSignalTrackQualityKaonPID);
    fOutputList->Add(fHistMomTOFBetaTrackQualityKaonPID);
    fOutputList->Add(fHistMomTPCSignalTrackQualityProtonPID);
    fOutputList->Add(fHistMomTOFBetaTrackQualityProtonPID);
  }
  
  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{  

  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;
  fEventCounter->Fill(fUtils->getCentClass());
  if(fOnK0sAna) V0TrackPairAnalysis();
  if(fOnKaonAna || fOnPionAna) PrimeTrackPairAnalysis();
}

bool AliAnalysisTaskAODTrackPair::TrackQA(AliAODTrack* track)
{
  
  if(!track) return false;
  
  float sigTPC = track->GetTPCsignal();
  float sigTOF =track->GetTOFsignal();
  float length =track->GetIntegratedLength();
  float beta =(sigTOF>0)?(double)length/ (2.99792457999999984e-02 * sigTOF):0;

  float dca_xy=9999;
  float dca_z=9999;
  track->GetImpactParameters(dca_xy,dca_z);

  if( fUtils->isAcceptPrimeTrackQuality(track) ){
      
    fHistPrimaryTrackNSPDClusterQuality->Fill(track->Pt(),track->GetITSNcls());
    fHistPrimaryTrackNITSClusterQuality->Fill(track->Pt(),track->GetITSNcls());
    fHistPrimaryTrackNTPCClusterQuality->Fill(track->Pt(),track->GetTPCNcls());
    fHistPrimaryTrackEtaQuality->Fill(track->Pt(),track->Eta());
    fHistPrimaryTrackChi2perNDFQuality->Fill(track->Pt(),track->Chi2perNDF());
    fHistPrimaryTrackITSChi2perNDFQuality->Fill(track->Pt(),track->GetITSNcls()>0 ? track->GetITSchi2()/track->GetITSNcls() : 0);
    fHistPrimaryTrackTPCChi2perNDFQuality->Fill(track->Pt(),track->GetTPCNcls()>0 ? track->GetTPCchi2()/track->GetTPCNcls() : 0);
    fHistPrimaryTrackTPCCrossRawsQuality->Fill(track->Pt(),track->GetTPCNCrossedRows());
    fHistPrimaryTrackTPCFindableCrossRawsQuality->Fill(track->Pt(),track->GetTPCNclsF());
    fHistPrimaryTrackTPCCrossRawsOverFindableQuality->Fill(track->Pt(),track->GetTPCNclsF()>0 ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0);
    fHistPrimaryTrackTOFBunchCrossingQuality->Fill(track->Pt(),track->GetTOFBunchCrossing());
    fHistPrimaryTrackKinkIndexQuality->Fill(track->Pt(),track->GetKinkIndex(0));
    fHistPrimaryTrackDCAxyQuality->Fill(track->Pt(),dca_xy);
    fHistPrimaryTrackDCAzQuality->Fill(track->Pt(),dca_z);
            
    fHistMomTPCSignalTrackQuality->Fill(track->P(),sigTPC);      
      
    if(track->GetTOFBunchCrossing()==0) {
      fHistMomTOFBetaTrackQuality->Fill(track->P(),beta);
    }

    if( fUtils->isAcceptPid(track,AliPID::kPion) ){
      fHistMomTPCSignalTrackQualityPionPID->Fill(track->P(),sigTPC);		
      if(track->GetTOFBunchCrossing()==0) {
	fHistMomTOFBetaTrackQualityPionPID->Fill(track->P(),beta);
      }
    }

    if( fUtils->isAcceptPid(track,AliPID::kKaon) ){
      fHistMomTPCSignalTrackQualityKaonPID->Fill(track->P(),sigTPC);		
      if(track->GetTOFBunchCrossing()==0) {
	fHistMomTOFBetaTrackQualityKaonPID->Fill(track->P(),beta);
      }
    }
      
    if( fUtils->isAcceptPid(track,AliPID::kProton) ){
      fHistMomTPCSignalTrackQualityProtonPID->Fill(track->P(),sigTPC);		
      if(track->GetTOFBunchCrossing()==0) {
	fHistMomTOFBetaTrackQualityProtonPID->Fill(track->P(),beta);
      }
    }    
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::Initialize()
{

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());

  if( !fUtils->setEvent(fEvent,fInputHandler) )
    return false;
  
  if(fIsMC){
    AliAODMCHeader* mcHeader = (AliAODMCHeader*) fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) return false;  
    TList* headerList = mcHeader->GetCocktailHeaders();
    if(!headerList) return false;
    //Double_t impPar = mcHeader->GetImpactParameter();    
    TClonesArray *mcarray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!mcarray) return false;
    fUtils->setMCArray(mcarray);
  }
  
  return true;

}


AliEventPool* AliAnalysisTaskAODTrackPair::setPool(int id)
{
  /*
  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;
  
  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi();
  
  AliEventPool* pool = NULL;;
  
  if(id==0){
    pool = (AliEventPool*)fPoolKaonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  }
  else if(id==1){
    pool = (AliEventPool*)fPoolPionTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  }
  else{
    pool = (AliEventPool*)fPoolK0sMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  }

  return pool;
  */

  return 0;
}

bool AliAnalysisTaskAODTrackPair::PrimeTrackPairAnalysis(){

  TObjArray* fKaonTrackArray = NULL;
  TObjArray* fPionTrackArray = NULL;

  if(fOnKaonAna){
    fKaonTrackArray = new TObjArray();
    fKaonTrackArray -> SetOwner();
  }
  if(fOnPionAna){
    fPionTrackArray = new TObjArray();  
    fPionTrackArray -> SetOwner();
  }

  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi(); 

  AliEventPool* poolKaon = (AliEventPool*)fPoolKaonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  AliEventPool* poolPion = (AliEventPool*)fPoolPionTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  
  /*
  cout<<"   poolKaon->GetCurrentNEvents():  "<<poolKaon->GetCurrentNEvents()<<endl;
  TObjArray* poolKaonTracks1 = (TObjArray*)poolKaon->GetEvent(0);
  TObjArray* poolKaonTracks2 = (TObjArray*)poolKaon->GetEvent(1);
  if(poolKaonTracks1)cout<<"      poolKaonTracks1:  "<<poolKaonTracks1->GetEntriesFast()<<endl;
  if(poolKaonTracks2)cout<<"      poolKaonTracks2:  "<<poolKaonTracks2->GetEntriesFast()<<endl;
  cout<<"   poolPion->GetCurrentNEvents():  "<<poolPion->GetCurrentNEvents()<<endl;
  TObjArray* poolPionTracks1 = (TObjArray*)poolPion->GetEvent(0);
  TObjArray* poolPionTracks2 = (TObjArray*)poolPion->GetEvent(1);
  if(poolPionTracks1)cout<<"      poolPionTracks1:  "<<poolPionTracks1->GetEntriesFast()<<endl;
  if(poolPionTracks2)cout<<"      poolPionTracks2:  "<<poolPionTracks2->GetEntriesFast()<<endl;
  */
  
  Int_t nTrack = fEvent->GetNumberOfTracks();
  
  AliAODTrack *track1 = NULL;
  AliAODTrack *track2 = NULL;

  TObjArray* poolTracks = NULL;

  TLorentzVector t1;
  TLorentzVector t2;
  TLorentzVector t12;

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){
    
    track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);        
    if(!track1) continue;
    
    if(fOnTrackQA) TrackQA(track1);
        
    for(Int_t iTrack2=iTrack1+1; iTrack2<nTrack; ++iTrack2){
      
      track2 = (AliAODTrack*)fEvent->GetTrack(iTrack2);              
      if(!track2) continue;

      if(fOnKaonAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon) && fUtils->isAcceptPrimeTrack(track2,AliPID::kKaon)){
	
	t1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),fUtils->fMassKaon);
	t2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),fUtils->fMassKaon);      
	t12 = t1 + t2;  
	
	if(!fUtils->isAcceptTrackPair(t12.Rapidity(),t12.Pt())) continue;

	if( track1->Charge() != track2->Charge() ){
	  fHistULSKKPairMassPt->Fill(t12.M(),t12.Pt());	
	}
	else {
	  if( track1->Charge()>0 ){
	    fHistLSppKKPairMassPt->Fill(t12.M(),t12.Pt());
	  }
	  else{
	    fHistLSmmKKPairMassPt->Fill(t12.M(),t12.Pt());
	  }
	}
      }

      if(fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion) && fUtils->isAcceptPrimeTrack(track2,AliPID::kPion)){
	
	t1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),fUtils->fMassPion);
	t2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),fUtils->fMassPion);      
	t12 = t1 + t2;  
	
	if(fUtils->isAcceptTrackPair(t12.Rapidity(),t12.Pt())){
	  if( track1->Charge() != track2->Charge() ){
	    //fHistULSPiPiPairMassPt->Fill(t12.M(),t12.Pt());	
	  }
	  else {
	    if( track1->Charge()>0 ){
	      //fHistLSppPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	    }
	    else{
	      //fHistLSmmPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	    }
	  }
	}       
      }      
    }//end of loop track2
    
    if(fOnKaonAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon)){
      
      if (poolKaon->IsReady()){
      
	for (Int_t iMixEvt=0; iMixEvt<poolKaon->GetCurrentNEvents(); iMixEvt++){

	  poolTracks = (TObjArray*)poolKaon->GetEvent(iMixEvt);	
	
	  for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){
	    
	    track2 = (AliAODTrack*)poolTracks->At(iTrack2);
	    
	    t1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),fUtils->fMassKaon);
	    t2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),fUtils->fMassKaon);      
	    t12 = t1 + t2;  
	
	    if(!fUtils->isAcceptTrackPair(t12.Rapidity(),t12.Pt())) continue;
	    
	    if( track1->Charge() != track2->Charge() ){
	      fHistMixULSKKPairMassPt->Fill(t12.M(),t12.Pt());	
	    }
	    else {
	      if( track1->Charge()>0 ){
		fHistMixLSppKKPairMassPt->Fill(t12.M(),t12.Pt());
	      }
	      else{
		fHistMixLSmmKKPairMassPt->Fill(t12.M(),t12.Pt());
	      }
	    }

	  }//end of loop track2
	}// end of loop iMixEvt
	
      }//poolKaonon->IsReady()    
    }//if(fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon))
    /*
    if(fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion)){

      if (poolPion->IsReady()){
      
	for (Int_t iMixEvt=0; iMixEvt<poolPion->GetCurrentNEvents(); iMixEvt++){

	  TObjArray* poolTracks = (TObjArray*)poolPion->GetEvent(iMixEvt);	
	
	  for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){
	  
	    AliAODTrack* __track2__ = (AliAODTrack*)poolTracks->At(iTrack2);
	    AliAODTrack* track2 = (AliAODTrack*)__track2__->Clone();

	    TLorentzVector t1;
	    TLorentzVector t2;
	    TLorentzVector t12;
	    
	    t1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),fUtils->fMassPion);
	    t2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),fUtils->fMassPion);      
	    t12 = t1 + t2;  
	
	    if(!fUtils->isAcceptTrackPair(t12.Rapidity(),t12.Pt())) continue;
	    
	    if( track1->Charge() != track2->Charge() ){
	      fHistMixULSPiPiPairMassPt->Fill(t12.M(),t12.Pt());	
	    }
	    else {
	      if( track1->Charge()>0 ){
		fHistMixLSppPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	      }
	      else{
		fHistMixLSmmPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	      }
	    }
	    
	    delete track2;

	  }//end of loop track2
	}// end of loop iMixEvt
	
      }//poolPionon->IsReady()    
    }//if(fUtils->isAcceptPrimeTrack(track1,AliPID::kPion))
    */
    if(fOnKaonAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon)){      
      fKaonTrackArray->Add(track1);
    }

    if(fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion)){            
      fPionTrackArray->Add(track1);
    }
  }//end of loop track1
  
  TObjArray* fKaonTrackArrayClone = NULL;
  TObjArray* fPionTrackArrayClone = NULL;

  if(fOnKaonAna){
    fKaonTrackArrayClone = (TObjArray*)fKaonTrackArray->Clone();
    fKaonTrackArrayClone->SetOwner();
    //cout<<"       This event stores  "<<fKaonTrackArrayClone->GetEntriesFast()<<endl;
    if(fKaonTrackArrayClone->GetEntriesFast()>0){
      poolKaon->UpdatePool(fKaonTrackArrayClone);
    }
  }

  if(fOnPionAna){
    fPionTrackArrayClone = (TObjArray*)fPionTrackArray->Clone();
    fPionTrackArrayClone->SetOwner();
    //cout<<"       This event stores  "<<fPionTrackArrayClone->GetEntriesFast()<<endl;
    if(fPionTrackArrayClone->GetEntriesFast()>0){
      poolPion->UpdatePool(fPionTrackArrayClone);
    }
  }
  
  //delete fKaonTrackArray;
  //delete fPionTrackArray;
 
  return true;
}

bool AliAnalysisTaskAODTrackPair::V0TrackPairAnalysis()
{ 
  /*
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  AliEventPool* pool = setPool(2);  

  int nV0Pair = fEvent->GetNumberOfV0s();

  for(int iV0=0; iV0<nV0Pair; ++iV0){
    
    AliAODv0 *v0_1 = fEvent->GetV0(iV0);

    if( !fUtils->isAcceptDefaultV0(v0_1) ) continue;
    
    fHistArmenterosPodolanskiPlot->Fill(v0_1->AlphaV0(),v0_1->PtArmV0());
    
    if( !fUtils->isK0sV0(v0_1) ) continue;

    fHistArmenterosPodolanskiPlot_K0s->Fill(v0_1->AlphaV0(),v0_1->PtArmV0());
    
    double fill_1[]={v0_1->MassK0Short(),v0_1->Pt(),fUtils->getCentClass()};

    fHistK0sMassPt->Fill(v0_1->MassK0Short(),v0_1->Pt());

    for(int jV0=iV0+1; jV0<nV0Pair; ++jV0){
    
      AliAODv0 *v0_2 = fEvent->GetV0(jV0);
    
      if( !fUtils->isAcceptDefaultV0(v0_2) ) continue;
      if( !fUtils->isK0sV0(v0_2) ) continue;
    
      if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;
      if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;
      
      TLorentzVector k1;
      TLorentzVector k2;      
      k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
      k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);
      
      TLorentzVector k12 = k1 + k2;

      double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};
      fHistK0sK0sMassPt->Fill(k12.M(),k12.Pt());

    }

    if (pool->IsReady()){

      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){

	TObjArray* poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);	

	for(Int_t iV02=0; iV02<poolTracks->GetEntriesFast(); ++iV02){

	  AliAODv0* __v0_2__ = (AliAODv0*)poolTracks->At(iV02);
	  AliAODv0* v0_2 = (AliAODv0*)__v0_2__->Clone();

	  if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;
	  if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;

	  TLorentzVector k1;
	  TLorentzVector k2;      
	  k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
	  k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);
      
	  TLorentzVector k12 = k1 + k2;

	  double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};
	  fHistMixK0sK0sMassPt->Fill(k12.M(),k12.Pt());
	  
	  delete v0_2;
	  
	}//end of loop track2
      }// end of loop iMixEvt
    
    }//poolPion->IsReady()    
    
    fTrackArray->Add(v0_1);
  }

  TObjArray* fTrackArrayClone = (TObjArray*)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if(fTrackArrayClone->GetEntriesFast()>0){
    pool->UpdatePool(fTrackArrayClone);
  }
  
  delete fTrackArray;
  */
  return true;
}
