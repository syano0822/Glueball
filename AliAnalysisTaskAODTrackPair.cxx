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
  fPoolTrackMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),  
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
  fSparseK0s(NULL),
  fSparseK0sK0s(NULL),
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
  fPoolTrackMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),
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
  fSparseK0s(NULL),
  fSparseK0sK0s(NULL),
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

  fPoolTrackMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,
					      fNCentBins,(double*)fCentBins,
					      fNVtxZBins,(double*)fVtxBins,
					      fNPsiBins,(double*)fPsiBins);  
  fPoolTrackMgr->SetTargetValues(fTrackDepth,(double)fReadyFraction,fPoolSize);  
  
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

  const Int_t binnum_K0s            = 3;
  Int_t bins_K0s[binnum_K0s]       = {200,100,binnum_cent_hist};
  double minbins_K0s[binnum_K0s] = {0.4,0,0};
  double maxbins_K0s[binnum_K0s] = {0.6,10,100};
  TString  namebins_K0s[]        = {"mass","pt","centrality"};

  fSparseK0s = new THnSparseF("fSparseK0s","",binnum_K0s,bins_K0s,minbins_K0s,maxbins_K0s);
  fOutputList->Add(fSparseK0s);

  const Int_t binnum_K0sK0s            = 3;
  Int_t bins_K0sK0s[binnum_K0sK0s]       = {2000,100,binnum_cent_hist};
  double minbins_K0sK0s[binnum_K0sK0s] = {0.5,0,0};
  double maxbins_K0sK0s[binnum_K0sK0s] = {2.5,10,100};
  TString  namebins_K0sK0s[]        = {"mass","pt","centrality"};
  fSparseK0sK0s = new THnSparseF("fSparseK0sK0s","",binnum_K0sK0s,bins_K0sK0s,minbins_K0sK0s,maxbins_K0sK0s);
  fOutputList->Add(fSparseK0sK0s);

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

  fHistMomTPCSigmaTrackQuality = new TH2F("fHistMomTPCSigmaTrackQuality","",1000,0,10,20,-5,5);
  fHistMomTOFSigmaTrackQuality = new TH2F("fHistMomTOFSigmaTrackQuality","",1000,0,10,20,-5,5);
  fHistMomTPCSigmaTrackQualityPionPID = new TH2F("fHistMomTPCSigmaTrackQualityPionPID","",1000,0,10,20,-5,5);
  fHistMomTOFSigmaTrackQualityPionPID = new TH2F("fHistMomTOFSigmaTrackQualityPionPID","",1000,0,10,20,-5,5);
  fHistMomTPCSigmaTrackQualityKaonPID = new TH2F("fHistMomTPCSigmaTrackQualityKaonPID","",1000,0,10,20,-5,5);
  fHistMomTOFSigmaTrackQualityKaonPID = new TH2F("fHistMomTOFSigmaTrackQualityKaonPID","",1000,0,10,20,-5,5);
  fHistMomTPCSigmaTrackQualityProtonPID = new TH2F("fHistMomTPCSigmaTrackQualityProtonPID","",1000,0,20,-5,5);
  fHistMomTOFSigmaTrackQualityProtonPID = new TH2F("fHistMomTOFSigmaTrackQualityProtonPID","",1000,0,20,-5,5);
  fOutputList->Add(fHistMomTPCSigmaTrackQuality);
  fOutputList->Add(fHistMomTOFSigmaTrackQuality);
  fOutputList->Add(fHistMomTPCSigmaTrackQualityPionPID);
  fOutputList->Add(fHistMomTOFSigmaTrackQualityPionPID);
  fOutputList->Add(fHistMomTPCSigmaTrackQualityKaonPID);
  fOutputList->Add(fHistMomTOFSigmaTrackQualityKaonPID);
  fOutputList->Add(fHistMomTPCSigmaTrackQualityProtonPID);
  fOutputList->Add(fHistMomTOFSigmaTrackQualityProtonPID);
  
  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{  
  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;
  
  fEventCounter->Fill(fUtils->getCentClass());

  TrackQA();  
  V0TrackPairAnalysis();
}

bool AliAnalysisTaskAODTrackPair::TrackQA()
{
  
  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track = (AliAODTrack*)fEvent->GetTrack(iTrack1);        
    
    float sigTPC = track->GetTPCsignal();
    float sigTOF =track->GetTOFsignal();
    float length =track->GetIntegratedLength();
    float beta =(sigTOF>0)?(double)length/ (2.99792457999999984e-02 * sigTOF):0;
    
    if( fUtils->isAcceptTrackQuality(track) ){

      fHistMomTPCSignalTrackQuality->Fill(track->P(),sigTPC);      
      if(track->GetTOFBunchCrossing()==0) {
	fHistMomTOFBetaTrackQuality->Fill(track->P(),beta);
      }
      if( fUtils->isAcceptPid(track,AliPID::kPion) ){
	fHistMomTPCSignalTrackQualityPionPID->Fill(track->P(),sigTPC);	
	fHistMomTPCSigmaTrackQualityPionPID->Fill(track->P(),fUtils->getTPCsigma(track,AliPID::kPion));

	if(track->GetTOFBunchCrossing()==0) {
	  fHistMomTOFBetaTrackQualityPionPID->Fill(track->P(),beta);
	  fHistMomTOFSigmaTrackQualityPionPID->Fill(track->P(),fUtils->getTOFsigma(track,AliPID::kPion));
	}
      }
      else if( fUtils->isAcceptPid(track,AliPID::kKaon) ){
	fHistMomTPCSignalTrackQualityKaonPID->Fill(track->P(),sigTPC);
	fHistMomTPCSigmaTrackQualityKaonPID->Fill(track->P(),fUtils->getTPCsigma(track,AliPID::kKaon));
	if(track->GetTOFBunchCrossing()==0) {
	  fHistMomTOFBetaTrackQualityKaonPID->Fill(track->P(),beta);
	  fHistMomTOFSigmaTrackQualityKaonPID->Fill(track->P(),fUtils->getTOFsigma(track,AliPID::kKaon));
	}
      }
      else if( fUtils->isAcceptPid(track,AliPID::kProton) ){
	fHistMomTPCSignalTrackQualityProtonPID->Fill(track->P(),sigTPC);
	fHistMomTPCSigmaTrackQualityProtonPID->Fill(track->P(),fUtils->getTPCsigma(track,AliPID::kProton));
	if(track->GetTOFBunchCrossing()==0) {
	  fHistMomTOFBetaTrackQualityProtonPID->Fill(track->P(),beta);
	  fHistMomTOFSigmaTrackQualityProtonPID->Fill(track->P(),fUtils->getTOFsigma(track,AliPID::kProton));
	}
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
  
  return true;
}


AliEventPool* AliAnalysisTaskAODTrackPair::setPool()
{
  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi();

  AliEventPool* pool = (AliEventPool*)fPoolTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  
  return pool;
}

bool AliAnalysisTaskAODTrackPair::TrackPairAnalysisEveMixing(){
  
  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForMixing) ) return false;
  
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  AliEventPool* pool = setPool();  
  
  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);        
    //if(!fUtils->isAcceptTrack(track1)) continue;
    
    if (pool->IsReady()){
      
      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){
	
	TObjArray* poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);	
	
	for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){

	  AliAODTrack* __track2__ = (AliAODTrack*)poolTracks->At(iTrack2);
	  
	  AliAODTrack* track2 = (AliAODTrack*)__track2__->Clone();
	  
	  //if(!fUtils->isAcceptTrack(track2)) continue;
	  
	  delete track2;

	}//end of loop track2
      }// end of loop iMixEvt
    
    }//poolPion->IsReady()    
    
    fTrackArray->Add(track1);
    
  }//end of loop track1

  TObjArray* fTrackArrayClone = (TObjArray*)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if(fTrackArrayClone->GetEntriesFast()>0){
    pool->UpdatePool(fTrackArrayClone);
  }
  
}

bool AliAnalysisTaskAODTrackPair::V0TrackPairAnalysis()
{ 
  
  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForSame) ) return false;
  
  int nV0Pair = fEvent->GetNumberOfV0s();

  for(int iV0=0; iV0<nV0Pair; ++iV0){
    
    AliAODv0 *v0_1 = fEvent->GetV0(iV0);

    if(!fUtils->isAcceptV0asK0s(v0_1)) continue;
    
    double fill_1[]={v0_1->MassK0Short(),v0_1->Pt(),fUtils->getCentClass()};
    
    fSparseK0s->Fill(fill_1);
    
    if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;
    
    for(int jV0=iV0+1; jV0<nV0Pair; ++jV0){
    
      AliAODv0 *v0_2 = fEvent->GetV0(jV0);
    
      if( !fUtils->isAcceptV0asK0s(v0_2) ) continue;
      
      if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;
      
      //fUtils->isAcceptV0Pair(v0_1,v0_2);

      TLorentzVector k1;
      TLorentzVector k2;
      
      k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
      k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);
      
      TLorentzVector k12 = k1 + k2;

      double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};
      fSparseK0sK0s->Fill(fill_12);
  
    }

    
  }


}

