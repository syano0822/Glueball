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
  fSparseK0sK0s(NULL)
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
  fSparseK0sK0s(NULL)
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
  
  double bins_event_hist[]={0,1,2,3,4,5,6,7,8,9,10};
  int binnum_event_hist = sizeof(bins_event_hist)/sizeof(double) - 1;

  std::string event_label[]={"CMUL7","CMLL7","CMUL7orCMLL7","CMUL7andCMLL7","CMUL7withDS","CMLL7withDS","CMUL7orCMLL7withDS","CMUL7andCMLL7withDS"};
  fEventCounter = new TH2F("fEventCounter","",binnum_event_hist,bins_event_hist,binnum_cent_hist,bins_cent_hist);
  for(int iname=0; iname<sizeof(event_label)/sizeof(std::string); ++iname)
    fEventCounter->GetXaxis()->SetBinLabel(iname+1,event_label[iname].c_str());
  fOutputList->Add(fEventCounter);

  const Int_t binnum_Cent            = 3;
  Int_t bins_Cent[binnum_Cent]       = {200,100,binnum_cent_hist};
  double minbins_Cent[binnum_Cent] = {0.4,0,0};
  double maxbins_Cent[binnum_Cent] = {0.6,10,100};
  TString  namebins_Cent[]        = {"mass","rap","pt","centrality"};

  fSparseK0s = new THnSparseF("fSparseK0s","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseK0sK0s = new THnSparseF("fSparseK0sK0s","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  
  fOutputList->Add(fSparseK0s);
  fOutputList->Add(fSparseK0sK0s);

  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{  
  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;
  V0TrackPairAnalysis();
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

    if(!fUtils->isAcceptK0sTrackPair(v0_1)) continue;
    
    double fill_1[]={v0_1->MassK0Short(),v0_1->Pt(),fUtils->getCentClass()};
    
    fSparseK0s->Fill(fill_1);
    
    if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;

    for(int jV0=iV0+1; jV0<nV0Pair; ++jV0){
    
      AliAODv0 *v0_2 = fEvent->GetV0(jV0);

      if(!fUtils->isAcceptK0sTrackPair(v0_2)) continue;
      
      if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;
            
      TLorentzVector k1;
      TLorentzVector k2;
      
      k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
      k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);
      
      TLorentzVector k12 = k1 + k2;

      double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};

      fSparseK0s->Fill(fill_12);
  
    }

    
  }


}

