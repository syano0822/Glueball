#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliEventplane.h"
#include "AliVEventHandler.h"
#include "AliAODVertex.h"
#include "AliMultSelection.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"

#include "AliESDtrack.h"


#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include "iostream"

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPairUtils)       

AliAnalysisTaskAODTrackPairUtils::AliAnalysisTaskAODTrackPairUtils() : TNamed(),  
  fEvent(NULL),
  fMultSelection(NULL),
  fPIDResponse(NULL),
  fInputHandler(NULL),
  fRunNumber(-99999),
  fIsMC(false),
  fIsVtxZcut(true),
  fMaxVertexCutZ(-10),
  fMinVertexCutZ(10),
  fIsPairRapCut(true),
  fMinPairRapCut(-4.0),
  fMaxPairRapCut(-2.5),
  fIsPUcut(true),

  fMinTrackPt(0.3),
  fMinTrackEta(-0.8),
  fMaxTrackEta(0.8),

  fMinMassK0sCand(0.49),
  fMaxMassK0sCand(0.51),
  
  fMinV0Radius(0.5),
  fMinPosDCA(0.06),
  fMinNegDCA(0.06),
  fMinCosPoint(0.97),
  fMinRejLambdaMass(1.105),
  fMaxRejLambdaMass(1.1255),

  fMinPionSigmaTPC(-3.),
  fMaxPionSigmaTPC(3),
  fMinPionSigmaTOF(-3),
  fMaxPionSigmaTOF(3),

  fMinKaonSigmaTPC(-3),
  fMaxKaonSigmaTPC(3),
  fMinKaonSigmaTOF(-3),
  fMaxKaonSigmaTOF(3),

  fMinProtonSigmaTPC(-3),
  fMaxProtonSigmaTPC(3),
  fMinProtonSigmaTOF(-3),
  fMaxProtonSigmaTOF(3),

  fMultiMethod("SPDtracklet"),

  fVtxXYZ(),
  fVtxZ(0),
  fCent(0),
  fPsi(0),
  
  fCentSPDTrk(0.),
  fCentV0A(0.),
  fCentV0C(0.),
  fCentV0M(0.),

  fNSPDTrk05(0),
  fNSPDTrk10(0),
  fNSPDTrk15(0),
  fNSPDTrk20(0),
  fNSPDTrkAll(0),

  fNClustSPD1(0),
  fNClustSPD2(0),

  fChV0A(0.),
  fChV0C(0.),
  fChV0M(0.),
  fTimeV0A(0.),
  fTimeV0C(0.)
{

}

AliAnalysisTaskAODTrackPairUtils::~AliAnalysisTaskAODTrackPairUtils()
{
 
  
 
}

void AliAnalysisTaskAODTrackPairUtils::setInit()
{

  fEvent=NULL;
  fMultSelection=NULL;

  fVtxXYZ[0] = -999;
  fVtxXYZ[1] = -999;
  fVtxXYZ[2] = -999;

  fVtxZ=-999;
  fCent=-999;
  fPsi=-999;
  
  fCentSPDTrk=-999;
  fCentV0A=-999;
  fCentV0C=-999;
  fCentV0M=-999;
  
  fNSPDTrk05=0;
  fNSPDTrk10=0;
  fNSPDTrk15=0;
  fNSPDTrk20=0;
  fNSPDTrkAll=0;

  fNClustSPD1=0;
  fNClustSPD2=0;

  fChV0A=0.;
  fChV0C=0.;
  fChV0M=0.;

  fTimeV0A=0.;
  fTimeV0C=0.;
}

bool AliAnalysisTaskAODTrackPairUtils::setEvent(AliAODEvent* event, AliVEventHandler* handler)
{  
  
  setInit();
  
  fEvent = event;  
  if(!fEvent) return false;
  
  fInputHandler = handler;
  if(!fInputHandler) return false;
  
  if(!fPIDResponse){
    fPIDResponse=fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return false;
  }  

  fMultSelection = (AliMultSelection *)fEvent->FindListObject("MultSelection");  
  if(!fMultSelection) return false;

  fRunNumber = fEvent->GetRunNumber();
  
  if(!setVtxZCentPsi())
    return false;  
  if(!setSPDTrk())
    return false;
  if(!setSPDClust())
    return false;
  if(!setVZERO())
    return false;
  
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isSameRunnumber()
{
  if(fRunNumber == fEvent->GetRunNumber()) return true;
  else return false;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptEvent()
{  
  
  if(!fEvent) return false;

  AliAODVertex * vtx  = fEvent->GetPrimaryVertex();

  if(fIsVtxZcut && (fMinVertexCutZ>vtx->GetZ() || vtx->GetZ()>fMaxVertexCutZ))
    return false;
  if(fIsVtxZcut && vtx->GetNContributors()<2)
    return false;
  if(fIsPUcut && fEvent->IsPileupFromSPDInMultBins())
    return false;
  if(fMultSelection && fMultSelection->GetMultiplicityPercentile(fMultiMethod,false)<0 && fMultSelection->GetMultiplicityPercentile(fMultiMethod,false)>100.){
    return false;
  }
  
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptTrack(AliAODTrack* track,AliPID::EParticleType pid){  
  
  if(!isAcceptTrackQuality(track)) return false;
  if(!isAcceptPid(track,pid)) return false;

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptTrackQuality(AliAODTrack* track){  
  
  if(track->Pt()<fMinTrackPt) return false;
  if(track->Eta()<fMinTrackEta && fMaxTrackEta<track->Eta()) return false;

  float crossedRows = track->GetTPCClusterInfo(2,1);

  if( !(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  if ( crossedRows < 70 ) return false;
  
  if( track->GetTPCNclsF()<=0 ) return false;

  float crossedRowsOverFindable = crossedRows / ((double)(track->GetTPCNclsF()));

  if( crossedRowsOverFindable< 0.8 ) return false;

  bool hasITShit = false;  
  if( track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1) ||
      track->HasPointOnITSLayer(4) || track->HasPointOnITSLayer(5))
    hasITShit = true;
  
  bool hasTOFhit = false;
  if( track->GetTOFBunchCrossing()==0 ) hasTOFhit=true;
  if(!hasITShit && !hasTOFhit) return false;
  

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptPid(AliAODTrack* track, AliPID::EParticleType pid){
  
  
  float fPionSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, pid);
  float fPionSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, pid);
  
  bool hasTOFhit = false;
  if( track->GetTOFBunchCrossing()==0 ) hasTOFhit=true;
  
  if(hasTOFhit){
    if(fMinPionSigmaTPC>fPionSigmaTPC || fPionSigmaTPC>fMaxPionSigmaTPC) return false;  
    if(fMinPionSigmaTOF>fPionSigmaTOF || fPionSigmaTPC>fMaxPionSigmaTOF) return false;
  }
  else{
    if(fMinPionSigmaTPC>fPionSigmaTPC || fPionSigmaTPC>fMaxPionSigmaTPC) return false;  
  }  

  return true;
}


bool AliAnalysisTaskAODTrackPairUtils::isK0sCandidate(float mass){
  
  if(fMinMassK0sCand<mass && mass<fMaxMassK0sCand){
    return true;
  }

  return false;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptTrackPair(AliAODTrack* track1, AliAODTrack* track2)
{      
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptV0Pair(AliAODv0* v0_1,AliAODv0* v0_2)
{      
  
  AliAODTrack *pTrack1=(AliAODTrack *)v0_1->GetDaughter(0); //0->Positive Daughter
  AliAODTrack *nTrack1=(AliAODTrack *)v0_1->GetDaughter(1); //1->Negative Daughter

  AliAODTrack *pTrack2=(AliAODTrack *)v0_2->GetDaughter(0); //0->Positive Daughter
  AliAODTrack *nTrack2=(AliAODTrack *)v0_2->GetDaughter(1); //1->Negative Daughter
  
  cout<<pTrack1->Pt()<<"   "<<nTrack1->Pt()<<"   "<<pTrack2->Pt()<<"   "<<nTrack2->Pt()<<endl;
  

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptK0sTrackPair(AliAODv0* v0)
{      

  AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
  AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
  
  float DcaPosToPrimVertex = v0->DcaPosToPrimVertex();
  float DcaNegToPrimVertex = v0->DcaNegToPrimVertex();
  
  if(!isAcceptTrack(pTrack,AliPID::kPion)) return false;
  if(!isAcceptTrack(nTrack,AliPID::kPion)) return false;
  
  int OnFlyStatus = v0->GetOnFlyStatus();

  float chi2 = v0->Chi2V0();
  float DcaV0Daughters = v0->DcaV0Daughters();
  float DcaV0ToPrimVertex = v0->DcaV0ToPrimVertex();
  float V0CosineOfPointingAngle = v0->CosPointingAngle((double*)fVtxXYZ);
  
  double decayV[3]={};
  v0->GetXYZ(decayV);

  float radius = sqrt(pow(decayV[0],2) + pow(decayV[1],2));
  
  if(radius < fMinV0Radius) return false;
  if(DcaPosToPrimVertex < fMinPosDCA) return false;
  if(DcaNegToPrimVertex < fMinNegDCA) return false;
  if(V0CosineOfPointingAngle<fMinCosPoint) return false;
  if(DcaV0Daughters>1.0) return false;

  double rap = v0->RapK0Short();
  
  if(fIsPairRapCut & !(fMinPairRapCut<rap && rap<fMaxPairRapCut)) return false;

  float massLambda = v0->MassLambda();
  if(fMinRejLambdaMass<massLambda && massLambda<fMaxRejLambdaMass) return false;

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::setVtxZCentPsi()
{
  
  if(!fEvent) return false;
  if(!fMultSelection) return false;

  AliAODVertex* vtx = (AliAODVertex*)fEvent->GetPrimaryVertex();
  
  if(!vtx) return false;

  fVtxZ = vtx->GetZ();
  fCent = fMultSelection->GetMultiplicityPercentile(fMultiMethod,false);
  fPsi  = 0;
  
  fVtxXYZ[0] = vtx->GetX();
  fVtxXYZ[1] = vtx->GetY();
  fVtxXYZ[2] = vtx->GetZ();

  fCentSPDTrk = fMultSelection->GetMultiplicityPercentile("SPDTracklets", false);
  fCentV0M    = fMultSelection->GetMultiplicityPercentile("V0M", false);
  fCentV0A    = fMultSelection->GetMultiplicityPercentile("V0A", false);
  fCentV0C    = fMultSelection->GetMultiplicityPercentile("V0C", false);

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::setSPDTrk(){

  if(!fEvent) return false;

  AliAODTracklets *tracklet = (AliAODTracklets*)fEvent->GetTracklets();

  if(!tracklet) return false;

  if(tracklet){
    for(int iTrk=0; iTrk < tracklet->GetNumberOfTracklets(); ++iTrk){
      Double_t eta = -1*TMath::Log(TMath::Tan(tracklet->GetTheta(iTrk)/2.));
      Double_t phi = tracklet->GetPhi(iTrk);
      if(fabs(eta)<0.5) ++fNSPDTrk05;
      if(fabs(eta)<1.0) ++fNSPDTrk10;
      if(fabs(eta)<1.5) ++fNSPDTrk15;
      if(fabs(eta)<2.0) ++fNSPDTrk20;
      ++fNSPDTrkAll;
    }
  }
  
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::setSPDClust()
{
  AliAODHeader *header = (AliAODHeader*)fEvent->GetHeader();
  if(!header) return false;
  fNClustSPD1 = header->GetNumberOfITSClusters(0);
  fNClustSPD2 = header->GetNumberOfITSClusters(1);
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::setVZERO(){
  
  AliAODVZERO * aodV0 = fEvent->GetVZEROData();  
  if(!aodV0) return false;

  fChV0A = aodV0->GetMTotV0A();
  fChV0C = aodV0->GetMTotV0C();
  fChV0M = fChV0A + fChV0C;
  fTimeV0A = aodV0->GetV0ATime();
  fTimeV0C = aodV0->GetV0CTime();
  
  return true;
}
