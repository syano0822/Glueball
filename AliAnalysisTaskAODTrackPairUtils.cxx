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
  fMCArray(NULL),
  fRunNumber(-99999),
  fIsMC(false),
  fIsVtxZcut(true),
  fMaxVertexCutZ(-10),
  fMinVertexCutZ(10),
  fIsPairRapCut(true),
  fMinPairRapCut(0.5),
  fMaxPairRapCut(-0.5),
  fIsPUcut(true),

  fPrimeTrackDCAFunc(NULL),
  fPrimeTrackFilterBit(5),
  fMaxTrackChi2perNDF(4.),
  fMaxITSTrackChi2perNDF(36.),
  fMaxTPCTrackChi2perNDF(5.),
  fMinFindableTPCclsFraction(0.8),
  fMinCrossRows(70),
  fMinNClustTPC(70),
  fMinTrackPt(0.15),
  fMaxTrackPt(999.),
  fMinTrackEta(-0.8),
  fMaxTrackEta(0.8),

  fMinMassK0sCand(0.49),
  fMaxMassK0sCand(0.51),
								       
  fMinPropLifetime(0.0),
  fMaxPropLifetime(20.0),
  fMinV0Radius(0.5),
  fMaxV0Radius(999.),
  fMinV0TrackDCA(0.06),
  fMaxV0TrackDCA(999.),
  fMinPosDCA(0.06),
  fMaxPosDCA(999.),
  fMinNegDCA(0.06),
  fMaxNegDCA(999.),
  fMinCosPoint(0.97),
  fMaxCosPoint(1.00),
  fMinRejLambdaMass(1.105),
  fMaxRejLambdaMass(1.1255),
  
  fRequireTOFhit(false),
  fOnTPCTOFCombPID(true),
  fPionSigmaCombSigma(2.),
  fKaonSigmaCombSigma(2.),
  fProtonSigmaCombSigma(2.),

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
  fTimeV0C(0.),
  
  fArmenterosBandWidth(0.05),
  fMinArmenterosLine(NULL),
  fMaxArmenterosLine(NULL)
{
  
}

AliAnalysisTaskAODTrackPairUtils::~AliAnalysisTaskAODTrackPairUtils()
{
 
  
 
}

void AliAnalysisTaskAODTrackPairUtils::Initialization()
{
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"                               CONFIGURATION PARAMS                             "<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"fMaxVertexCutZ:        "<<fMaxVertexCutZ<<endl;
  cout<<"fMinVertexCutZ:        "<<fMinVertexCutZ<<endl;
  cout<<"fMinPairRapCut:        "<<fMinPairRapCut<<endl;
  cout<<"fMaxPairRapCut:        "<<fMaxPairRapCut<<endl;  
  cout<<"fMaxTrackChi2perNDF:   "<<fMaxTrackChi2perNDF<<endl;
  cout<<"fMinCrossRows:         "<<fMinCrossRows<<endl;
  cout<<"fMinTrackPt:           "<<fMinTrackPt<<endl;
  cout<<"fMaxTrackPt:           "<<fMaxTrackPt<<endl;
  cout<<"fMinTrackEta:          "<<fMinTrackEta<<endl;
  cout<<"fMaxTrackEta:          "<<fMaxTrackEta<<endl;  
  cout<<"fMinMassK0sCand:       "<<fMinMassK0sCand<<endl;
  cout<<"fMaxMassK0sCand:       "<<fMaxMassK0sCand<<endl;  
  cout<<"fMinPropLifetime:      "<<fMinPropLifetime<<endl;
  cout<<"fMaxPropLifetime:      "<<fMaxPropLifetime<<endl;
  cout<<"fMinV0Radius:          "<<fMinV0Radius<<endl;
  cout<<"fMaxV0Radius:          "<<fMaxV0Radius<<endl;
  cout<<"fMinPosDCA:            "<<fMinPosDCA<<endl;
  cout<<"fMaxPosDCA:            "<<fMaxPosDCA<<endl;
  cout<<"fMinNegDCA:            "<<fMinNegDCA<<endl;
  cout<<"fMaxNegDCA:            "<<fMaxNegDCA<<endl;
  cout<<"fMinCosPoint:          "<<fMinCosPoint<<endl;
  cout<<"fMaxCosPoint:          "<<fMaxCosPoint<<endl;
  cout<<"fMinRejLambdaMass:     "<<fMinRejLambdaMass<<endl;
  cout<<"fMaxRejLambdaMass:     "<<fMaxRejLambdaMass<<endl;  
  cout<<"fPionSigmaCombSigma :  "<<fPionSigmaCombSigma<<endl;
  cout<<"fKaonSigmaCombSigma:   "<<fKaonSigmaCombSigma<<endl;
  cout<<"fProtonSigmaCombSigma: "<<fProtonSigmaCombSigma<<endl;  
  cout<<"fMultiMethod:          "<<fMultiMethod<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;
  cout<<"================================================================================"<<endl;  
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

bool AliAnalysisTaskAODTrackPairUtils::isAcceptPrimeTrack(AliAODTrack* track,AliPID::EParticleType pid){  
  
  if(!isAcceptPrimeTrackQuality(track)) return false;
  if(!isAcceptPid(track,pid)) return false;

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptV0Track(AliAODTrack* track,AliPID::EParticleType pid){  
  
  if(!isAcceptV0TrackQuality(track)) return false;
  if(!isAcceptPid(track,pid)) return false;

  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptTrackKinematics(AliAODTrack* track){  
  if(track->Pt()<fMinTrackPt || fMaxTrackPt<track->Pt()) return false;
  if(track->Eta()<fMinTrackEta || fMaxTrackEta<track->Eta()) return false;  
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptPrimeTrackQuality(AliAODTrack* track){  
  
  return track->TestFilterBit(AliAODTrack::kTrkGlobal) ? true : false;
  /*
  if( !track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) ) return false;    
  if( track->GetTPCNcls() < 1 ) return false;
  if( !isAcceptTrackKinematics(track) ) return false;  
  if( track->Chi2perNDF() > fMaxTrackChi2perNDF ) return false;
  if( track->GetITSchi2()/track->GetITSNcls() > fMaxITSTrackChi2perNDF ) return false;
  if( track->GetTPCchi2()/track->GetTPCNcls() > fMaxTPCTrackChi2perNDF ) return false;  
  if( !(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  if( !(track->GetStatus() & AliESDtrack::kITSrefit)) return false;  
  if( track->GetTPCNcls() < fMinNClustTPC ) return false;  
  if( track->GetTPCNclsF()<=0 ) return false;
  float crossedRows = track->GetTPCNCrossedRows();
  float crossedRowsOverFindable = crossedRows / ((double)(track->GetTPCNclsF()));
  if( crossedRowsOverFindable < fMinFindableTPCclsFraction ) return false;

  if( fRequireTOFhit && track->GetTOFBunchCrossing()!=0 ) return false;
  if( track->GetKinkIndex(0) > 0 ) return false;  

  float max_dca = fPrimeTrackDCAFunc->Eval(track->Pt());  
  float dca_xy=9999;
  float dca_z=9999;
  track->GetImpactParameters(dca_xy,dca_z);
  float dca = sqrt(dca_xy*dca_xy + dca_z*dca_z);
  
  if( fabs(dca_xy)>max_dca ) return false;
  if( fabs(dca_z) > 2.0 ) return false;
  */
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptV0TrackQuality(AliAODTrack* track){  
  
  if( track->GetTPCNcls() < 1 ) return false;
  if( !isAcceptTrackKinematics(track) ) return false; 
  if( !(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  if( track->GetTPCClusterInfo(2,1) < fMinCrossRows ) return false;
  if( track->GetTPCNclsF()<=0 ) return false;
  float crossedRows = track->GetTPCNCrossedRows();
  float crossedRowsOverFindable = crossedRows / ((double)(track->GetTPCNclsF()));
  if( crossedRowsOverFindable < fMinFindableTPCclsFraction ) return false;

  if( track->GetKinkIndex(0) > 0 ) return false;  

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
  
  float sigmaRadiusComb=0;
  float minSigmaRangeTPC=0;
  float maxSigmaRangeTPC=0;
  float minSigmaRangeTOF=0;
  float maxSigmaRangeTOF=0;

  if(pid == AliPID::kPion){
    sigmaRadiusComb = fPionSigmaCombSigma;
    minSigmaRangeTPC = fMinPionSigmaTPC;
    maxSigmaRangeTPC = fMaxPionSigmaTPC;
    minSigmaRangeTOF = fMinPionSigmaTOF;
    maxSigmaRangeTOF = fMaxPionSigmaTOF;
  }
  else if(pid == AliPID::kKaon){
    sigmaRadiusComb = fKaonSigmaCombSigma;
    minSigmaRangeTPC = fMinKaonSigmaTPC;
    maxSigmaRangeTPC = fMaxKaonSigmaTPC;
    minSigmaRangeTOF = fMinKaonSigmaTOF;
    maxSigmaRangeTOF = fMaxKaonSigmaTOF;    
  }
  else if(pid == AliPID::kProton){
    sigmaRadiusComb = fProtonSigmaCombSigma;
    minSigmaRangeTPC = fMinProtonSigmaTPC;
    maxSigmaRangeTPC = fMaxProtonSigmaTPC;
    minSigmaRangeTOF = fMinProtonSigmaTOF;
    maxSigmaRangeTOF = fMaxProtonSigmaTOF;
  }
  else{
    return false;
  }

  float fSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, pid);
  float fSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, pid);
  
  float fSigmaTPCTOF = sqrt(fSigmaTPC*fSigmaTPC + fSigmaTOF*fSigmaTOF);
  
  bool hasTOFhit = true;
  
  if( fabs(fSigmaTOF)>990 ) hasTOFhit = false;
  
  if(hasTOFhit){    
    if(minSigmaRangeTPC<fSigmaTPC && fSigmaTPC<maxSigmaRangeTPC && minSigmaRangeTOF<fSigmaTOF && fSigmaTOF<maxSigmaRangeTOF) return true;
    else return false;
  }
  else{
    if(minSigmaRangeTPC<fSigmaTPC && fSigmaTPC<maxSigmaRangeTPC) return true;
    else return false;
  }

}

bool AliAnalysisTaskAODTrackPairUtils::isK0sCandidate(float mass){
  
  if(fMinMassK0sCand<mass && mass<fMaxMassK0sCand){
    return true;
  }

  return false;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptTrackPair(float rap, float pt)
{      
  if(fMinPairRapCut>rap || rap<fMaxPairRapCut) return false;
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptV0Pair(AliAODv0* v0_1,AliAODv0* v0_2)
{      
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isAcceptDefaultV0(AliAODv0* v0)
{      
  
  if(!v0) return false;
  
  if(v0->GetNDaughters() != 2) return false;
  if(v0->GetNProngs() != 2) return false;
  if(v0->GetCharge() != 0) return false;
  if(v0->GetOnFlyStatus()) return false;
  
  float V0CosineOfPointingAngle = v0->CosPointingAngle((double*)fVtxXYZ);
  if(V0CosineOfPointingAngle<fMinCosPoint) return false;
  
  float DcaV0Daughters = v0->DcaV0Daughters();
  if(DcaV0Daughters>1.5) return false;
  
  double decayV[3]={};
  v0->GetXYZ(decayV);
  
  float radius = sqrt(pow(decayV[0]-fVtxXYZ[0],2) + pow(decayV[1]-fVtxXYZ[1],2));
  if(radius < fMinV0Radius || radius>fMaxV0Radius) return false;
  
  float length = sqrt(pow(decayV[0]-fVtxXYZ[0],2) + pow(decayV[1]-fVtxXYZ[1],2) + pow(decayV[2]-fVtxXYZ[2],2));
  float proper_length = v0->MassK0Short() * length/v0->P();   
  if(fMinPropLifetime>proper_length || proper_length>fMaxPropLifetime) return false;
  
  return true;
}

bool AliAnalysisTaskAODTrackPairUtils::isK0sV0(AliAODv0 *v0)
{

  AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0);
  AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1);  
  if(!pTrack) return false;
  if(!nTrack) return false;
  
  if(!isAcceptV0Track(pTrack,AliPID::kPion)) return false;
  if(!isAcceptV0Track(nTrack,AliPID::kPion)) return false;

  double rap = v0->RapK0Short();
  if( fIsPairRapCut & !(fMinPairRapCut<rap && rap<fMaxPairRapCut) ) return false;
  
  double alpha = v0->AlphaV0();
  double pt = v0->PtArmV0();
  
  if(fabs(alpha)>0.84) return false;

  double min_pt = fMinArmenterosLine->Eval(alpha);
  double max_pt = fMaxArmenterosLine->Eval(alpha);
  
  if(min_pt>pt || pt>max_pt){
    return false;
  }
  
  return true;
}

int AliAnalysisTaskAODTrackPairUtils::isSameMother(AliAODTrack* track1, AliAODTrack* track2)
{
  
  if(!track1) return -999;
  if(!track2) return -999;

  int label1 = track1->GetLabel();
  int label2 = track2->GetLabel();

  if(label1<1) return -999;
  if(label2<1) return -999;

  AliAODMCParticle *particle1 = (AliAODMCParticle*)fMCArray->At(label1);
  AliAODMCParticle *particle2 = (AliAODMCParticle*)fMCArray->At(label2);
  if(!particle1) return -999;
  if(!particle2) return -999;
  
  int mom_label1 = particle1->GetMother();
  int mom_label2 = particle2->GetMother();

  if(mom_label1<1) return -999;
  if(mom_label2<1) return -999;

  if(mom_label1 != mom_label2) return -999;
  
  AliAODMCParticle *mom_particle = (AliAODMCParticle*)fMCArray->At(mom_label1);
  if(!mom_particle) return -999;

  return mom_particle->GetPdgCode();
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

bool AliAnalysisTaskAODTrackPairUtils::hasTOF(AliAODTrack* track){
  
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  float len = track->GetIntegratedLength();

  if( !hasTOFout || !hasTOFtime || (len < 350.) ){
    return false;
  }

  return true;  
}
