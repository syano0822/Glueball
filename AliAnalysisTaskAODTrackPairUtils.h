#ifndef AliAnalysisTaskAODTrackPairUtils_cxx
#define AliAnalysisTaskAODTrackPairUtils_cxx

#include "TNamed.h"
#include "TFile.h"
#include "AliMuonTrackCuts.h"
#include "AliVEventHandler.h"
#include "AliAnalysisMuonUtility.h"
#include <iostream>

class AliAnalysisTaskAODTrackPairUtils : public TNamed {
  
 public:
  
  AliAnalysisTaskAODTrackPairUtils();
  ~AliAnalysisTaskAODTrackPairUtils();

  
  bool setEvent(AliAODEvent* event, AliVEventHandler* handler);
  bool setMCArray(TClonesArray* array){fMCArray = array;}
  
  void setMC(bool isMC){fIsMC = isMC;}
  
  bool isAcceptEvent();
  bool isAcceptPrimeTrack(AliAODTrack* track, AliPID::EParticleType pid);
  bool isAcceptV0Track(AliAODTrack* track, AliPID::EParticleType pid);
  bool isAcceptTrackKinematics(AliAODTrack* track);
  bool isAcceptPrimeTrackQuality(AliAODTrack* track);
  bool isAcceptV0TrackQuality(AliAODTrack* track);
  bool isAcceptPid(AliAODTrack* track, AliPID::EParticleType pid);
  bool isAcceptTrackPair(float rap, float pt);
  bool isAcceptV0Pair(AliAODv0* v0_1,AliAODv0* v0_2);
  bool isAcceptDefaultV0(AliAODv0* v0);
  bool isK0sCandidate(float mass);
  bool isK0sV0(AliAODv0 *v0);
  bool isSameMother(AliAODTrack* track1,AliAODTrack* track2);
  
  bool hasTOF(AliAODTrack* track);

  bool isMC(){
    return fIsMC;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //Set analysis cut flags
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  void setVertexCut(float min, float max)
  {
    fMinVertexCutZ = min;
    fMaxVertexCutZ = max;
    fIsVtxZcut = true;
  }  
  void setPairRapidityCut(float min, float max)
  {
    fMinPairRapCut = min;
    fMaxPairRapCut = max;
    fIsPairRapCut = true;
  }  
  void setPileupRejectionCut(bool flag)
  {
    fIsPUcut = flag;
  }  
  void setPionSelectSigmaTPC(float min, float max)
  {
    fMinPionSigmaTPC = min;
    fMaxPionSigmaTPC = max;
  }
  void setPionSelectSigmaTOF(float min, float max)
  {
    fMinPionSigmaTOF = min;
    fMaxPionSigmaTOF = max;
  }
  void setMinNClustTPC(float min)
  {
    fMinNClustTPC = min;
  }
  void setMinTPCCrossRows(float cross)
  {
    fMinCrossRows = cross;
  }
  void setMinFindableTPCclsFrac(float min)
  {
    fMinFindableTPCclsFraction = min;
  }
  void setTrackChi2perNDF(float max){
    fMaxTrackChi2perNDF = max;
  }
  void setTPCTrackChi2perNDF(float max){
    fMaxTPCTrackChi2perNDF = max;
  }
  void setITSTrackChi2perNDF(float max){
    fMaxITSTrackChi2perNDF = max;
  }

  void setPrimeTrackDCAFunction(std::string func)
  {
    fPrimeTrackDCAFunc = new TF1("fPrimeTrackDCAFunc",func.c_str(),0,1000);
  }
  void setPrimeTrackFilterBit(int bit)
  {
    fPrimeTrackFilterBit = bit;
  }
  void setPromenterosLimit(double pcm, double r0, double width){
    double e = 0.0001;
    fMinArmenterosLine = new TF1("fMinArmenterosLine","[0]*pow(1-pow(x/[1],2),1./2.)*sqrt(1.-[2])",-r0+e,r0-e);
    fMinArmenterosLine->SetParameters(pcm,r0,width);
    fMaxArmenterosLine = new TF1("fMaxArmenterosLine","[0]*pow(1-pow(x/[1],2),1./2.)*sqrt(1.+[2])",-r0+e,r0-e);
    fMaxArmenterosLine->SetParameters(pcm,r0,width);
  }
  void setTrackKinematicRange(float minpt,float maxpt,float mineta, float maxeta)
  {
    fMinTrackPt = minpt;
    fMaxTrackPt = maxpt;
    fMinTrackEta = mineta;
    fMaxTrackEta = maxeta;
  }
  void setKaonSelectSigmaTPC(float min, float max)
  {
    fMinKaonSigmaTPC = min;
    fMaxKaonSigmaTPC = max;
  }
  void setKaonSelectSigmaTOF(float min, float max)
  {
    fMinKaonSigmaTOF = min;
    fMaxKaonSigmaTOF = max;
  }  
  void setProtonSelectSigmaTPC(float min, float max)
  {
    fMinProtonSigmaTPC = min;
    fMaxProtonSigmaTPC = max;
  }
  void setProtonSelectSigmaTOF(float min, float max)
  {
    fMinProtonSigmaTOF = min;
    fMaxProtonSigmaTOF = max;
  }
  void setPropLifeTime(float min, float max)
  {
    fMinPropLifetime = min;
    fMaxPropLifetime = max;
  }
  void setMinV0Radius(float min,float max)
  {
    fMinV0Radius = min;
    fMaxV0Radius = max;
  }
  void setV0TrackDCA(float min, float max)
  {
    fMinV0TrackDCA = min;
    fMaxV0TrackDCA = max;
  }

  void setMinV0PosTrackDCA(float min,float max)
  {
    fMinPosDCA = min;
    fMaxPosDCA = max;
  }
  void setMinV0NegTrackDCA(float min,float max)
  {
    fMinNegDCA = min;
    fMaxNegDCA = max;
  }
  void setMinV0CosPoint(float min,float max)
  {
    fMinCosPoint = min;
    fMaxCosPoint = max;
  }
  void setRejLambdaMassRange(float min, float max)
  {
    fMinRejLambdaMass=min;
    fMaxRejLambdaMass=max;
  }
  void setTPCTOFCombPID(float pion, float kaon, float proton)
  {
    fOnTPCTOFCombPID = true;
    fPionSigmaCombSigma = pion;
    fKaonSigmaCombSigma = kaon;
    fProtonSigmaCombSigma = proton;
  }
  void setRequireTOFhit(bool flag)
  {
    fRequireTOFhit = flag;
  } 

  //////////////////////////////////////////////////////////////////////////////////////////////
  //Set analysis object
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  void setMultiEstimateMethod(std::string method)
  {
    fMultiMethod = method;
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////
  //Get the analysis variables
  //////////////////////////////////////////////////////////////////////////////////////////////

  float getVtxZ()
  {
    return fVtxZ;
  }
  float getCentClass(int spec)
  {
    if(spec==0)
      return fCentSPDTrk;
    else if(spec==1)
      return fCentV0M;
    else if(spec==2)
      return fCentV0A;
    else if(spec==3)
      return fCentV0C;
    else
      return -999;
  }
  float getCentClass()
  {
    return fCent;
  }
  float getPsi()
  {
    return fPsi;
  }
  int getRunnumber(){
    return fRunNumber;
  }  
  int getNSPDTrkInfo(int spec){
    if(spec==0)
      return fNSPDTrk05;
    else if(spec==1)
      return fNSPDTrk10;
    else if(spec==2)
      return fNSPDTrk15;
    else if(spec==3)
      return fNSPDTrk20;
    else if(spec==4)
      return fNSPDTrkAll;
    else 
      return 0;
  }
  int getNSPDClustInfo(int spec){
    if(spec==0)
      return fNClustSPD1;
    else if(spec==1)
      return fNClustSPD2;
    else 
      return 0;
  }

  float getVzeroInfo(int spec)
  {
    if(spec==0)
      return fChV0A;
    else if(spec==1)
      return fChV0C;
    else if(spec==2)
      return fChV0M;
    else if(spec==3)
      return fTimeV0A-fTimeV0C;
    else if(spec==4)
      return fTimeV0A+fTimeV0C;
    else
      return 0;
  }
  
  float getTPCsigma(AliAODTrack* track, AliPID::EParticleType pid)
  {
    if(!track || !fPIDResponse) return -999;
    return fPIDResponse->NumberOfSigmasTPC(track, pid);
  }
  float getTOFsigma(AliAODTrack* track, AliPID::EParticleType pid)
  {
    if(!track || !fPIDResponse) return -999;
    return fPIDResponse->NumberOfSigmasTOF(track, pid);
  }
  
  void Initialization();
  void setInit();
  bool isSameRunnumber();
  bool setVtxZCentPsi();

  bool setSPDTrk();
  bool setSPDClust();
  bool setVZERO();

  AliAODEvent* fEvent;
  AliMultSelection* fMultSelection;  
  AliVEventHandler* fInputHandler;
  AliPIDResponse* fPIDResponse;
  TClonesArray* fMCArray;

  int fRunNumber;
  
  bool fIsMC;
 
  bool fIsVtxZcut;
  float fMaxVertexCutZ;
  float fMinVertexCutZ;

  bool fIsPairRapCut;
  float fMinPairRapCut;
  float fMaxPairRapCut;
  
  bool fIsPUcut;  
  
  int fPrimeTrackFilterBit;
  TF1* fPrimeTrackDCAFunc;
  float fMaxTrackChi2perNDF;
  float fMaxTPCTrackChi2perNDF;
  float fMaxITSTrackChi2perNDF;
  float fMinCrossRows;
  float fMinNClustTPC;
  float fMinFindableTPCclsFraction;
  float fMinTrackPt;
  float fMaxTrackPt;
  float fMinTrackEta;
  float fMaxTrackEta;

  float fMinMassK0sCand;
  float fMaxMassK0sCand;
  
  float fMinPropLifetime;
  float fMaxPropLifetime;  
  float fMinV0Radius;
  float fMaxV0Radius;
  float fMinV0TrackDCA;
  float fMaxV0TrackDCA;
  float fMinPosDCA;
  float fMaxPosDCA;
  float fMinNegDCA;
  float fMaxNegDCA;
  float fMinCosPoint;
  float fMaxCosPoint;
  float fMaxRejLambdaMass;
  float fMinRejLambdaMass;

  bool fRequireTOFhit;
  bool fOnTPCTOFCombPID;
  float fPionSigmaCombSigma;
  float fKaonSigmaCombSigma;
  float fProtonSigmaCombSigma;

  float fMinPionSigmaTPC;
  float fMaxPionSigmaTPC;
  float fMinPionSigmaTOF;
  float fMaxPionSigmaTOF;

  float fMinKaonSigmaTPC;
  float fMaxKaonSigmaTPC;
  float fMinKaonSigmaTOF;
  float fMaxKaonSigmaTOF;

  float fMinProtonSigmaTPC;
  float fMaxProtonSigmaTPC;
  float fMinProtonSigmaTOF;
  float fMaxProtonSigmaTOF;

  std::string fMultiMethod;
  
  float fVtxXYZ[3];
  float fVtxZ;
  float fCent;
  float fPsi;

  float fCentSPDTrk;
  float fCentV0A;
  float fCentV0C;
  float fCentV0M;

  int fNSPDTrk05;
  int fNSPDTrk10;
  int fNSPDTrk15;
  int fNSPDTrk20;
  int fNSPDTrkAll;
  
  int fNClustSPD1;
  int fNClustSPD2;

  float fChV0A;
  float fChV0C;
  float fChV0M;
  float fTimeV0A;
  float fTimeV0C;

  float fMassK0s = 0.497614;
  float fMassKaon = 0.493677;
  float fMassPion = 0.13957061;

  TF1* fMinArmenterosLine;
  TF1* fMaxArmenterosLine;
  double fArmenterosBandWidth;

  ClassDef(AliAnalysisTaskAODTrackPairUtils, 1); // example of analysis
};

#endif
