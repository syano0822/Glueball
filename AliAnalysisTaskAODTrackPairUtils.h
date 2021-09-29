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

  void setMC(bool isMC){fIsMC = isMC;}
  
  bool isAcceptEvent();
  bool isAcceptTrack(AliAODTrack* track, AliPID::EParticleType pid);
  bool isAcceptTrackQuality(AliAODTrack* track);
  bool isAcceptPid(AliAODTrack* track, AliPID::EParticleType pid);
  bool isAcceptTrackPair(AliAODTrack* track1, AliAODTrack* track2);
  bool isAcceptV0Pair(AliAODv0* v0_1,AliAODv0* v0_2);
  bool isAcceptK0sTrackPair(AliAODv0* v0);
  bool isK0sCandidate(float mass);

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
  
  void setMinV0Radius(float min)
  {
    fMinV0Radius = min;
  }
  void setMinV0PosTrackDCA(float min)
  {
    fMinPosDCA = min;
  }
  void setMinV0NegTrackDCA(float min)
  {
    fMinNegDCA = min;
  }
  void setMinV0CosPoint(float min)
  {
    fMinCosPoint = min;
  }
  void setRejLambdaMassRange(float min, float max)
  {
    fMinRejLambdaMass=min;
    fMaxRejLambdaMass=max;
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

  int fRunNumber;
  
  bool fIsMC;
 
  bool fIsVtxZcut;
  float fMaxVertexCutZ;
  float fMinVertexCutZ;

  bool fIsPairRapCut;
  float fMinPairRapCut;
  float fMaxPairRapCut;
  
  bool fIsPUcut;  
  
  float fMinTrackPt;
  float fMinTrackEta;
  float fMaxTrackEta;

  float fMinMassK0sCand;
  float fMaxMassK0sCand;

  float fMinV0Radius;
  float fMinPosDCA;
  float fMinNegDCA;
  float fMinCosPoint;
  float fMaxRejLambdaMass;
  float fMinRejLambdaMass;

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

  ClassDef(AliAnalysisTaskAODTrackPairUtils, 1); // example of analysis
};

#endif
