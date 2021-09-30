#ifndef AliAnalysisTaskAODTrackPair_cxx
#define AliAnalysisTaskAODTrackPair_cxx

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliAnalysisTaskAODTrackPairUtils.h"

class AliAnalysisTaskAODTrackPair : public AliAnalysisTaskSE {
  
 public:
  
  AliAnalysisTaskAODTrackPair();
  AliAnalysisTaskAODTrackPair(const char* name);
  virtual ~AliAnalysisTaskAODTrackPair();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
  void setMC(bool isMC){fIsMC = isMC;}
  
  void setUtils(AliAnalysisTaskAODTrackPairUtils *utils)
  {
    fUtils = utils;
  }
  void setEvtMixingTrackDepth(int depth)
  {
    fTrackDepth=depth;
  }
  void setEvtMixingPoolSize(int size)
  {
    fPoolSize=size;
  }
  void setEvtMixingReadyFraction(float frac)
  {
    fReadyFraction=frac;
  }
  void setEvtMixingPoolVtxZ(bool flag)
  {
    onEvtMixingPoolVtxZ=flag;
  }
  void setEvtMixingPoolCent(bool flag)
  {
    onEvtMixingPoolCent=flag;
  }
  void setEvtMixingPoolPsi(bool flag)
  {
    onEvtMixingPoolPsi=flag;
  }
  void setMixingEventTrigger(unsigned int mask)
  {
    fTriggerMaskForMixing = mask;
  }
  void setSameEventTrigger(unsigned int mask)
  {
    fTriggerMaskForSame = mask;
  }

 private:
  
  AliAnalysisTaskAODTrackPair(const AliAnalysisTaskAODTrackPair&); // not implemented
  AliAnalysisTaskAODTrackPair& operator=(const AliAnalysisTaskAODTrackPair&); // not implemented
  
  bool Initialize();
  bool TrackQA();
  bool V0TrackPairAnalysis();  
  bool TrackPairAnalysisEveMixing();  
  
  AliEventPool* setPool();
  
  AliAODEvent    *fEvent;  
  AliEventPoolManager *fPoolTrackMgr;
  AliAnalysisTaskAODTrackPairUtils* fUtils;
  TClonesArray   *fMCTrackArray;
  
  bool fIsMC;  
  int fRunNumber;
  int fTrackDepth;
  int fPoolSize;
  
  unsigned int fTriggerMaskForSame;
  unsigned int fTriggerMaskForMixing;
  
  float fReadyFraction;
  
  bool onEvtMixingPoolVtxZ;
  bool onEvtMixingPoolCent;
  bool onEvtMixingPoolPsi;
  
  ////////////////////////////////////////////////
  //Output histos
  ////////////////////////////////////////////////

  TList* fOutputList;
  
  TH1F* fEventCounter;
  
  THnSparse* fSparseK0s;
  THnSparse* fSparseK0sK0s;
  
  TH2F* fHistMomTPCSignalTrackQuality;
  TH2F* fHistMomTOFBetaTrackQuality;
  TH2F* fHistMomTPCSignalTrackQualityPionPID;
  TH2F* fHistMomTOFBetaTrackQualityPionPID;
  TH2F* fHistMomTPCSignalTrackQualityKaonPID;
  TH2F* fHistMomTOFBetaTrackQualityKaonPID;
  TH2F* fHistMomTPCSignalTrackQualityProtonPID;
  TH2F* fHistMomTOFBetaTrackQualityProtonPID;

  TH2F* fHistMomTPCSigmaTrackQuality;
  TH2F* fHistMomTOFSigmaTrackQuality;
  TH2F* fHistMomTPCSigmaTrackQualityPionPID;
  TH2F* fHistMomTOFSigmaTrackQualityPionPID;
  TH2F* fHistMomTPCSigmaTrackQualityKaonPID;
  TH2F* fHistMomTOFSigmaTrackQualityKaonPID;
  TH2F* fHistMomTPCSigmaTrackQualityProtonPID;
  TH2F* fHistMomTOFSigmaTrackQualityProtonPID;

  
  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
