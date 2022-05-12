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
  void setME(bool flag){fOnME = flag;}
  void setKaonAnalysis(bool flag){fOnKaonAna = flag;}
  void setPionAnalysis(bool flag){fOnPionAna = flag;}
  void setK0sAnalysis(bool flag){fOnK0sAna = flag;}
  void setTrackQA(bool flag){fOnTrackQA = flag;}
  void setMLSample(bool flag){fOnMLTrainSampleout = flag;}
  
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
  bool TrackQA(AliAODTrack *track);
  bool V0TrackPairAnalysis();  
  bool primeTrackPairAnalysis();  
  bool processMC();
  bool RegisterTreeForTrackPID();
  
  AliAODEvent    *fEvent;  
  
  AliEventPoolManager *fPoolKaonTrackMgr;
  AliEventPoolManager *fPoolPionTrackMgr;
  AliEventPoolManager *fPoolK0sMgr;
  
  AliAnalysisTaskAODTrackPairUtils* fUtils;
  TClonesArray   *fMCTrackArray;
  
  bool fIsMC;  
  bool fOnME;
  bool fOnKaonAna;
  bool fOnPionAna;
  bool fOnK0sAna;
  bool fOnTrackQA;
  bool fOnMLTrainSampleout;
  bool fOnTrackMLTrainSampleout;

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
  
  TTree *fTreeML;
  float fTrackPt;
  float fTrackP;
  float fTrackTheta;
  float fTrackPhi;
  float fTrackLength;
  float fTrackBeta;
  float fTrackTrackChi2perNDF;
  float fTrackTrackITSNcls;
  float fTrackTrackTPCNcls;
  float fTrackTrackTOFNcls;
  float fTrackTrackTPCChi2;
  float fTrackTrackITSChi2;
  float fTrackTPCCrossedRows;
  float fTrackTPCFindableNcls;
  float fTrackTOFBCTime;
  float fTrackTOFKinkIndex;
  float fTrackDCAxy;
  float fTrackDCAz;
  float fTrackTPCsigmaK;
  float fTrackTOFsigmaK;
  bool fTrackTrueKaonLabel;

  TH1F* fEventCounter;

  TH2F* fHistArmenterosPodolanskiPlot; 
  TH2F* fHistArmenterosPodolanskiPlot_K0s; 
  
  TH2F* fHistULSKKPairMassPt;
  TH2F* fHistLSppKKPairMassPt;
  TH2F* fHistLSmmKKPairMassPt;

  TH2F* fHistMixULSKKPairMassPt;
  TH2F* fHistMixLSppKKPairMassPt;
  TH2F* fHistMixLSmmKKPairMassPt;

  TH2F* fHistULSPiPiPairMassPt;
  TH2F* fHistLSppPiPiPairMassPt;
  TH2F* fHistLSmmPiPiPairMassPt;

  TH2F* fHistMixULSPiPiPairMassPt;
  TH2F* fHistMixLSppPiPiPairMassPt;
  TH2F* fHistMixLSmmPiPiPairMassPt;

  TH2F* fHistULSKKPairMassPt_PionID;
  TH2F* fHistLSppKKPairMassPt_PionID;
  TH2F* fHistLSmmKKPairMassPt_PionID;

  TH2F* fHistMixULSKKPairMassPt_PionID;
  TH2F* fHistMixLSppKKPairMassPt_PionID;
  TH2F* fHistMixLSmmKKPairMassPt_PionID;
  
  TH2F* fHistK0sMassPt;
  TH2F* fHistK0sK0sMassPt;
  TH2F* fHistMixK0sK0sMassPt;
  
  TH2F* fHistPrimaryTrackNSPDCluster;
  TH2F* fHistPrimaryTrackNITSCluster;
  TH2F* fHistPrimaryTrackNTPCCluster;
  TH2F* fHistPrimaryTrackEta;
  TH2F* fHistPrimaryTrackChi2perNDF;
  TH2F* fHistPrimaryTrackITSChi2perNDF;
  TH2F* fHistPrimaryTrackTPCChi2perNDF;
  TH2F* fHistPrimaryTrackTPCCrossRaws;
  TH2F* fHistPrimaryTrackTPCFindableCrossRaws;
  TH2F* fHistPrimaryTrackTPCCrossRawsOverFindable;
  TH2F* fHistPrimaryTrackTOFBunchCrossing;
  TH2F* fHistPrimaryTrackKinkIndex;
  TH2F* fHistPrimaryTrackDCAxy;
  TH2F* fHistPrimaryTrackDCAz;

  TH2F* fHistPrimaryTrackNSPDClusterQuality;
  TH2F* fHistPrimaryTrackNITSClusterQuality;
  TH2F* fHistPrimaryTrackNTPCClusterQuality;
  TH2F* fHistPrimaryTrackEtaQuality;
  TH2F* fHistPrimaryTrackChi2perNDFQuality;
  TH2F* fHistPrimaryTrackITSChi2perNDFQuality;
  TH2F* fHistPrimaryTrackTPCChi2perNDFQuality;
  TH2F* fHistPrimaryTrackTPCCrossRawsQuality;
  TH2F* fHistPrimaryTrackTPCFindableCrossRawsQuality;
  TH2F* fHistPrimaryTrackTPCCrossRawsOverFindableQuality;
  TH2F* fHistPrimaryTrackTOFBunchCrossingQuality;
  TH2F* fHistPrimaryTrackKinkIndexQuality;
  TH2F* fHistPrimaryTrackDCAxyQuality;
  TH2F* fHistPrimaryTrackDCAzQuality;

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

  TH2F* fHistULSKKPairMassPt_FromSameSource;
  TH2F* fHistLSppKKPairMassPt_FromSameSource;
  TH2F* fHistLSmmKKPairMassPt_FromSameSource;
  
  TH2F* fHistULSKKPairMassPt_FromPhi;
  TH2F* fHistULSKKPairMassPt_FromK0s;
  TH2F* fHistULSKKPairMassPt_FromEta;
  TH2F* fHistULSKKPairMassPt_FromEtaPrime;
  TH2F* fHistULSKKPairMassPt_FromRho;
  TH2F* fHistULSKKPairMassPt_FromOmega;
  TH2F* fHistULSKKPairMassPt_FromK0star;
  TH2F* fHistULSKKPairMassPt_FromF980;
  TH2F* fHistULSKKPairMassPt_FromF1270;
  TH2F* fHistULSKKPairMassPt_FromF1370;
  TH2F* fHistULSKKPairMassPt_FromF1500;
  TH2F* fHistULSKKPairMassPt_FromF1525;
  TH2F* fHistULSKKPairMassPt_FromF1710;
  TH2F* fHistULSKKPairMassPt_FromPhi_MisID;
  TH2F* fHistULSKKPairMassPt_FromK0s_MisID;
  TH2F* fHistULSKKPairMassPt_FromEta_MisID;
  TH2F* fHistULSKKPairMassPt_FromEtaPrime_MisID;
  TH2F* fHistULSKKPairMassPt_FromRho_MisID;
  TH2F* fHistULSKKPairMassPt_FromOmega_MisID;
  TH2F* fHistULSKKPairMassPt_FromK0star_MisID;
  TH2F* fHistULSKKPairMassPt_FromF980_MisID;
  TH2F* fHistULSKKPairMassPt_FromF1270_MisID;
  TH2F* fHistULSKKPairMassPt_FromF1370_MisID;
  TH2F* fHistULSKKPairMassPt_FromF1500_MisID;
  TH2F* fHistULSKKPairMassPt_FromF1525_MisID;
  TH2F* fHistULSKKPairMassPt_FromF1710_MisID;
     
  TH2F* fHistULSKKPairMassPt_FromLambda1520;
  TH2F* fHistULSKKPairMassPt_FromK01430;
  TH2F* fHistULSKKPairMassPt_FromKp1430;
  TH2F* fHistULSKKPairMassPt_FromXi0star1530;
  TH2F* fHistULSKKPairMassPt_FromXi01820;
  TH2F* fHistULSKKPairMassPt_FromXip1820;
  TH2F* fHistULSKKPairMassPt_FromSigmastarp1385;
  TH2F* fHistULSKKPairMassPt_FromSigmastarm1385;
  TH2F* fHistULSKKPairMassPt_FromLambda1520_MisID;
  TH2F* fHistULSKKPairMassPt_FromK01430_MisID;
  TH2F* fHistULSKKPairMassPt_FromKp1430_MisID;
  TH2F* fHistULSKKPairMassPt_FromXi0star1530_MisID;
  TH2F* fHistULSKKPairMassPt_FromXi01820_MisID;
  TH2F* fHistULSKKPairMassPt_FromXip1820_MisID;
  TH2F* fHistULSKKPairMassPt_FromSigmastarp1385_MisID;
  TH2F* fHistULSKKPairMassPt_FromSigmastarm1385_MisID;

  TH2F* fHistULSPiPiPairMassPt_FromSameSource;
  TH2F* fHistLSppPiPiPairMassPt_FromSameSource;
  TH2F* fHistLSmmPiPiPairMassPt_FromSameSource;

  TH2F* fHistULSPiPiPairMassPt_FromPhi;
  TH2F* fHistULSPiPiPairMassPt_FromK0s;
  TH2F* fHistULSPiPiPairMassPt_FromEta;
  TH2F* fHistULSPiPiPairMassPt_FromEtaPrime;
  TH2F* fHistULSPiPiPairMassPt_FromRho;
  TH2F* fHistULSPiPiPairMassPt_FromOmega;
  TH2F* fHistULSPiPiPairMassPt_FromK0star;
  TH2F* fHistULSPiPiPairMassPt_FromF980;
  TH2F* fHistULSPiPiPairMassPt_FromF1270;
  TH2F* fHistULSPiPiPairMassPt_FromF1370;
  TH2F* fHistULSPiPiPairMassPt_FromF1500;
  TH2F* fHistULSPiPiPairMassPt_FromF1525;
  TH2F* fHistULSPiPiPairMassPt_FromF1710;
  TH2F* fHistULSPiPiPairMassPt_FromPhi_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromK0s_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromEta_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromEtaPrime_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromRho_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromOmega_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromK0star_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF980_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF1270_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF1370_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF1500_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF1525_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromF1710_MisID;
  
  TH2F* fHistULSPiPiPairMassPt_FromLambda1520;
  TH2F* fHistULSPiPiPairMassPt_FromK01430;
  TH2F* fHistULSPiPiPairMassPt_FromKp1430;
  TH2F* fHistULSPiPiPairMassPt_FromXi0star1530;
  TH2F* fHistULSPiPiPairMassPt_FromXi01820;
  TH2F* fHistULSPiPiPairMassPt_FromXip1820;
  TH2F* fHistULSPiPiPairMassPt_FromSigmastarp1385;
  TH2F* fHistULSPiPiPairMassPt_FromSigmastarm1385;
  TH2F* fHistULSPiPiPairMassPt_FromLambda1520_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromK01430_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromKp1430_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromXi0star1530_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromXi01820_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromXip1820_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID;
  TH2F* fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID;


  TH1F* fHistRhoMassTrue;
  TH1F* fHistPhiMassTrue;
  TH1F* fHistF980MassTrue;
  TH1F* fHistF1270MassTrue;
  TH1F* fHistK0sMassTrue;

  TH1F* fHistPhiPtTrue;
  TH1F* fHistEtaPtTrue;
  TH1F* fHistEtaPrimePtTrue;
  TH1F* fHistRhoPtTrue;
  TH1F* fHistOmegaPtTrue;
  TH1F* fHistK0sPtTrue;
  TH1F* fHistK0starPtTrue;
  TH1F* fHistF980PtTrue;
  TH1F* fHistF1270PtTrue;
  TH1F* fHistF1500PtTrue;
  TH1F* fHistF1525PtTrue;
  TH1F* fHistF1710PtTrue;

  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
