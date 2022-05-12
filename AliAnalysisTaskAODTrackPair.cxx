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
  fOnME(false),
  fOnKaonAna(false),
  fOnPionAna(false),
  fOnK0sAna(false),
  fOnTrackQA(false),
  fOnTrackMLTrainSampleout(false),
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
  
  fTreeML(NULL),
  fTrackPt(0.),
  fTrackP(0.),
  fTrackTheta(0.),
  fTrackPhi(0.),
  fTrackLength(0.),
  fTrackBeta(0.),
  fTrackTrackChi2perNDF(0.),
  fTrackTrackITSNcls(0.),
  fTrackTrackTPCNcls(0.),
  fTrackTrackTOFNcls(0.),
  fTrackTrackTPCChi2(0.),
  fTrackTrackITSChi2(0.),
  fTrackTPCCrossedRows(0.),
  fTrackTPCFindableNcls(0.),
  fTrackTOFBCTime(0.),
  fTrackTOFKinkIndex(0.),
  fTrackDCAxy(0.),
  fTrackDCAz(0.),
  fTrackTPCsigmaK(0.),
  fTrackTOFsigmaK(0.),
  fTrackTrueKaonLabel(false),

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
  fHistMomTOFSigmaTrackQualityProtonPID(NULL),

  fHistULSKKPairMassPt_FromSameSource(NULL),
  fHistLSppKKPairMassPt_FromSameSource(NULL),
  fHistLSmmKKPairMassPt_FromSameSource(NULL),  
  fHistULSKKPairMassPt_FromPhi(NULL),
  fHistULSKKPairMassPt_FromK0s(NULL),
  fHistULSKKPairMassPt_FromEta(NULL),
  fHistULSKKPairMassPt_FromEtaPrime(NULL),
  fHistULSKKPairMassPt_FromRho(NULL),
  fHistULSKKPairMassPt_FromOmega(NULL),
  fHistULSKKPairMassPt_FromK0star(NULL),
  fHistULSKKPairMassPt_FromF980(NULL),
  fHistULSKKPairMassPt_FromF1270(NULL),
  fHistULSKKPairMassPt_FromF1370(NULL),
  fHistULSKKPairMassPt_FromF1500(NULL),
  fHistULSKKPairMassPt_FromF1525(NULL),
  fHistULSKKPairMassPt_FromF1710(NULL),
  fHistULSKKPairMassPt_FromPhi_MisID(NULL),
  fHistULSKKPairMassPt_FromK0s_MisID(NULL),
  fHistULSKKPairMassPt_FromEta_MisID(NULL),
  fHistULSKKPairMassPt_FromEtaPrime_MisID(NULL),
  fHistULSKKPairMassPt_FromRho_MisID(NULL),
  fHistULSKKPairMassPt_FromOmega_MisID(NULL),
  fHistULSKKPairMassPt_FromK0star_MisID(NULL),
  fHistULSKKPairMassPt_FromF980_MisID(NULL),
  fHistULSKKPairMassPt_FromF1270_MisID(NULL),
  fHistULSKKPairMassPt_FromF1370_MisID(NULL),
  fHistULSKKPairMassPt_FromF1500_MisID(NULL),
  fHistULSKKPairMassPt_FromF1525_MisID(NULL),
  fHistULSKKPairMassPt_FromF1710_MisID(NULL),
  fHistULSKKPairMassPt_FromLambda1520(NULL),
  fHistULSKKPairMassPt_FromK01430(NULL),
  fHistULSKKPairMassPt_FromKp1430(NULL),
  fHistULSKKPairMassPt_FromXi0star1530(NULL),
  fHistULSKKPairMassPt_FromXi01820(NULL),
  fHistULSKKPairMassPt_FromXip1820(NULL),
  fHistULSKKPairMassPt_FromSigmastarp1385(NULL),
  fHistULSKKPairMassPt_FromSigmastarm1385(NULL),
  fHistULSKKPairMassPt_FromLambda1520_MisID(NULL),
  fHistULSKKPairMassPt_FromK01430_MisID(NULL),
  fHistULSKKPairMassPt_FromKp1430_MisID(NULL),
  fHistULSKKPairMassPt_FromXi0star1530_MisID(NULL),
  fHistULSKKPairMassPt_FromXi01820_MisID(NULL),
  fHistULSKKPairMassPt_FromXip1820_MisID(NULL),
  fHistULSKKPairMassPt_FromSigmastarp1385_MisID(NULL),
  fHistULSKKPairMassPt_FromSigmastarm1385_MisID(NULL),

  fHistULSPiPiPairMassPt_FromSameSource(NULL),
  fHistLSppPiPiPairMassPt_FromSameSource(NULL),
  fHistLSmmPiPiPairMassPt_FromSameSource(NULL),
  fHistULSPiPiPairMassPt_FromPhi(NULL),
  fHistULSPiPiPairMassPt_FromK0s(NULL),
  fHistULSPiPiPairMassPt_FromEta(NULL),
  fHistULSPiPiPairMassPt_FromEtaPrime(NULL),
  fHistULSPiPiPairMassPt_FromRho(NULL),
  fHistULSPiPiPairMassPt_FromOmega(NULL),
  fHistULSPiPiPairMassPt_FromK0star(NULL),
  fHistULSPiPiPairMassPt_FromF980(NULL),
  fHistULSPiPiPairMassPt_FromF1270(NULL),
  fHistULSPiPiPairMassPt_FromF1370(NULL),
  fHistULSPiPiPairMassPt_FromF1500(NULL),
  fHistULSPiPiPairMassPt_FromF1525(NULL),
  fHistULSPiPiPairMassPt_FromF1710(NULL),
  fHistULSPiPiPairMassPt_FromPhi_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK0s_MisID(NULL),
  fHistULSPiPiPairMassPt_FromEta_MisID(NULL),
  fHistULSPiPiPairMassPt_FromEtaPrime_MisID(NULL),
  fHistULSPiPiPairMassPt_FromRho_MisID(NULL),
  fHistULSPiPiPairMassPt_FromOmega_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK0star_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF980_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1270_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1370_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1500_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1525_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1710_MisID(NULL),
  fHistULSPiPiPairMassPt_FromLambda1520(NULL),
  fHistULSPiPiPairMassPt_FromK01430(NULL),
  fHistULSPiPiPairMassPt_FromKp1430(NULL),
  fHistULSPiPiPairMassPt_FromXi0star1530(NULL),
  fHistULSPiPiPairMassPt_FromXi01820(NULL),
  fHistULSPiPiPairMassPt_FromXip1820(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarp1385(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarm1385(NULL),
  fHistULSPiPiPairMassPt_FromLambda1520_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK01430_MisID(NULL),
  fHistULSPiPiPairMassPt_FromKp1430_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXi0star1530_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXi01820_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXip1820_MisID(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID(NULL),

  fHistRhoMassTrue(NULL),
  fHistPhiMassTrue(NULL),
  fHistF980MassTrue(NULL),
  fHistF1270MassTrue(NULL),
  fHistK0sMassTrue(NULL),

  fHistPhiPtTrue(NULL),
  fHistEtaPtTrue(NULL),
  fHistEtaPrimePtTrue(NULL),
  fHistRhoPtTrue(NULL),
  fHistOmegaPtTrue(NULL),
  fHistK0sPtTrue(NULL),
  fHistK0starPtTrue(NULL),
  fHistF980PtTrue(NULL),
  fHistF1270PtTrue(NULL),
  fHistF1500PtTrue(NULL),
  fHistF1525PtTrue(NULL),
  fHistF1710PtTrue(NULL)


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
  fOnME(false),
  fOnKaonAna(false),
  fOnPionAna(false),
  fOnK0sAna(false),
  fOnTrackQA(false),
  fOnTrackMLTrainSampleout(false),
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

  fTreeML(NULL),
  fTrackPt(0.),
  fTrackP(0.),
  fTrackTheta(0.),
  fTrackPhi(0.),
  fTrackLength(0.),
  fTrackBeta(0.),
  fTrackTrackChi2perNDF(0.),
  fTrackTrackITSNcls(0.),
  fTrackTrackTPCNcls(0.),
  fTrackTrackTOFNcls(0.),
  fTrackTrackTPCChi2(0.),
  fTrackTrackITSChi2(0.),
  fTrackTPCCrossedRows(0.),
  fTrackTPCFindableNcls(0.),
  fTrackTOFBCTime(0.),
  fTrackTOFKinkIndex(0.),
  fTrackDCAxy(0.),
  fTrackDCAz(0.),
  fTrackTPCsigmaK(0.),
  fTrackTOFsigmaK(0.),
  fTrackTrueKaonLabel(false),

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
  fHistMomTOFSigmaTrackQualityProtonPID(NULL),

  fHistULSKKPairMassPt_FromSameSource(NULL),
  fHistLSppKKPairMassPt_FromSameSource(NULL),
  fHistLSmmKKPairMassPt_FromSameSource(NULL),
  fHistULSKKPairMassPt_FromPhi(NULL),
  fHistULSKKPairMassPt_FromK0s(NULL),
  fHistULSKKPairMassPt_FromEta(NULL),
  fHistULSKKPairMassPt_FromEtaPrime(NULL),
  fHistULSKKPairMassPt_FromRho(NULL),
  fHistULSKKPairMassPt_FromOmega(NULL),
  fHistULSKKPairMassPt_FromK0star(NULL),
  fHistULSKKPairMassPt_FromF980(NULL),
  fHistULSKKPairMassPt_FromF1270(NULL),
  fHistULSKKPairMassPt_FromF1370(NULL),
  fHistULSKKPairMassPt_FromF1500(NULL),
  fHistULSKKPairMassPt_FromF1525(NULL),
  fHistULSKKPairMassPt_FromF1710(NULL),
  fHistULSKKPairMassPt_FromPhi_MisID(NULL),
  fHistULSKKPairMassPt_FromK0s_MisID(NULL),
  fHistULSKKPairMassPt_FromEta_MisID(NULL),
  fHistULSKKPairMassPt_FromEtaPrime_MisID(NULL),
  fHistULSKKPairMassPt_FromRho_MisID(NULL),
  fHistULSKKPairMassPt_FromOmega_MisID(NULL),
  fHistULSKKPairMassPt_FromK0star_MisID(NULL),
  fHistULSKKPairMassPt_FromF980_MisID(NULL),
  fHistULSKKPairMassPt_FromF1270_MisID(NULL),
  fHistULSKKPairMassPt_FromF1370_MisID(NULL),
  fHistULSKKPairMassPt_FromF1500_MisID(NULL),
  fHistULSKKPairMassPt_FromF1525_MisID(NULL),
  fHistULSKKPairMassPt_FromF1710_MisID(NULL),
  fHistULSKKPairMassPt_FromLambda1520(NULL),
  fHistULSKKPairMassPt_FromK01430(NULL),
  fHistULSKKPairMassPt_FromKp1430(NULL),
  fHistULSKKPairMassPt_FromXi0star1530(NULL),
  fHistULSKKPairMassPt_FromXi01820(NULL),
  fHistULSKKPairMassPt_FromXip1820(NULL),
  fHistULSKKPairMassPt_FromSigmastarp1385(NULL),
  fHistULSKKPairMassPt_FromSigmastarm1385(NULL),
  fHistULSKKPairMassPt_FromLambda1520_MisID(NULL),
  fHistULSKKPairMassPt_FromK01430_MisID(NULL),
  fHistULSKKPairMassPt_FromKp1430_MisID(NULL),
  fHistULSKKPairMassPt_FromXi0star1530_MisID(NULL),
  fHistULSKKPairMassPt_FromXi01820_MisID(NULL),
  fHistULSKKPairMassPt_FromXip1820_MisID(NULL),
  fHistULSKKPairMassPt_FromSigmastarp1385_MisID(NULL),
  fHistULSKKPairMassPt_FromSigmastarm1385_MisID(NULL),

  fHistULSPiPiPairMassPt_FromSameSource(NULL),
  fHistLSppPiPiPairMassPt_FromSameSource(NULL),
  fHistLSmmPiPiPairMassPt_FromSameSource(NULL),
  fHistULSPiPiPairMassPt_FromPhi(NULL),
  fHistULSPiPiPairMassPt_FromK0s(NULL),
  fHistULSPiPiPairMassPt_FromEta(NULL),
  fHistULSPiPiPairMassPt_FromEtaPrime(NULL),
  fHistULSPiPiPairMassPt_FromRho(NULL),
  fHistULSPiPiPairMassPt_FromOmega(NULL),
  fHistULSPiPiPairMassPt_FromK0star(NULL),
  fHistULSPiPiPairMassPt_FromF980(NULL),
  fHistULSPiPiPairMassPt_FromF1270(NULL),
  fHistULSPiPiPairMassPt_FromF1370(NULL),
  fHistULSPiPiPairMassPt_FromF1500(NULL),
  fHistULSPiPiPairMassPt_FromF1525(NULL),
  fHistULSPiPiPairMassPt_FromF1710(NULL),
  fHistULSPiPiPairMassPt_FromPhi_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK0s_MisID(NULL),
  fHistULSPiPiPairMassPt_FromEta_MisID(NULL),
  fHistULSPiPiPairMassPt_FromEtaPrime_MisID(NULL),
  fHistULSPiPiPairMassPt_FromRho_MisID(NULL),
  fHistULSPiPiPairMassPt_FromOmega_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK0star_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF980_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1270_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1370_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1500_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1525_MisID(NULL),
  fHistULSPiPiPairMassPt_FromF1710_MisID(NULL),
  fHistULSPiPiPairMassPt_FromLambda1520(NULL),
  fHistULSPiPiPairMassPt_FromK01430(NULL),
  fHistULSPiPiPairMassPt_FromKp1430(NULL),
  fHistULSPiPiPairMassPt_FromXi0star1530(NULL),
  fHistULSPiPiPairMassPt_FromXi01820(NULL),
  fHistULSPiPiPairMassPt_FromXip1820(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarp1385(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarm1385(NULL),
  fHistULSPiPiPairMassPt_FromLambda1520_MisID(NULL),
  fHistULSPiPiPairMassPt_FromK01430_MisID(NULL),
  fHistULSPiPiPairMassPt_FromKp1430_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXi0star1530_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXi01820_MisID(NULL),
  fHistULSPiPiPairMassPt_FromXip1820_MisID(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID(NULL),
  fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID(NULL),

  fHistRhoMassTrue(NULL),
  fHistPhiMassTrue(NULL),
  fHistF980MassTrue(NULL),
  fHistF1270MassTrue(NULL),
  fHistK0sMassTrue(NULL),

  fHistPhiPtTrue(NULL),
  fHistEtaPtTrue(NULL),
  fHistEtaPrimePtTrue(NULL),
  fHistRhoPtTrue(NULL),
  fHistOmegaPtTrue(NULL),
  fHistK0sPtTrue(NULL),
  fHistK0starPtTrue(NULL),
  fHistF980PtTrue(NULL),
  fHistF1270PtTrue(NULL),
  fHistF1500PtTrue(NULL),
  fHistF1525PtTrue(NULL),
  fHistF1710PtTrue(NULL)
  
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
  

  if(fOnMLTrainSampleout){
    fTreeML = new TTree("fTreeML","Tree for machine leraning for PID");
    fTreeML->Branch("fTrackPt", &fTrackPt, "fTrackPt/F");
    fTreeML->Branch("fTrackP", &fTrackP, "fTrackP/F");
    fTreeML->Branch("fTrackTheta", &fTrackTheta, "fTrackTheta/F");
    fTreeML->Branch("fTrackPhi", &fTrackPhi, "fTrackPhi/F");
    fTreeML->Branch("fTrackLength", &fTrackLength, "fTrackLength/F");
    fTreeML->Branch("fTrackBeta", &fTrackBeta, "fTrackBeta/F");
    fTreeML->Branch("fTrackTrackChi2perNDF", &fTrackTrackChi2perNDF, "fTrackTrackChi2perNDF/F");
    fTreeML->Branch("fTrackTrackITSNcls", &fTrackTrackITSNcls, "fTrackTrackITSNcls/F");
    fTreeML->Branch("fTrackTrackTPCNcls", &fTrackTrackTPCNcls, "fTrackTrackTPCNcls/F");
    fTreeML->Branch("fTrackTrackTOFNcls", &fTrackTrackTOFNcls, "fTrackTrackTOFNcls/F");
    fTreeML->Branch("fTrackTrackTPCChi2", &fTrackTrackTPCChi2, "fTrackTrackTPCChi2/F");
    fTreeML->Branch("fTrackTrackITSChi2", &fTrackTrackITSChi2, "fTrackTrackITSChi2/F");
    fTreeML->Branch("fTrackTPCCrossedRows", &fTrackTPCCrossedRows, "fTrackTPCCrossedRows/F");
    fTreeML->Branch("fTrackTPCFindableNcls", &fTrackTPCFindableNcls, "fTrackTPCFindableNcls/F");
    fTreeML->Branch("fTrackTOFBCTime", &fTrackTOFBCTime, "fTrackTOFBCTime/F");
    fTreeML->Branch("fTrackTOFKinkIndex", &fTrackTOFKinkIndex, "fTrackTOFKinkIndex/F");
    fTreeML->Branch("fTrackDCAxy", &fTrackDCAxy, "fTrackDCAxy/F");
    fTreeML->Branch("fTrackDCAz", &fTrackDCAz, "fTrackDCAz/F");
    fTreeML->Branch("fTrackTPCsigmaK", &fTrackTPCsigmaK, "fTrackTPCsigmaK/F");
    fTreeML->Branch("fTrackTOFsigmaK", &fTrackTOFsigmaK, "fTrackTOFsigmaK/F");
    fTreeML->Branch("fTrackTrueKaonLabel", &fTrackTrueKaonLabel, "fTrackTrueKaonLabel/B");
    fOutputList->Add(fTreeML);
  }

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
    if(fOnME){
      fHistMixULSKKPairMassPt = new TH2F("fHistMixULSKKPairMassPt","",700,0.5,4.,100,0,10);
      fHistMixLSppKKPairMassPt = new TH2F("fHistMixLSppKKPairMassPt","",700,0.5,4.,100,0,10);
      fHistMixLSmmKKPairMassPt = new TH2F("fHistMixLSmmKKPairMassPt","",700,0.5,4.,100,0,10);
      fOutputList->Add(fHistMixULSKKPairMassPt);
      fOutputList->Add(fHistMixLSppKKPairMassPt);
      fOutputList->Add(fHistMixLSmmKKPairMassPt);
    }
  }

  if(fOnPionAna){
    fHistULSPiPiPairMassPt = new TH2F("fHistULSPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistLSppPiPiPairMassPt = new TH2F("fHistLSppPiPiPairMassPt","",800,0.,4.,100,0,10);
    fHistLSmmPiPiPairMassPt = new TH2F("fHistLSmmPiPiPairMassPt","",800,0.,4.,100,0,10);
    fOutputList->Add(fHistULSPiPiPairMassPt);
    fOutputList->Add(fHistLSppPiPiPairMassPt);
    fOutputList->Add(fHistLSmmPiPiPairMassPt);
    if(fOnME){
      fHistMixULSPiPiPairMassPt = new TH2F("fHistMixULSPiPiPairMassPt","",800,0.,4.,100,0,10);
      fHistMixLSppPiPiPairMassPt = new TH2F("fHistMixLSppPiPiPairMassPt","",800,0.,4.,100,0,10);
      fHistMixLSmmPiPiPairMassPt = new TH2F("fHistMixLSmmPiPiPairMassPt","",800,0.,4.,100,0,10);
      fOutputList->Add(fHistMixULSPiPiPairMassPt);
      fOutputList->Add(fHistMixLSppPiPiPairMassPt);
      fOutputList->Add(fHistMixLSmmPiPiPairMassPt);
    }
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

  if(fIsMC){    

    double addbin = 5;
    
    fHistRhoMassTrue = new TH1F("fHistRhoMassTrue","",700*addbin,0.,3.5);
    fHistPhiMassTrue = new TH1F("fHistPhiMassTrue","",700*addbin,0.,3.5);
    fHistF980MassTrue = new TH1F("fHistF980MassTrue","",700*addbin,0.,3.5);
    fHistF1270MassTrue = new TH1F("fHistF1270MassTrue","",700*addbin,0.0,3.5);
    fHistK0sMassTrue = new TH1F("fHistK0sMassTrue","",700*addbin,0.,3.5);

    fOutputList->Add(fHistRhoMassTrue);
    fOutputList->Add(fHistPhiMassTrue);
    fOutputList->Add(fHistF980MassTrue);
    fOutputList->Add(fHistF1270MassTrue);
    fOutputList->Add(fHistK0sMassTrue);

    fHistPhiPtTrue = new TH1F("fHistPhiPtTrue","",300,0,30);
    fHistEtaPtTrue = new TH1F("fHistEtaPtTrue","",300,0,30);
    fHistEtaPrimePtTrue = new TH1F("fHistEtaPrimePtTrue","",300,0,30);
    fHistRhoPtTrue = new TH1F("fHistRhoPtTrue","",300,0,30);
    fHistOmegaPtTrue = new TH1F("fHistOmegaPtTrue","",300,0,30);
    fHistK0sPtTrue = new TH1F("fHistK0sPtTrue","",300,0,30);
    fHistK0starPtTrue = new TH1F("fHistK0starPtTrue","",300,0,30);
    fHistF980PtTrue = new TH1F("fHistF980PtTrue","",300,0,30);
    fHistF1270PtTrue = new TH1F("fHistF1270PtTrue","",300,0,30);
    fHistF1500PtTrue = new TH1F("fHistF1500PtTrue","",300,0,30);
    fHistF1525PtTrue = new TH1F("fHistF1525PtTrue","",300,0,30);
    fHistF1710PtTrue = new TH1F("fHistF1710PtTrue","",300,0,30);
	
    fOutputList->Add(fHistPhiPtTrue);
    fOutputList->Add(fHistEtaPtTrue);
    fOutputList->Add(fHistEtaPrimePtTrue);
    fOutputList->Add(fHistRhoPtTrue);
    fOutputList->Add(fHistOmegaPtTrue);
    fOutputList->Add(fHistK0sPtTrue);
    fOutputList->Add(fHistK0starPtTrue);
    fOutputList->Add(fHistF980PtTrue);
    fOutputList->Add(fHistF1270PtTrue);
    fOutputList->Add(fHistF1500PtTrue);
    fOutputList->Add(fHistF1525PtTrue);
    fOutputList->Add(fHistF1710PtTrue);

    if(fOnKaonAna){      
      fHistULSKKPairMassPt_FromSameSource = new TH2F("fHistULSKKPairMassPt_FromSameSource","",700*addbin,0.5,4.,100,0,10);
      fHistLSppKKPairMassPt_FromSameSource = new TH2F("fHistLSppKKPairMassPt_FromSameSource","",700*addbin,0.5,4.,100,0,10);
      fHistLSmmKKPairMassPt_FromSameSource = new TH2F("fHistLSmmKKPairMassPt_FromSameSource","",700*addbin,0.5,4.,100,0,10);
      fOutputList->Add(fHistULSKKPairMassPt_FromSameSource);
      fOutputList->Add(fHistLSppKKPairMassPt_FromSameSource);
      fOutputList->Add(fHistLSmmKKPairMassPt_FromSameSource);
      
      fHistULSKKPairMassPt_FromPhi = new TH2F("fHistULSKKPairMassPt_FromPhi","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK0s = new TH2F("fHistULSKKPairMassPt_FromK0s","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromEta = new TH2F("fHistULSKKPairMassPt_FromEta","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromEtaPrime = new TH2F("fHistULSKKPairMassPt_FromEtaPrime","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromRho = new TH2F("fHistULSKKPairMassPt_FromRho","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromOmega = new TH2F("fHistULSKKPairMassPt_FromOmega","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK0star = new TH2F("fHistULSKKPairMassPt_FromK0star","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF980 = new TH2F("fHistULSKKPairMassPt_FromF980","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1270 = new TH2F("fHistULSKKPairMassPt_FromF1270","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1370 = new TH2F("fHistULSKKPairMassPt_FromF1370","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1500 = new TH2F("fHistULSKKPairMassPt_FromF1500","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1525 = new TH2F("fHistULSKKPairMassPt_FromF1525","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1710 = new TH2F("fHistULSKKPairMassPt_FromF1710","",700*addbin,0.5,4.,100,0,10);      
      fHistULSKKPairMassPt_FromLambda1520 = new TH2F("fHistULSKKPairMassPt_FromLambda1520","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK01430 = new TH2F("fHistULSKKPairMassPt_FromK01430","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromKp1430 = new TH2F("fHistULSKKPairMassPt_FromKp1430","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXi0star1530 = new TH2F("fHistULSKKPairMassPt_FromXi0star1530","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXi01820 = new TH2F("fHistULSKKPairMassPt_FromXi1820","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXip1820 = new TH2F("fHistULSKKPairMassPt_FromXip1820","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromSigmastarp1385 = new TH2F("fHistULSKKPairMassPt_FromSigmastarp1385","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromSigmastarm1385 = new TH2F("fHistULSKKPairMassPt_FromSigmastarm1385","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromPhi_MisID = new TH2F("fHistULSKKPairMassPt_FromPhi_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK0s_MisID = new TH2F("fHistULSKKPairMassPt_FromK0s_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromEta_MisID = new TH2F("fHistULSKKPairMassPt_FromEta_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromEtaPrime_MisID = new TH2F("fHistULSKKPairMassPt_FromEtaPrime_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromRho_MisID = new TH2F("fHistULSKKPairMassPt_FromRho_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromOmega_MisID = new TH2F("fHistULSKKPairMassPt_FromOmega_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK0star_MisID = new TH2F("fHistULSKKPairMassPt_FromK0star_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF980_MisID = new TH2F("fHistULSKKPairMassPt_FromF980_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1270_MisID = new TH2F("fHistULSKKPairMassPt_FromF1270_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1370_MisID = new TH2F("fHistULSKKPairMassPt_FromF1370_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1500_MisID = new TH2F("fHistULSKKPairMassPt_FromF1500_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1525_MisID = new TH2F("fHistULSKKPairMassPt_FromF1525_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromF1710_MisID = new TH2F("fHistULSKKPairMassPt_FromF1710_MisID","",700*addbin,0.5,4.,100,0,10);      
      fHistULSKKPairMassPt_FromLambda1520_MisID = new TH2F("fHistULSKKPairMassPt_FromLambda1520_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromK01430_MisID = new TH2F("fHistULSKKPairMassPt_FromK01430_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromKp1430_MisID = new TH2F("fHistULSKKPairMassPt_FromKp1430_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXi0star1530_MisID = new TH2F("fHistULSKKPairMassPt_FromXi0star1530_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXi01820_MisID = new TH2F("fHistULSKKPairMassPt_FromXi1820_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromXip1820_MisID = new TH2F("fHistULSKKPairMassPt_FromXip1820_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromSigmastarp1385_MisID = new TH2F("fHistULSKKPairMassPt_FromSigmastarp1385_MisID","",700*addbin,0.5,4.,100,0,10);
      fHistULSKKPairMassPt_FromSigmastarm1385_MisID = new TH2F("fHistULSKKPairMassPt_FromSigmastarm1385_MisID","",700*addbin,0.5,4.,100,0,10);      
      fOutputList->Add(fHistULSKKPairMassPt_FromPhi);
      fOutputList->Add(fHistULSKKPairMassPt_FromK0s);
      fOutputList->Add(fHistULSKKPairMassPt_FromEta);
      fOutputList->Add(fHistULSKKPairMassPt_FromEtaPrime);
      fOutputList->Add(fHistULSKKPairMassPt_FromRho);
      fOutputList->Add(fHistULSKKPairMassPt_FromOmega);
      fOutputList->Add(fHistULSKKPairMassPt_FromK0star);
      fOutputList->Add(fHistULSKKPairMassPt_FromF980);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1270);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1370);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1500);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1525);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1710);  
      fOutputList->Add(fHistULSKKPairMassPt_FromLambda1520);
      fOutputList->Add(fHistULSKKPairMassPt_FromK01430);
      fOutputList->Add(fHistULSKKPairMassPt_FromKp1430);
      fOutputList->Add(fHistULSKKPairMassPt_FromXi0star1530);
      fOutputList->Add(fHistULSKKPairMassPt_FromXi01820);
      fOutputList->Add(fHistULSKKPairMassPt_FromXip1820);
      fOutputList->Add(fHistULSKKPairMassPt_FromSigmastarp1385);
      fOutputList->Add(fHistULSKKPairMassPt_FromSigmastarm1385);
      fOutputList->Add(fHistULSKKPairMassPt_FromPhi_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromK0s_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromEta_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromEtaPrime_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromRho_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromOmega_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromK0star_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF980_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1270_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1370_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1500_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1525_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromF1710_MisID);  
      fOutputList->Add(fHistULSKKPairMassPt_FromLambda1520_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromK01430_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromKp1430_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromXi0star1530_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromXi01820_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromXip1820_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromSigmastarp1385_MisID);
      fOutputList->Add(fHistULSKKPairMassPt_FromSigmastarm1385_MisID);

      
    }
    if(fOnPionAna){
      fHistULSPiPiPairMassPt_FromSameSource = new TH2F("fHistULSPiPiPairMassPt_FromSameSource","",800*addbin,0.,4.,100,0,10);
      fHistLSppPiPiPairMassPt_FromSameSource = new TH2F("fHistLSppPiPiPairMassPt_FromSameSource","",800*addbin,0.,4.,100,0,10);
      fHistLSmmPiPiPairMassPt_FromSameSource = new TH2F("fHistLSmmPiPiPairMassPt_FromSameSource","",800*addbin,0.,4.,100,0,10);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromSameSource);
      fOutputList->Add(fHistLSppPiPiPairMassPt_FromSameSource);
      fOutputList->Add(fHistLSmmPiPiPairMassPt_FromSameSource);

      
      fHistULSPiPiPairMassPt_FromPhi = new TH2F("fHistULSPiPiPairMassPt_FromPhi","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK0s = new TH2F("fHistULSPiPiPairMassPt_FromK0s","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromEta = new TH2F("fHistULSPiPiPairMassPt_FromEta","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromEtaPrime = new TH2F("fHistULSPiPiPairMassPt_FromEtaPrime","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromRho = new TH2F("fHistULSPiPiPairMassPt_FromRho","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromOmega = new TH2F("fHistULSPiPiPairMassPt_FromOmega","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK0star = new TH2F("fHistULSPiPiPairMassPt_FromK0star","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF980 = new TH2F("fHistULSPiPiPairMassPt_FromF980","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1270 = new TH2F("fHistULSPiPiPairMassPt_FromF1270","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1370 = new TH2F("fHistULSPiPiPairMassPt_FromF1370","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1500 = new TH2F("fHistULSPiPiPairMassPt_FromF1500","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1525 = new TH2F("fHistULSPiPiPairMassPt_FromF1525","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1710 = new TH2F("fHistULSPiPiPairMassPt_FromF1710","",800*addbin,0.,4.,100,0,10);      
      fHistULSPiPiPairMassPt_FromLambda1520 = new TH2F("fHistULSPiPiPairMassPt_FromLambda1520","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK01430 = new TH2F("fHistULSPiPiPairMassPt_FromK01430","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromKp1430 = new TH2F("fHistULSPiPiPairMassPt_FromKp1430","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXi0star1530 = new TH2F("fHistULSPiPiPairMassPt_FromXi0star1530","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXi01820 = new TH2F("fHistULSPiPiPairMassPt_FromXi1820","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXip1820 = new TH2F("fHistULSPiPiPairMassPt_FromXip1820","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromSigmastarp1385 = new TH2F("fHistULSPiPiPairMassPt_FromSigmastarp1385","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromSigmastarm1385 = new TH2F("fHistULSPiPiPairMassPt_FromSigmastarm1385","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromPhi_MisID = new TH2F("fHistULSPiPiPairMassPt_FromPhi_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK0s_MisID = new TH2F("fHistULSPiPiPairMassPt_FromK0s_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromEta_MisID = new TH2F("fHistULSPiPiPairMassPt_FromEta_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromEtaPrime_MisID = new TH2F("fHistULSPiPiPairMassPt_FromEtaPrime_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromRho_MisID = new TH2F("fHistULSPiPiPairMassPt_FromRho_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromOmega_MisID = new TH2F("fHistULSPiPiPairMassPt_FromOmega_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK0star_MisID = new TH2F("fHistULSPiPiPairMassPt_FromK0star_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF980_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF980_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1270_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF1270_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1370_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF1370_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1500_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF1500_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1525_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF1525_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromF1710_MisID = new TH2F("fHistULSPiPiPairMassPt_FromF1710_MisID","",800*addbin,0.,4.,100,0,10);      
      fHistULSPiPiPairMassPt_FromLambda1520_MisID = new TH2F("fHistULSPiPiPairMassPt_FromLambda1520_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromK01430_MisID = new TH2F("fHistULSPiPiPairMassPt_FromK01430_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromKp1430_MisID = new TH2F("fHistULSPiPiPairMassPt_FromKp1430_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXi0star1530_MisID = new TH2F("fHistULSPiPiPairMassPt_FromXi0star1530_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXi01820_MisID = new TH2F("fHistULSPiPiPairMassPt_FromXi1820_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromXip1820_MisID = new TH2F("fHistULSPiPiPairMassPt_FromXip1820_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID = new TH2F("fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID","",800*addbin,0.,4.,100,0,10);
      fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID = new TH2F("fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID","",800*addbin,0.,4.,100,0,10);      
      fOutputList->Add(fHistULSPiPiPairMassPt_FromPhi);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK0s);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromEta);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromEtaPrime);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromRho);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromOmega);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK0star);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF980);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1270);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1370);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1500);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1525);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1710);  
      fOutputList->Add(fHistULSPiPiPairMassPt_FromLambda1520);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK01430);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromKp1430);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXi0star1530);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXi01820);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXip1820);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromSigmastarp1385);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromSigmastarm1385);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromPhi_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK0s_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromEta_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromEtaPrime_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromRho_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromOmega_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK0star_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF980_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1270_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1370_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1500_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1525_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromF1710_MisID);  
      fOutputList->Add(fHistULSPiPiPairMassPt_FromLambda1520_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromK01430_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromKp1430_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXi0star1530_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXi01820_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromXip1820_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID);
      fOutputList->Add(fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID);

    }
  }

  
  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{  

  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;
  fEventCounter->Fill(fUtils->getCentClass());
  if(fIsMC)processMC();
  if(fOnK0sAna) V0TrackPairAnalysis();
  if(fOnKaonAna || fOnPionAna) primeTrackPairAnalysis();
  RegisterTreeForTrackPID();
}

bool AliAnalysisTaskAODTrackPair::RegisterTreeForTrackPID()
{

  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track = NULL;
  AliAODMCParticle* particle = NULL;

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    track = (AliAODTrack*)fEvent->GetTrack(iTrack1);
    
    if(!track) continue;

    if( !fUtils->isAcceptPrimeTrackQuality(track) ) continue;

    float sigTOF =track->GetTOFsignal();
    float length =track->GetIntegratedLength();
    float beta =(sigTOF>0)?(double)length/ (2.99792457999999984e-02 * sigTOF):0;

    float dca_xy=9999;
    float dca_z=9999;
    track->GetImpactParameters(dca_xy,dca_z);
    
    AliTOFHeader * tofHeader = (AliTOFHeader*)track->GetTOFHeader();

    fTrackPt = track->Pt();
    fTrackP = track->P();
    fTrackTheta = track->Theta();
    fTrackPhi = track->Phi();
    fTrackLength = track->GetIntegratedLength();
    fTrackBeta = beta;
    fTrackTrackChi2perNDF = track->Chi2perNDF();
    fTrackTrackITSNcls = track->GetITSNcls();
    fTrackTrackTPCNcls = track->GetTPCNcls();
    fTrackTrackTOFNcls = tofHeader->GetNumberOfTOFclusters();
    fTrackTrackTPCChi2 = track->GetTPCchi2();
    fTrackTrackITSChi2 = track->GetITSchi2();
    fTrackTPCCrossedRows = track->GetTPCCrossedRows();
    fTrackTPCFindableNcls = track->GetTPCNclsF();
    fTrackTOFBCTime = track->GetTOFBunchCrossing();
    fTrackTOFKinkIndex = track->GetKinkIndex(0);
    fTrackDCAxy = dca_xy;
    fTrackDCAz = dca_z;
    fTrackTPCsigmaK = fUtils->getTPCSigma(track,AliPID::kKaon);
    fTrackTOFsigmaK = fUtils->getTOFSigma(track,AliPID::kKaon);
    
    fTrackTrueKaonLabel = false;

    if(fIsMC && fMCTrackArray && track->GetLabel()>0){
      particle = (AliAODMCParticle*)fMCTrackArray->At(track->GetLabel());
      if(fabs(particle->GetPdgCode()) == fUtils->fPdgCodeKaon){
	fTrackTrueKaonLabel = true;
      }
    }

    fTreeML->Fill();
    
  }

  return true;
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
    fMCTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCTrackArray) return false;
    fUtils->setMCArray(fMCTrackArray);
  }
  
  return true;

}

bool AliAnalysisTaskAODTrackPair::processMC(){
  
  if(!fMCTrackArray) return false;  

  AliAODMCParticle *particle1;
  AliAODMCParticle *dauparticle1;
  AliAODMCParticle *dauparticle2;

  TLorentzVector t1;
  TLorentzVector t2;
  TLorentzVector t12;
     
  for(Int_t iTrack1=0; iTrack1<fMCTrackArray->GetSize(); ++iTrack1){

    particle1 = (AliAODMCParticle*)fMCTrackArray->At(iTrack1);

    if(!particle1) continue;
    if(!TDatabasePDG::Instance()->GetParticle(particle1->GetPdgCode())) continue;

    int pdgcode = particle1->GetPdgCode();
    
    if(fabs(particle1->Y())<0.5){
      if(pdgcode == fUtils->fPdgCodePhi){
	fHistPhiPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeEta){
	fHistEtaPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeEtaPrime){
	fHistEtaPrimePtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeRho){
	fHistRhoPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeOmega){
	fHistOmegaPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeK0s){
	fHistK0sPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeK0star){
	fHistK0starPtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeF980){
	fHistF980PtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeF1270){
	fHistF1270PtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeF1500){
	fHistF1500PtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeF1525){
	fHistF1525PtTrue->Fill(particle1->Pt());
      }
      else if(pdgcode == fUtils->fPdgCodeF1710){
	fHistF1710PtTrue->Fill(particle1->Pt());
      }


    }

    int ndau = particle1->GetDaughterLast() - particle1->GetDaughterFirst() +1;

    if(ndau<2) continue;

    dauparticle1 = (AliAODMCParticle*)fMCTrackArray->At(particle1->GetDaughterLabel(0)); 
    dauparticle2 = (AliAODMCParticle*)fMCTrackArray->At(particle1->GetDaughterLabel(1)); 
    if(!dauparticle1) continue;
    if(!dauparticle2) continue;

    if(fabs(dauparticle1->GetPdgCode())==fUtils->fPdgCodePion && fabs(dauparticle2->GetPdgCode())==fUtils->fPdgCodePion){
      t1.SetPtEtaPhiM(dauparticle1->Pt(),dauparticle1->Eta(),dauparticle1->Phi(),fUtils->fMassPion);
      t2.SetPtEtaPhiM(dauparticle2->Pt(),dauparticle2->Eta(),dauparticle2->Phi(),fUtils->fMassPion);      
    }
    else if(fabs(dauparticle1->GetPdgCode())==fUtils->fPdgCodeKaon && fabs(dauparticle2->GetPdgCode())==fUtils->fPdgCodeKaon){
      t1.SetPtEtaPhiM(dauparticle1->Pt(),dauparticle1->Eta(),dauparticle1->Phi(),fUtils->fMassKaon);
      t2.SetPtEtaPhiM(dauparticle2->Pt(),dauparticle2->Eta(),dauparticle2->Phi(),fUtils->fMassKaon);      
    }
    else {
      continue;
    }

    t12 = t1 + t2;  

    if(pdgcode == fUtils->fPdgCodePhi && ndau==2){
      fHistPhiMassTrue -> Fill(t12.M());
    }
    else if(pdgcode == fUtils->fPdgCodeRho && ndau==2){
      fHistRhoMassTrue -> Fill(t12.M());
    }
    else if(pdgcode == fUtils->fPdgCodeF980 && ndau==2){
      fHistF980MassTrue -> Fill(t12.M());
    }
    else if(pdgcode == fUtils->fPdgCodeF1270 && ndau==2){
      fHistF1270MassTrue -> Fill(t12.M());
    }
    else if(pdgcode == fUtils->fPdgCodeK0s && ndau==2){
      fHistK0sMassTrue -> Fill(t12.M());
    }
    else{
      continue;
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::primeTrackPairAnalysis(){

  TObjArray* fKaonTrackArray = NULL;
  TObjArray* fPionTrackArray = NULL;

  if(fOnME && fOnKaonAna){
    fKaonTrackArray = new TObjArray();
    fKaonTrackArray -> SetOwner();
  }
  if(fOnME && fOnPionAna){
    fPionTrackArray = new TObjArray();  
    fPionTrackArray -> SetOwner();
  }

  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;
  if(fOnME){
    if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
    if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
    if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi(); 
  }    

  AliEventPool* poolKaon = NULL;
  AliEventPool* poolPion = NULL;
  
  if(fOnME && fOnKaonAna){
    poolKaon = (AliEventPool*)fPoolKaonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  }
  if(fOnME && fOnPionAna){
    poolPion = (AliEventPool*)fPoolPionTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  }
  
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

	if(fIsMC){
	  if(fabs(fUtils->isSameMother(track1,track2))>100 && fUtils->isDiquarkOrigin(fUtils->isSameMother(track1,track2))==false ){	    
	    
	    if(fabs(fUtils->getTrackTruePID(track1)) == fUtils->fPdgCodeKaon && 
	       fabs(fUtils->getTrackTruePID(track2)) == fUtils->fPdgCodeKaon){	      
	      if( track1->Charge() != track2->Charge() ){
		fHistULSKKPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());			
		
		if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodePhi){
		  fHistULSKKPairMassPt_FromPhi->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0s){
		  fHistULSKKPairMassPt_FromK0s->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEta){
		  fHistULSKKPairMassPt_FromEta->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEtaPrime){
		  fHistULSKKPairMassPt_FromEtaPrime->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeRho){
		  fHistULSKKPairMassPt_FromRho->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeOmega){
		  fHistULSKKPairMassPt_FromOmega->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0star){
		  fHistULSKKPairMassPt_FromK0star->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF980){
		  fHistULSKKPairMassPt_FromF980->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1270){
		  fHistULSKKPairMassPt_FromF1270->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1370){
		  fHistULSKKPairMassPt_FromF1370->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1525){
		  fHistULSKKPairMassPt_FromF1525->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1500){
		  fHistULSKKPairMassPt_FromF1500->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1710){
		  fHistULSKKPairMassPt_FromF1710->Fill(t12.M(),t12.Pt());
		}		
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeLambda1520){
		  fHistULSKKPairMassPt_FromLambda1520->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK01430){
		  fHistULSKKPairMassPt_FromK01430->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeKp1430){
		  fHistULSKKPairMassPt_FromKp1430->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi0star1530){
		  fHistULSKKPairMassPt_FromXi0star1530->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi01820){
		  fHistULSKKPairMassPt_FromXi01820->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXip1820){
		  fHistULSKKPairMassPt_FromXip1820->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarp1385){
		  fHistULSKKPairMassPt_FromSigmastarp1385->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarm1385){
		  fHistULSKKPairMassPt_FromSigmastarm1385->Fill(t12.M(),t12.Pt());
		}


		else {

		}
	      }	      
	      else {
		if( track1->Charge()>0 ){
		  fHistLSppKKPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());
		}
		else{
		  fHistLSmmKKPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());
		}
	      }	      
	    }//fabs(fUtils->getTrackTruePID(track1)) == fUtils->fPdgCodeKaon
	    else{	      
	      if( track1->Charge() != track2->Charge() ){		
		if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodePhi){
		  fHistULSKKPairMassPt_FromPhi_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0s){
		  fHistULSKKPairMassPt_FromK0s_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEta){
		  fHistULSKKPairMassPt_FromEta_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEtaPrime){
		  fHistULSKKPairMassPt_FromEtaPrime_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeRho){
		  fHistULSKKPairMassPt_FromRho_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeOmega){
		  fHistULSKKPairMassPt_FromOmega_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0star){
		  fHistULSKKPairMassPt_FromK0star_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF980){
		  fHistULSKKPairMassPt_FromF980_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1270){
		  fHistULSKKPairMassPt_FromF1270_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1370){
		  fHistULSKKPairMassPt_FromF1370_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1525){
		  fHistULSKKPairMassPt_FromF1525_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1500){
		  fHistULSKKPairMassPt_FromF1500_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1710){
		  fHistULSKKPairMassPt_FromF1710_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeLambda1520){
		  fHistULSKKPairMassPt_FromLambda1520_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK01430){
		  fHistULSKKPairMassPt_FromK01430_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeKp1430){
		  fHistULSKKPairMassPt_FromKp1430_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi0star1530){
		  fHistULSKKPairMassPt_FromXi0star1530_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi01820){
		  fHistULSKKPairMassPt_FromXi01820_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXip1820){
		  fHistULSKKPairMassPt_FromXip1820_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarp1385){
		  fHistULSKKPairMassPt_FromSigmastarp1385_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarm1385){
		  fHistULSKKPairMassPt_FromSigmastarm1385_MisID->Fill(t12.M(),t12.Pt());
		}
		else {

		}
	      }
	    }
	    
	  }//fUtils->isSameMother(track1,track2)
	}//fIsMC

      }
      
      if(fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion) && fUtils->isAcceptPrimeTrack(track2,AliPID::kPion)){
	
	t1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),fUtils->fMassPion);
	t2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),fUtils->fMassPion);      
	t12 = t1 + t2;  
	
	if(fUtils->isAcceptTrackPair(t12.Rapidity(),t12.Pt())){
	  if( track1->Charge() != track2->Charge() ){
	    fHistULSPiPiPairMassPt->Fill(t12.M(),t12.Pt());	
	  }
	  else {
	    if( track1->Charge()>0 ){
	      fHistLSppPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	    }
	    else{
	      fHistLSmmPiPiPairMassPt->Fill(t12.M(),t12.Pt());
	    }
	  }
	}       

	if(fIsMC){
	  if(fabs(fUtils->isSameMother(track1,track2))>100 && fUtils->isDiquarkOrigin(fUtils->isSameMother(track1,track2))==false ){	    
	    
	    if(fabs(fUtils->getTrackTruePID(track1)) == fUtils->fPdgCodePion && 
	       fabs(fUtils->getTrackTruePID(track2)) == fUtils->fPdgCodePion){	      
	      if( track1->Charge() != track2->Charge() ){
		fHistULSPiPiPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());			
		
		if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodePhi){
		  fHistULSPiPiPairMassPt_FromPhi->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0s){
		  fHistULSPiPiPairMassPt_FromK0s->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEta){
		  fHistULSPiPiPairMassPt_FromEta->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEtaPrime){
		  fHistULSPiPiPairMassPt_FromEtaPrime->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeRho){
		  fHistULSPiPiPairMassPt_FromRho->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeOmega){
		  fHistULSPiPiPairMassPt_FromOmega->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0star){
		  fHistULSPiPiPairMassPt_FromK0star->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF980){
		  fHistULSPiPiPairMassPt_FromF980->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1270){
		  fHistULSPiPiPairMassPt_FromF1270->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1370){
		  fHistULSPiPiPairMassPt_FromF1370->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1525){
		  fHistULSPiPiPairMassPt_FromF1525->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1500){
		  fHistULSPiPiPairMassPt_FromF1500->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1710){
		  fHistULSPiPiPairMassPt_FromF1710->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeLambda1520){
		  fHistULSPiPiPairMassPt_FromLambda1520->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK01430){
		  fHistULSPiPiPairMassPt_FromK01430->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeKp1430){
		  fHistULSPiPiPairMassPt_FromKp1430->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi0star1530){
		  fHistULSPiPiPairMassPt_FromXi0star1530->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi01820){
		  fHistULSPiPiPairMassPt_FromXi01820->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXip1820){
		  fHistULSPiPiPairMassPt_FromXip1820->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarp1385){
		  fHistULSPiPiPairMassPt_FromSigmastarp1385->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarm1385){
		  fHistULSPiPiPairMassPt_FromSigmastarm1385->Fill(t12.M(),t12.Pt());
		}
		else {

		}
	      }	      
	      else {
		if( track1->Charge()>0 ){
		  fHistLSppPiPiPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());
		}
		else{
		  fHistLSmmPiPiPairMassPt_FromSameSource->Fill(t12.M(),t12.Pt());
		}
	      }	      
	    }//fabs(fUtils->getTrackTruePID(track1)) == fUtils->fPdgCodeKaon
	    else{	      
	      if( track1->Charge() != track2->Charge() ){		
		if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodePhi){
		  fHistULSPiPiPairMassPt_FromPhi_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0s){
		  fHistULSPiPiPairMassPt_FromK0s_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEta){
		  fHistULSPiPiPairMassPt_FromEta_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeEtaPrime){
		  fHistULSPiPiPairMassPt_FromEtaPrime_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeRho){
		  fHistULSPiPiPairMassPt_FromRho_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeOmega){
		  fHistULSPiPiPairMassPt_FromOmega_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK0star){
		  fHistULSPiPiPairMassPt_FromK0star_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF980){
		  fHistULSPiPiPairMassPt_FromF980_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1270){
		  fHistULSPiPiPairMassPt_FromF1270_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1370){
		  fHistULSPiPiPairMassPt_FromF1370_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1525){
		  fHistULSPiPiPairMassPt_FromF1525_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1500){
		  fHistULSPiPiPairMassPt_FromF1500_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeF1710){
		  fHistULSPiPiPairMassPt_FromF1710_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeLambda1520){
		  fHistULSPiPiPairMassPt_FromLambda1520_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeK01430){
		  fHistULSPiPiPairMassPt_FromK01430_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeKp1430){
		  fHistULSPiPiPairMassPt_FromKp1430_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi0star1530){
		  fHistULSPiPiPairMassPt_FromXi0star1530_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXi01820){
		  fHistULSPiPiPairMassPt_FromXi01820_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeXip1820){
		  fHistULSPiPiPairMassPt_FromXip1820_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarp1385){
		  fHistULSPiPiPairMassPt_FromSigmastarp1385_MisID->Fill(t12.M(),t12.Pt());
		}
		else if(fabs(fUtils->isSameMother(track1,track2)) == fUtils->fPdgCodeSigmastarm1385){
		  fHistULSPiPiPairMassPt_FromSigmastarm1385_MisID->Fill(t12.M(),t12.Pt());
		}		
		else {

		}
	      }
	    }
	    
	  }//fUtils->isSameMother(track1,track2)
	}//fIsMC
      }
      
    }//end of loop track2
    
    if(fOnKaonAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon)){
      
      if (fOnME && poolKaon->IsReady()){
      
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
    
    if(fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion)){

      if (fOnME && poolPion->IsReady()){
      
	for (Int_t iMixEvt=0; iMixEvt<poolPion->GetCurrentNEvents(); iMixEvt++){

	  poolTracks = (TObjArray*)poolPion->GetEvent(iMixEvt);	
	  
	  for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){
	  
	    track2 = (AliAODTrack*)poolTracks->At(iTrack2);
	    
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
	    
	  }//end of loop track2
	}// end of loop iMixEvt
	
      }//poolPionon->IsReady()    
    }//if(fUtils->isAcceptPrimeTrack(track1,AliPID::kPion))
    
    if(fOnME && fOnKaonAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kKaon)){      
      fKaonTrackArray->Add(track1);
    }

    if(fOnME && fOnPionAna && fUtils->isAcceptPrimeTrack(track1,AliPID::kPion)){            
      fPionTrackArray->Add(track1);
    }
  }//end of loop track1
  
  TObjArray* fKaonTrackArrayClone = NULL;
  TObjArray* fPionTrackArrayClone = NULL;
  
  if(fOnME && fOnKaonAna){
    fKaonTrackArrayClone = (TObjArray*)fKaonTrackArray->Clone();
    fKaonTrackArrayClone->SetOwner();
    //cout<<"       This event stores  "<<fKaonTrackArrayClone->GetEntriesFast()<<endl;
    if(fKaonTrackArrayClone->GetEntriesFast()>0){
      poolKaon->UpdatePool(fKaonTrackArrayClone);
    }
  }
  
  if(fOnME && fOnPionAna){
    fPionTrackArrayClone = (TObjArray*)fPionTrackArray->Clone();
    fPionTrackArrayClone->SetOwner();
    //cout<<"       This event stores  "<<fPionTrackArrayClone->GetEntriesFast()<<endl;
    if(fPionTrackArrayClone->GetEntriesFast()>0){
      poolPion->UpdatePool(fPionTrackArrayClone);
    }  
  }
 
  return true;
}

bool AliAnalysisTaskAODTrackPair::V0TrackPairAnalysis()
{ 

  
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi(); 

  AliEventPool* pool = NULL;  
  pool = (AliEventPool*)fPoolK0sMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  

  AliAODv0 *v0_1 = NULL;
  AliAODv0 *v0_2 = NULL;

  TLorentzVector k1;
  TLorentzVector k2;
  TLorentzVector k12;      

  TObjArray* poolTracks = NULL;

  int nV0Pair = fEvent->GetNumberOfV0s();

  for(int iV0=0; iV0<nV0Pair; ++iV0){
    
    v0_1 = fEvent->GetV0(iV0);

    if( !fUtils->isAcceptDefaultV0(v0_1) ) continue;
    
    fHistArmenterosPodolanskiPlot->Fill(v0_1->AlphaV0(),v0_1->PtArmV0());
    
    if( !fUtils->isK0sV0(v0_1) ) continue;

    fHistArmenterosPodolanskiPlot_K0s->Fill(v0_1->AlphaV0(),v0_1->PtArmV0());
    
    double fill_1[]={v0_1->MassK0Short(),v0_1->Pt(),fUtils->getCentClass()};

    fHistK0sMassPt->Fill(v0_1->MassK0Short(),v0_1->Pt());

    for(int jV0=iV0+1; jV0<nV0Pair; ++jV0){
    
      v0_2 = fEvent->GetV0(jV0);
    
      if( !fUtils->isAcceptDefaultV0(v0_2) ) continue;
      if( !fUtils->isK0sV0(v0_2) ) continue;
    
      if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;
      if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;
      
      k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
      k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);
      
      k1 + k2;

      double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};
      fHistK0sK0sMassPt->Fill(k12.M(),k12.Pt());
      
    }
    
    if (pool->IsReady()){

      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){

	poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);	

	for(Int_t iV02=0; iV02<poolTracks->GetEntriesFast(); ++iV02){

	  v0_2 = (AliAODv0*)poolTracks->At(iV02);

	  if( !fUtils->isK0sCandidate(v0_1->MassK0Short()) ) continue;
	  if( !fUtils->isK0sCandidate(v0_2->MassK0Short()) ) continue;

	  k1.SetPtEtaPhiM(v0_1->Pt(),v0_1->Eta(),v0_1->Phi(),fUtils->fMassK0s);
	  k2.SetPtEtaPhiM(v0_2->Pt(),v0_2->Eta(),v0_2->Phi(),fUtils->fMassK0s);

	  k12 = k1 + k2;

	  double fill_12[]={k12.M(),k12.Pt(),fUtils->getCentClass()};
	  fHistMixK0sK0sMassPt->Fill(k12.M(),k12.Pt());

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
  
  return true;
}
