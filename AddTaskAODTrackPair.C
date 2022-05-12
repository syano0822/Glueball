AliAnalysisTaskAODTrackPair* AddTaskAODTrackPair(bool isMC=false,UInt_t offlineTriggerMask = AliVEvent::kINT7,
						 float combPISigmaPion = 2.0,
						 float combPISigmaKaon = 2.0,
						 string multi_method="V0M",
						 bool reqTOFhit = false,
						 bool onPURej = true)

						 
{
  

  float min_vtxz =-10;
  float max_vtxz = 10;
  float min_pair_rap = -0.5;
  float max_pair_rap =  0.5;
  float maxTrackChi2perNDF = 4;
  float maxTPCTrackChi2perNDF = 4;
  float maxITSTrackChi2perNDF = 36;
  float minTPCNClust = 70;
  float minCrossRows = 70;
  float minFindTPCclsFrac=0.8;
  float minTrackPt = 0.15;
  float maxTrackPt = 999.0;
  float minTrackEta = -0.8;
  float maxTrackEta = 0.8;
  float minPionSigmaTPC = -2;
  float maxPionSigmaTPC =  2;
  float minKaonSigmaTPC = -2;
  float maxKaonSigmaTPC =  2;
  float minProtonSigmaTPC = -2;
  float maxProtonSigmaTPC =  2;
  float minPionSigmaTOF = -2;
  float maxPionSigmaTOF =  2;
  float minKaonSigmaTOF = -2;
  float maxKaonSigmaTOF =  2;
  float minProtonSigmaTOF = -2;
  float maxProtonSigmaTOF =  2;
  float minV0Radius = 0.5;
  float maxV0Radius = 999.;
  float minV0TrackDCA = 0.06;
  float maxV0TrackDCA = 999.;
  float minV0CosPoint = 0.97;
  float maxV0CosPoint = 1.00;
  float minPropLifetime = 0.;
  float maxPropLifetime = 20.;
  float minRejLambdaMass = 1.105;
  float maxRejLambdaMass = 1.255;
  float combPISigmaProton = 2.0;
  float pcm = 0.207;
  float r0 = 0.85;
  float width = 0.05;
  string primeDcaFunc="0.0105+0.0350/pow(x;1.1)";


  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setVertexCut(min_vtxz,max_vtxz);
  utils->setPromenterosLimit(pcm,r0,width);
  utils->setPairRapidityCut(min_pair_rap,max_pair_rap);
  utils->setMinTPCCrossRows(minCrossRows);
  utils->setPileupRejectionCut(onPURej);
  utils->setMultiEstimateMethod(multi_method);
  utils->setMinNClustTPC(minTPCNClust);
  utils->setTrackChi2perNDF(maxTrackChi2perNDF);
  utils->setTPCTrackChi2perNDF(maxTPCTrackChi2perNDF);
  utils->setITSTrackChi2perNDF(maxITSTrackChi2perNDF);
  utils->setTrackKinematicRange(minTrackPt,maxTrackPt,minTrackEta,maxTrackEta);
  utils->setMinFindableTPCclsFrac(minFindTPCclsFrac);
  utils->setPionSelectSigmaTPC(minPionSigmaTPC, maxPionSigmaTPC);
  utils->setPionSelectSigmaTOF(minPionSigmaTOF, maxPionSigmaTOF);
  utils->setKaonSelectSigmaTPC(minKaonSigmaTPC, maxKaonSigmaTPC);
  utils->setKaonSelectSigmaTOF(minKaonSigmaTOF, maxKaonSigmaTOF);
  utils->setProtonSelectSigmaTPC(minProtonSigmaTPC, maxProtonSigmaTPC);
  utils->setProtonSelectSigmaTOF(minProtonSigmaTOF, maxProtonSigmaTOF);
  utils->setPropLifeTime(minPropLifetime,maxPropLifetime);
  utils->setMinV0Radius(minV0Radius,maxV0Radius);  
  utils->setV0TrackDCA(minV0TrackDCA,maxV0TrackDCA);
  utils->setMinV0CosPoint(minV0CosPoint,maxV0CosPoint);
  utils->setRejLambdaMassRange(minRejLambdaMass,maxRejLambdaMass);
  utils->setRequireTOFhit(reqTOFhit);
  utils->setTPCTOFCombPID(combPISigmaPion,combPISigmaKaon,combPISigmaKaon);
  utils->setPrimeTrackDCAFunction(primeDcaFunc);
  utils->setMultiEstimateMethod("V0M");
  utils->Initialization();
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAODMuonEventSelection", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAODMuonEventSelection", "This task requires an input event handler");
    return NULL;
  }
    
  AliAnalysisTaskAODTrackPair* task = new AliAnalysisTaskAODTrackPair("glueball");
  task->SelectCollisionCandidates(offlineTriggerMask);
  task->setSameEventTrigger(offlineTriggerMask);
  task->setMixingEventTrigger(offlineTriggerMask);
  task->setUtils(utils);
  task->setEvtMixingTrackDepth(1);
  task->setEvtMixingPoolSize(1);
  task->setEvtMixingReadyFraction(0.1);
  task->setEvtMixingPoolVtxZ(true);
  task->setEvtMixingPoolCent(true);
  task->setEvtMixingPoolPsi(true);
  task->setKaonAnalysis(false);
  task->setPionAnalysis(false);
  task->setK0sAnalysis(false);
  task->setTrackQA(false);
  task->setMLSample(true);
  task->setMC(isMC);
  task->setME(false);
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("glueball",TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput1);

  return task;

}
