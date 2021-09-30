AliAnalysisTaskAODTrackPair* AddTaskAODTrackPair(UInt_t offlineTriggerMask = AliVEvent::kINT7,
						 float min_vtxz =-10,
						 float max_vtxz = 10,
						 float min_pair_rap = -0.5,
						 float max_pair_rap =  0.5,
						 float maxTrackChi2perNDF = 4,
						 float minCrossRows = 70,
						 float minTrackPt = 0.15,
						 float maxTrackPt = 999.,
						 float minTrackEta = -0.8,
						 float maxTrackEta = 0.8,
						 float minPionSigmaTPC = -2,
						 float maxPionSigmaTPC =  2,
						 float minKaonSigmaTPC = -2,
						 float maxKaonSigmaTPC =  2,
						 float minProtonSigmaTPC = -2,
						 float maxProtonSigmaTPC =  2,
						 float minPionSigmaTOF = -2,
						 float maxPionSigmaTOF =  2,
						 float minKaonSigmaTOF = -2,
						 float maxKaonSigmaTOF =  2,
						 float minProtonSigmaTOF = -2,
						 float maxProtonSigmaTOF =  2,
						 float minV0Radius = 0.5,
						 float maxV0Radius = 999.,
						 float minV0PosTrackDCA = 0.06,
						 float maxV0PosTrackDCA = 999.,
						 float minV0NegTrackDCA = 0.06,
						 float maxV0NegTrackDCA = 999.,
						 float minV0CosPoint = 0.97,
						 float maxV0CosPoint = 1.00,
						 float minPropLifetime = 0.,
						 float maxPropLifetime = 20.,
						 float minRejLambdaMass = 1.105,
						 float maxRejLambdaMass = 1.255,
						 float combPISigmaPion = 3.0,
						 float combPISigmaKaon = 3.0,
						 float combPISigmaProton = 3.0,
						 string multi_method="SPDTracklets",
						 bool onPURej = true,
						 bool isMC=false)
{
  
  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setVertexCut(min_vtxz,max_vtxz);
  utils->setPairRapidityCut(min_pair_rap,max_pair_rap);
  utils->setMinTPCCrossRows(minCrossRows);
  utils->setPileupRejectionCut(onPURej);
  utils->setMultiEstimateMethod(multi_method);
  utils->setTrackChi2perNDF(maxTrackChi2perNDF);
  utils->setTrackKinematicRange(minTrackPt,maxTrackPt,minTrackEta,maxTrackEta);
  utils->setPionSelectSigmaTPC(minPionSigmaTPC, maxPionSigmaTPC);
  utils->setPionSelectSigmaTOF(minPionSigmaTOF, maxPionSigmaTOF);
  utils->setKaonSelectSigmaTPC(minKaonSigmaTPC, maxKaonSigmaTPC);
  utils->setKaonSelectSigmaTOF(minKaonSigmaTOF, maxKaonSigmaTOF);
  utils->setProtonSelectSigmaTPC(minProtonSigmaTPC, maxProtonSigmaTPC);
  utils->setProtonSelectSigmaTOF(minProtonSigmaTOF, maxProtonSigmaTOF);
  utils->setPropLifeTime(minPropLifetime,maxPropLifetime);
  utils->setMinV0Radius(minV0Radius,maxV0Radius);
  utils->setMinV0PosTrackDCA(minV0PosTrackDCA,maxV0PosTrackDCA);
  utils->setMinV0NegTrackDCA(minV0NegTrackDCA,maxV0NegTrackDCA);
  utils->setMinV0CosPoint(minV0CosPoint,maxV0CosPoint);
  utils->setRejLambdaMassRange(minRejLambdaMass,maxRejLambdaMass);
  utils->setTPCTOFCombPID(combPISigmaPion,combPISigmaKaon,combPISigmaKaon);
  
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
  task->setEvtMixingTrackDepth(100);
  task->setEvtMixingPoolSize(100);
  task->setEvtMixingReadyFraction(0.1);
  task->setEvtMixingPoolVtxZ(true);
  task->setEvtMixingPoolCent(true);
  task->setEvtMixingPoolPsi(true);
  task->setMC(isMC);
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("glueball",TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput1);

  return task;

}
