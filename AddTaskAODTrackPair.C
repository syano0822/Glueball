AliAnalysisTaskAODTrackPair* AddTaskAODTrackPair(UInt_t offlineTriggerMask = AliVEvent::kINT7,float min_vtxz =-10,float max_vtxz = 10,float min_pair_rap = -0.5,float max_pair_rap =  0.5,string multi_method="SPDTracklets",bool onPURej = true,bool isMC=false)
{
  
  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setVertexCut(min_vtxz,max_vtxz);
  utils->setPairRapidityCut(min_pair_rap,max_pair_rap);
  utils->setPileupRejectionCut(onPURej);
  utils->setMultiEstimateMethod(multi_method);
  
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
