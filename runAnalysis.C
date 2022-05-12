#if !defined(__CINT__) || defined(__CLING__)
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/AddTaskPIDResponse.C>
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <PWGLF/RESONANCES/macros/mini/AddTaskF1500.C>
//#include <./AddTaskF1500.C>
#include <./AddTaskAODEventQA.C>
#endif

AliAnalysisGrid* CreateAlienHandler(TString period, TString run_mode, Bool_t isJDL, TString type);

void runAnalysis(TString runPeriod = "LHC16k",
		 TString run_mode  = "test",
		 Bool_t isJDL      = true,
		 TString type      = "data",
		 Bool_t  local     = false)
{
  
  // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
  gInterpreter->ProcessLine("./");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include ./");
#endif
  
  Bool_t isMC = false;
  
  if(type == "data"){

  }
  else{
    isMC = true;
  }

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  if(local) {
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("AliAOD.root");
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
  }
  else {
    // also specify the include (header) paths on grid        
    AliAnalysisGrid *alienHandler = CreateAlienHandler(runPeriod,run_mode,isJDL,type);
    if (!alienHandler) return;
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
  }
  
  // -----------------------------------------
  //            PHYSICS SELECTION
  // -----------------------------------------
#if !defined (__CINT__) || defined (__CLING__)
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);//isMC,PileupCuts
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  if(run_mode == "full") AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);
#endif
    
  // -----------------------------------------
  //               MULT SELECTION
  // -----------------------------------------
#if !defined (__CINT__) || defined (__CLING__)
  AliMultSelectionTask *multSelTask=reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C ()"));
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *multSelTask = AddTaskMultSelection(false);
  multSelTask->SetUseDefaultCalib(true);
  if(isMC) multSelTask->SetUseDefaultMCCalib(kTRUE);
  else     multSelTask->SetUseDefaultCalib(kTRUE);
#endif
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* taskPID = AddTaskPIDResponse(isMC,false,true,1);
  
  gROOT->LoadMacro("./AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gROOT->LoadMacro("./AliAnalysisTaskAODTrackPair.cxx++g");
  gROOT->LoadMacro("./AddTaskAODTrackPair.C");

  
  /*
  AliAnalysisTaskAODTrackPair* trackpair = AddTaskAODTrackPair(AliVEvent::kINT7,
							       -10,
							       10,
							       -0.5,
							       0.5,
							       70,
							       0.15,
							       999.
							       -0.8,
							       0.8,
							       -2,2,
							       -2,2,
							       -2,2,
							       -2,2,
							       -2,2,
							       -2,2,
							       0.5,200.,
							       0.06,999.,
							       0.06,999.,
							       0.97,1.00,
							       1.109,1.121,
							       2.0,2.0,2.0,
							       "V0M",
							       1,
							       0);
  */

  UInt_t offlineTriggerMask = AliVEvent::kINT7;
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
  float combPISigmaPion = 2.0;
  float combPISigmaKaon = 2.0;
  float combPISigmaProton = 2.0;
  float pcm = 0.207;
  float r0 = 0.85;
  float width = 0.05;
  string multi_method="V0M";
  string primeDcaFunc="0.0105+0.0350/pow(x;1.1)";
  bool reqTOFhit = false;
  bool onPURej = true;

  AliAnalysisTaskAODTrackPair* trackpair = AddTaskAODTrackPair(isMC);

  // -----------------------------------------
  //               Add Task V0Reader
  // -----------------------------------------  
#if !defined (__CINT__) || defined (__CLING__)
  //AddTaskGlueball();
#else
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddTaskGlueball.C");
  
  //gROOT->LoadMacro("./AddTaskGlueball.C");
  //AliRsnMiniAnalysisTask *task = AddTaskGlueball("gluball",isMC,1,AliPIDResponse::kPP,AliVEvent::kINT7,0,0,0,32,1,AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,1.0,1.0,0.2,2.5,2300,0.0,10.0,20,0.001,10.0,1,1);
  //AddTaskGlueball("glueball",isMC,1,AliPIDResponse::kPP,AliVEvent::kINT7,0,0,0,5,1,AliRsnCutSetDaughterParticle::kTPCpidTOFveto,2.0,2.0,0.2,2.5,2300,0.0,10.0,20,0.001,10.0,1,1);
#endif
  if(!mgr->InitAnalysis()) return;
  
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);

  mgr->StartAnalysis("grid");
  
}

AliAnalysisGrid* CreateAlienHandler(TString runPeriod, TString run_mode, Bool_t isJDL, TString type){
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  plugin->SetMergeViaJDL(isJDL);
  if(run_mode=="terminate")
    plugin->SetRunMode("terminate");
  if(run_mode=="full")
    plugin->SetRunMode("full");
  if(run_mode=="test")
    plugin->SetRunMode("test");
  if(run_mode=="offline")
    plugin->SetRunMode("offline");
  
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20210822-1");
  plugin->SetDefaultOutputs(kTRUE);
  
  plugin->SetGridWorkingDir("PWGLF/AOD/"+runPeriod+"/GlueballK0sK0s/"+type);
  plugin->SetGridOutputDir("output");
    
  if(type == "data"){
    if (runPeriod.Contains("LHC17p")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17p");    
      if(runPeriod.Contains("woSDD")) plugin->SetDataPattern("/pass1_CENT_woSDD/AOD234/ AliAOD.root");
      else                            plugin->SetDataPattern("/pass1_FAST/AOD234/ AliAOD.root");      
      plugin->SetSplitMaxInputFileNumber(20);
      //plugin->SetSplitMaxInputFileNumber(1);
    }
    else if (runPeriod.Contains("LHC17q")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17q");    
      if(runPeriod.Contains("woSDD")) plugin->SetDataPattern("/pass1_CENT_woSDD/AOD234/ AliAOD.root");
      else                         plugin->SetDataPattern("/pass1_FAST/AOD234/ AliAOD.root");      
      plugin->SetSplitMaxInputFileNumber(20);
      //plugin->SetSplitMaxInputFileNumber(1);
    }
    else if (runPeriod.Contains("LHC16q")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16q");    
      if(runPeriod.Contains("woSDD")) plugin->SetDataPattern("/pass2_CENT_woSDD/AOD244/ AliAOD.root");
      else                         plugin->SetDataPattern("/pass2_FAST/AOD244/ AliAOD.root");      
      plugin->SetSplitMaxInputFileNumber(20);
    }
    else if (runPeriod.Contains("LHC16t")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16t");
      if(runPeriod.Contains("woSDD")) plugin->SetDataPattern("/pass2_CENT_woSDD/AOD244/ AliAOD.root");
      else                         plugin->SetDataPattern("/pass2_FAST/AOD244/ AliAOD.root");      
      plugin->SetSplitMaxInputFileNumber(20);
    }
    
    plugin->SetRunPrefix("000");    
  }

  else if(type == "LHC21j1a_cent_woSDD"){
    plugin->SetGridDataDir("/alice/sim/2021/LHC21j1a_cent_woSDD");
    plugin->SetDataPattern("/AOD/ AliAOD.root");
  }
  else if(type == "LHC21j1a_fast"){
    plugin->SetGridDataDir("/alice/sim/2021/LHC21j1a_fast");
    plugin->SetDataPattern("/AOD/ AliAOD.root");
  }
  
  else if(type == "LHC20g14a"){
    plugin->SetGridDataDir("/alice/sim/2020/LHC20g14a");
    plugin->SetDataPattern("/AOD243/ AliAOD.root");
  }

  else if(type == "LHC18c8b_fast"){
    plugin->SetGridDataDir("/alice/sim/2018/LHC18c8b_fast");
    plugin->SetDataPattern("/AOD/ AliAOD.root");
  }

  else if(type == "LHC18c8b_cent_woSDD"){
    plugin->SetGridDataDir("/alice/sim/2018/LHC18c8b_cent_woSDD");
    plugin->SetDataPattern("/AOD/ AliAOD.root");
  }

  else if(type == "LHC18c8b_cent"){
    plugin->SetGridDataDir("/alice/sim/2018/LHC18c8b_cent");
    plugin->SetDataPattern("/AOD/ AliAOD.root");
  }
  
  else if(type == "LHC19h6a"){
    plugin->SetGridDataDir("/alice/sim/2020/LHC19h6a");
    plugin->SetDataPattern("/AOD243/ AliAOD.root");
  }
   
  if(runPeriod.Contains("LHC17q")){//p-p 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282365);
    plugin->AddRunNumber(282366);
    plugin->AddRunNumber(282367);
  }
  
  if(runPeriod.Contains("LHC17p_FAST") || runPeriod.Contains("LHC17l4")){//p-p 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282343); plugin->AddRunNumber(282342); plugin->AddRunNumber(282341); plugin->AddRunNumber(282340); plugin->AddRunNumber(282314); plugin->AddRunNumber(282313); 
    plugin->AddRunNumber(282312); plugin->AddRunNumber(282309); plugin->AddRunNumber(282307); plugin->AddRunNumber(282306); plugin->AddRunNumber(282305); plugin->AddRunNumber(282304);
    plugin->AddRunNumber(282303); plugin->AddRunNumber(282302); plugin->AddRunNumber(282247); plugin->AddRunNumber(282230); plugin->AddRunNumber(282229); plugin->AddRunNumber(282227); 
    plugin->AddRunNumber(282224); plugin->AddRunNumber(282206); plugin->AddRunNumber(282189); plugin->AddRunNumber(282147); plugin->AddRunNumber(282146); plugin->AddRunNumber(282127); 
    plugin->AddRunNumber(282126); plugin->AddRunNumber(282125); plugin->AddRunNumber(282123); plugin->AddRunNumber(282122); plugin->AddRunNumber(282120); plugin->AddRunNumber(282119); 
    plugin->AddRunNumber(282118); plugin->AddRunNumber(282099); plugin->AddRunNumber(282098); plugin->AddRunNumber(282078); plugin->AddRunNumber(282051); plugin->AddRunNumber(282050);
    plugin->AddRunNumber(282031); plugin->AddRunNumber(282030); plugin->AddRunNumber(282025); plugin->AddRunNumber(282021); plugin->AddRunNumber(282016); plugin->AddRunNumber(282008);
 }

  if(runPeriod.Contains("LHC17p_woSDD_No1")){//p-p 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282343); plugin->AddRunNumber(282342); plugin->AddRunNumber(282341); plugin->AddRunNumber(282340); plugin->AddRunNumber(282314); plugin->AddRunNumber(282313); 
    plugin->AddRunNumber(282312); plugin->AddRunNumber(282309); plugin->AddRunNumber(282307); plugin->AddRunNumber(282306); plugin->AddRunNumber(282305); plugin->AddRunNumber(282304); 
    plugin->AddRunNumber(282303); plugin->AddRunNumber(282302); plugin->AddRunNumber(282247); plugin->AddRunNumber(282230); plugin->AddRunNumber(282229); plugin->AddRunNumber(282227); 
    plugin->AddRunNumber(282224); plugin->AddRunNumber(282206); plugin->AddRunNumber(282189); plugin->AddRunNumber(282147); plugin->AddRunNumber(282146); plugin->AddRunNumber(282127); 
    plugin->AddRunNumber(282126); plugin->AddRunNumber(282125); plugin->AddRunNumber(282123); plugin->AddRunNumber(282122); plugin->AddRunNumber(282120); plugin->AddRunNumber(282119); 
    plugin->AddRunNumber(282118); plugin->AddRunNumber(282098); plugin->AddRunNumber(282078); plugin->AddRunNumber(282051); plugin->AddRunNumber(282050);
    plugin->AddRunNumber(282031); plugin->AddRunNumber(282030); plugin->AddRunNumber(282025); plugin->AddRunNumber(282021); plugin->AddRunNumber(282016); plugin->AddRunNumber(282008);
  }
  if(runPeriod.Contains("LHC17p_woSDD_No2")){//p-p 5 TeV
    plugin->AddRunNumber(282099);
  }


  if(runPeriod.Contains("LHC16q")){//p-Pb 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(265525); plugin->AddRunNumber(265521); plugin->AddRunNumber(265501); plugin->AddRunNumber(265500); plugin->AddRunNumber(265499); plugin->AddRunNumber(265435);
    plugin->AddRunNumber(265427); plugin->AddRunNumber(265426); plugin->AddRunNumber(265425); plugin->AddRunNumber(265424); plugin->AddRunNumber(265422); plugin->AddRunNumber(265421);
    plugin->AddRunNumber(265420); plugin->AddRunNumber(265419); plugin->AddRunNumber(265388); plugin->AddRunNumber(265387); plugin->AddRunNumber(265385); plugin->AddRunNumber(265384);
    plugin->AddRunNumber(265383); plugin->AddRunNumber(265381); plugin->AddRunNumber(265378); plugin->AddRunNumber(265377); plugin->AddRunNumber(265344); plugin->AddRunNumber(265343);
    plugin->AddRunNumber(265342); plugin->AddRunNumber(265339); plugin->AddRunNumber(265338); plugin->AddRunNumber(265336); plugin->AddRunNumber(265334); plugin->AddRunNumber(265332);
    plugin->AddRunNumber(265309);
  }

  if(runPeriod.Contains("LHC16t")){//p-Pb 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(267166); plugin->AddRunNumber(267165); plugin->AddRunNumber(267164); plugin->AddRunNumber(267163);
  }

  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  plugin->SetAnalysisSource("AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODTrackPair.cxx");
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libANALYSISalice.so libANALYSISaliceBase.so libCORRFW.so libOADB.so libCore.so libTree.so libGeom.so libVMC.so libPhysics.so "
			    "AliAnalysisTaskAODTrackPairUtils.h AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODTrackPair.h AliAnalysisTaskAODTrackPair.cxx");
  
  //Set Job
  plugin->SetExecutableCommand("aliroot -b -q");
  plugin->SetAnalysisMacro("analysis_syano_"+type+"_"+runPeriod+".C");
  plugin->SetExecutable("analysis_syano_"+type+"_"+runPeriod+".sh");
  
  if(type == "data") plugin->SetNtestFiles(2);
  else               plugin->SetNtestFiles(5);
  
  
  
  plugin->SetOutputToRunNo();
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName("analysis_syano_MC.jdl");
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");
  
  return plugin;
}
