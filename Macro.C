#include "StRoot/StV0Maker/StV0Type.h"
//#include "St_base/Stypes.h"
#include "TString.h"
void Macro(TString inputFileList, Int_t nFiles, TString outputDir, TString jobID){

  //TString inputFile("test.27GeV.MuDst.root");
  TString outPicoFileName = outputDir + jobID + "2013Pico.root";
  TString outHistoFileName = outputDir + jobID + "2013Histo.root";
  Int_t nEvents = 50000;

  StV0Type V0Type = kLambda;

  gROOT -> Macro("loadMuDst.C");
  gSystem -> Load("StV0Maker.so");
  gSystem -> Load("StXiMaker.so");

  StChain* chain = new StChain; 
  StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",inputFileList,"MuDst", nFiles );// nFiles is the max file number to chain together, just a upper bound.

  muDstMaker -> SetStatus("*",0);  
  muDstMaker -> SetStatus("MuEvent",1);
  muDstMaker -> SetStatus("GlobalTracks",1);
  //muDstMaker->SetSatus("BTofHeader",1);
  muDstMaker -> SetDebug(0);

 
  StV0Maker* preSignal = new StV0Maker(muDstMaker, outPicoFileName, outHistoFileName, V0Type);
  preSignal -> SetDebug(0);

  TString xiHistoName = outputDir + jobID + "CascadeHisto.root";
  TString xiPicoName = outputDir + jobID + "CascadePico.root";
  StXiMaker* xiSignal = new StXiMaker(muDstMaker, preSignal, xiPicoName, xiHistoName);
  xiSignal -> SetDebug(0);

  chain -> Init();

  Int_t i = 0;
// EreturnCodes fileStatus=0;
  Int_t fileStatus = 0;
//printf("Event Loop begins\n");
  while(i < nEvents && fileStatus != 2 ){
       fileStatus = chain -> Make(i);
       i++;
       //printf("envent %f", i);
  }

  chain -> Finish();

  printf("%d Events Loop finished\n", i);
  
  delete chain;
}
