//Xianglei Zhu, created on Apr. 9, 2008
//topologically reconstruct V0 (Ks or Lambda) from the globaltracks in MuDst.
//
#include "StRoot/StV0Maker/StV0Type.h"
#include "StRoot/StXiMaker/StXiType.h"

//void reconstruct_v0(TString InputFileList, Int_t nFiles = 1000, Int_t nEvents = 0, TString OutputDir = "output/", TString JobIdName = "testrun" );

void reconstruct_v0(TString InputFileList, Int_t nFiles, Int_t nEvents, TString OutputDir, TString JobIdName ) 
{
 
  // Load libraries
  gROOT   -> Macro("loadMuDst.C");

  //gSystem->Load("StDbBroker");
  //gSystem->Load("StDbUtilities");
  //gSystem->Load("St_db_Maker");  //They are used for Run9 to calculate the Primary Vertex

  gSystem -> Load("StV0Maker.so") ;
  gSystem -> Load("StXiMaker.so") ;

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;

  StMuDstMaker*          muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;
  // Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
  muDstMaker -> SetStatus("*",0) ;               // Turn off all branches
  muDstMaker -> SetStatus("MuEvent",1) ;         // Turn on the Event data (esp. Event number)
  muDstMaker -> SetStatus("GlobalTracks",1) ;    // Turn on the global track data
  muDstMaker -> SetStatus("BTofHeader",1) ;    // Turn on the global track data
  muDstMaker -> SetDebug(0) ;                    // Turn off Debug information
//muDstMaker -> SetDebug(1) ;                  //Turn on Debug information
  //cout<<endl<<"============  Data Base ========="<<endl;
  //dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
  //dbMk->SetMaxEntryTime(20100130,0);  //They are used for Run9 to calculate the Primary Vertex

  //first pass to get the raw signal of Lambda
  StV0Maker* rawsig  =  new StV0Maker(muDstMaker,"v0makerfirst") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".lambda.histo.root") ;
  rawsig -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawsig -> setV0TreeFileName(OutputDir+JobIdName+".lambda.picodst.root"); // V0 candidate tree file for further cuts. 
  rawsig -> setV0Type(kLambda);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  rawsig -> setRotate(false); 
  rawsig -> setDumpNull(true);
  rawsig -> SetDebug(0);
  
  //second pass to get the raw signal of AntiLambda
  StV0Maker* antirawsig  =  new StV0Maker(muDstMaker,"v0makersecond") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antilambda.histo.root") ;
  antirawsig -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  antirawsig -> setV0TreeFileName(OutputDir+JobIdName+".antilambda.picodst.root"); // V0 candidate tree file for further cuts. 
  antirawsig -> setV0Type(kAntiLambda);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  antirawsig -> setRotate(false); 
  antirawsig -> SetDebug(0);

  
  //third pass to get the raw signal of Kshort
  StV0Maker* ksrawsig  =  new StV0Maker(muDstMaker,"v0makerthird") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".ks.histo.root") ;
  ksrawsig -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  ksrawsig -> setV0TreeFileName(OutputDir+JobIdName+".ks.picodst.root"); // V0 candidate tree file for further cuts. 
  ksrawsig -> setV0Type(kKs);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  ksrawsig -> setRotate(false); 
//  ksrawsig -> setDumpNull(true); //FIXME: essential option to get correct normalization!!!!
  ksrawsig -> SetDebug(0);
 
 /*
  //fourth pass to estimate the background by rotating the transverse coordinate and momentum of second daughter.
  StV0Maker* backgr  =  new StV0Maker(muDstMaker,"v0makerfourth") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".lambdarot.histo.root") ;
  backgr -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  backgr -> setV0TreeFileName(OutputDir+JobIdName+".lambdarot.picodst.root"); // V0 candidate tree file for further cuts. 
  backgr -> setV0Type(kLambda);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  backgr -> setRotate(true); 
  backgr -> SetDebug(0);

  //fifth pass to estimate the background by rotating the transverse coordinate and momentum of second daughter.
  StV0Maker* backgr2  =  new StV0Maker(muDstMaker,"v0makerfifth") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antilambdarot.histo.root") ;
  backgr2 -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  backgr2 -> setV0TreeFileName(OutputDir+JobIdName+".antilambdarot.picodst.root"); // V0 candidate tree file for further cuts. 
  backgr2 -> setV0Type(kAntiLambda);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  backgr2 -> setRotate(true); 
  backgr2 -> SetDebug(0);

  
  //sixth pass to estimate the background by rotating the transverse coordinate and momentum of second daughter.
  StV0Maker* backgr3  =  new StV0Maker(muDstMaker,"v0makersixth") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".ksrot.histo.root") ;
  backgr3 -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  backgr3 -> setV0TreeFileName(OutputDir+JobIdName+".ksrot.picodst.root"); // V0 candidate tree file for further cuts. 
  backgr3 -> setV0Type(kKs);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! 
  backgr3 -> setRotate(true); 
  backgr3 -> SetDebug(0);
*/  

  
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawxi  =  new StXiMaker(muDstMaker, rawsig, "StXiMakerfirst") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".xi.histo.root") ;
  rawxi -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawxi -> setXiTreeFileName(OutputDir+JobIdName+".xi.picodst.root"); // V0 candidate tree file for further cuts. 
  rawxi -> setXiType(kXi);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawxi -> setRotate(false); 
  rawxi -> SetDebug(0);
/*
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawxirotate  =  new StXiMaker(muDstMaker, rawsig, "StXiMakersecond") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".xirot.histo.root") ;
  rawxirotate -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawxirotate -> setXiTreeFileName(OutputDir+JobIdName+".xirot.picodst.root"); // V0 candidate tree file for further cuts. 
  rawxirotate -> setXiType(kXi);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawxirotate -> setRotate(true); 
  rawxirotate -> SetDebug(0);
 */
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawantixi  =  new StXiMaker(muDstMaker, antirawsig, "StXiMakerAntiXi") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antixi.histo.root") ;
  rawantixi -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawantixi -> setXiTreeFileName(OutputDir+JobIdName+".antixi.picodst.root"); // V0 candidate tree file for further cuts. 
  rawantixi -> setXiType(kAntiXi);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawantixi -> setRotate(false); 
  rawantixi -> SetDebug(0);
/*
   // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawantixirotate  =  new StXiMaker(muDstMaker, antirawsig, "StXiMakerAntiXiRot") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antixirot.histo.root") ;
  rawantixirotate -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawantixirotate -> setXiTreeFileName(OutputDir+JobIdName+".antixirot.picodst.root"); // V0 candidate tree file for further cuts. 
  rawantixirotate -> setXiType(kAntiXi);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawantixirotate -> setRotate(true); 
  rawantixirotate -> SetDebug(0);
 */
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawomega  =  new StXiMaker(muDstMaker, rawsig, "StXiMakerthird") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".omega.histo.root") ;
  rawomega -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawomega -> setXiTreeFileName(OutputDir+JobIdName+".omega.picodst.root"); // V0 candidate tree file for further cuts. 
  rawomega -> setXiType(kOmega);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawomega -> setRotate(false); 
  rawomega -> SetDebug(0);
/*
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawomegarotate  =  new StXiMaker(muDstMaker, rawsig, "StXiMakerfourth") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".omegarot.histo.root") ;
  rawomegarotate -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawomegarotate -> setXiTreeFileName(OutputDir+JobIdName+".omegarot.picodst.root"); // V0 candidate tree file for further cuts. 
  rawomegarotate -> setXiType(kOmega);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawomegarotate -> setRotate(true); 
  rawomegarotate -> SetDebug(0);
 */
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawantiomega  =  new StXiMaker(muDstMaker, antirawsig, "StXiMakerAntiOmega") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antiomega.histo.root") ;
  rawantiomega -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawantiomega -> setXiTreeFileName(OutputDir+JobIdName+".antiomega.picodst.root"); // V0 candidate tree file for further cuts. 
  rawantiomega -> setXiType(kAntiOmega);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawantiomega -> setRotate(false); 
  rawantiomega -> SetDebug(0);
/*
  // use StXiMaker to reconstruct Xi or Omega
  StXiMaker* rawantiomegarotate  =  new StXiMaker(muDstMaker, antirawsig, "StXiMakerAntiOmegaRot") ;
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".antiomegarot.histo.root") ;
  rawantiomegarotate -> setHistoFileName(OutputDir+Name) ; // Name the output file for histograms
  rawantiomegarotate -> setXiTreeFileName(OutputDir+JobIdName+".antiomegarot.picodst.root"); // V0 candidate tree file for further cuts. 
  rawantiomegarotate -> setXiType(kAntiOmega);   //set V0 type. kKs, or kLambda, kAntiLambda. once a time! do not try to mess things up for the time being.
  rawantiomegarotate -> setRotate(true); 
  rawantiomegarotate -> SetDebug(0);
*/
  if ( nEvents == 0 )  nEvents = 10000000 ;       // Take all events in nFiles if nEvents = 0

  // Loop over the links in the chain
  Int_t iInit = chain -> Init() ;
  if (iInit) chain->Fatal(iInit,"on init");
  
  // chain -> EventLoop(1,nEvents) ;  //will output lots of useless debugging info.
  Int_t istat = 0, i = 1;
  while (i <= nEvents && istat != 2) {
     //if(i==740){i++; continue;}
     if(i%10000==0)cout << endl << "== Event " << i << " start ==" << endl;
     chain->Clear();
     istat = chain->Make(i);
     //if(i%1000==0)cout << endl << "== Event " << i << " finish =="<< endl;
     if (istat == 2)
	  cout << "Last  event processed. Status = " << istat << endl;
     if (istat == 3)
	  cout << "Error event processed. Status = " << istat << endl;
     i++;
  }

  if (nEvents > 1) chain -> Finish() ;

  // Cleanup
  delete chain ;
}
