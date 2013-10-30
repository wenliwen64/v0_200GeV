#include "StV0Maker.h"
#include "StDcaService.h"
#include "StarClassLibrary/SystemOfUnits.h"

#include <iostream>

#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "tables/St_vertexSeed_Table.h"
#include "StBTofHeader.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"

ClassImp(StV0Maker)                   // Macro for CINT compatibility

StV0Maker::StV0Maker( StMuDstMaker* maker, const char * name) : StMaker(name)
{ // Initialize and/or zero all public/private data members here.

  //for ( Int_t i = 0 ; i < kMaxNumberOfTH1F ; i++ )  // Zero the histogram pointers, not necessary. it is NULL naturaly.
  //  {
  //    histogram[i] = NULL ;
  //  }

  mMuDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
  mV0Type = kLambda;	//Lambda as default!

  mRotate = false;
  mSameSignPlus = false;
  mSameSignMinus = false;
  
  mDcaAlgoLong = true;

  mDumpNull = false;

  histogram_output = NULL  ;                    // Zero the Pointer to histogram output file
  v0tree_output = NULL  ;                    // Zero the Pointer to v0 tree output file
  mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro
  mV0TreeOutputFileName = "" ;               // V0 Output File Name will be set inside the "analysis".C macro

  mV0Tree = NULL ;

  mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker 
  mTestNTrack = 0;
  mTestVZ = 0;

  mBeamHelix = NULL;

}

StV0Maker::~StV0Maker() 
{ // Destroy and/or zero out all public/private data members here.
}

void StV0Maker::initConst(){
   //initialize the constant for different V0 types.
   
   if(mV0Type == kLambda || mV0Type == kAntiLambda){
	// for Lambda and AntiLambda
	mMass1      = 0.93827; // mass of proton
	mMass2      = 0.13957; // mass of pion
	mMassV0     = 1.115684;// mass of Lambda

	if(mV0Type == kLambda) { 
	   mCharge1	= 1;
	   mCharge2	= -1;
	}
	else {
	   mCharge1	= -1;
	   mCharge2	= 1;
	}
	//do not setup the cut values here. those are parameters
	
	//parameters for StDcaService.cxx
	kShiftConnect = 0.3;
	kShiftContain = 0.3;
   }
   else if(mV0Type == kKs){
	// for Ks
	mMass1     = 0.13957;     // mass of pion+
	mMass2     = 0.13957;     // mass of pion-
	mMassV0    = 0.49768;     // mass of K0s
	
	mCharge1	= 1;
	mCharge2	= -1;
	//do not setup the cut values here. those are parameters
	
	//parameters for StDcaService.cxx
	kShiftConnect = 0.3;
	kShiftContain = 0.3;
   }
   else {
	// for photon
	assert(mV0Type == kPhoton); 
	mMass1	= 0.51099907e-3;
	mMass2	= 0.51099907e-3;
	mMassV0	= 0;

	mCharge1	= 1;
	mCharge2	= -1;
	if(mSameSignPlus){
	   mCharge1    = 1;
	   mCharge2    = 1;
	}
	if(mSameSignMinus){
	   mCharge1    = -1;
	   mCharge2    = -1;
	}
	//parameters for StDcaService.cxx
	kShiftConnect = 2;
	kShiftContain = 2;
   }

   return;
}

void StV0Maker::initParam(){
   //setup the cut values here. do not hard-code them in ::Make()

   cutAbsVertexZLeEq  = 50.;
   //cutTriggerIdEq  = 210000; //dau run8 zdce
//   cutTriggerIdEq  = 260001; //run10 test
   //cutTriggerIdEq  = 200020; //auau run7
   //cutTriggerIdEq  = 66007;	 //cucu run5
   cutTriggerIdEq  = 360001;	//auau Run11 27GeV
   
   if(mV0Type == kPhoton){
	cutNHitsGr = 10;
	cutPtGrEq = 0.006;

	cutAbsNSigma1Le = 4.;
	cutAbsNSigma2Le = 4.;
	cutDca1GrEq  = 0;
	cutDca2GrEq  = 0;
	//cutDca1LeEq  = 1.5;
	//cutDca2LeEq  = 1.5;

	cutDca1to2LeEq = 1.5;
	cutV0MassWidthLeEq = 0.2;
	//cutDauPtArmLeEq = 0.3;
	//cutAbsDausPtShoulderDiffLeEq = 1.2;
	//cutDau1DecAngGr = 0.;
	//cutDau2DecAngGr = 0.;
	//cutV0rdotpGr  = 0.;
	cutDcaV0Le    = 1.5;
	cutV0DecLenGrEq = 0;
	cutDau1Dau2Ang3DLe = 0.2;
	cutDau1Dau2DipAngDiffLe = 0.05;
   }
   else if(mV0Type == kLambda || mV0Type == kAntiLambda){
	cutNHitsGr = 15;
	cutPtGrEq = 0.15;

	cutAbsNSigma1Le = 4.;
	cutAbsNSigma2Le = 4.;
	cutDca1GrEq  = 0; //0.5
	cutDca2GrEq  = 0; //1.5
	//cutDca1LeEq  = 3.0;
	//cutDca2LeEq  = 1.5;

	cutDca1to2LeEq = 1.0; //0.8
	cutV0MassWidthLeEq = 0.07;
	//cutDauPtArmLeEq = 0.3;
	//cutAbsDausPtShoulderDiffLeEq = 1.2;
	//cutDau1DecAngGr = 0.;
	//cutDau2DecAngGr = 0.;
	cutV0rdotpGr  = 0.;
	cutDcaV0Le    = 5.0; //3.5
	cutV0DecLenGrEq =3.0;
   }
   else {
	cutNHitsGr = 15;
	cutPtGrEq = 0.15;

	cutAbsNSigma1Le = 4.;
	cutAbsNSigma2Le = 4.;
	cutDca1GrEq  = 0;
	cutDca2GrEq  = 0;
	//cutDca1LeEq  = 1.5;
	//cutDca2LeEq  = 1.5;

	cutDca1to2LeEq = 1.0;
	cutV0MassWidthLeEq = 0.07;
	//cutDauPtArmLeEq = 0.3;
	//cutAbsDausPtShoulderDiffLeEq = 1.2;
	//cutDau1DecAngGr = 0.;
	//cutDau2DecAngGr = 0.;
	cutV0rdotpGr  = 0.;
	cutDcaV0Le    = 3.0;
	cutV0DecLenGrEq =2.0;
   }


   return;
}

void StV0Maker::initHisto()
{
   // Create Histograms
   // there is no better way to set QA histograms. do not use histogram arrays or vectors.
   // it is not useful. there are no need to operate all the histograms at the same time.

   const Int_t    nbins    =  100   ;

   //QA for events
   hNPrimVertex  = new TH1F( "PrimVertex", "Number of Primary Vertex", 10, 0.0, 10.0 ) ;
   hVertexZ  = new TH1F( "VertexZ", "Event Vertex Z Position", nbins*4, -100.0, 100.0 ) ; 
   hNRefMult  = new TH1F( "RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
   hSelectNRefMult  = new TH1F( "SelectRefMult", "Reference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
   hSelectNRefMultM1  = new TH1F( "SelectRefMultM1", "Reference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
   hSelectNRefMultM2  = new TH1F( "SelectRefMultM2", "Reference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
   hSelectNRefMultFtpcE  = new TH1F( "SelectRefMultFtpcE", "Reference Multiplicity of selected events (FTPC East)", 1000, 0.0, 1000.0 ) ;
   hSelectBbcEAdcSum = new TH1F("SelectBbcEAdcSum", "BBC East Adc Sum", 5000, 0.0, 5000.0);

   //QA for global tracks
   hPtRaw  = new TH1F( "PtRaw", "Transverse Momentum for all particles", nbins*6, 0.0, 30.0 ) ;
   hEtaRaw  = new TH1F( "EtaRaw", "Eta for all particles", nbins*5, -2, 2 ) ;
   hPhiRaw  = new TH1F( "PhiRaw", "Phi for all particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
   hPt  = new TH1F( "Pt", "Transverse Momentum for selected particles", nbins*6, 0.0, 30.0 ) ;
   hEta  = new TH1F( "Eta", "Eta for selected particles", nbins*5, -2, 2 ) ;
   hPhi  = new TH1F( "Phi", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
   hPhiLowPt  = new TH1F( "PhiLowPt", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
   hPhiHighPt  = new TH1F( "PhiHighPt", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
   hDedxP  = new TH2F( "DedxP", "dEdx for selected particles", nbins*3, -3, 3, nbins*3, 0, 3e-5) ;
   hNSigmaPion  = new TH1F( "nSigmaPion", "nSigmaPion for selected particles", nbins*2, -10, 10 ) ;
   hNSigmaProton  = new TH1F( "nSigmaProton", "nSigmaProton for selected particles", nbins*2, -10, 10 ) ;
   hNSigmaKaon  = new TH1F( "nSigmaKaon", "nSigmaKaon for selected particles", nbins*2, -10, 10 ) ;
   hNHitsFit  = new TH1F( "nHitsFit", "nHitsFit for all particles", 80, 0, 80 ) ;
   hNHits  = new TH1F( "nHits", "nHits for all particles", 80, 0, 80 ) ;
   
   hPtDiff  = new TH1F( "ptdiff", "Reference Multiplicity", 100, 0.0, 0.02 ) ;
   hOrDiff  = new TH1F( "ordiff", "Reference Multiplicity", 100, 0.0, 0.2 ) ;
   hPDiff  = new TH1F( "pdiff", "Reference Multiplicity", 100, 0.0, 0.4 ) ;
   hDcaDiff  = new TH1F( "dcadiff", "Reference Multiplicity", 100, -0.4, 0.4 ) ;

   hmDau1Dca = new TH1F("hmDau1Dca","Dca to PV of dau1", 500,0,25);
   hmDau2Dca = new TH1F("hmDau2Dca","Dca to PV of dau2", 500,0,25);
   hmDau1Dca2 = new TH1F("hmDau1Dca2","Dca to PV of dau1", 500,0,25);
   hmDau2Dca2 = new TH1F("hmDau2Dca2","Dca to PV of dau2", 500,0,25);

   //QA for V0's
   hInvMass  = new TH1F( "V0Mass", "V0 Inv. Mass", 200, mMassV0-cutV0MassWidthLeEq, mMassV0+cutV0MassWidthLeEq ) ;

   return;
}

void StV0Maker::initTree()
{
   //initialize the TTree for StV0Dst
   mV0Tree = new TTree("V0PicoDst","V0PicoDst from StV0Maker");

   mV0Tree->SetDirectory(v0tree_output);

   mV0Tree->Branch("runnumber",&mV0Dst.runnumber,"runnumber/I");
   mV0Tree->Branch("evtnumber",&mV0Dst.evtnumber,"evtnumber/I");
   mV0Tree->Branch("trgmode",&mV0Dst.trgmode,"trgmode/I");
   mV0Tree->Branch("nrefmult",&mV0Dst.nrefmult,"nrefmult/I");
   mV0Tree->Branch("nrefmultTOF",&mV0Dst.nrefmultTOF,"nrefmultTOF/I");
   mV0Tree->Branch("bbcadcsumeast",&mV0Dst.bbcadcsumeast,"bbcadcsumeast/I");
   mV0Tree->Branch("zdcadc0",&mV0Dst.zdcadc0,"zdcadc0/F");
   mV0Tree->Branch("bbccirate",&mV0Dst.bbccirate,"bbccirate/D");
   mV0Tree->Branch("primvertexX",&mV0Dst.primvertexX,"primvertexX/F");
   mV0Tree->Branch("primvertexY",&mV0Dst.primvertexY,"primvertexY/F");
   mV0Tree->Branch("primvertexZ",&mV0Dst.primvertexZ,"primvertexZ/F");
   mV0Tree->Branch("magn",&mV0Dst.magn,"magn/F");
   mV0Tree->Branch("nv0",&mV0Dst.nv0,"nv0/I");

   mV0Tree->Branch("v0mass",mV0Dst.v0mass,"v0mass[nv0]/F");
   mV0Tree->Branch("v0pt",mV0Dst.v0pt,"v0pt[nv0]/F");
   mV0Tree->Branch("v0rapidity",mV0Dst.v0rapidity,"v0rapidity[nv0]/F");
   mV0Tree->Branch("v0eta",mV0Dst.v0eta,"v0eta[nv0]/F");
   mV0Tree->Branch("v0x",mV0Dst.v0x,"v0x[nv0]/F");
   mV0Tree->Branch("v0y",mV0Dst.v0y,"v0y[nv0]/F");
   mV0Tree->Branch("v0z",mV0Dst.v0z,"v0z[nv0]/F");
   mV0Tree->Branch("v0px",mV0Dst.v0px,"v0px[nv0]/F");
   mV0Tree->Branch("v0py",mV0Dst.v0py,"v0py[nv0]/F");
   mV0Tree->Branch("v0pz",mV0Dst.v0pz,"v0pz[nv0]/F");
   mV0Tree->Branch("v0declen",mV0Dst.v0declen,"v0declen[nv0]/F");
   mV0Tree->Branch("v0dca",mV0Dst.v0dca,"v0dca[nv0]/F");
   mV0Tree->Branch("v0dca2d",mV0Dst.v0dca2d,"v0dca2d[nv0]/F");
   mV0Tree->Branch("v0pathlen",mV0Dst.v0pathlen,"v0pathlen[nv0]/F");
   mV0Tree->Branch("dau1id",mV0Dst.dau1id,"dau1id[nv0]/I");
   mV0Tree->Branch("dau2id",mV0Dst.dau2id,"dau2id[nv0]/I");
   mV0Tree->Branch("dau1dca",mV0Dst.dau1dca,"dau1dca[nv0]/F");
   mV0Tree->Branch("dau1dca2d",mV0Dst.dau1dca2d,"dau1dca2d[nv0]/F");
   mV0Tree->Branch("dau1nhits",mV0Dst.dau1nhits,"dau1nhits[nv0]/I");
   mV0Tree->Branch("dau1dedx",mV0Dst.dau1dedx,"dau1dedx[nv0]/F");
   mV0Tree->Branch("dau1nsigma",mV0Dst.dau1nsigma,"dau1nsigma[nv0]/F");
   mV0Tree->Branch("dau1eta",mV0Dst.dau1eta,"dau1eta[nv0]/F");
   mV0Tree->Branch("dau1pt",mV0Dst.dau1pt,"dau1pt[nv0]/F");
   mV0Tree->Branch("dau1px",mV0Dst.dau1px,"dau1px[nv0]/F");
   mV0Tree->Branch("dau1py",mV0Dst.dau1py,"dau1py[nv0]/F");
   mV0Tree->Branch("dau1pz",mV0Dst.dau1pz,"dau1pz[nv0]/F");
   mV0Tree->Branch("dau1tpc",mV0Dst.dau1tpc,"dau1tpc[nv0]/I");
   mV0Tree->Branch("dau1ssd",mV0Dst.dau1ssd,"dau1ssd[nv0]/I");
   mV0Tree->Branch("dau1svt",mV0Dst.dau1svt,"dau1svt[nv0]/I");
   mV0Tree->Branch("dau1tofflag",mV0Dst.dau1tofflag,"dau1tofflag[nv0]/I");
   mV0Tree->Branch("dau1tof",mV0Dst.dau1tof,"dau1tof[nv0]/F");
   mV0Tree->Branch("dau1pathlen",mV0Dst.dau1pathlen,"dau1pathlen[nv0]/F");
   mV0Tree->Branch("dau2dca",mV0Dst.dau2dca,"dau2dca[nv0]/F");
   mV0Tree->Branch("dau2dca2d",mV0Dst.dau2dca2d,"dau2dca2d[nv0]/F");
   mV0Tree->Branch("dau2nhits",mV0Dst.dau2nhits,"dau2nhits[nv0]/I");
   mV0Tree->Branch("dau2dedx",mV0Dst.dau2dedx,"dau2dedx[nv0]/F");
   mV0Tree->Branch("dau2nsigma",mV0Dst.dau2nsigma,"dau2nsigma[nv0]/F");
   mV0Tree->Branch("dau2eta",mV0Dst.dau2eta,"dau2eta[nv0]/F");
   mV0Tree->Branch("dau2pt",mV0Dst.dau2pt,"dau2pt[nv0]/F");
   mV0Tree->Branch("dau2px",mV0Dst.dau2px,"dau2px[nv0]/F");
   mV0Tree->Branch("dau2py",mV0Dst.dau2py,"dau2py[nv0]/F");
   mV0Tree->Branch("dau2pz",mV0Dst.dau2pz,"dau2pz[nv0]/F");
   mV0Tree->Branch("dau2tpc",mV0Dst.dau2tpc,"dau2tpc[nv0]/I");
   mV0Tree->Branch("dau2ssd",mV0Dst.dau2ssd,"dau2ssd[nv0]/I");
   mV0Tree->Branch("dau2svt",mV0Dst.dau2svt,"dau2svt[nv0]/I");
   mV0Tree->Branch("dau2tofflag",mV0Dst.dau2tofflag,"dau2tofflag[nv0]/I");
   mV0Tree->Branch("dau2tof",mV0Dst.dau2tof,"dau2tof[nv0]/F");
   mV0Tree->Branch("dau2pathlen",mV0Dst.dau2pathlen,"dau2pathlen[nv0]/F");
   mV0Tree->Branch("dca1to2",mV0Dst.dca1to2,"dca1to2[nv0]/F");

   return;
}

Int_t StV0Maker::Init( )
{
  // setup the constants according to mV0Type
  initConst();

  // initialize parameters (cuts)
  initParam();
  
  // Create Histogram output file
  if(mHistogramOutputFileName == "") { 
     //CAUTION: ALWAYS USE { } HERE!!! LOG_XXX is a if()xxx macro!!!
     LOG_ERROR << "StV0Maker: Please specify the histrogram output file" <<endm;
     exit(-1);
  }
  else {
     histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;  
  }
  // Book histograms
  initHisto();
  
  // Create V0 Tree output file
  if(mV0TreeOutputFileName == "") {
     LOG_WARN << "StV0Maker: The V0 tree output file is not specified! output is smeared!" <<endm;
  }
  else {
     v0tree_output = new TFile( mV0TreeOutputFileName, "recreate" ) ;
     // Create V0 Tree
     initTree();
  }

  // Clear daughter vectors.
  mDauDcaVec1.clear();
  mDauDcaVec2.clear();
  mDauVec1.clear();
  mDauVec2.clear();
  return kStOK ; 
}
/*
Int_t StV0Maker::InitRun(int runnumber) {

  //========== Set Beam Line =====================
  double x0 = 0.;
  double y0 = 0.;
  double dxdz = 0.;
  double dydz = 0.;

  //Get Current Beam Line Constraint from database
  TDataSet* dbDataSet = this->GetDataBase("Calibrations/rhic");

  if (dbDataSet) {
    vertexSeed_st* vSeed = ((St_vertexSeed*) (dbDataSet->FindObject("vertexSeed")))->GetTable();

    x0 = vSeed->x0;
    y0 = vSeed->y0;
    dxdz = vSeed->dxdz;
    dydz = vSeed->dydz;
  }
  else {
    LOG_INFO << "StV0Maker -- No Database for beamline" << endm;
  }

  LOG_INFO << "BeamLine Constraint (StV0Maker): " << endm;
  LOG_INFO << "x(z) = " << x0 << " + " << dxdz << " * z" << endm;
  LOG_INFO << "y(z) = " << y0 << " + " << dydz << " * z" << endm;


  //set by hand
//  x0=0.;
//  y0=0.;
  StThreeVectorD origin(x0,y0,0.0);
  double pt = 88889999;
  double nxy=::sqrt(dxdz*dxdz +  dydz*dydz);
  if(nxy<1.e-5){ // beam line _MUST_ be tilted
    LOG_WARN << "StTofrNtupleMaker:: Beam line must be tilted!" << endm;
    nxy=dxdz=1.e-5;
  }
  double p0=pt/nxy;
  double px   = p0*dxdz;
  double py   = p0*dydz;
  double pz   = p0; // approximation: nx,ny<<0
  StThreeVectorD MomFstPt(px*GeV, py*GeV, pz*GeV);
  //delete mBeamHelix;
  mBeamHelix = new StPhysicalHelixD(MomFstPt,origin,0.5*tesla,1.);

  return kStOK;
}

Int_t StV0Maker::FinishRun(int runnumber)
{
   if(mBeamHelix) delete mBeamHelix;
   return kStOK;
}
*/
Int_t StV0Maker::Make( )
{ // Do each event
   
  if(GetDebug()) LOG_QA<<"in StV0Maker::Make"<<endm;
  // Do some cleaning here, used for StXiMaker or other subsequent makers
  mV0Dst.nv0 = 0;
  mPassEventCut = false;

  // Get 'event' data 
  StMuEvent* muEvent      =  mMuDstMaker->muDst()->event() ;

  // further selection on events. find qualified events!!
  // triggerid, skip events that do not comply with the required trigger id.
  if (!muEvent) return kStOK;
  if ( !muEvent->triggerIdCollection().nominal().isTrigger(cutTriggerIdEq) ) return kStOK ; 
  //if ( !muEvent->triggerIdCollection().nominal().isTrigger(2001) && !muEvent->triggerIdCollection().nominal().isTrigger(2003)) return kStOK ; //run3 dAu
  mV0Dst.trgmode=0; //dummy for run 9
  /*
  if( muEvent->triggerIdCollection().nominal().isTrigger(210020) && !muEvent->triggerIdCollection().nominal().isTrigger(210000)) mV0Dst.trgmode = 0; //vpd-zdce
  if( !muEvent->triggerIdCollection().nominal().isTrigger(210020) && muEvent->triggerIdCollection().nominal().isTrigger(210000)) mV0Dst.trgmode = 1; //zdce
  if( muEvent->triggerIdCollection().nominal().isTrigger(210020) && muEvent->triggerIdCollection().nominal().isTrigger(210000)) mV0Dst.trgmode = 2;  //mix
  */

  hNPrimVertex -> Fill( mMuDstMaker->muDst()->numberOfPrimaryVertices());

  mV0Dst.runnumber = muEvent->runNumber();
  mV0Dst.evtnumber = muEvent->eventNumber();

  
  // Cut on the number of vertices in the event.  On old tapes, no-vertex gets reported as VtxPosition=(0,0,0).
  // Skip events that do not have a primary vertex. the small value '1e-5' is ok for hard coding.
  if ( fabs(muEvent->primaryVertexPosition().x()) < 1e-5 && fabs(muEvent->primaryVertexPosition().y()) < 1e-5 && fabs(muEvent->primaryVertexPosition().z()) < 1e-5 )  return kStOK ;  
  //if(GetDebug()) LOG_QA<<"in StV0Maker::Make : has pv"<<endm;
  

  // possible duplicate events.
  if ( mMuDstMaker->muDst()->numberOfPrimaryTracks() == mTestNTrack && mEventsProcessed !=0 && mTestVZ !=0 &&  muEvent->primaryVertexPosition().z() == mTestVZ ) {
     LOG_WARN << mEventsProcessed <<" "<<"seems a duplicated event!"<<endm;
     return kStOK ;
  }
  mTestVZ = muEvent->primaryVertexPosition().z();
  mTestNTrack = mMuDstMaker->muDst()->numberOfPrimaryTracks();
/*
  // find primary vertex from VPD vertex and TOF matched global tracks beam projection.
  // if there is no VpdVz, skip events; if there are no TOF primary tracks, skip events.
  // condition for TOF primary tracks, |DCA_xy| < 1cm; |Zproj - vpzvz|<6 cm;
  // the xy coordinates of pv will be beam line, z is from dca piont of closest primary tracks.
  StBTofHeader *tofHeader = mMuDstMaker->muDst()->btofHeader();
  //if(tofHeader)cout<<tofHeader->vpdVz()<<endl;
  if(!tofHeader) return kStOK;
  double vpdVz = tofHeader->vpdVz();
  if(vpdVz < -990.) return kStOK;

  double dca_2d_smallvalue = 100.0;
  StThreeVectorD tofpv(0,0,0);
  TObjArray* tracks = mMuDstMaker->muDst()->globalTracks() ;    // Create a TObject array containing the global tracks
  TObjArrayIter GetTracks(tracks) ;  */                            // Create an iterator to step through the tracks
  StMuTrack* track ;                                             // Pointer to a track
/*  while ( ( track = (StMuTrack*)GetTracks.Next() ) )             // Main loop for Iterating over tracks
  {
     if(track->btofPidTraits().matchFlag())
     {
	  StThreeVectorD tofPos =  track->helix().at(track->helix().pathLengths(*mBeamHelix).first);
	  StThreeVectorD beamPos = mBeamHelix->at(track->helix().pathLengths(*mBeamHelix).second);
	  StThreeVectorD dcatof = tofPos - beamPos;
	  if(Debug()) {
	     LOG_INFO<<" tofPos(x,y,z) = "<<tofPos.x()<<","<<tofPos.y()<<","<<tofPos.z()<<endm;
	     LOG_INFO<<" beamPos(x,y,z) = "<<beamPos.x()<<","<<beamPos.y()<<","<<beamPos.z()<<endm;
	     LOG_INFO<<"  dca  (x,y,z) = "<<dcatof.x()<<","<<dcatof.y()<<","<<dcatof.z()<<endm;
	     LOG_INFO<<" 2D dca        = "<<sqrt(pow(dcatof.x(),2)+pow(dcatof.y(),2))<<endm;
	     LOG_INFO<<" 2D signed dca = "<<track->helix().geometricSignedDistance(beamPos.x(),beamPos.y())<<endm;
	  }
	  double dca_2d = sqrt(pow(dcatof.x(),2)+pow(dcatof.y(),2));
	  if(dca_2d < dca_2d_smallvalue) 
	  {
	     dca_2d_smallvalue = dca_2d;
	     tofpv.setX(beamPos.x());
	     tofpv.setY(beamPos.y());
	     tofpv.setZ(tofPos.z());
	  }
     }
  }
  if(dca_2d_smallvalue>1.) return kStOK; //no tof primary tracks. therefore no p.v. from tof.
*/
 
  // Fill some QA plots
  hVertexZ -> Fill(muEvent->primaryVertexPosition().z()) ; // Make histogram of the vertex Z distribution
  hNRefMult -> Fill( muEvent->refMult() );		

  // cut on vertexZ
  if (fabs(muEvent->primaryVertexPosition().z()) > cutAbsVertexZLeEq ) return kStOK ;
  mV0Dst.primvertexX = muEvent->primaryVertexPosition().x();
  mV0Dst.primvertexY = muEvent->primaryVertexPosition().y();
  mV0Dst.primvertexZ = muEvent->primaryVertexPosition().z();

  // cut on centrality or reference multiplicity.
  if ( muEvent->refMult() ) {}   //TODO: need to check whether this is the same as in old code. the old code might ignore the case of pile-up.
  mV0Dst.nrefmult = muEvent->refMult();
  mV0Dst.nrefmultTOF = muEvent->btofTrayMultiplicity();
  mV0Dst.bbcadcsumeast = muEvent->bbcTriggerDetector().adcSumEast();
  mV0Dst.zdcadc0 = muEvent->zdcTriggerDetector().adc(0);
  mV0Dst.bbccirate = muEvent->runInfo().bbcCoincidenceRate();

  mPassEventCut = true;
  // Do 'event' analysis based on event data 

  // Record some information...
  hSelectNRefMult -> Fill( muEvent->refMult() ); //this is an ESSENTIAL histogram to record the total number of events for certain centrality. always make sure it is filled AFTER event selection!
  if(mV0Dst.trgmode==1)hSelectNRefMultM1 -> Fill( muEvent->refMult() ); //this is an ESSENTIAL histogram to record the total number of events for certain centrality. always make sure it is filled AFTER event selection!
  if(mV0Dst.trgmode==2)hSelectNRefMultM2 -> Fill( muEvent->refMult() ); //this is an ESSENTIAL histogram to record the total number of events for certain centrality. always make sure it is filled AFTER event selection!
  hSelectNRefMultFtpcE -> Fill( muEvent->refMultFtpcEast() );
  hSelectBbcEAdcSum -> Fill( muEvent->bbcTriggerDetector().adcSumEast() );
  
  Double_t magn = muEvent->runInfo().magneticField();	
  //Double_t magn = muEvent->magneticField();	//checked! the same as above
  mV0Dst.magn = magn;
  
  // Get 'track' data, make cuts on tracks, do physics analysis, histogram results.
  TObjArray* tracks = mMuDstMaker->muDst()->globalTracks() ;    // Create a TObject array containing the global tracks
  TObjArrayIter GetTracks(tracks) ;                              // Create an iterator to step through the tracks

  //StMuTrack* track ;                                             // Pointer to a track
  //GetTracks.Reset();
  while ( ( track = (StMuTrack*)GetTracks.Next() ) )             // Main loop for Iterating over tracks
  {
     hPtRaw -> Fill( track->pt() ) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.
     hEtaRaw -> Fill( track->eta() ) ;		  //at dca to PV
     hPhiRaw -> Fill( track->phi() ) ;		  //at dca to PV
     short flag = track->flag();
     unsigned short nHits = track->nHits();	//total # of hits in all available detectors
     unsigned short nHitsFit = track->nHitsFit(kTpcId);
     short charge = track->charge();
     //StThreeVectorF p = track->p();	  //at dca to PV
     //StThreeVectorF origin = track->firstPoint();  //? NOT SURE whether firstPoint from detectorinfo is the same as helix.orgin()!!!
     double nsigmapion = track->nSigmaPion();
     double nsigmaproton = track->nSigmaProton();
     double nsigmakaon = track->nSigmaKaon();
     double nsigmaelectron = track->nSigmaElectron();
     double dedx = track->dEdx();

     StPhysicalHelixD helix = track->helix(); //inner helix. good for dca to PV.
     StThreeVectorD p = helix.momentum(magn*kilogauss);    //momentum at origin
     StThreeVectorD origin = helix.origin();  //origin of helix
     double pt = p.perp();

     //some checks.
     hNHitsFit -> Fill( nHitsFit ) ;
     hNHits -> Fill( nHits ) ;
     hPtDiff -> Fill( track->pt() - p.perp() ) ; 
     hOrDiff -> Fill( (track->firstPoint() - origin).mag() ) ; 
     hPDiff -> Fill( fabs(track->p().mag() - p.mag()) ) ; 
     //comments: there are difference between the values above. But they seem to be acceptably small!

     int hrot;		//helicity of helix, sign of -charge*magn
     if (-charge*magn > 0) hrot = 1;
     else hrot = -1;

     double nsigma;
     if(charge == mCharge1) {
	  if(mV0Type == kLambda || mV0Type == kAntiLambda) nsigma = nsigmaproton;
	  else if(mV0Type == kKs) nsigma = nsigmapion;
	  else nsigma = nsigmaelectron;
     }
     else {
	  if(mV0Type != kPhoton)  nsigma = nsigmapion;
	  else nsigma = nsigmaelectron;
     }

     //if(track->vertexIndex()!=StMuDst::currentVertexIndex())continue;
     //if you want to use track->dca(), turn this on. if it is not turned on, that function crashes.
     //OK. let's cut tracks
     //if(nHitsFit<=cutNHitsGr)continue;
     if(nHits<=cutNHitsGr)continue;
     if(flag <=0  )continue; //or <=0 ?
     //if(track->bad() )continue; 
     if(abs(charge)!=1) continue;
    // if(!(track->btofPidTraits().matchFlag())&&(track->pt())<1.5) continue;//make sure we have the TOF information of track for pt<1.5
     hPt -> Fill( track->pt() ) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.
     hEta -> Fill( track->eta() ) ;		  //at dca to PV
     hPhi -> Fill( track->phi() ) ;		  //at dca to PV
     if(pt<0.5)hPhiLowPt->Fill(track->phi());
     else hPhiHighPt->Fill(track->phi());
     hNSigmaPion->Fill(nsigmapion);
     hNSigmaProton->Fill(nsigmaproton);
     hNSigmaKaon->Fill(nsigmakaon);
     hDedxP->Fill(p.mag()*charge,dedx);

     if(pt<cutPtGrEq)continue; //should be larger. like 0.15 or 0.2

     if(charge == mCharge1 && fabs(nsigma)<cutAbsNSigma1Le){
	  //record the first daughter 
	  //fill the vector
	  StPhysicalHelixD helix = track->helix(); //inner helix. good for dca to PV.
	  double pathlength = helix.pathLength(muEvent->primaryVertexPosition(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
	  StThreeVectorF dca = helix.at(pathlength)-muEvent->primaryVertexPosition();
	  
	  //some tests on dca functions.
	  //StThreeVectorF dca1 = track->dca(muEvent->primaryVertexPosition()); //it simply crash!
	  //StThreeVectorF dca11 = track->dcaGlobal(); //it doesn't crash, but give some zero values
	  //hDcaDiff -> Fill((dca1-dca11).mag());
	  //LOG_QA<<"KK "<<(dca1-dca11).mag()<<" "<<dca1.mag()<<endm;
	  
	  //if(getDcaToPV(track, muEvent->primaryVertexPosition()) - dca.mag() > 0.4)
	  //LOG_QA << getDcaToPV(track, muEvent->primaryVertexPosition()) << " "<< dca.mag()<<endm;
	  //if(getDcaToPV(track, muEvent->primaryVertexPosition()) < 2)
	  //hDcaDiff -> Fill(getDcaToPV(track, muEvent->primaryVertexPosition())-dca.mag());
	  
	  //double dca = getDcaToPV(track, muEvent->primaryVertexPosition());
	  
	  //CONCLUSION: StHelix method seems more strict, always use it instead of getDcaToPV!!!
//	  if(!(track->btofPidTraits().matchFlag())&&(track->pt())<1.5) continue;
          hmDau1Dca->Fill((track->dcaGlobal()).mag());
          hmDau1Dca2->Fill(dca.mag());

          if(dca.mag()<cutDca1GrEq)continue;
	  //if((mV0Type==kLambda||mV0Type==kAntiLambda) && dca.mag()>cutDca1LeEq)continue;
	  mDauDcaVec1.push_back(dca.mag());
	  mDauVec1.push_back(track);
     }

     if(charge == mCharge2 && fabs(nsigma)<cutAbsNSigma2Le){
	  //record the second daughter
	  StPhysicalHelixD helix = track->helix(); //inner helix. good for dca to PV.

	  StThreeVectorF pv = muEvent->primaryVertexPosition();
	  //rotate transverse coordinates and momenta for background estimation
	  if(mRotate){
	     StThreeVectorD p1 = helix.momentum(magn*kilogauss);    //momentum at origin
	     StThreeVectorD x1 = helix.origin();    //origin
	     p1.setX(-p1.x());
	     p1.setY(-p1.y());
	     x1.setX(-(x1.x()-pv.x())+pv.x());
	     x1.setY(-(x1.y()-pv.y())+pv.y());
	     StPhysicalHelixD helixtmp(p1, x1, magn*kilogauss, track->charge());
	     helix = helixtmp;
	  }
	  
	  double pathlength = helix.pathLength(muEvent->primaryVertexPosition(), false); // do scan periods. NOTE: the default is false.
	  StThreeVectorF dca = helix.at(pathlength)-muEvent->primaryVertexPosition();

	  //LOG_QA << getDcaToPV(track, muEvent->primaryVertexPosition()) << " "<< dca.mag()<<endm;
	  //if(getDcaToPV(track, muEvent->primaryVertexPosition()) < 2)
	  //hDcaDiff -> Fill(getDcaToPV(track, muEvent->primaryVertexPosition())-dca.mag());

	  //double dca = getDcaToPV(track, muEvent->primaryVertexPosition());
          hmDau2Dca->Fill((track->dcaGlobal()).mag());
          hmDau2Dca2->Fill(dca.mag());

	  if(dca.mag()<cutDca2GrEq)continue;
	  //if(dca.mag()>cutDca2LeEq)continue;
	  mDauDcaVec2.push_back(dca.mag());
	  mDauVec2.push_back(track);
     }
  }

  assert(mDauVec1.size() == mDauDcaVec1.size() && mDauVec2.size() == mDauDcaVec2.size());
  
  //reconstruct V0
  int nV0 = 0;
  //cout<< "daughter vector 1 size = "<<mDauVec1.size()<<endl;
  //cout<< "daughter vector 2 size = "<<mDauVec2.size()<<endl;
  for(unsigned int i=0;i<mDauVec1.size();i++){
    // cout<< "daughter vector1 number = "<<i<<endl;
     //get proton track info here
     //cut them before in track selection
     StMuTrack * track1 = mDauVec1[i];
     StPhysicalHelixD helix1 = track1->helix(); //inner helix. good for dca to PV.
     StThreeVectorD p1 = helix1.momentum(magn*kilogauss);    //momentum at origin
     //StThreeVectorD origin1 = helix1.origin();  //origin of helix
     double pt1 = p1.perp();
     double dca1 = mDauDcaVec1[i];

     //record TOF information...
     int tofflag1 = track1->btofPidTraits().matchFlag(); 
     double tof1 = track1->btofPidTraits().timeOfFlight();
     StThreeVectorF tofpos1 = track1->btofPidTraits().position();

     for(unsigned int j=0;j<mDauVec2.size();j++){
      //    cout<< "daughter vector2 number = "<<j<<endl;
	  //get pion track info here
	  StMuTrack * track2 = mDauVec2[j];
	  if(track2->id() == track1->id())continue;  //for same sign
	  if((mSameSignPlus || mSameSignMinus) && j<=i)continue; //avoid double counting in s.s. REQURING vec1 and vec2 are the same!
	  StPhysicalHelixD helix2 = track2->helix(); //inner helix. good for dca to PV.
	  
	  StThreeVectorF pv = muEvent->primaryVertexPosition();
	  if(mRotate){
	     StThreeVectorD tp1 = helix2.momentum(magn*kilogauss);    //momentum at origin
	     StThreeVectorD tx1 = helix2.origin();    //origin
	     tp1.setX(-tp1.x());
	     tp1.setY(-tp1.y());
	     tx1.setX(-(tx1.x()-pv.x())+pv.x());
	     tx1.setY(-(tx1.y()-pv.y())+pv.y());
	     StPhysicalHelixD helixtmp(tp1, tx1, magn*kilogauss, track2->charge());
	     helix2 = helixtmp;
	  }
	  
	  StThreeVectorD p2 = helix2.momentum(magn*kilogauss);    //momentum at origin
	  double pt2 = p2.perp();
	  double dca2 = mDauDcaVec2[j];
	  //cut them before in track selection
	  
	  StThreeVectorF xv0, op1, op2;
	  double dca1to2;
	  //double stdca1to2;
	  StThreeVectorF ox1,ox2;
	  if(!mDcaAlgoLong){
	     pair<double,double> tmps = helix1.pathLengths(helix2); 
	     //StThreeVectorD ox1 = helix1.at(tmps.first);
	     //StThreeVectorD ox2 = helix2.at(tmps.second);
	     ox1 = helix1.at(tmps.first);
	     ox2 = helix2.at(tmps.second);
	     //stdca1to2 = (ox1-ox2).mag();
	     dca1to2 = (ox1-ox2).mag();
	     xv0 = (ox1 + ox2)/2.;
	     op1 = helix1.momentumAt(tmps.first, magn*kilogauss);
	     op2 = helix2.momentumAt(tmps.second, magn*kilogauss);
	  }

	  //kStHelixDca = true;
	  //kMinimize = true;
	  //StHelix method above is VERY SLOW. use it only for checking the consistency of 
	  //long's code
	  if(mDcaAlgoLong) dca1to2 = closestDistance(helix1, helix2, magn, pv, xv0, op1, op2);
	  //cut on dca1to2
	  //if(stdca1to2 <= cutDca1to2LeEq)
	  //   LOG_QA << stdca1to2<<" "<<dca1to2<<endm;
	  //LOG_QA << dca1to2<<endm;
	  if(dca1to2 > cutDca1to2LeEq) continue;

	  double oe1 = sqrt(op1.mag2() + mMass1*mMass1);
	  double oe2 = sqrt(op2.mag2() + mMass2*mMass2);
	  double v0mass = sqrt(mMass1*mMass1 + mMass2*mMass2 + 2.*oe1*oe2 - 2.*op1.dot(op2));
	  //cut on v0mass
	  if(fabs(v0mass-mMassV0) > cutV0MassWidthLeEq)continue;

	  StThreeVectorD pv0 = op1 + op2;
	  StThreeVectorD xv0toPV = xv0 - pv;

	  //helix of v0: straight line
	  StPhysicalHelixD helixv0(pv0,xv0,0,0);

	  //pthead, ptarm cut
	  double pthead1 = op1.dot(pv0)/pv0.mag();
	  double pthead2 = op2.dot(pv0)/pv0.mag();
	  double ptarm   = sqrt(op1.mag2()-pthead1*pthead1);
	  if(mV0Type!=kKs && mV0Type!=kPhoton ){
	     //if(ptarm> cutDauPtArmLeEq || fabs((pthead1-pthead2)/(pthead1+pthead2))> cutAbsDausPtShoulderDiffLeEq )continue;
	  }

	  //forward decay cut
	  double ang1 = op1.x()*xv0toPV.x() + op1.y()*xv0toPV.y();
	  double ang2 = op2.x()*xv0toPV.x() + op2.y()*xv0toPV.y();
	  if(mV0Type!=kKs && mV0Type!=kPhoton ){
	     //if(ang1<=cutDau1DecAngGr || ang2<=cutDau2DecAngGr)continue;
	  }
	  
	  //r dot p for v0. cut on it. should be larger than 0. 
	  double rdotp = xv0toPV.dot(pv0) ;
	  if(mV0Type!=kPhoton && rdotp<=cutV0rdotpGr)continue;
	  
	  //calculate v0 to PV dca. v0 carry no charge. straight line. cut on dca
	  double dcav0toPV = rdotp*rdotp/pv0.mag2();
	  dcav0toPV = sqrt( xv0toPV.mag2() - dcav0toPV);
	  if(dcav0toPV>=cutDcaV0Le)continue;

	  //cut on decay length
	  double v0decaylength = xv0toPV.mag();
	  if(v0decaylength < cutV0DecLenGrEq )continue;

	  if(mV0Type==kPhoton && acos(op1.dot(op2)/op1.mag()/op2.mag()) > cutDau1Dau2Ang3DLe)continue;

	  if(mV0Type==kPhoton && fabs(helix1.dipAngle()-helix2.dipAngle()) > cutDau1Dau2DipAngDiffLe)continue;

	  //cut on sinth, or theta
	  double sinth = (xv0toPV.cross(pv0)).mag()/xv0toPV.mag()/pv0.mag();
	  double theta = atan2(sinth, rdotp/xv0toPV.mag()/pv0.mag()); //theta range: from 0 to pi

	  //record TOF information
	  int tofflag2 = track2->btofPidTraits().matchFlag(); 
	  double tof2 = track2->btofPidTraits().timeOfFlight();
	  StThreeVectorF tofpos2 = track2->btofPidTraits().position();

	  pair<double,double> pthpair = helix2.pathLength(215.0);
	 // cout<<tofflag2<<" "<<tof2<<" "<<tofpos2<<" "<<tofpos2.perp()<<" "<<op1.pseudoRapidity()<<endl;
	 // cout<<pthpair.first<<" "<<pthpair.second<<endl;
	 // cout<<helix2.at(pthpair.second)<<endl;

	  double tofpathlen1 = -999.;
	  double tofpathlen2 = -999.;
	  if(tofflag1>0) tofpathlen1 = helix1.pathLength(tofpos1) - helix1.pathLength(xv0); 
	  if(tofflag2>0) tofpathlen2 = helix2.pathLength(tofpos2) - helix2.pathLength(xv0); 

	  double v0pathlen = helixv0.pathLength(xv0)-helixv0.pathLength(pv);

	  mV0Dst.v0mass[nV0] = v0mass;
	  mV0Dst.v0pt[nV0]   = pv0.perp();
	  mV0Dst.v0rapidity[nV0]    = log( (sqrt(v0mass*v0mass+pv0.mag2()) + pv0.z())/sqrt(v0mass*v0mass+pv0.perp2()));
	  mV0Dst.v0eta[nV0]    = 0.5*log( (pv0.mag() + pv0.z())/(pv0.mag() - pv0.z()) );
	  mV0Dst.v0x[nV0]	     = xv0.x();
	  mV0Dst.v0y[nV0]	     = xv0.y();
	  mV0Dst.v0z[nV0]	     = xv0.z();
	  mV0Dst.v0px[nV0]     = pv0.x();
	  mV0Dst.v0py[nV0]     = pv0.y();
	  mV0Dst.v0pz[nV0]     = pv0.z();
	  mV0Dst.v0declen[nV0] = v0decaylength;
	  mV0Dst.v0dca[nV0]    = dcav0toPV;
	  mV0Dst.v0dca2d[nV0]  = helixv0.geometricSignedDistance(pv.x(),pv.y());
	  mV0Dst.v0pathlen[nV0]= v0pathlen; 
	  mV0Dst.dau1id[nV0]   = track1->id();
	  mV0Dst.dau2id[nV0]   = track2->id();
	  mV0Dst.dau1dca[nV0]  = dca1;
	  mV0Dst.dau1dca2d[nV0]  = helix1.geometricSignedDistance(pv.x(),pv.y());
	  mV0Dst.dau1nhits[nV0]= track1->nHits();
	  mV0Dst.dau1dedx[nV0] = track1->dEdx();
	  mV0Dst.dau1nsigma[nV0] = (mV0Type!=kKs)? ((mV0Type!=kPhoton)?track1->nSigmaProton():track1->nSigmaElectron()) : track1->nSigmaPion();
	  mV0Dst.dau1eta[nV0]   = op1.pseudoRapidity();
	  mV0Dst.dau1pt[nV0]   = op1.perp();
	  mV0Dst.dau1px[nV0]   = op1.x();
	  mV0Dst.dau1py[nV0]   = op1.y();
	  mV0Dst.dau1pz[nV0]   = op1.z();
	  mV0Dst.dau1tpc[nV0]  = track1->nHitsFit(kTpcId);
	  mV0Dst.dau1ssd[nV0]  = track1->nHitsFit(kSsdId);
	  mV0Dst.dau1svt[nV0]  = track1->nHitsFit(kSvtId);
	  mV0Dst.dau1tofflag[nV0] = tofflag1;
	  mV0Dst.dau1tof[nV0] = tof1;
	  mV0Dst.dau1pathlen[nV0] = tofpathlen1;

	  mV0Dst.dau2dca[nV0]  = dca2;
	  mV0Dst.dau2dca2d[nV0]  = helix2.geometricSignedDistance(pv.x(),pv.y());
	  mV0Dst.dau2nhits[nV0]= track2->nHits();
	  mV0Dst.dau2dedx[nV0] = track2->dEdx();
	  mV0Dst.dau2nsigma[nV0] = (mV0Type!=kPhoton)? track2->nSigmaPion():track1->nSigmaElectron();
	  mV0Dst.dau2eta[nV0]   = op2.pseudoRapidity();
	  mV0Dst.dau2pt[nV0]   = op2.perp();
	  mV0Dst.dau2px[nV0]   = op2.x();
	  mV0Dst.dau2py[nV0]   = op2.y();
	  mV0Dst.dau2pz[nV0]   = op2.z();
	  mV0Dst.dau2tpc[nV0]  = track2->nHitsFit(kTpcId);
	  mV0Dst.dau2ssd[nV0]  = track2->nHitsFit(kSsdId);
	  mV0Dst.dau2svt[nV0]  = track2->nHitsFit(kSvtId);
	  mV0Dst.dau2tofflag[nV0] = tofflag2;
	  mV0Dst.dau2tof[nV0] = tof2;
	  mV0Dst.dau2pathlen[nV0] = tofpathlen2;

	  mV0Dst.dca1to2[nV0]  = dca1to2;

	  nV0 ++;

	  hInvMass->Fill(v0mass);
	  //mV0Vec.push_back(v0data);

     }
  }

  mV0Dst.nv0 = nV0;

  if(nV0 ==0 && mDumpNull && mV0Tree) mV0Tree->Fill();
  if(nV0 > 0 && mV0Tree) mV0Tree->Fill();
  //dump v0 vector into a TTree
  
  mEventsProcessed++ ;
  //LOG_QA << mDauVec1.size() <<" "<< mDauVec2.size() <<" "<< mV0Vec.size()<<endm;
  mDauDcaVec1.clear();
  mDauDcaVec2.clear();
  mDauVec1.clear();
  mDauVec2.clear();
  //mV0Vec.clear();

  return kStOK ;
  
}

Int_t StV0Maker::Finish( )
{ // Do once at the end the analysis

  // Write histograms to disk, output miscellaneous other information
  if(histogram_output!=NULL) histogram_output -> Write() ;   // Write all histograms to disk 
  if(v0tree_output!=NULL) v0tree_output -> Write() ;   // Write all histograms to disk 

  cout << "Total Events Processed in StV0Maker " << mEventsProcessed << endl ;

  return kStOk ;  

}


