#include "StXiMaker.h"
#include "StarClassLibrary/SystemOfUnits.h"

#include <iostream>

#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StV0Maker/StV0Maker.h"
#include "StV0Maker/StDcaService.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"

ClassImp(StXiMaker)                   // Macro for CINT compatibility

StXiMaker::StXiMaker( StMuDstMaker* maker, StV0Maker* v0maker, const char * name) : StMaker(name)
{ // Initialize and/or zero all public/private data members here.

  //for ( Int_t i = 0 ; i < kMaxNumberOfTH1F ; i++ )  // Zero the histogram pointers, not necessary. it is NULL naturaly.
  //  {
  //    histogram[i] = NULL ;
  //  }

  mMuDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
  mV0Maker      	 = v0maker ;                  // Pass V0Maker pointer to DstAnlysisMaker Class member functions

  mXiType = kXi;	//Lambda as default!

  mRotate = false;
  mDcaAlgoLong = true; //use LONG's dca method as default

  histogram_output = NULL  ;                    // Zero the Pointer to histogram output file
  xitree_output = NULL  ;                    // Zero the Pointer to v0 tree output file
  mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro
  mXiTreeOutputFileName = "" ;               // Xi Output File Name will be set inside the "analysis".C macro

  mXiTree = NULL ;

  mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker 
  mTestNTrack = 0;
  mTestVZ = 0;

}

StXiMaker::~StXiMaker() 
{ // Destroy and/or zero out all public/private data members here.
}

void StXiMaker::initConst(){
   //initialize the constant for different Xi types.
   
   if(mXiType == kXi || mXiType == kAntiXi){
	// for Xi or AntiXi
	mMass1      = 0.93827; // mass of proton or pbar
	mMass2      = 0.13957; // mass of pion- or pion+
	mMassV0     = 1.115684;// mass of Lambda or AntiLambda
	mMassBachelor = 0.13957; //mass of pion- or pion+
	mMassXi	= 1.32131;   //mass of Xi or AntiXi

	if(mXiType == kXi) { 
	   mCharge1	= 1;
	   mCharge2	= -1;
	   mChargeBachelor = -1;
	   mChargeXi = -1;
	}
	else {
	   mCharge1	= -1;
	   mCharge2	= 1;
	   mChargeBachelor = 1;
	   mChargeXi = 1;
	}
	//do not setup the cut values here. those are parameters
	
	//parameters for StDcaService.cxx
	kShiftConnect = 0.3;
	kShiftContain = 0.3;
   }
   else {
	assert(mXiType == kOmega || mXiType == kAntiOmega); 
	// for Omega or AntiOmega 
	mMass1      = 0.93827; // mass of proton or pbar
	mMass2      = 0.13957; // mass of pion- or pion+
	mMassV0     = 1.115684;// mass of Lambda or AntiLambda
	mMassBachelor = 0.493677; //mass of K- or K+
	mMassXi	= 1.67245;   //mass of Omega or AntiOmega
	
	if(mXiType == kOmega) { 
	   mCharge1	= 1;
	   mCharge2	= -1;
	   mChargeBachelor = -1;
	   mChargeXi = -1;
	}
	else {
	   mCharge1	= -1;
	   mCharge2	= 1;
	   mChargeBachelor = 1;
	   mChargeXi = 1;
	}
	//do not setup the cut values here. those are parameters
	
	//parameters for StDcaService.cxx
	kShiftConnect = 0.3;
	kShiftContain = 0.3;
   }

   return;
}

void StXiMaker::initParam(){
   //setup the cut values here. do not hard-code them in ::Make()

   cutNHitsGr = 15;
   cutPtGrEq = 0.15;

   cutAbsNSigma1Le = 4.;
   cutDca1GrEq  = 0.5; //1.0

   cutV0DcaGrEq = 0; //0.3
   cutV0MassWidthLeEq = 0.02;
   
   cutDca1to2LeEq = 1.0; //0.8
   cutXiMassWidthLeEq = 0.05;
   //cutDauPtArmLeEq = 0.3;
   //cutAbsDausPtShoulderDiffLeEq = 1.2;
   //cutDau1DecAngGr = 0.;
   cutXirdotpGr  = 0;
   cutXircrosspLeEq = 0.5; //0.2
   //cutDcaXiLe    = 3.3;
   cutXiDecLenGrEq = 2.0; //2.0

   return;
}

void StXiMaker::initHisto()
{
   // Create Histograms
   // there is no better way to set QA histograms. do not use histogram arrays or vectors.
   // it is not useful. there are no need to operate all the histograms at the same time.

   const Int_t    nbins    =  100   ;

   hVertexZ  = new TH1F( "Vertex", "Event Vertex Z Position", nbins, -25.0, 25.0 ) ; 
   hPt  = new TH1F( "Pt", "Transverse Momentum for all particles", nbins, 0.0, 10.0 ) ;
   hNPrimVertex  = new TH1F( "PrimVertex", "Number of Primary Vertex", 10, 0.0, 10.0 ) ;
   hNRefMult  = new TH1F( "RefMult", "Reference Multiplicity", 200, 0.0, 1000.0 ) ;
   hPtDiff  = new TH1F( "ptdiff", "Reference Multiplicity", 100, 0.0, 0.02 ) ;
   hOrDiff  = new TH1F( "ordiff", "Reference Multiplicity", 100, 0.0, 0.2 ) ;
   hPDiff  = new TH1F( "pdiff", "Reference Multiplicity", 100, 0.0, 0.4 ) ;
   hDcaDiff  = new TH1F( "dcadiff", "Reference Multiplicity", 100, -0.4, 0.4 ) ;
   hInvMass  = new TH1F( "invmass", "Xi Inv. Mass", 100, mMassXi-cutXiMassWidthLeEq, mMassXi+cutXiMassWidthLeEq ) ;
   hInvMassSTD  = new TH1F( "invmassSTD", "STD Xi Inv. Mass", 100, mMassXi-cutXiMassWidthLeEq, mMassXi+cutXiMassWidthLeEq ) ;

   return;
}

void StXiMaker::initTree()
{
   //initialize the TTree for StXiDst
   mXiTree = new TTree("XiPicoDst","XiPicoDst from StXiMaker");

   mXiTree->SetDirectory(xitree_output);

   mXiTree->Branch("runnumber",&mXiDst.runnumber,"runnumber/I");
   mXiTree->Branch("evtnumber",&mXiDst.evtnumber,"evtnumber/I");
   mXiTree->Branch("trgmode",&mXiDst.trgmode,"trgmode/I");
   mXiTree->Branch("nrefmult",&mXiDst.nrefmult,"nrefmult/I");
   mXiTree->Branch("nrefmultTOF",&mXiDst.nrefmultTOF,"nrefmultTOF/I");
   mXiTree->Branch("bbcadcsumeast",&mXiDst.bbcadcsumeast,"bbcadcsumeast/I");
   mXiTree->Branch("zdcadc0",&mXiDst.zdcadc0,"zdcadc0/F");
   mXiTree->Branch("bbccirate",&mXiDst.bbccirate,"bbccirate/D");
   mXiTree->Branch("primvertexX",&mXiDst.primvertexX,"primvertexX/F");
   mXiTree->Branch("primvertexY",&mXiDst.primvertexY,"primvertexY/F");
   mXiTree->Branch("primvertexZ",&mXiDst.primvertexZ,"primvertexZ/F");   
   mXiTree->Branch("magn",&mXiDst.magn,"magn/F");
   mXiTree->Branch("nxi",&mXiDst.nxi,"nxi/I");

   mXiTree->Branch("v0mass",mXiDst.v0mass,"v0mass[nxi]/F");
   mXiTree->Branch("v0pt",mXiDst.v0pt,"v0pt[nxi]/F");
   mXiTree->Branch("v0rapidity",mXiDst.v0rapidity,"v0rapidity[nxi]/F");
   mXiTree->Branch("v0eta",mXiDst.v0eta,"v0eta[nxi]/F");
   mXiTree->Branch("v0x",mXiDst.v0x,"v0x[nxi]/F");
   mXiTree->Branch("v0y",mXiDst.v0y,"v0y[nxi]/F");
   mXiTree->Branch("v0z",mXiDst.v0z,"v0z[nxi]/F");
   mXiTree->Branch("v0px",mXiDst.v0px,"v0px[nxi]/F");
   mXiTree->Branch("v0py",mXiDst.v0py,"v0py[nxi]/F");
   mXiTree->Branch("v0pz",mXiDst.v0pz,"v0pz[nxi]/F");
   mXiTree->Branch("v0declen",mXiDst.v0declen,"v0declen[nxi]/F");
   mXiTree->Branch("v0dca",mXiDst.v0dca,"v0dca[nxi]/F");
   mXiTree->Branch("v0dca2d",mXiDst.v0dca2d,"v0dca2d[nxi]/F");
   mXiTree->Branch("v0pathlen",mXiDst.v0pathlen,"v0pathlen[nxi]/F");

   mXiTree->Branch("dau1id",mXiDst.dau1id,"dau1id[nxi]/I");
   mXiTree->Branch("dau1dca",mXiDst.dau1dca,"dau1dca[nxi]/F");
   mXiTree->Branch("dau1dca2d",mXiDst.dau1dca2d,"dau1dca2d[nxi]/F");
   mXiTree->Branch("dau1nhits",mXiDst.dau1nhits,"dau1nhits[nxi]/I");
   mXiTree->Branch("dau1dedx",mXiDst.dau1dedx,"dau1dedx[nxi]/F");
   mXiTree->Branch("dau1nsigma",mXiDst.dau1nsigma,"dau1nsigma[nxi]/F");
   mXiTree->Branch("dau1eta",mXiDst.dau1eta,"dau1eta[nxi]/F");
   mXiTree->Branch("dau1pt",mXiDst.dau1pt,"dau1pt[nxi]/F");
   mXiTree->Branch("dau1px",mXiDst.dau1px,"dau1px[nxi]/F");
   mXiTree->Branch("dau1py",mXiDst.dau1py,"dau1py[nxi]/F");
   mXiTree->Branch("dau1pz",mXiDst.dau1pz,"dau1pz[nxi]/F");
   mXiTree->Branch("dau1tpc",mXiDst.dau1tpc,"dau1tpc[nxi]/I");
   mXiTree->Branch("dau1ssd",mXiDst.dau1ssd,"dau1ssd[nxi]/I");
   mXiTree->Branch("dau1svt",mXiDst.dau1svt,"dau1svt[nxi]/I");
   mXiTree->Branch("dau1tofflag",mXiDst.dau1tofflag,"dau1tofflag[nxi]/I");
   mXiTree->Branch("dau1tof",mXiDst.dau1tof,"dau1tof[nxi]/F");
   mXiTree->Branch("dau1pathlen",mXiDst.dau1pathlen,"dau1pathlen[nxi]/F");
 
   mXiTree->Branch("dau2id",mXiDst.dau2id,"dau2id[nxi]/I");
   mXiTree->Branch("dau2dca",mXiDst.dau2dca,"dau2dca[nxi]/F");
   mXiTree->Branch("dau2dca2d",mXiDst.dau2dca2d,"dau2dca2d[nxi]/F");
   mXiTree->Branch("dau2nhits",mXiDst.dau2nhits,"dau2nhits[nxi]/I");
   mXiTree->Branch("dau2dedx",mXiDst.dau2dedx,"dau2dedx[nxi]/F");
   mXiTree->Branch("dau2nsigma",mXiDst.dau2nsigma,"dau2nsigma[nxi]/F");
   mXiTree->Branch("dau2eta",mXiDst.dau2eta,"dau2eta[nxi]/F");
   mXiTree->Branch("dau2pt",mXiDst.dau2pt,"dau2pt[nxi]/F");
   mXiTree->Branch("dau2px",mXiDst.dau2px,"dau2px[nxi]/F");
   mXiTree->Branch("dau2py",mXiDst.dau2py,"dau2py[nxi]/F");
   mXiTree->Branch("dau2pz",mXiDst.dau2pz,"dau2pz[nxi]/F");
   mXiTree->Branch("dau2tpc",mXiDst.dau2tpc,"dau2tpc[nxi]/I");
   mXiTree->Branch("dau2ssd",mXiDst.dau2ssd,"dau2ssd[nxi]/I");
   mXiTree->Branch("dau2svt",mXiDst.dau2svt,"dau2svt[nxi]/I");
   mXiTree->Branch("dau2tofflag",mXiDst.dau2tofflag,"dau2tofflag[nxi]/I");
   mXiTree->Branch("dau2tof",mXiDst.dau2tof,"dau2tof[nxi]/F");
   mXiTree->Branch("dau2pathlen",mXiDst.dau2pathlen,"dau2pathlen[nxi]/F");

   mXiTree->Branch("dca1to2",mXiDst.dca1to2,"dca1to2[nxi]/F");

   mXiTree->Branch("bachid",mXiDst.bachid,"bachid[nxi]/I");
   mXiTree->Branch("bachdca",mXiDst.bachdca,"bachdca[nxi]/F");
   mXiTree->Branch("bachdca2d",mXiDst.bachdca2d,"bachdca2d[nxi]/F");
   mXiTree->Branch("bachnhits",mXiDst.bachnhits,"bachnhits[nxi]/I");
   mXiTree->Branch("bachdedx",mXiDst.bachdedx,"bachdedx[nxi]/F");
   mXiTree->Branch("bachnsigma",mXiDst.bachnsigma,"bachnsigma[nxi]/F");
   mXiTree->Branch("bacheta",mXiDst.bacheta,"bacheta[nxi]/F");
   mXiTree->Branch("bachpt",mXiDst.bachpt,"bachpt[nxi]/F");
   mXiTree->Branch("bachpx",mXiDst.bachpx,"bachpx[nxi]/F");
   mXiTree->Branch("bachpy",mXiDst.bachpy,"bachpy[nxi]/F");
   mXiTree->Branch("bachpz",mXiDst.bachpz,"bachpz[nxi]/F");
   mXiTree->Branch("bachtpc",mXiDst.bachtpc,"bachtpc[nxi]/I");
   mXiTree->Branch("bachssd",mXiDst.bachssd,"bachssd[nxi]/I");
   mXiTree->Branch("bachsvt",mXiDst.bachsvt,"bachsvt[nxi]/I");
   mXiTree->Branch("bachtofflag",mXiDst.bachtofflag,"bachtofflag[nxi]/I");
   mXiTree->Branch("bachtof",mXiDst.bachtof,"bachtof[nxi]/F");
   mXiTree->Branch("bachpathlen",mXiDst.bachpathlen,"bachpathlen[nxi]/F");
   
   mXiTree->Branch("dcav0tobach",mXiDst.dcav0tobach,"dcav0tobach[nxi]/F");
   //mXiTree->Branch("stdcav0tobach",mXiDst.stdcav0tobach,"stdcav0tobach[nxi]/F");
   
   mXiTree->Branch("ximass",mXiDst.ximass,"ximass[nxi]/F");
   mXiTree->Branch("xipt",mXiDst.xipt,"xipt[nxi]/F");
   mXiTree->Branch("xirapidity",mXiDst.xirapidity,"xirapidity[nxi]/F");
   mXiTree->Branch("xieta",mXiDst.xieta,"xieta[nxi]/F");
   mXiTree->Branch("xix",mXiDst.xix,"xix[nxi]/F");
   mXiTree->Branch("xiy",mXiDst.xiy,"xiy[nxi]/F");
   mXiTree->Branch("xiz",mXiDst.xiz,"xiz[nxi]/F");
   mXiTree->Branch("xipx",mXiDst.xipx,"xipx[nxi]/F");
   mXiTree->Branch("xipy",mXiDst.xipy,"xipy[nxi]/F");
   mXiTree->Branch("xipz",mXiDst.xipz,"xipz[nxi]/F");
   mXiTree->Branch("xideclen",mXiDst.xideclen,"xideclen[nxi]/F");
   mXiTree->Branch("xidca",mXiDst.xidca,"xidca[nxi]/F");
   mXiTree->Branch("xidca2d",mXiDst.xidca2d,"xidca2d[nxi]/F");
   mXiTree->Branch("xisinth",mXiDst.xisinth,"xisinth[nxi]/F");
   mXiTree->Branch("xipathlen",mXiDst.xipathlen,"xipathlen[nxi]/F");

   return;
}

Int_t StXiMaker::Init( )
{
  // setup the constants according to mXiType
  initConst();

  // initialize parameters (cuts)
  initParam();
  
  // Create Histogram output file
  if(mHistogramOutputFileName == "") { 
     //CAUTION: ALWAYS USE { } HERE!!! LOG_XXX is a if()xxx macro!!!
     LOG_ERROR << "StXiMaker: Please specify the histrogram output file" <<endm;
     exit(-1);
  }
  else {
     histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;  
  }
  // Book histograms
  initHisto();
  
  // Create Xi Tree output file
  if(mXiTreeOutputFileName == "") {
     LOG_WARN << "StXiMaker: The Xi tree output file is not specified! output is smeared!" <<endm;
  }
  else {
     xitree_output = new TFile( mXiTreeOutputFileName, "recreate" ) ;
     // Create Xi Tree
     initTree();
  }

  // Clear daughter vectors.
  mDauDcaVec1.clear();
  mDauDcaVec2.clear();
  mDauVec1.clear();
  mDauVec2.clear();
  //mXiVec.clear();
  return kStOK ; 
}

Int_t StXiMaker::Make( )
{ // Do each event
   
  if(GetDebug()) LOG_QA<<"in StXiMaker::Make"<<endm;
  // Do some cleaning here
  mXiDst.nxi = 0;

  //skip the event if it did not passed the event cut in StV0Maker
  if(!mV0Maker->passEventCut()) return kStOK;

  // Get 'event' data 
  StMuEvent* muEvent      =  mMuDstMaker->muDst()->event() ;
  // Get 'v0' data 
  const StV0Dst& v0Dst    =  mV0Maker->getV0Dst();

  //if no v0 candidates were found in StV0Maker, do not search Xi.
  if(v0Dst.nv0==0) return kStOK;
  
  mXiDst.runnumber = v0Dst.runnumber;
  mXiDst.evtnumber = v0Dst.evtnumber;
  mXiDst.trgmode = v0Dst.trgmode;
  mXiDst.nrefmult = v0Dst.nrefmult;
  mXiDst.nrefmultTOF = v0Dst.nrefmultTOF;
  mXiDst.bbcadcsumeast = v0Dst.bbcadcsumeast;
  mXiDst.zdcadc0 = v0Dst.zdcadc0;
  mXiDst.bbccirate = v0Dst.bbccirate;
  mXiDst.primvertexX = v0Dst.primvertexX;
  mXiDst.primvertexY = v0Dst.primvertexY;
  mXiDst.primvertexZ = v0Dst.primvertexZ;
  // Do 'event' analysis based on event data 
  
  StThreeVectorF pv;
  pv.setX(mXiDst.primvertexX);
  pv.setY(mXiDst.primvertexY);
  pv.setZ(mXiDst.primvertexZ);

  // Record some information...
  Double_t magn = muEvent->runInfo().magneticField();	
  //Double_t magn = muEvent->magneticField();	//checked! the same as above
  mXiDst.magn = magn;
  
  // Select Bachelor tracks ...
  // Get 'track' data, make cuts on tracks, do physics analysis, histogram results.
  TObjArray* tracks = mMuDstMaker->muDst()->globalTracks() ;    // Create a TObject array containing the global tracks
  TObjArrayIter GetTracks(tracks) ;                              // Create an iterator to step through the tracks

  StMuTrack* track ;                                             // Pointer to a track
  while ( ( track = (StMuTrack*)GetTracks.Next() ) )             // Main loop for Iterating over tracks
  {
     hPt -> Fill( track->pt() ) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.
     //double eta = track->eta();		  //at dca to PV
     //double phi = track->phi();		  //at dca to PV
     short flag = track->flag();
     unsigned short nHits = track->nHits();
     unsigned short nHitsFit = track->nHitsFit(kTpcId);
     short charge = track->charge();
     //StThreeVectorF p = track->p();	  //at dca to PV
     //StThreeVectorF origin = track->firstPoint();  //? NOT SURE whether firstPoint from detectorinfo is the same as helix.orgin()!!!
     double nsigmapion = track->nSigmaPion();
     double nsigmaproton = track->nSigmaProton();
     double nsigmakaon = track->nSigmaKaon();
     double dedx = track->dEdx();

     StPhysicalHelixD helix = track->helix(); //inner helix. good for dca to PV.
     StThreeVectorD p = helix.momentum(magn*kilogauss);    //momentum at origin
     StThreeVectorD origin = helix.origin();  //origin of helix
     double pt = p.perp();

     //some checks.
     hPtDiff -> Fill( track->pt() - p.perp() ) ; 
     hOrDiff -> Fill( (track->firstPoint() - origin).mag() ) ; 
     hPDiff -> Fill( fabs(track->p().mag() - p.mag()) ) ; 
     //comments: there are difference between the values above. But they seem to be acceptably small!

     int hrot;		//helicity of helix, sign of -charge*magn
     if (-charge*magn > 0) hrot = 1;
     else hrot = -1;

     double nsigma;
     if(mXiType == kXi || mXiType == kAntiXi) nsigma = nsigmapion;
     else nsigma = nsigmakaon;

     //if(track->vertexIndex()!=StMuDst::currentVertexIndex())continue;
     //if you want to use track->dca(), turn this on. if it is not turned on, that function crashes.
     //OK. let's cut tracks
     //if(nHitsFit<=cutNHitsGr)continue;
     if(nHits<=cutNHitsGr)continue;
     if(flag <=0  )continue; //or <=0 ?
     //if(track->bad() )continue; 
     if(abs(charge)!=1) continue;
     if(pt<cutPtGrEq)continue; //should be larger. like 0.15 or 0.2

     if(charge == mChargeBachelor && fabs(nsigma)<cutAbsNSigma1Le){
	  //record the bachelor 
	  //fill the vector
	  StPhysicalHelixD helix = track->helix(); //inner helix. good for dca to PV.

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
	  if(dca.mag()<cutDca1GrEq)continue;
	  mDauDcaVec1.push_back(dca.mag());
	  mDauVec1.push_back(track);
     }
  }

  assert(mDauVec1.size() == mDauDcaVec1.size());
 // cout<<mDauVec1.size()<<" == "<<mDauDcaVec1.size()<<endl; 
  //reconstruct Xi
  int nXi = 0;
 // cout<< "there are "<<v0Dst.nv0<<" V0"<<endl;
  for(int i=0;i<v0Dst.nv0;i++){
     //get v0 info here
     //cut them before in StV0Maker

     if(v0Dst.v0dca[i]<cutV0DcaGrEq)continue;
     if(fabs(v0Dst.v0mass[i]-mMassV0)>cutV0MassWidthLeEq)continue;
     
     StThreeVectorF xv0(v0Dst.v0x[i],v0Dst.v0y[i],v0Dst.v0z[i]), 
			  pv0(v0Dst.v0px[i],v0Dst.v0py[i],v0Dst.v0pz[i]);
   //  cout<<"get "<<i<<"th v0 fine"<<endl;
     
     for(unsigned int j=0;j<mDauVec1.size();j++){
	  //get pion track info here
	  StMuTrack * track1 = mDauVec1[j];
	  if(track1->id() == v0Dst.dau2id[i])continue;  
	  StPhysicalHelixD helix1 = track1->helix(); //inner helix. good for dca to PV.
	  
	  if(mRotate){
	     StThreeVectorD tp1 = helix1.momentum(magn*kilogauss);    //momentum at origin
	     StThreeVectorD tx1 = helix1.origin();    //origin
	     tp1.setX(-tp1.x());
	     tp1.setY(-tp1.y());
	     tx1.setX(-(tx1.x()-pv.x())+pv.x());
	     tx1.setY(-(tx1.y()-pv.y())+pv.y());
	     StPhysicalHelixD helixtmp(tp1, tx1, magn*kilogauss, track1->charge());
	     helix1 = helixtmp;
	  }
	
	  double dca1 = mDauDcaVec1[j];
	  StThreeVectorD p1 = helix1.momentum(magn*kilogauss);    //momentum at origin
	  //cut them before in track selection
	  
	  StThreeVectorF xxi, op1;
	  
	  StPhysicalHelixD helixv0(pv0/pv0.perp()*100, xv0, magn*kilogauss, 1);
	  //StHelix dca algorithm is FAST in case of helix to staight line. 
	  //We will use it as default algorithm, 
	  //if you want to try out Long's method, switch on 
	  //mDcaAlgoLong in the constructor.
	  double dca1tov0;
	  if(!mDcaAlgoLong){
	     //v0 is a straight line. to simulate it, we boost its transverse momentum to 100GeV/c.
	     //about 700 meters transverse radii. 
	     //ALWAYS call helix1's pathLengths method, treat helixv0 as a parameter.
	     //Otherwise, bachelor helix1's periods will be scaned, the background will be HUGE!
	     //StPhysicalHelixD helixv0(pv0/pv0.perp()*100, xv0, magn*kilogauss, 1);
	     //LOG_QA << 1.0/helixv0.curvature()<<endm;
	     pair<double,double> tmps = helix1.pathLengths(helixv0); 
	     StThreeVectorD ox1 = helix1.at(tmps.first);
	     StThreeVectorD ox2 = helixv0.at(tmps.second);
	     
	     dca1tov0 = (ox1-ox2).mag();
	     xxi = (ox1+ox2)/2.;
	     op1 = helix1.momentumAt(tmps.first, magn*kilogauss);
	  }

	  //LONG's two helices dca method, should also work. also need to construct v0 helix first.
	  //this method work well for smaller v0 pt (less than 100GeV). but has some difference 
	  //in some rare cases. should not be used for serious Xi analysis. 
	  //the following lines are for debug only!
	  if(mDcaAlgoLong){
	     StThreeVectorF op2;
	     //StPhysicalHelixD helixv0(pv0/pv0.perp()*100., xv0, magn*kilogauss, 1);
	     dca1tov0 = closestDistance(helix1, helixv0, magn, pv, xxi, op1, op2);
	  }
	  
	  //if(mDcaAlgoLong) dca1tov0 = closestDistance(xv0, pv0, helix1, magn, xxi, op1);
	  //cut on dca1to2
	  if(dca1tov0 > cutDca1to2LeEq) continue;
	  //LOG_QA<<stdca1tov0<<" "<<longdca1tov0<<" "<<dca1tov0<<endm;
	  
	  double oe1 = sqrt(op1.mag2() + mMassBachelor*mMassBachelor);
	  double oev0 = sqrt(pv0.mag2() + mMassV0*mMassV0);
	  double ximass = sqrt(mMassBachelor*mMassBachelor + mMassV0*mMassV0 + 2.*oe1*oev0 - 2.*op1.dot(pv0));
	  //cut on ximass
	  if(fabs(ximass-mMassXi) > cutXiMassWidthLeEq)continue;

	  //StThreeVectorD xv0 = (ox1 + ox2)/2.;
	  StThreeVectorD pxi = op1 + pv0;
	  StThreeVectorD xxitoPV = xxi - pv;

	  /*
	  //pthead, ptarm cut
	  double pthead1 = op1.dot(pv0)/pv0.mag();
	  double pthead2 = op2.dot(pv0)/pv0.mag();
	  double ptarm   = sqrt(op1.mag2()-pthead1*pthead1);
	  if(mXiType!=kKs){
	     if(ptarm> cutDauPtArmLeEq || fabs((pthead1-pthead2)/(pthead1+pthead2))> cutAbsDausPtShoulderDiffLeEq )continue;
	  }
	  */

	  //forward decay cut
	  double ang1 = (pxi.x()*xxitoPV.x() + pxi.y()*xxitoPV.y())/pxi.perp()/xxitoPV.perp();
	  //if(ang1<=cutDau1DecAngGr)continue;
	  
	  //r dot p for xi. cut on it. should be larger than 0. 
	  double rdotp = xxitoPV.dot(pxi);
	  if(rdotp/xxitoPV.mag()/pxi.mag() <= cutXirdotpGr)continue;
	  if(pv0.dot(xv0-xxi) <= cutXirdotpGr)continue;
	  
	  //calculate xi to PV dca. cut on dca
	  StPhysicalHelixD helixxi(pxi, xxi, magn*kilogauss, mChargeXi);
	  double pathlength = helixxi.pathLength(pv, false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
	  StThreeVectorF dca = helixxi.at(pathlength)-pv;
	  double dcaxitoPV = dca.mag();
	  //if(dcaxitoPV>=cutDcaXiLe)continue;

	  //cut on decay length
	  double xidecaylength = xxitoPV.mag();
	  //if(xidecaylength > v0Dst.v0declen[i] )continue;
	  if(xidecaylength < cutXiDecLenGrEq) continue;

	  //cut on sinth, or theta
	  double sinth = (xxitoPV.cross(pxi)).mag()/xxitoPV.mag()/pxi.mag();
	  double theta = atan2(sinth, rdotp/xxitoPV.mag()/pxi.mag()); //theta range: from 0 to pi
	  if(sinth > cutXircrosspLeEq) continue;

	  //record TOF information...
	  int tofflag1 = track1->btofPidTraits().matchFlag();
	  double tof1 = track1->btofPidTraits().timeOfFlight();
	  StThreeVectorF tofpos1 = track1->btofPidTraits().position();
	  double tofpathlen1 = -999.;
	  if(tofflag1>0) tofpathlen1 = helix1.pathLength(tofpos1) - helix1.pathLength(xxi);

	  double v0pathlen = helixv0.pathLength(xv0)-helixv0.pathLength(xxi);
	  double xipathlen = helixxi.pathLength(xxi)-helixxi.pathLength(pv);

	  //save the xi information (including cut variables) into StXiDst. 
	  mXiDst.v0mass[nXi] = v0Dst.v0mass[i];
	  mXiDst.v0pt[nXi]   = v0Dst.v0pt[i];
	  mXiDst.v0rapidity[nXi]    = v0Dst.v0rapidity[i];
	  mXiDst.v0eta[nXi]    = v0Dst.v0eta[i];
	  mXiDst.v0x[nXi]    = v0Dst.v0x[i];
	  mXiDst.v0y[nXi]    = v0Dst.v0y[i];
	  mXiDst.v0z[nXi]    = v0Dst.v0z[i];
	  mXiDst.v0px[nXi]    = v0Dst.v0px[i];
	  mXiDst.v0py[nXi]    = v0Dst.v0py[i];
	  mXiDst.v0pz[nXi]    = v0Dst.v0pz[i];
	  mXiDst.v0declen[nXi] = v0Dst.v0declen[i];
	  mXiDst.v0dca[nXi]    = v0Dst.v0dca[i];
	  mXiDst.v0dca2d[nXi]    = v0Dst.v0dca2d[i];
	  mXiDst.v0pathlen[nXi]  = v0pathlen;

	  mXiDst.dau1id[nXi]    = v0Dst.dau1id[i];
	  mXiDst.dau1dca[nXi]    = v0Dst.dau1dca[i];
	  mXiDst.dau1dca2d[nXi]    = v0Dst.dau1dca2d[i];
	  mXiDst.dau1nhits[nXi]    = v0Dst.dau1nhits[i];
	  mXiDst.dau1dedx[nXi]    = v0Dst.dau1dedx[i];
	  mXiDst.dau1nsigma[nXi]    = v0Dst.dau1nsigma[i];
	  mXiDst.dau1eta[nXi]    = v0Dst.dau1eta[i];
	  mXiDst.dau1pt[nXi]    = v0Dst.dau1pt[i];
	  mXiDst.dau1px[nXi]    = v0Dst.dau1px[i];
	  mXiDst.dau1py[nXi]    = v0Dst.dau1py[i];
	  mXiDst.dau1pz[nXi]    = v0Dst.dau1pz[i];
	  mXiDst.dau1tpc[nXi]    = v0Dst.dau1tpc[i];
	  mXiDst.dau1ssd[nXi]    = v0Dst.dau1ssd[i];
	  mXiDst.dau1svt[nXi]    = v0Dst.dau1svt[i];
	  mXiDst.dau1tofflag[nXi]    = v0Dst.dau1tofflag[i];
	  mXiDst.dau1tof[nXi]    = v0Dst.dau1tof[i];
	  mXiDst.dau1pathlen[nXi]    = v0Dst.dau1pathlen[i];

	  mXiDst.dau2id[nXi]    = v0Dst.dau2id[i];
	  mXiDst.dau2dca[nXi]    = v0Dst.dau2dca[i];
	  mXiDst.dau2dca2d[nXi]    = v0Dst.dau2dca2d[i];
	  mXiDst.dau2nhits[nXi]    = v0Dst.dau2nhits[i];
	  mXiDst.dau2dedx[nXi]    = v0Dst.dau2dedx[i];
	  mXiDst.dau2nsigma[nXi]    = v0Dst.dau2nsigma[i];
	  mXiDst.dau2eta[nXi]    = v0Dst.dau2eta[i];
	  mXiDst.dau2pt[nXi]    = v0Dst.dau2pt[i];
	  mXiDst.dau2px[nXi]    = v0Dst.dau2px[i];
	  mXiDst.dau2py[nXi]    = v0Dst.dau2py[i];
	  mXiDst.dau2pz[nXi]    = v0Dst.dau2pz[i];
	  mXiDst.dau2tpc[nXi]    = v0Dst.dau2tpc[i];
	  mXiDst.dau2ssd[nXi]    = v0Dst.dau2ssd[i];
	  mXiDst.dau2svt[nXi]    = v0Dst.dau2svt[i];
	  mXiDst.dau2tofflag[nXi]    = v0Dst.dau2tofflag[i];
	  mXiDst.dau2tof[nXi]    = v0Dst.dau2tof[i];
	  mXiDst.dau2pathlen[nXi]    = v0Dst.dau2pathlen[i];
	  
	  mXiDst.dca1to2[nXi]  = v0Dst.dca1to2[i];

	  mXiDst.bachid[nXi] = track1->id();
        mXiDst.bachdca[nXi]  = dca1;
        mXiDst.bachdca2d[nXi]  = helix1.geometricSignedDistance(pv.x(),pv.y());
        mXiDst.bachnhits[nXi]= track1->nHits();
        mXiDst.bachdedx[nXi] = track1->dEdx();
        mXiDst.bachnsigma[nXi] = (mXiType==kXi||mXiType==kAntiXi)? track1->nSigmaPion() : track1->nSigmaKaon();
	  mXiDst.bacheta[nXi]  = op1.pseudoRapidity();
        mXiDst.bachpt[nXi]   = op1.perp();
        mXiDst.bachpx[nXi]   = op1.x();
        mXiDst.bachpy[nXi]   = op1.y();
        mXiDst.bachpz[nXi]   = op1.z();
        mXiDst.bachtpc[nXi]  = track1->nHitsFit(kTpcId);
        mXiDst.bachssd[nXi]  = track1->nHitsFit(kSsdId);
        mXiDst.bachsvt[nXi]  = track1->nHitsFit(kSvtId);
	  mXiDst.bachtofflag[nXi] = tofflag1;
	  mXiDst.bachtof[nXi] = tof1;
	  mXiDst.bachpathlen[nXi] = tofpathlen1;

	  
	  mXiDst.dcav0tobach[nXi] = dca1tov0;
	  //mXiDst.stdcav0tobach[nXi] = stdca1tov0;
	  
	  mXiDst.ximass[nXi]   = ximass;
        mXiDst.xidca[nXi]    = dcaxitoPV;
        mXiDst.xidca2d[nXi]    = helixxi.geometricSignedDistance(pv.x(),pv.y());
        mXiDst.xipt[nXi]   = pxi.perp();
        mXiDst.xirapidity[nXi]    = log( (sqrt(ximass*ximass+pxi.mag2()) + pxi.z())/sqrt(ximass*ximass+pxi.perp2()));
        mXiDst.xieta[nXi]    = 0.5*log( (pxi.mag() + pxi.z())/(pxi.mag() - pxi.z()) );
        mXiDst.xix[nXi]      = xxi.x();
        mXiDst.xiy[nXi]      = xxi.y();
        mXiDst.xiz[nXi]      = xxi.z();
        mXiDst.xipx[nXi]     = pxi.x();
        mXiDst.xipy[nXi]     = pxi.y();
        mXiDst.xipz[nXi]     = pxi.z();
        mXiDst.xideclen[nXi] = xidecaylength;
	  mXiDst.xisinth[nXi]  = sinth;
	  mXiDst.xipathlen[nXi]  = xipathlen;

	  nXi ++;

	  
	  hInvMass->Fill(ximass);
	  //mXiVec.push_back(v0data);

     }
  }

  mXiDst.nxi = nXi;

  if(GetDebug()) {
     LOG_QA<<"in StV0Maker::Make : has "<<v0Dst.nv0<<" V0"<<endm;
     LOG_QA<<"in StXiMaker::Make : has "<<nXi<<" Xi"<<endm;
     LOG_QA<<".."<<mEventsProcessed<<endm;
  }
  if(nXi > 0 && mXiTree) mXiTree->Fill();
  //dump Xi vector into a TTree, record those events with Xi candidate only.
  //refer to hSelectNRefMult generated by StV0Maker for total number of events
  
  mEventsProcessed++ ;
  //LOG_QA << mDauVec1.size() <<" "<< mDauVec2.size() <<" "<< mXiVec.size()<<endm;
  mDauDcaVec1.clear();
  mDauDcaVec2.clear();
  mDauVec1.clear();
  mDauVec2.clear();
  //mXiVec.clear();
  
  return kStOK ;
  
}


Int_t StXiMaker::Finish( )
{ // Do once at the end the analysis

  // Write histograms to disk, output miscellaneous other information
  if(histogram_output!=NULL) histogram_output -> Write() ;   // Write all histograms to disk 
  if(xitree_output!=NULL) xitree_output -> Write() ;   // Write all histograms to disk 

  cout << "Total Events Processed in StXiMaker " << mEventsProcessed << endl ;

  return kStOk ;  

}


