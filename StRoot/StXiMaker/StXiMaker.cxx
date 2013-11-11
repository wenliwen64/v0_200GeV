#define DEBUGXI
//#define DEBUGXI_HISTO 1
#include "TTree.h"
#include "TFile.h"
#include "StXiMaker.h"
#include "StXiDst.h"
//#include "StV0Maker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "TString.h"
//#include "StV0Type.h"
//#include "StV0Dst.h"
#include "TH1.h"
#include "TH2.h"
#include <stdio.h>
#include "StMuDSTMaker/COMMON/StMuEvent.h"
//#include "StTypeDefs.h"
#include "TClonesArray.h"
#include <vector>
#include <cmath>
#include <utility>
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StTriggerId.h"
#include "TMath.h"

#include "StDcaService.h"
ClassImp(StXiMaker)


StXiMaker::StXiMaker(StMuDstMaker* muDstMaker, 
StV0Maker* V0Maker,
TString outDSTFileName,
TString outHistoFileName
): 
mOutHistoFile(NULL),  
mXiPicoDstFile(NULL),
mOutHistoFileName(outHistoFileName), 
mXiPicoDstFileName(outDSTFileName), 
//mV0ParticleName(""),
mRotation(0),
mSameSign(0),
mXiZeroMass(1.314),
mXiMinusMass(1.321),
mXiMassWidth(0.2), // The default unit of energy or mass in STAR is GeV.
mProtonMass(0.938),//???
mAntiProtonMass(0.938),//???
mPionMinusMass(0.140),//???
mPionPlusMass(0.140),//???
mBuffer(""), 
mHInvariantMassName(""), 
//mV0Type(V0Type),
//mGlobalTracks(0),// ??
//mPrimaryTracks(0),// ??
mMuDstMaker(muDstMaker),
mV0Maker(V0Maker),
mEventFailed(0),
mEventPassed(0),
//mV0Dst(),
mXiTree(NULL),
//mCutTriggerIdEq(360001),
mCutAbsPrimaryVertexZLeEq(50.0),
mCutNHitsGrEq(15),
mCutPtGrEq(0.15),
mCutNSigmaPionLe(2.0),
mCutNSigmaProtonLe(1.5),
mCutNSigmaElectronLe(4.0),
mCutNSigmaKaonLe(4.0),
mCutDau1Dca2PVGrEq(3),
mCutDau2Dca2PVGrEq(3),
mCutTwoTracksDcaLeEq(1.5),
mCutDecLengthGrEq(3.0),
mCutDcaXi2PVLe(5.0){

//    std::vector<std::pair<StMuTrack*, StMuTrack*> > pairTracksVector;
//    std::vector<std::pair<StMuTrack*, StMuTrack*> > antiPairTracksVector;

//    std::vector<std::pair<Double_t, Double_t> > momentumPairVector;
//    std::vector<std::pair<Double_t, Double_t> > antiMomentumPairVector;

#ifdef DEBUGXI_HISTO
    mOutHistoFile = NULL;
    mXiPicoDstFile = NULL;
#endif
} 
  

Int_t  StXiMaker::Init(){

    //-------------------Initialize output files----------------------------------------------------------------------------------------
    mXiPicoDstFile = new TFile(mXiPicoDstFileName, "recreate");       // mXiPicoDstFileName is given by user code
    mOutHistoFile = new TFile(mOutHistoFileName, "recreate");         // mOutHistoFileName is given by user code
#ifdef DEBUGXI_HISTO 
    printf("Initialization happens right now\n");
#endif

    //-------------------Book Histograms----------------------------------------------------------------------------------------
//    hNPV = new TH1F();
    hVz = new TH1F("VertexZ","Event Vertex Z position", 400, -100, 100);
          //-------------Initialize invariant-mass spectrum histogram 
    TString mParticleName = "CascadeMinus";
    mHInvariantMassName.Append(mParticleName);
    mHInvariantMassName.Append("InvariantMass");
    //TString antiHInvariantMassName = mHInvariantMassName + "Anti";
         //------------- 
    //hNPV = new TH1F("numberOfPV", "number of primary vertices", 10, 0, 10);
    hRefMult = new TH1F("refMultiplicity", "from event->refmult()", 1000, 0, 1000);
    //hSelectNRefMult = new TH1F("SelectRefMult", "Reference Multiplicity of selected events", 1000, 0, 1000);         //-------------?
    //hSelectNRefMultM1 = new TH1F("SelectRefMult1", "Reference Multiplicity of selected events", 1000, 0, 1000);       //-------------?
    //hSelectNRefMultM2 = new TH1F("SelectRefMult2", "Reference Multiplicity of selected events", 1000, 0, 1000);       //-------------?
    //hSelectNRefMultFtpcE = new TH1F("SelectRefMultFtpcE", "Reference Multiplicity of selected events(FTPC East)", 1000, 0, 1000);
    //hSelectBbcEAdcSum = new TH1F("SelectBbcAdcSum", "BBC East Adc Sum", 5000, 0, 5000);
         //---------------QA for global tracks-?
    hPtRaw = new TH1F("PtRaw", "Transverse Momentum for all particles", 600, 0, 30);   
    hEtaRaw = new TH1F("EtaRaw", "Eta for all particles", 500, -2, 2);
    hPhiRaw = new TH1F("PhiRaw", "Phi for all particles", 1000, -TMath::Pi(), TMath::Pi());
    hPt = new TH1F("Pt", "tracks pt after some initial cuts", 600, 0, 30);
    hEta = new TH1F("Eta", "tracks eta after some initial cuts", 500, -2, 2);
    hPhi = new TH1F("Phi", "tracks phi after some initial cuts", 1000, -TMath::Pi(), TMath::Pi());
    //hPhiLowPt = new TH1F("PhiLowPt", "Phi for selected particles of low Pt", 1000, -Pi, Pi); 
    //hPhiHighPt = new TH1F("PhiHighPt", "Phi for selected particles of high Pt", 1000, -Pi, Pi); 
    hDeDxP = new TH2F("DeDxP", "dEdx for selected particles", 300, -3, 3, 300, 0, 3e-5);    //-------------dedxP what does the P stand for, how to determine the range of histograms
    hNHitsFit = new TH1F("nHitsFit", "nHitsFit for all particles", 80, 0, 80);
    hNHits = new TH1F("nHits", "nHitsfor all particles", 80, 0, 80);
    hNSigmaPion = new TH1F("nSigmaPion", "nSigmaPion for selected particles", 200, -10, 10);
    hNSigmaProton = new TH1F("nSigmaProton", "nSigmaProton for selected particles", 200, -10, 10);
//    hNSigmaKaon = new TH1F("nSigmaKaon", "nSigmaKaon for selected particles", 200, -10, 10);
  
//    hPtDiff = new TH1F("ptdiff", "Reference Multiplicity", 100, 0, 0.02); //---------what the hell is this
//    hOrientationDiff = new TH1F("orientationdiff", "Reference Multiplicity", 100, 0, 0.2);
//    hPDiff = new TH1F("pdiff", "Reference Multiplicity", 100, 0, 4);
//    hDcaDiff = new TH1F("dcadiff", "Reference Multiplicity",100, -0.4, 0.4);

//    hDau1Dca = new TH1F("hmDau1Dca", "Dca to PV of dau1 using dcaGlobal()", 500, 0, 25);
//    hDau2Dca = new TH1F("hmDau2Dca", "Dca to PV of dau2 using dcaGlobal()", 500, 0, 25);
//    hDau1Dca2 = new TH1F("hmDau1Dca2", "Dca to PV of dau1 using another method", 500, 0, 25);
//    hDau2Dca2 = new TH1F("hmDau2Dca2", "Dca to PV of dau2 using another method", 500, 0, 25);

    //hXiZeroInvariantMass = new TH1F(mHInvariantMassName, "Invariant Mass", 300, mXiZeroMass - mXiMassWidth, mXiZeroMass + mXiMassWidth);
    hXiMinusInvariantMass = new TH1F(mHInvariantMassName, "Anti Particle Invariant Mass", 300, mXiMinusMass - mXiMassWidth, mXiMinusMass + mXiMassWidth);
#ifdef DEBUGXI
    printf("histogram OK!\n");
#endif

    //------------Initialize tree-----------------------------------------------------------------------------------

    //mXiPicoDstFileName.Append(mXiParticleName);    
    //mXiPicoDstFile = new TFile(mXiPicoDstFileName, "recreate");       // mXiPicoDstFileName is given by user code
    //mOutHistoFile = new TFile(mOutHistoFileName, "recreate");      // mOutHistoFileName is given by user code
#ifdef DEBUGXI_HISTO
    gDirectory -> pwd();
    if(mOutHistoFile) printf("histoFile open OK\n");
    if(mXiPicoDstFile) printf("picoDataFile open OK\n");
#endif
    mXiTree = new TTree("XiPicoDst", "From MuDst");

    mXiTree -> SetDirectory(mXiPicoDstFile ); 

    mXiTree -> Branch("runNumber", &mXiDst.runNumber, "runNumber/I");
    mXiTree -> Branch("eventId", &mXiDst.eventId, "eventId/I");
    //mXiTree -> Branch("triggerMode", &mXiDst.triggerMode, "triggerMode/I");
    mXiTree -> Branch("nRefMult", &mXiDst.nRefMult, "nRefMult/I");
    mXiTree -> Branch("nRefMultTof", &mXiDst.nRefMultTof, "nRefMultTof/I");
    mXiTree -> Branch("magneticField", &mXiDst.magneticField, "magneticField/F");
    mXiTree -> Branch("PVX", &mXiDst.PVX, "PVX/F");
    mXiTree -> Branch("PVY", &mXiDst.PVY, "PVY/F");
    mXiTree -> Branch("PVZ", &mXiDst.PVZ, "PVZ/F");
    mXiTree -> Branch("numberOfXi", &mXiDst.numberOfXi, "numberOfXi/I");

    mXiTree -> Branch("xiType", mXiDst.xiType, "xiType[numberOfXi]/I");
    mXiTree -> Branch("xiMass", mXiDst.xiMass, "xiMass[numberOfXi]/F");
    mXiTree -> Branch("xiPt", mXiDst.xiPt, "xiPt[numberOfXi]/F");
    mXiTree -> Branch("xiRapidity", mXiDst.xiRapidity, "xiRapidity[numberOfXi]/F");
    mXiTree -> Branch("xiEta", mXiDst.xiEta, "xiEta[numberOfXi]/F");
    mXiTree -> Branch("xiX", mXiDst.xiX, "xiX[numberOfXi]/F");
    mXiTree -> Branch("xiY", mXiDst.xiY, "xiY[numberOfXi]/F");
    mXiTree -> Branch("xiZ", mXiDst.xiZ, "xiZ[numberOfXi]/F");
    mXiTree -> Branch("xiPx", mXiDst.xiPx, "xiPx[numberOfXi]/F");
    mXiTree -> Branch("xiPy", mXiDst.xiPy, "xiPy[numberOfXi]/F");
    mXiTree -> Branch("xiPz", mXiDst.xiPz, "xiPz[numberOfXi]/F");
    mXiTree -> Branch("xiDecLength", mXiDst.xiDecLength, "xiDecLength[numberOfXi]/F");
    mXiTree -> Branch("xiDca2PV", mXiDst.xiDca2PV, "xiDca2PV[numberOfXi]/F");
    mXiTree -> Branch("xiDca2PV2d", mXiDst.xiDca2PV2d, "xiDca2PV2d[numberOfXi]/F");
    mXiTree -> Branch("xiDca2Dau", mXiDst.xiDca2Dau, "xiDca2Dau[numberOfXi]/F");

    //mXiTree -> Branch("dau1TrackId", mXiDst.dau1TrackId, "dau1TrackId[numberOfXi]/I");
    //mXiTree -> Branch("dau1Charge", mXiDst.dau1Charge, "dau1Charge[numberOfXi]/I");
    //mXiTree -> Branch("v0Dca2PV", mXiDst.v0Dca2PV, "v0Dca2PV[numberOfXi]/F");
    //mXiTree -> Branch("dau1Dca2PV2d", mXiDst.dau1Dca2PV2d, "dau1Dca2PV2d[numberOfXi]/F");
    //mXiTree -> Branch("dau1NHits", mXiDst.dau1NHits, "dau1NHits[numberOfXi]/I");
      //mXiTree -> Branch("dau1DeDx", mXiDst.dau1DeDx, "dau1DeDx[numberOfXi]/F");
      //mXiTree -> Branch("dau1DeDx", mXiDst.dau1DeDx, "dau1DeDx[numberOfXi]/F");
    //mXiTree -> Branch("dau1NSigmaPion", mXiDst.dau1NSigmaPion, "dau1NSigmaPion[numberOfXi]/F");
    //mXiTree -> Branch("dau1NSigmaProton", mXiDst.dau1NSigmaProton, "dau1NSigmaProton[numberOfXi]/F");
    mXiTree -> Branch("v0Pt", mXiDst.v0Pt, "v0Pt[numberOfXi]/F");
    mXiTree -> Branch("v0Px", mXiDst.v0Px, "v0Px[numberOfXi]/F");
    mXiTree -> Branch("v0Py", mXiDst.v0Py, "v0Py[numberOfXi]/F");
    mXiTree -> Branch("v0Pz", mXiDst.v0Pz, "v0Pz[numberOfXi]/F");
    //mXiTree -> Branch("dau1NHitsTpc", mXiDst.dau1NHitsTpc, "dau1NHitsTpc[numberOfXi]/I");
    //mXiTree -> Branch("dau1NHitsSsd", mXiDst.dau1NHitsSsd, "dau1NHitsSsd[numberOfXi]/I");
    //mXiTree -> Branch("dau1NHitsSvt", mXiDst.dau1NHitsSvt, "dau1NHitsSvt[numberOfXi]/I");
    //mXiTree -> Branch("dau1TofFlag", mXiDst.dau1TofFlag, "dau1TofFlag[numberOfXi]/I");
    //mXiTree -> Branch("dau1Tof", mXiDst.dau1Tof, "dau1Tof[numberOfXi]/F");

    mXiTree -> Branch("pionMinusTrackId", mXiDst.pionMinusTrackId, "pionMinusTrackId[numberOfXi]/I");
    mXiTree -> Branch("pionMinusCharge", mXiDst.pionMinusCharge, "pionMinusCharge[numberOfXi]/I");
    mXiTree -> Branch("pionMinusDca2PV", mXiDst.pionMinusDca2PV, "pionMinusDca2PV[numberOfXi]/F");
    mXiTree -> Branch("pionMinusDca2PV2d", mXiDst.pionMinusDca2PV2d, "pionMinusDca2PV2d[numberOfXi]/F");
    mXiTree -> Branch("pionMinusNHits", mXiDst.pionMinusNHits, "pionMinusNHits[numberOfXi]/I");
      //mXiTree -> Branch("pionMinusDeDx", mXiDst.dau1DeDx, "pionMinusDeDx[numberOfXi]/F");
      //mXiTree -> Branch("pionMinusDeDx", mXiDst.dau1DeDx, "pionMinusDeDx[numberOfXi]/F");
    mXiTree -> Branch("pionMinusNSigmaPion", mXiDst.pionMinusNSigmaPion, "pionMinusNSigmaPion[numberOfXi]/F");
    //mXiTree -> Branch("dau2NSigmaProton", mXiDst.dau2NSigmaProton, "dau2NSigmaProton[numberOfXi]/F");
    mXiTree -> Branch("pionMinusPt", mXiDst.pionMinusPt, "pionMinusPt[numberOfXi]/F");
    mXiTree -> Branch("pionMinusPx", mXiDst.pionMinusPx, "pionMinusPx[numberOfXi]/F");
    mXiTree -> Branch("pionMinusPy", mXiDst.pionMinusPy, "pionMinusPy[numberOfXi]/F");
    mXiTree -> Branch("pionMinusPz", mXiDst.pionMinusPz, "pionMinusPz[numberOfXi]/F");
    mXiTree -> Branch("pionMinusNHitsTpc", mXiDst.pionMinusNHitsTpc, "pionMinusNHitsTpc[numberOfXi]/I");
    mXiTree -> Branch("pionMinusNHitsSsd", mXiDst.pionMinusNHitsSsd, "pionMinusNHitsSsd[numberOfXi]/I");
    mXiTree -> Branch("pionMinusNHitsSvt", mXiDst.pionMinusNHitsSvt, "pionMinusNHitsSvt[numberOfXi]/I");
    mXiTree -> Branch("pionMinusTofFlag", mXiDst.pionMinusTofFlag, "pionMinusTofFlag[numberOfXi]/I");
    mXiTree -> Branch("pionMinusTof", mXiDst.pionMinusTof, "pionMinusTof[numberOfXi]/F");

    return kStOK;
}

Int_t StXiMaker::Make() {
    const StV0Dst& v0Dst = mV0Maker -> getV0Dst(); 
#ifdef DEBUGXI_HISTO
    gDirectory -> pwd();
#endif
#ifdef DEBUGXI
    printf("New Event in Cascade!\n");
#endif
#ifdef DEBUGXI
    printf("Now I am making sth Cascade\n");
#endif
    StMuEvent* muEvent = mMuDstMaker -> muDst() -> event(); // I think mMuDstMaker would be automatically move to next event as well as collection of tracks
    
    Double_t magneticField = muEvent -> runInfo().magneticField();

    if(!pCutGoodEvents(muEvent)) {  mEventFailed++; return kStOK; }
#ifdef DEBUGXI
    printf("event good!\n");
#endif
    if(!pCutEventsTrigger(muEvent)) return kStOK;
#ifdef DEBUGXI
    printf("trigger good! in Cascade\n");
#endif
    if(!pCutEventsPVZ(muEvent)) return kStOK;
#ifdef DEBUGXI
    printf("PVZ good!\n");
#endif

    //---------Event-wise information---------------------------------------------------------------------------------------------------------------
    mXiDst.eventId = muEvent -> eventId();
    mXiDst.runNumber = muEvent -> runId();
    //mXiDst.triggerMode = mCutTriggerIdEq; 
    mXiDst.nRefMult = muEvent -> refMult();
    mXiDst.nRefMultTof = muEvent -> btofTrayMultiplicity();
    mXiDst.magneticField = muEvent -> runInfo().magneticField(); 
    mXiDst.PVX = v0Dst.PVX;//TODO: check this item. mXiDst.PVX = muEvent -> primaryVertexPosition().x();  
    mXiDst.PVY = v0Dst.PVY;//muEvent -> primaryVertexPosition().y();  
    mXiDst.PVZ = v0Dst.PVZ;//muEvent -> primaryVertexPosition().z();  
    //-----------Fill event-wise information histograms------------------------------------------------------------------------------- 
    hRefMult -> Fill(muEvent -> refMult());//I dont get what does refMult mean
    //------------------------------------------------------------------------------------------------------------------------------------------------
    StThreeVectorF pv(mXiDst.PVX, mXiDst.PVY, mXiDst.PVZ);//TODO: check
    int numberOfXi = 0;
    //#ifdef DEBUGXI
    //    else printf("Event OK!\n");  
    //#endif
    //if(!muEvent){
    //  mNEventFailed++; //
    //   return kStOk; // What if I return other variables like kStErr?? I guess maybe if I return kStErr, the chain will  not be able to continue working.
    //}

    //Fill some histograms
    hVz -> Fill(muEvent -> primaryVertexPosition().z());
    // the following code will finish i) trigger ID selection; ii) eliminate duplicate vertex( only existing in subsequent events?); iii) eliminate events without vertex.
    // Apply cuts: i)vertexZ; ii) centrality or reference multiplicity; iii) track cuts---Pt, eta, phi, flag(??), nHits, nHitsFit(??), charge, nsigmapion, proton, kaon, electron, dEdx(??);

    TClonesArray* mGlobalTracks = (TClonesArray*) mMuDstMaker -> muDst() -> globalTracks();
    TClonesArray* mPrimaryTracks = (TClonesArray*) mMuDstMaker -> muDst() -> primaryTracks();

    //Int_t globalTrackNumber=mGlobalTracks->GetEntries();
#ifdef DEBUGXI
    //printf("the number of total global tracks is %d\n", globalTrackNumber);
#endif
    //if(!StXiMaker.isGlobalTracksOk(muEvent)) return kStOk;
    //if(!StXiMaker.isPrimaryTracksOk(muEvent)) return kStOk;

    TIter globalTrackIter(mGlobalTracks);
    TIter primaryTrackIter(mPrimaryTracks);

    StMuTrack* singleTrack = new StMuTrack;

    //vector<StMuTrack*> dau1TracksVector;// for Lambda reconstruction------proton particles
    //vector<StMuTrack*> antiDau1TracksVector;// for AntiLambda reconstruction------anti pion particles
    //vector<StMuTrack*> dau2TracksVector;// for Lambda reconstruction------pion minus particles
    //vector<StMuTrack*> antiDau2TracksVector;// for AntiLambda reconstruction------anti proton particles
    vector<StMuTrack*> pionMinusTracksVector;

    singleTrack=(StMuTrack*) globalTrackIter.Next();
#ifdef DEBUGXI
    if(singleTrack) printf("I get first track!\n");
#endif
    Int_t charge = 0; 
    // For test
    Int_t trackPairNumber = 0;
    Int_t trackNumber = 0;

    StThreeVectorF primaryVertexPosition = muEvent -> primaryVertexPosition();

    //---------------Get the tracks-----------------------------------------------------------------------------------------------------------------
    while(singleTrack=(StMuTrack*) globalTrackIter.Next()){

	trackNumber++;
        //-------------Get the momentum------------
        StThreeVectorF p = singleTrack -> helix().momentum(magneticField*kilogauss);
        //-----------Fill some QA histograms------------------------------------
        hPtRaw -> Fill(singleTrack -> pt());
        hEtaRaw -> Fill(singleTrack -> eta());
        hPhiRaw -> Fill(singleTrack -> phi());
        hNHitsFit -> Fill(singleTrack -> nHitsFit());
        hNHits -> Fill(singleTrack -> nHits());
        hNSigmaProton -> Fill(singleTrack -> nSigmaProton());
        hNSigmaPion -> Fill(singleTrack -> nSigmaPion());
        if (abs(singleTrack -> charge()) == 1 ) hDeDxP -> Fill(p.mag() * (singleTrack -> charge()), singleTrack -> dEdx());
        //----------------------------------------------------------------------
	unsigned short nHits = singleTrack -> nHits();
	if (!pCutNHits(nHits)) continue;
#ifdef DEBUGXI
        printf("pass nHits cut!\n");
#endif

	Double_t pt = singleTrack -> pt();
	if (!pCutPt(pt)) continue;
#ifdef DEBUGXI
        printf("pass pt cut!\n");
#endif
	int trackFlag = singleTrack -> flag();  
	if (trackFlag <= 0) continue; //escape if the tof info is not available
#ifdef DEBUGXI
        printf("pass trackflag cut!\n");
#endif

#ifdef DEBUGXI
	printf("trackNumber=%d\n", trackNumber);
#endif

	//  if (trackNumber<400) continue;
#ifdef DEBUGXI
	printf("Another Track!\n"); 
#endif 


	charge = singleTrack -> charge();
	if( abs(charge) != 1 ){
#ifdef DEBUGXI
	    printf("Charge Not +-1!\n"); 
#endif
	    continue;
	}


#ifdef DEBUGXI
	else printf("Charge +-1!\n");
#endif
	//StThreeVectorF primaryVertexPosition = muEvent->primaryVertexPosition();
	if(1){ //TODO: Include Lambda and AntiLambda
#ifdef DEBUGXI
	    printf("kLambda begins\n");
#endif
	    StPhysicalHelixD helix = singleTrack -> helix();
	    Double_t pathlength = helix.pathLength( primaryVertexPosition, false );
	    StThreeVectorF singleHelixDca = helix.at( pathlength ) - primaryVertexPosition;

	    // Select pionMinus tracks    
	    if(charge == -1){
		StThreeVectorD singleTrackDca2PV = GetTrackDca2PV( helix, primaryVertexPosition);
		//applyCutsSingleHelix();
		if(!pCutPionMinusNSigma(singleTrack)) {
#ifdef DEBUGXI
		    printf("Not a Pion Track!\n"); 
#endif
		    continue;
		}

		// To eliminate tracks directly from primary vertex
		if(!pCutDau2Dca2PV(singleTrackDca2PV.mag())) {
#ifdef DEBUGXI
		    printf("Primary Pion Track!\n"); 
#endif
		    continue;
		}

		    pionMinusTracksVector.push_back(singleTrack);

		}

	    } 
	    //singleTrack=(StMuTrack*) globalTrackIter.Next();
#ifdef DEBUGXI
	    printf("Get new track!\n ");
#endif
	    // For test
	    trackPairNumber++;
#ifdef DEBUGXI
	    printf("cumulative track pair number is %d \n", trackPairNumber);
#endif

	}
#ifdef DEBUGXI
    printf("Now reconstruction begins\n");
#endif

    // Get excited! Let's  reconstruct Xi particles!
    // Pair any two tracks from dau1TracksVector and dau2TracksVector
    // Get size of proton  and  pionminus tracks vector 

//----------------------------------Lambda Reconstruction-------------------------------------------
    //Int_t sizeOfProtonTrVec = dau1TracksVector.size();
    Int_t sizeOfPionMinusTrVec = pionMinusTracksVector.size();
    Int_t numberOfV0 = v0Dst.numberOfV0;

    for(Int_t i = 0; i < numberOfV0; i++){
        //eliminate anti-Lambda particles
        if(v0Dst.v0Type[i] == 0) continue;
        //reconstruct v0(Lambda) particle's pv and momentum
        StThreeVectorF  pv0(v0Dst.v0Px[i], v0Dst.v0Py[i], v0Dst.v0Pz[i]);
        StThreeVectorF  xv0(v0Dst.v0X[i], v0Dst.v0Y[i], v0Dst.v0Z[i]); 
        //reconstruct the v0 straight helix using a large momentum
        StPhysicalHelixD v0Helix(pv0/pv0.perp()*100, xv0, magneticField*kilogauss, 1);
         
	//StPhysicalHelixD dau1Helix = dau1TracksVector[i] -> helix();

	for(Int_t j = 0; j < sizeOfPionMinusTrVec; j++) {

            if(pionMinusTracksVector[j] -> id() == v0Dst.dau2TrackId[i]) continue;
	    StPhysicalHelixD pionMinusHelix = pionMinusTracksVector[j] -> helix();
     

	    //std::pair<StMuTrack*, StMuTrack*> pairTracks(dau1TracksVector[i], dau2TracksVector[j]);
	    //std::pair<Double_t, Double_t> pathlengthPair=dau1Helix.pathLengths(dau2Helix);// return value in centimenters

	    // Apply Long's code: helix1, helix2, magnetic field, pv------event's primary vertex,xXi------xi's position, op1------dau1's momentum at xi position, op2------......

	    StThreeVectorF xXi, op1, op2;
	    Double_t dcaTwoTracks = closestDistance(v0Helix, pionMinusHelix, magneticField, primaryVertexPosition, xXi, op1, op2); 

	    StThreeVectorF xXi2PV = xXi - primaryVertexPosition;
	    StThreeVectorD pXi = op1 + op2;

	    Double_t rDotP = xXi2PV.dot(pXi);
	    Double_t dcaXi2PV = sqrt(xXi2PV.mag2() -(rDotP * rDotP)/(pXi * pXi));
            StPhysicalHelixD helixXi(pXi, xXi, 0, 0);//TODO: check the initialization.

	    if(!pCutDcaXi2PV(dcaXi2PV)) continue;
#ifdef DEBUGXI
            printf("pass xi dca cut!\n");
#endif
	    if(!pCutDecLength(xXi2PV.mag())) continue;
#ifdef DEBUGXI
            printf("pass xi decay length cut!\n");
#endif
	    if(pCutTwoTracksDca(dcaTwoTracks)) {
#ifdef DEBUGXI
		printf("good pair!\n");
#endif
                //-----------Record some basic information---------------------------------------------
                StThreeVectorF v0Dca2PV = GetTrackDca2PV(v0Helix, primaryVertexPosition);
                StThreeVectorF pionMinusDca2PV = GetTrackDca2PV(pionMinusHelix, primaryVertexPosition);
                //-------------------------------------------------
		//pairTracksVector.push_back(pairTracks);

		std::pair<StThreeVectorF, StThreeVectorF> momentumPair(op1, op2);
		//momentumPairVector.push_back(momentumPair);

		Double_t v0E = sqrt(v0Dst.v0Mass[i]*v0Dst.v0Mass[i] + momentumPair.first.mag2());
		Double_t pionMinusE = sqrt(mPionMinusMass * mPionMinusMass + momentumPair.second.mag2());

		StThreeVectorF momentumSumTwoTracks = momentumPair.first + momentumPair.second;// what is the difference between helix's and track's momenta
		Double_t xiMinusInvariantMass = sqrt((v0E + pionMinusE) * (v0E + pionMinusE) - momentumSumTwoTracks.mag2());
#ifdef DEBUGXI
		printf("Mass Calculation OK!\n");
#endif

                //---------Fill the invariant mass spectrum
                //hXiPt -> Fill(pXi.perp());
                hXiMinusInvariantMass -> Fill(xiMinusInvariantMass);
#ifdef DEBUGXI
		printf("Histo Filling OK!\n");
#endif

		//------------------- Assign the values to the mXiDst, 1 -> Lambda---------------------------------------
		mXiDst.xiType[numberOfXi] = 1;
#ifdef DEBUGXI
		printf("xiType OK!\n");
#endif

		mXiDst.xiMass[numberOfXi] = xiMinusInvariantMass;
                mXiDst.xiRapidity[numberOfXi] = TMath::Log((sqrt(xiMinusInvariantMass * xiMinusInvariantMass + pXi.mag2()))/(xiMinusInvariantMass * xiMinusInvariantMass + pXi.perp2()));
                mXiDst.xiEta[numberOfXi] = 0.5 * TMath::Log((pXi.mag()+pXi.z())/(pXi.mag()-pXi.z()));
		mXiDst.xiPt[numberOfXi] = pXi.perp();
		mXiDst.xiX[numberOfXi] = xXi.x();
		mXiDst.xiY[numberOfXi] = xXi.y();
		mXiDst.xiZ[numberOfXi] = xXi.z();
		mXiDst.xiPx[numberOfXi] = pXi.x();
		mXiDst.xiPy[numberOfXi] = pXi.y();
		mXiDst.xiPz[numberOfXi] = pXi.z();
#ifdef DEBUGXI
		printf("Dst xi filling OK!\n");
#endif

		mXiDst.xiDecLength[numberOfXi] = xXi2PV.mag();
		mXiDst.xiDca2PV[numberOfXi] = dcaXi2PV;
                mXiDst.xiDca2PV2d[numberOfXi] = helixXi.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
		mXiDst.xiDca2Dau[numberOfXi] = dcaTwoTracks;
               

                //mXiDst.v0TrackId[numberOfXi] = v0TracksVector[i] -> id(); 
                //mXiDst.v0Charge[numberOfXi] = v0TracksVector[i] -> charge();
                //mXiDst.v0Dca2PV[numberOfXi] = v0dca2PV.mag();
                //mXiDst.v0Dca2PV2d[numberOfXi] = v0Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                //mXiDst.v0NHits[numberOfXi] = v0TracksVector[i] -> nHits();
                //mXiDst.v0NSigmaPion[numberOfXi] = v0TracksVector[i] -> nSigmaPion();
                //mXiDst.v0NSigmaProton[numberOfXi] = v0TracksVector[i] -> nSigmaProton();
                mXiDst.v0Pt[numberOfXi] = op1.perp();
                mXiDst.v0Px[numberOfXi] = op1.x();
                mXiDst.v0Py[numberOfXi] = op1.y();
                mXiDst.v0Pz[numberOfXi] = op1.z();
                //mXiDst.v0NHitsTpc[numberOfXi] = v0Dst.v0[i] -> nHitsFit(kTpcId);
                //mXiDst.v0NHitsSsd[numberOfXi] = v0TracksVector[i] -> nHitsFit(kSsdId);
                //mXiDst.v0NHitsSvt[numberOfXi] = v0TracksVector[i] -> nHitsFit(kSvtId);
                //mXiDst.v0TofFlag[numberOfXi] = v0TracksVector[i] -> btofPidTraits().matchFlag();
                //mXiDst.v0Tof[numberOfXi] = v0TracksVector[i] -> btofPidTraits().timeOfFlight();
#ifdef DEBUGXI
		printf("Dst xi v0 filling OK!\n");
#endif
             
                mXiDst.pionMinusTrackId[numberOfXi] = pionMinusTracksVector[j] -> id(); 
                mXiDst.pionMinusCharge[numberOfXi] = pionMinusTracksVector[j] -> charge();
                mXiDst.pionMinusDca2PV[numberOfXi] = pionMinusDca2PV.mag();
                mXiDst.pionMinusDca2PV2d[numberOfXi] = pionMinusHelix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mXiDst.pionMinusNHits[numberOfXi] = pionMinusTracksVector[j] -> nHits();
                mXiDst.pionMinusNSigmaPion[numberOfXi] = pionMinusTracksVector[j] -> nSigmaPion();
                //mXiDst.dau2NSigmaProton[numberOfXi] = dau2TracksVector[j] -> nSigmaProton();
                mXiDst.pionMinusPt[numberOfXi] = op2.perp();
                mXiDst.pionMinusPx[numberOfXi] = op2.x();
                mXiDst.pionMinusPy[numberOfXi] = op2.y();
                mXiDst.pionMinusPz[numberOfXi] = op2.z();
                mXiDst.pionMinusNHitsTpc[numberOfXi] = pionMinusTracksVector[j] -> nHitsFit(kTpcId);
                mXiDst.pionMinusNHitsSsd[numberOfXi] = pionMinusTracksVector[j] -> nHitsFit(kSsdId);
                mXiDst.pionMinusNHitsSvt[numberOfXi] = pionMinusTracksVector[j] -> nHitsFit(kSvtId);
                mXiDst.pionMinusTofFlag[numberOfXi] = pionMinusTracksVector[j] -> btofPidTraits().matchFlag();
                mXiDst.pionMinusTof[numberOfXi] = pionMinusTracksVector[j] -> btofPidTraits().timeOfFlight();

		numberOfXi++;
#ifdef DEBUGXI
                printf("xi picodst write ok!\n");
#endif
              
	    }

	}

    }

/*-------------------------AntiLambda Reconstruction--------------------------------------------
    Double_t sizeOfPionPlusTrVec = antiDau1TracksVector.size();
    Double_t sizeOfAntiProtonTrVec = antiDau2TracksVector.size();

    for(Int_t i = 0; i < sizeOfPionPlusTrVec; i++){

	StPhysicalHelixD antiDau1Helix=antiDau1TracksVector[i]->helix();

	for(Int_t j = 0; j < sizeOfAntiProtonTrVec; j++) {

	    StPhysicalHelixD antiDau2Helix = antiDau2TracksVector[j] -> helix();

	    std::pair<StMuTrack*, StMuTrack*> pairTracks(antiDau1TracksVector[i], antiDau2TracksVector[j]);
	    std::pair<Double_t, Double_t> pathlengthPair = antiDau1Helix.pathLengths(antiDau2Helix);// return value in centimenters

	    //   StThreeVectorD dcaTwoTracks=v0Helix.at(pathlengthPair.first)-dau2Helix.at(pathlengthPair.second);

	    //      std::vector<std::pair<StMuTrack*, StMuTrack*> > pairTracksVector;
	    //      std::vector<std::pair<Double_t, Double_t> > pathlengthPairVector;
	    // Apply Long's code: helix1, helix2, magnetic field, pv------event's primary vertex,xxi------xi's position, op1------v0's momentum at xi position, op2------......

	    StThreeVectorF xxi, op1, op2;
	    Double_t dcaTwoTracks = closestDistance(antiDau1Helix, antiDau2Helix, mXiDst.magneticField, primaryVertexPosition, xxi, op1, op2); 

	    StThreeVectorF xXi2PV = xxi - primaryVertexPosition;
	    StThreeVectorD pXi = op1 + op2;

	    Double_t rDotP = xXi2PV.dot(pXi);
	    Double_t dcaXi2PV = sqrt(xXi2PV.mag2() -(rDotP * rDotP)/(pXi * pXi));
            StPhysicalHelixD helixxi(pXi, xxi, 0, 0);//TODO: check the initialization.

	    if(!pCutDcaXi2PV(dcaXi2PV)) continue;
	    if(!pCutDecLength(xXi2PV.mag())) continue;
	    if(pCutTwoTracksDca(dcaTwoTracks)) {
#ifdef DEBUGXI
		printf("good anti pair!\n");
#endif 
		antiPairTracksVector.push_back(pairTracks);

                //-----------Record some basic information---------------------------------------------
                StThreeVectorF dau1dca2PV = GetTrackDca2PV(antiDau1Helix, primaryVertexPosition);
                StThreeVectorF dau2dca2PV = GetTrackDca2PV(antiDau2Helix, primaryVertexPosition);
                //-------------------------------------------------

		antiPairTracksVector.push_back(pairTracks);

		std::pair<StThreeVectorF, StThreeVectorF> antiMomentumPair(op1, op2);
		antiMomentumPairVector.push_back(antiMomentumPair);

		Double_t protonE = sqrt(mProtonMass * mProtonMass + antiMomentumPair.second.mag2());
		Double_t pionE = sqrt(mPionPlusMass * mPionPlusMass + antiMomentumPair.first.mag2());

		StThreeVectorF momentumSumTwoTracks = antiMomentumPair.first + antiMomentumPair.second;// what is the difference between helix's and track's momenta
		Double_t antiXiInvariantMass = sqrt((protonE + pionE) * (protonE + pionE) - momentumSumTwoTracks.mag2());
                // Fill the antiparticle invariant mass spectrum
                //hXiPt -> Fill(pXi.perp());
                hAntiXiInvariantMass -> Fill(antiXiInvariantMass);
		// Assign the values to the mXiDst, 0 -> AntiLambda
		mXiDst.xiType[numberOfXi] = 0;
		mXiDst.xiMass[numberOfXi] = antiXiInvariantMass;
                mXiDst.xiRapidity[numberOfXi] = TMath::Log((sqrt(antiXiInvariantMass * antiXiInvariantMass + pXi.mag2()))/(antiXiInvariantMass * antiXiInvariantMass + pXi.perp2()));
                mXiDst.xiEta[numberOfXi] = 0.5 * TMath::Log((pXi.mag()+pXi.z())/(pXi.mag()-pXi.z()));
		mXiDst.xiPt[numberOfXi] = pXi.perp();
		mXiDst.xiX[numberOfXi] = xxi.x();
		mXiDst.xiY[numberOfXi] = xxi.y();
		mXiDst.xiZ[numberOfXi] = xxi.z();
		mXiDst.xiPx[numberOfXi] = pXi.x();
		mXiDst.xiPy[numberOfXi] = pXi.y();
		mXiDst.xiPz[numberOfXi] = pXi.z();

		mXiDst.xiDecLength[numberOfXi] = xXi2PV.mag();
		mXiDst.xiDca2PV[numberOfXi] = dcaXi2PV;
                mXiDst.xiDca2PV2d[numberOfXi] = helixxi.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
		mXiDst.xiDca2Dau[numberOfXi] = dcaTwoTracks;
               

                mXiDst.dau1TrackId[numberOfXi] = antiDau1TracksVector[i] -> id(); 
                mXiDst.dau1Charge[numberOfXi] = antiDau1TracksVector[i] -> charge();
                mXiDst.dau1Dca2PV[numberOfXi] = dau1dca2PV.mag();
                mXiDst.dau1Dca2PV2d[numberOfXi] = antiDau1Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mXiDst.dau1NHits[numberOfXi] = antiDau1TracksVector[i] -> nHits();
                mXiDst.dau1NSigmaPion[numberOfXi] = antiDau1TracksVector[i] -> nSigmaPion();
                mXiDst.dau1NSigmaProton[numberOfXi] = antiDau1TracksVector[i] -> nSigmaProton();
                mXiDst.dau1Pt[numberOfXi] = op1.perp();
                mXiDst.dau1Px[numberOfXi] = op1.x();
                mXiDst.dau1Py[numberOfXi] = op1.y();
                mXiDst.dau1Pz[numberOfXi] = op1.z();
                mXiDst.dau1NHitsTpc[numberOfXi] = antiDau1TracksVector[i] -> nHitsFit(kTpcId);
                mXiDst.dau1NHitsSsd[numberOfXi] = antiDau1TracksVector[i] -> nHitsFit(kSsdId);
                mXiDst.dau1NHitsSvt[numberOfXi] = antiDau1TracksVector[i] -> nHitsFit(kSvtId);
                mXiDst.dau1TofFlag[numberOfXi] = antiDau1TracksVector[i] -> btofPidTraits().matchFlag();
                mXiDst.dau1Tof[numberOfXi] = antiDau1TracksVector[i] -> btofPidTraits().timeOfFlight();
              
                mXiDst.dau2TrackId[numberOfXi] = antiDau2TracksVector[j] -> id(); 
                mXiDst.dau2Charge[numberOfXi] = antiDau2TracksVector[j] -> charge();
                mXiDst.dau2Dca2PV[numberOfXi] = dau2dca2PV.mag();
                mXiDst.dau2Dca2PV2d[numberOfXi] = antiDau2Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mXiDst.dau2NHits[numberOfXi] = antiDau2TracksVector[j] -> nHits();
                mXiDst.dau2NSigmaPion[numberOfXi] = antiDau2TracksVector[j] -> nSigmaPion();
                mXiDst.dau2NSigmaProton[numberOfXi] = antiDau2TracksVector[j] -> nSigmaProton();
                mXiDst.dau2Pt[numberOfXi] = op2.perp();
                mXiDst.dau2Px[numberOfXi] = op2.x();
                mXiDst.dau2Py[numberOfXi] = op2.y();
                mXiDst.dau2Pz[numberOfXi] = op2.z();
                mXiDst.dau2NHitsTpc[numberOfXi] = antiDau2TracksVector[j] -> nHitsFit(kTpcId);
                mXiDst.dau2NHitsSsd[numberOfXi] = antiDau2TracksVector[j] -> nHitsFit(kSsdId);
                mXiDst.dau2NHitsSvt[numberOfXi] = antiDau2TracksVector[j] -> nHitsFit(kSvtId);
                mXiDst.dau2TofFlag[numberOfXi] = antiDau2TracksVector[j] -> btofPidTraits().matchFlag();
                mXiDst.dau2Tof[numberOfXi] = antiDau2TracksVector[j] -> btofPidTraits().timeOfFlight();

		numberOfXi++;
	    }
	}
    }


    //  Int_t sizeOfPairTracksVector = pairTracksVector.size();
    //  Int_t sizeOfAntiPairTracksVector = antiPairTracksVector.size();
#ifdef DEBUGXI
    //printf("the number of tracks passing dca cut is %d\n", sizeOfPairTracksVector); 
#endif

      for(Int_t i = 0; i < sizeOfPairTracksVector; i++){
    //StThreeVectorF momentumProton=pairTracksVector[i].first->helix().momentumAt(pathlengthPairVector[i].first, magneticField);
    //StThreeVectorF momentumPion=pairTracksVector[i].second->helix().momentumAt(pathlengthPairVector[i].second, magneticField);


    //Double_t protonE=sqrt(mProtonMass*mProtonMass+momentumProton.dot(momentumProton));
    Double_t protonE = sqrt(mProtonMass * mProtonMass + (momentumPairVector[i].first).mag2());
    Double_t pionE = sqrt(mPionMass * mPionMass + (momentumPairVector[i].second).mag2());

    StThreeVectorF momentumSumTwoTracks = momentumPairVector[i].first + momentumPairVector[i].second;// what is the difference between helix's and track's momenta
    //StThreeVectorF momentumSumTwoTracks=momentumProton+momentumPion;// what is the difference between helix's and track's momenta
    Double_t invXiMass = sqrt((protonE + pionE) * (protonE + pionE) - momentumSumTwoTracks.mag2());
    hInvariantXiMass -> Fill(invXiMass);

#ifdef DEBUGXI
printf("Mass is %f\n", invXiMass); 
#endif
}

    / for(Int_t i = 0; i < sizeOfAntiPairTracksVector; i++){
    //StThreeVectorF momentumProton=pairTracksVector[i].first->helix().momentumAt(pathlengthPairVector[i].first, pionMinusmagneticField);
    //StThreeVectorF momentumPion=pairTracksVector[i].second->helix().momentumAt(pathlengthPairVector[i].second, magneticField);


    //Double_t protonE=sqrt(mAntiProtonMass*mAntiProtonMass+momentumProton.dot(momentumProton));
    Double_t antiProtonE=sqrt(mAntiProtonMass*mAntiProtonMass+(antiMomentumPairVector[i].first).mag2());
    Double_t pionPlusE=sqrt(mPionPlusMass*mPionPlusMass+(antiMomentumPairVector[i].second).mag2());

    StThreeVectorF momentumSumTwoTracks=momentumPairVector[i].first+momentumPairVector[i].second;// what is the difference between helix's and track's momenta
    //StThreeVectorF momentumSumTwoTracks=momentumProton+momentumPion;// what is the difference between helix's and track's momenta
    Double_t invXiMass=sqrt((protonE+pionE)*(protonE+pionE)-momentumSumTwoTracks.mag2());
    hInvariantXiMass->Fill(invXiMass);

#ifdef DEBUGXI
printf("Mass is %f\n", invXiMass); 
#endif
}
*/
    mXiDst.numberOfXi = numberOfXi;

    mXiTree -> Fill();
    return kStOK; 
}

Int_t StXiMaker::Finish(){

    mOutHistoFile -> Write();
#ifdef DEBUGXI_HISTO
//    gDirectory -> pwd();
//    if(mOutHistoFile -> Write()) printf("Histo write OK!\n");
#endif
#ifdef DEBUGXI_HISTO
//    if(mXiPicoDstFile -> Write()) printf("picoDst write OK!\n");
#endif
    mXiPicoDstFile -> Write();
    printf("File Closed!\n");
    return kStOK;

}

StThreeVectorF StXiMaker::GetTrackDca2PV(StPhysicalHelixD& helix, StThreeVectorF& primaryVertexPosition){

    Double_t pathLength = helix.pathLength( primaryVertexPosition, true); 
    StThreeVectorF dca = helix.at( pathLength) - primaryVertexPosition;
    return dca;
}

Bool_t StXiMaker::SetOutHistoFileName(TString outHistoFileName){ mOutHistoFileName = outHistoFileName; return 1;}

Bool_t StXiMaker::SetOutDSTFileName(TString outDSTFileName){ mXiPicoDstFileName = outDSTFileName; return 1;}

Bool_t StXiMaker::pCutProtonNSigma(StMuTrack* singleTrack){

    Double_t nSigmaProton = singleTrack -> nSigmaProton();   
    if (fabs(nSigmaProton) < mCutNSigmaProtonLe) return true;//kStOK has been defined in StTypeDefs.h
    else return false;
}

Bool_t StXiMaker::pCutAntiProtonNSigma(StMuTrack* singleTrack){
    Double_t nSigmaProton = singleTrack -> nSigmaProton();
    if(fabs(nSigmaProton) < mCutNSigmaProtonLe) return true;
    else return false;
}

Bool_t StXiMaker::pCutPionMinusNSigma(StMuTrack* singleTrack){
    Double_t nSigmaPion = singleTrack -> nSigmaPion();   
    if (fabs(nSigmaPion) < mCutNSigmaPionLe) return true;//kStOK has been defined in StTypeDefs.h
    else return false;
}

Bool_t StXiMaker::pCutPionPlusNSigma(StMuTrack* singleTrack){
    Double_t nSigmaPion = singleTrack -> nSigmaPion();
    if(fabs(nSigmaPion) < mCutNSigmaPionLe) return true;
    else return false;
}

Bool_t StXiMaker::pCutGoodEvents(StMuEvent* event){
    if(!event) return false;
    else return true;
}

Bool_t StXiMaker::pCutEventsTrigger(StMuEvent* event){
    StTriggerId triggerId =event -> triggerIdCollection().nominal();
    if(!triggerId.isTrigger(350003) && !triggerId.isTrigger(350013) && !triggerId.isTrigger(350023) && !triggerId.isTrigger(350033) && !triggerId.isTrigger(350043)) return false;
    else return true;
}

Bool_t StXiMaker::pCutEventsPVZ(StMuEvent* event){
    if(fabs(event -> primaryVertexPosition().z()) > mCutAbsPrimaryVertexZLeEq ) return false;
    else return true;
}

Bool_t StXiMaker::pCutNHits(Int_t nHits){
    if(nHits < mCutNHitsGrEq) return false;
    else return true;
}

Bool_t StXiMaker::pCutPt(Double_t pt){
    if(pt < mCutPtGrEq) return false;
    else return true;
}
 
Bool_t StXiMaker::pCutDau1Dca2PV(Double_t dau1dca2PV){
    if(dau1dca2PV < mCutDau1Dca2PVGrEq) return false;
    else return true;
}

Bool_t StXiMaker::pCutDau2Dca2PV(Double_t dau2dca2PV){
    if(dau2dca2PV < mCutDau2Dca2PVGrEq) return false;
    else return true;
}

Bool_t StXiMaker::pCutTwoTracksDca(Double_t twoTracksDca){
    if(twoTracksDca > mCutTwoTracksDcaLeEq) return false;
    else return true;
}

Bool_t StXiMaker::pCutDecLength(Double_t decayLength){
    if(decayLength < mCutDecLengthGrEq) return false;
    else return true;
}

Bool_t StXiMaker::pCutDcaXi2PV(Double_t xiDca2PV){
    if(xiDca2PV > mCutDcaXi2PVLe) return false;
    else return true;
}


