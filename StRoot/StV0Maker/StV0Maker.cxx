//#define DEBUG
//#define DEBUG_HISTO 1
#include "TTree.h"
#include "TFile.h"
#include "StV0Maker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "TString.h"
#include "StV0Type.h"
#include "StV0Dst.h"
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
#include "TMath.h"

#include "StDcaService.h"
ClassImp(StV0Maker)




StV0Maker::StV0Maker(StMuDstMaker* muDstMaker, 
TString outDSTFileName,
TString outHistoFileName,
StV0Type V0Type): 
mOutHistoFile(NULL),  
mV0PicoDstFile(NULL),
mOutHistoFileName(outHistoFileName), 
mV0PicoDstFileName(outDSTFileName), 
mV0ParticleName(""),
mRotation(0),
mSameSign(0),
mV0Mass(1.11),
mV0MassWidth(0.2), // The default unit of energy or mass in STAR is GeV.
mProtonMass(0.938),
mAntiProtonMass(0.938),
mPionMinusMass(0.140),
mPionPlusMass(0.140),
mBuffer(""), 
mHInvariantMassName(""), 
mV0Type(V0Type),
mGlobalTracks(0),// ??
mPrimaryTracks(0),// ??
mMuDstMaker(muDstMaker),
mEventFailed(0),
mEventPassed(0),
//mV0Dst(),
mV0Tree(NULL),
mCutTriggerIdEq(360001),
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
mCutDcaV02PVLe(5.0){

    std::vector<std::pair<StMuTrack*, StMuTrack*> > pairTracksVector;
    std::vector<std::pair<StMuTrack*, StMuTrack*> > antiPairTracksVector;

    std::vector<std::pair<Double_t, Double_t> > momentumPairVector;
    std::vector<std::pair<Double_t, Double_t> > antiMomentumPairVector;

#ifdef DEBUG_HISTO
    mOutHistoFile = NULL;
    mV0PicoDstFile = NULL;
#endif
} 
  

Int_t  StV0Maker::Init(){

    //-------------------Initialize output files----------------------------------------------------------------------------------------
    mV0PicoDstFile = new TFile(mV0PicoDstFileName, "recreate");       // mV0PicoDstFileName is given by user code
    mOutHistoFile = new TFile(mOutHistoFileName, "recreate");         // mOutHistoFileName is given by user code
#ifdef DEBUG_HISTO 
    printf("Initialization happens right now\n");
#endif

    //-------------------Book Histograms----------------------------------------------------------------------------------------
//    hNPV = new TH1F();
    hVz = new TH1F("VertexZ","Event Vertex Z position", 400, -100, 100);
          //-------------Initialize invariant-mass spectrum histogram 
    if(mV0Type == kLambda) mV0ParticleName = "Lambda";
    mHInvariantMassName.Append(mV0ParticleName);
    mHInvariantMassName.Append("InvariantMass");
    TString antiHInvariantMassName = mHInvariantMassName + "Anti";
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

    hV0InvariantMass = new TH1F(mHInvariantMassName, "Invariant Mass", 300, mV0Mass - mV0MassWidth, mV0Mass + mV0MassWidth);
    hAntiV0InvariantMass = new TH1F(antiHInvariantMassName, "Anti Particle Invariant Mass", 300, mV0Mass - mV0MassWidth, mV0Mass + mV0MassWidth);
#ifdef DEBUG
    printf("histogram OK!\n");
#endif

    //------------Initialize tree-----------------------------------------------------------------------------------

    //mV0PicoDstFileName.Append(mV0ParticleName);    
    //mV0PicoDstFile = new TFile(mV0PicoDstFileName, "recreate");       // mV0PicoDstFileName is given by user code
    //mOutHistoFile = new TFile(mOutHistoFileName, "recreate");      // mOutHistoFileName is given by user code
#ifdef DEBUG_HISTO
    gDirectory -> pwd();
    if(mOutHistoFile) printf("histoFile open OK\n");
    if(mV0PicoDstFile) printf("picoDataFile open OK\n");
#endif
    mV0Tree = new TTree("V0PicoDst", "From MuDst");

    mV0Tree -> SetDirectory(mV0PicoDstFile ); 

    mV0Tree -> Branch("runNumber", &mV0Dst.runNumber, "runNumber/I");
    mV0Tree -> Branch("eventId", &mV0Dst.eventId, "eventId/I");
    mV0Tree -> Branch("triggerMode", &mV0Dst.triggerMode, "triggerMode/I");
    mV0Tree -> Branch("nRefMult", &mV0Dst.nRefMult, "nRefMult/I");
    mV0Tree -> Branch("nRefMultTof", &mV0Dst.nRefMultTof, "nRefMultTof/I");
    mV0Tree -> Branch("magneticField", &mV0Dst.magneticField, "magneticField/F");
    mV0Tree -> Branch("PVX", &mV0Dst.PVX, "PVX/F");
    mV0Tree -> Branch("PVY", &mV0Dst.PVY, "PVY/F");
    mV0Tree -> Branch("PVZ", &mV0Dst.PVZ, "PVZ/F");
    mV0Tree -> Branch("numberOfV0", &mV0Dst.numberOfV0, "numberOfV0/I");

    mV0Tree -> Branch("v0type", mV0Dst.v0Type, "v0Type[numberOfV0]/I");
    mV0Tree -> Branch("v0Mass", mV0Dst.v0Mass, "v0Mass[numberOfV0]/F");
    mV0Tree -> Branch("v0Pt", mV0Dst.v0Pt, "v0Pt[numberOfV0]/F");
    mV0Tree -> Branch("v0Rapidity", mV0Dst.v0Rapidity, "v0Rapidity[numberOfV0]/F");
    mV0Tree -> Branch("v0Eta", mV0Dst.v0Eta, "v0Eta[numberOfV0]/F");
    mV0Tree -> Branch("v0X", mV0Dst.v0X, "v0X[numberOfV0]/F");
    mV0Tree -> Branch("v0Y", mV0Dst.v0Y, "v0Y[numberOfV0]/F");
    mV0Tree -> Branch("v0Z", mV0Dst.v0Z, "v0Z[numberOfV0]/F");
    mV0Tree -> Branch("v0Px", mV0Dst.v0Px, "v0Px[numberOfV0]/F");
    mV0Tree -> Branch("v0Py", mV0Dst.v0Py, "v0Py[numberOfV0]/F");
    mV0Tree -> Branch("v0Pz", mV0Dst.v0Pz, "v0Pz[numberOfV0]/F");
    mV0Tree -> Branch("v0DecLength", mV0Dst.v0DecLength, "v0DecLength[numberOfV0]/F");
    mV0Tree -> Branch("v0Dca2PV", mV0Dst.v0Dca2PV, "v0Dca2PV[numberOfV0]/F");
    mV0Tree -> Branch("v0Dca2PV2d", mV0Dst.v0Dca2PV2d, "v0Dca2PV2d[numberOfV0]/F");
    mV0Tree -> Branch("v0Dca2Dau", mV0Dst.v0Dca2Dau, "v0Dca2Dau[numberOfV0]/F");

    mV0Tree -> Branch("dau1TrackId", mV0Dst.dau1TrackId, "dau1TrackId[numberOfV0]/I");
    mV0Tree -> Branch("dau1Charge", mV0Dst.dau1Charge, "dau1Charge[numberOfV0]/I");
    mV0Tree -> Branch("dau1Dca2PV", mV0Dst.dau1Dca2PV, "dau1Dca2PV[numberOfV0]/F");
    mV0Tree -> Branch("dau1Dca2PV2d", mV0Dst.dau1Dca2PV2d, "dau1Dca2PV2d[numberOfV0]/F");
    mV0Tree -> Branch("dau1NHits", mV0Dst.dau1NHits, "dau1NHits[numberOfV0]/I");
      //mV0Tree -> Branch("dau1DeDx", mV0Dst.dau1DeDx, "dau1DeDx[numberOfV0]/F");
      //mV0Tree -> Branch("dau1DeDx", mV0Dst.dau1DeDx, "dau1DeDx[numberOfV0]/F");
    mV0Tree -> Branch("dau1NSigmaPion", mV0Dst.dau1NSigmaPion, "dau1NSigmaPion[numberOfV0]/F");
    mV0Tree -> Branch("dau1NSigmaProton", mV0Dst.dau1NSigmaProton, "dau1NSigmaProton[numberOfV0]/F");
    mV0Tree -> Branch("dau1Pt", mV0Dst.dau1Pt, "dau1Pt[numberOfV0]/F");
    mV0Tree -> Branch("dau1Px", mV0Dst.dau1Px, "dau1Px[numberOfV0]/F");
    mV0Tree -> Branch("dau1Py", mV0Dst.dau1Py, "dau1Py[numberOfV0]/F");
    mV0Tree -> Branch("dau1Pz", mV0Dst.dau1Pz, "dau1Pz[numberOfV0]/F");
    mV0Tree -> Branch("dau1NHitsTpc", mV0Dst.dau1NHitsTpc, "dau1NHitsTpc[numberOfV0]/I");
    mV0Tree -> Branch("dau1NHitsSsd", mV0Dst.dau1NHitsSsd, "dau1NHitsSsd[numberOfV0]/I");
    mV0Tree -> Branch("dau1NHitsSvt", mV0Dst.dau1NHitsSvt, "dau1NHitsSvt[numberOfV0]/I");
    mV0Tree -> Branch("dau1TofFlag", mV0Dst.dau1TofFlag, "dau1TofFlag[numberOfV0]/I");
    mV0Tree -> Branch("dau1Tof", mV0Dst.dau1Tof, "dau1Tof[numberOfV0]/F");

    mV0Tree -> Branch("dau2TrackId", mV0Dst.dau2TrackId, "dau2TrackId[numberOfV0]/I");
    mV0Tree -> Branch("dau2Charge", mV0Dst.dau2Charge, "dau2Charge[numberOfV0]/I");
    mV0Tree -> Branch("dau2Dca2PV", mV0Dst.dau2Dca2PV, "dau2Dca2PV[numberOfV0]/F");
    mV0Tree -> Branch("dau2Dca2PV2d", mV0Dst.dau2Dca2PV2d, "dau2Dca2PV2d[numberOfV0]/F");
    mV0Tree -> Branch("dau2NHits", mV0Dst.dau2NHits, "dau2NHits[numberOfV0]/I");
      //mV0Tree -> Branch("dau2DeDx", mV0Dst.dau1DeDx, "dau2DeDx[numberOfV0]/F");
      //mV0Tree -> Branch("dau2DeDx", mV0Dst.dau1DeDx, "dau2DeDx[numberOfV0]/F");
    mV0Tree -> Branch("dau2NSigmaPion", mV0Dst.dau2NSigmaPion, "dau2NSigmaPion[numberOfV0]/F");
    mV0Tree -> Branch("dau2NSigmaProton", mV0Dst.dau2NSigmaProton, "dau2NSigmaProton[numberOfV0]/F");
    mV0Tree -> Branch("dau2Pt", mV0Dst.dau2Pt, "dau2Pt[numberOfV0]/F");
    mV0Tree -> Branch("dau2Px", mV0Dst.dau2Px, "dau2Px[numberOfV0]/F");
    mV0Tree -> Branch("dau2Py", mV0Dst.dau2Py, "dau2Py[numberOfV0]/F");
    mV0Tree -> Branch("dau2Pz", mV0Dst.dau2Pz, "dau2Pz[numberOfV0]/F");
    mV0Tree -> Branch("dau2NHitsTpc", mV0Dst.dau2NHitsTpc, "dau2NHitsTpc[numberOfV0]/I");
    mV0Tree -> Branch("dau2NHitsSsd", mV0Dst.dau2NHitsSsd, "dau2NHitsSsd[numberOfV0]/I");
    mV0Tree -> Branch("dau2NHitsSvt", mV0Dst.dau2NHitsSvt, "dau2NHitsSvt[numberOfV0]/I");
    mV0Tree -> Branch("dau2TofFlag", mV0Dst.dau2TofFlag, "dau2TofFlag[numberOfV0]/I");
    mV0Tree -> Branch("dau2Tof", mV0Dst.dau2Tof, "dau2Tof[numberOfV0]/F");

    return kStOK;
}

Int_t StV0Maker::Make() {
#ifdef DEBUG_HISTO
    gDirectory -> pwd();
#endif
#ifdef DEBUG
    printf("New Event!\n");
#endif
#ifdef DEBUG
    printf("Now I am making sth\n");
#endif
    StMuEvent* muEvent = mMuDstMaker -> muDst() -> event(); // I think mMuDstMaker would be automatically move to next event as well as collection of tracks
    
    Double_t magneticField = muEvent -> runInfo().magneticField();

    if(!pCutGoodEvents(muEvent)) {  mEventFailed++; return kStOK; }
#ifdef DEBUG
    printf("event good!\n");
#endif
    if(!pCutEventsTrigger(muEvent)) return kStOK;
#ifdef DEBUG
    printf("trigger good!\n");
#endif
    if(!pCutEventsPVZ(muEvent)) return kStOK;
#ifdef DEBUG
    printf("PVZ good!\n");
#endif

    //---------Event-wise information---------------------------------------------------------------------------------------------------------------
    mV0Dst.eventId = muEvent -> eventId();
    mV0Dst.runNumber = muEvent -> runId();
    mV0Dst.triggerMode = mCutTriggerIdEq; 
    mV0Dst.nRefMult = muEvent -> refMult();
    mV0Dst.nRefMultTof = muEvent -> btofTrayMultiplicity();
    mV0Dst.magneticField = muEvent -> runInfo().magneticField(); 
    mV0Dst.PVX = muEvent -> primaryVertexPosition().x();  
    mV0Dst.PVY = muEvent -> primaryVertexPosition().y();  
    mV0Dst.PVZ = muEvent -> primaryVertexPosition().z();  
    //-----------Fill event-wise information histograms------------------------------------------------------------------------------- 
    hRefMult -> Fill(muEvent -> refMult());//I dont get what does refMult mean
    //------------------------------------------------------------------------------------------------------------------------------------------------
    int numberOfV0 = 0;
    //#ifdef DEBUG
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

    Int_t globalTrackNumber=mGlobalTracks->GetEntries();
#ifdef DEBUG
    printf("the number of total global tracks is %d\n", globalTrackNumber);
#endif
    //if(!StV0Maker.isGlobalTracksOk(muEvent)) return kStOk;
    //if(!StV0Maker.isPrimaryTracksOk(muEvent)) return kStOk;

    TIter globalTrackIter(mGlobalTracks);
    TIter primaryTrackIter(mPrimaryTracks);

    StMuTrack* singleTrack= new StMuTrack;

    vector<StMuTrack*> dau1TracksVector;// for Lambda reconstruction------proton particles
    vector<StMuTrack*> antiDau1TracksVector;// for AntiLambda reconstruction------anti pion particles
    vector<StMuTrack*> dau2TracksVector;// for Lambda reconstruction------pion minus particles
    vector<StMuTrack*> antiDau2TracksVector;// for AntiLambda reconstruction------anti proton particles

    singleTrack=(StMuTrack*) globalTrackIter.Next();
#ifdef DEBUG
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
#ifdef DEBUG
        printf("pass nHits cut!\n");
#endif

	Double_t pt = singleTrack -> pt();
	if (!pCutPt(pt)) continue;
#ifdef DEBUG
        printf("pass pt cut!\n");
#endif
	int trackFlag = singleTrack -> flag();  
	if (trackFlag <= 0) continue; //escape if the tof info is not available
#ifdef DEBUG
        printf("pass trackflag cut!\n");
#endif

#ifdef DEBUG
	printf("trackNumber=%d\n", trackNumber);
#endif

	//  if (trackNumber<400) continue;
#ifdef DEBUG
	printf("Another Track!\n"); 
#endif 


	charge = singleTrack -> charge();
	if( abs(charge) != 1 ){
#ifdef DEBUG
	    printf("Charge Not +-1!\n"); 
#endif
	    continue;
	}


#ifdef DEBUG
	else printf("Charge +-1!\n");
#endif
	//StThreeVectorF primaryVertexPostion=muEvent->primaryVertexPosition();
	if(mV0Type == kLambda ){ // Include Lambda and AntiLambda
#ifdef DEBUG
	    printf("kLambda begins\n");
#endif
	    StPhysicalHelixD helix = singleTrack -> helix();
	    Double_t pathlength = helix.pathLength( primaryVertexPosition, false );
	    StThreeVectorF singleHelixDca = helix.at( pathlength ) - primaryVertexPosition;
	    // Select proton tracks
	    if(charge == 1){
		StThreeVectorD singleTrackDca2PV = GetTrackDca2PV( helix, primaryVertexPosition );         

		if(!pCutProtonNSigma(singleTrack) && !pCutPionPlusNSigma(singleTrack)) { 
#ifdef DEBUG
		    printf("Not a proton!\n");
#endif
		    continue;
		}

		if(!pCutDau1Dca2PV(singleTrackDca2PV.mag())) {
#ifdef DEBUG
		    printf("Primary track!\n");
#endif
		    continue;
		}

		Double_t dau1PidDistance = abs(singleTrack -> nSigmaProton());
		Double_t antiDau1PidDistance = abs(singleTrack -> nSigmaPion()); 

		Double_t candPidDistance = dau1PidDistance - antiDau1PidDistance;

		if (pCutProtonNSigma(singleTrack) && !pCutPionPlusNSigma(singleTrack)) dau1TracksVector.push_back(singleTrack);
		if (pCutPionPlusNSigma(singleTrack) && !pCutProtonNSigma(singleTrack)) antiDau1TracksVector.push_back(singleTrack);
		if ((pCutProtonNSigma(singleTrack)) && (pCutPionPlusNSigma(singleTrack))) {

		    dau1TracksVector.push_back(singleTrack);
		    antiDau1TracksVector.push_back(singleTrack);

		}
	    }

	    // Select pionMinus tracks    
	    if(charge == -1){
		StThreeVectorD singleTrackDca2PV = GetTrackDca2PV( helix, primaryVertexPosition);
		//applyCutsSingleHelix();
		if(!pCutPionMinusNSigma(singleTrack) && !pCutAntiProtonNSigma(singleTrack)) {
#ifdef DEBUG
		    printf("Not a Pion Track!\n"); 
#endif
		    continue;
		}

		// To eliminate tracks directly from primary vertex
		if(!pCutDau2Dca2PV(singleTrackDca2PV.mag())) {
#ifdef DEBUG
		    printf("Primary Pion Track!\n"); 
#endif
		    continue;
		}


		Double_t dau2PidDistance = abs(singleTrack -> nSigmaPion()); 
		Double_t antiDau2PidDistance = abs(singleTrack -> nSigmaProton());

		Double_t candPidDistance = dau2PidDistance - antiDau2PidDistance; 

		if (pCutPionMinusNSigma(singleTrack) && !pCutAntiProtonNSigma(singleTrack)) dau2TracksVector.push_back(singleTrack);
		if (pCutAntiProtonNSigma(singleTrack) && !pCutPionMinusNSigma(singleTrack)) antiDau2TracksVector.push_back(singleTrack);

		if (pCutAntiProtonNSigma(singleTrack) && pCutPionMinusNSigma(singleTrack)) {

		    dau2TracksVector.push_back(singleTrack);
		    antiDau2TracksVector.push_back(singleTrack);

		}

	    } 
	    //singleTrack=(StMuTrack*) globalTrackIter.Next();
#ifdef DEBUG
	    printf("Get new track!\n ");
#endif
	    // For test
	    trackPairNumber++;
#ifdef DEBUG
	    printf("cumulative track pair number is %d \n", trackPairNumber);
#endif

	}
    }
#ifdef DEBUG
    printf("Now reconstruction begins\n");
#endif

    // Get excited! Let's  reconstruct V0 particles!
    // Pair any two tracks from dau1TracksVector and dau2TracksVector
    // Get size of proton  and  pionminus tracks vector 

//----------------------------------Lambda Reconstruction-------------------------------------------
    Double_t sizeOfProtonTrVec = dau1TracksVector.size();
    Double_t sizeOfPionMinusTrVec = dau2TracksVector.size();

    for(Int_t i = 0; i < sizeOfProtonTrVec; i++){

	StPhysicalHelixD dau1Helix = dau1TracksVector[i] -> helix();

	for(Int_t j = 0; j < sizeOfPionMinusTrVec; j++) {

	    StPhysicalHelixD dau2Helix = dau2TracksVector[j] -> helix();

	    std::pair<StMuTrack*, StMuTrack*> pairTracks(dau1TracksVector[i], dau2TracksVector[j]);
	    std::pair<Double_t, Double_t> pathlengthPair=dau1Helix.pathLengths(dau2Helix);// return value in centimenters

	    // Apply Long's code: helix1, helix2, magnetic field, pv------event's primary vertex,xv0------v0's position, op1------dau1's momentum at v0 position, op2------......

	    StThreeVectorF xv0, op1, op2;
	    Double_t dcaTwoTracks = closestDistance(dau1Helix, dau2Helix, mV0Dst.magneticField, primaryVertexPosition, xv0, op1, op2); 

	    StThreeVectorF xV02PV = xv0 - primaryVertexPosition;
	    StThreeVectorD pV0 = op1 + op2;

	    Double_t rDotP = xV02PV.dot(pV0);
	    Double_t dcaV02PV = sqrt(xV02PV.mag2() -(rDotP * rDotP)/(pV0 * pV0));
            StPhysicalHelixD helixv0(pV0, xv0, 0, 0);//TODO: check the initialization.

	    if(!pCutDcaV02PV(dcaV02PV)) continue;
#ifdef DEBUG
            printf("pass v0 dca cut!\n");
#endif
	    if(!pCutDecLength(xV02PV.mag())) continue;
#ifdef DEBUG
            printf("pass v0 decay length cut!\n");
#endif
	    if(pCutTwoTracksDca(dcaTwoTracks)) {
#ifdef DEBUG
		printf("good pair!\n");
#endif
                //-----------Record some basic information---------------------------------------------
                StThreeVectorF dau1dca2PV = GetTrackDca2PV(dau1Helix, primaryVertexPosition);
                StThreeVectorF dau2dca2PV = GetTrackDca2PV(dau2Helix, primaryVertexPosition);
                //-------------------------------------------------
		pairTracksVector.push_back(pairTracks);

		std::pair<StThreeVectorF, StThreeVectorF> momentumPair(op1, op2);
		momentumPairVector.push_back(momentumPair);

		Double_t protonE = sqrt(mProtonMass * mProtonMass + momentumPair.first.mag2());
		Double_t pionE = sqrt(mPionMinusMass * mPionMinusMass + momentumPair.second.mag2());

		StThreeVectorF momentumSumTwoTracks = momentumPair.first + momentumPair.second;// what is the difference between helix's and track's momenta
		Double_t v0InvariantMass = sqrt((protonE + pionE) * (protonE + pionE) - momentumSumTwoTracks.mag2());
#ifdef DEBUG
		printf("Mass Calculation OK!\n");
#endif

                //---------Fill the invariant mass spectrum
                //hV0Pt -> Fill(pV0.perp());
                hV0InvariantMass -> Fill(v0InvariantMass);
#ifdef DEBUG
		printf("Histo Filling OK!\n");
#endif

		//------------------- Assign the values to the mV0Dst, 1 -> Lambda---------------------------------------
		mV0Dst.v0Type[numberOfV0] = 1;
#ifdef DEBUG
		printf("v0Type OK!\n");
#endif

		mV0Dst.v0Mass[numberOfV0] = v0InvariantMass;
                mV0Dst.v0Rapidity[numberOfV0] = TMath::Log((sqrt(v0InvariantMass * v0InvariantMass + pV0.mag2()))/(v0InvariantMass * v0InvariantMass + pV0.perp2()));
                mV0Dst.v0Eta[numberOfV0] = 0.5 * TMath::Log((pV0.mag()+pV0.z())/(pV0.mag()-pV0.z()));
		mV0Dst.v0Pt[numberOfV0] = pV0.perp();
		mV0Dst.v0X[numberOfV0] = xv0.x();
		mV0Dst.v0Y[numberOfV0] = xv0.y();
		mV0Dst.v0Z[numberOfV0] = xv0.z();
		mV0Dst.v0Px[numberOfV0] = pV0.x();
		mV0Dst.v0Py[numberOfV0] = pV0.y();
		mV0Dst.v0Pz[numberOfV0] = pV0.z();
#ifdef DEBUG
		printf("Dst v0 filling OK!\n");
#endif

		mV0Dst.v0DecLength[numberOfV0] = xV02PV.mag();
		mV0Dst.v0Dca2PV[numberOfV0] = dcaV02PV;
                mV0Dst.v0Dca2PV2d[numberOfV0] = helixv0.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
		mV0Dst.v0Dca2Dau[numberOfV0] = dcaTwoTracks;
               

                mV0Dst.dau1TrackId[numberOfV0] = dau1TracksVector[i] -> id(); 
                mV0Dst.dau1Charge[numberOfV0] = dau1TracksVector[i] -> charge();
                mV0Dst.dau1Dca2PV[numberOfV0] = dau1dca2PV.mag();
                mV0Dst.dau1Dca2PV2d[numberOfV0] = dau1Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mV0Dst.dau1NHits[numberOfV0] = dau1TracksVector[i] -> nHits();
                mV0Dst.dau1NSigmaPion[numberOfV0] = dau1TracksVector[i] -> nSigmaPion();
                mV0Dst.dau1NSigmaProton[numberOfV0] = dau1TracksVector[i] -> nSigmaProton();
                mV0Dst.dau1Pt[numberOfV0] = op1.perp();
                mV0Dst.dau1Px[numberOfV0] = op1.x();
                mV0Dst.dau1Py[numberOfV0] = op1.y();
                mV0Dst.dau1Pz[numberOfV0] = op1.z();
                mV0Dst.dau1NHitsTpc[numberOfV0] = dau1TracksVector[i] -> nHitsFit(kTpcId);
                mV0Dst.dau1NHitsSsd[numberOfV0] = dau1TracksVector[i] -> nHitsFit(kSsdId);
                mV0Dst.dau1NHitsSvt[numberOfV0] = dau1TracksVector[i] -> nHitsFit(kSvtId);
                mV0Dst.dau1TofFlag[numberOfV0] = dau1TracksVector[i] -> btofPidTraits().matchFlag();
                mV0Dst.dau1Tof[numberOfV0] = dau1TracksVector[i] -> btofPidTraits().timeOfFlight();
#ifdef DEBUG
		printf("Dst v0 dau1 filling OK!\n");
#endif
             
                mV0Dst.dau2TrackId[numberOfV0] = dau2TracksVector[j] -> id(); 
                mV0Dst.dau2Charge[numberOfV0] = dau2TracksVector[j] -> charge();
                mV0Dst.dau2Dca2PV[numberOfV0] = dau2dca2PV.mag();
                mV0Dst.dau2Dca2PV2d[numberOfV0] = dau2Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mV0Dst.dau2NHits[numberOfV0] = dau2TracksVector[j] -> nHits();
                mV0Dst.dau2NSigmaPion[numberOfV0] = dau2TracksVector[j] -> nSigmaPion();
                mV0Dst.dau2NSigmaProton[numberOfV0] = dau2TracksVector[j] -> nSigmaProton();
                mV0Dst.dau2Pt[numberOfV0] = op2.perp();
                mV0Dst.dau2Px[numberOfV0] = op2.x();
                mV0Dst.dau2Py[numberOfV0] = op2.y();
                mV0Dst.dau2Pz[numberOfV0] = op2.z();
                mV0Dst.dau2NHitsTpc[numberOfV0] = dau2TracksVector[j] -> nHitsFit(kTpcId);
                mV0Dst.dau2NHitsSsd[numberOfV0] = dau2TracksVector[j] -> nHitsFit(kSsdId);
                mV0Dst.dau2NHitsSvt[numberOfV0] = dau2TracksVector[j] -> nHitsFit(kSvtId);
                mV0Dst.dau2TofFlag[numberOfV0] = dau2TracksVector[j] -> btofPidTraits().matchFlag();
                mV0Dst.dau2Tof[numberOfV0] = dau2TracksVector[j] -> btofPidTraits().timeOfFlight();

		numberOfV0++;
#ifdef DEBUG
                printf("v0 picodst write ok!\n");
#endif
              
	    }

	}

    }

//-------------------------AntiLambda Reconstruction--------------------------------------------
    Double_t sizeOfPionPlusTrVec = antiDau1TracksVector.size();
    Double_t sizeOfAntiProtonTrVec = antiDau2TracksVector.size();

    for(Int_t i = 0; i < sizeOfPionPlusTrVec; i++){

	StPhysicalHelixD antiDau1Helix=antiDau1TracksVector[i]->helix();

	for(Int_t j = 0; j < sizeOfAntiProtonTrVec; j++) {

	    StPhysicalHelixD antiDau2Helix = antiDau2TracksVector[j] -> helix();

	    std::pair<StMuTrack*, StMuTrack*> pairTracks(antiDau1TracksVector[i], antiDau2TracksVector[j]);
	    std::pair<Double_t, Double_t> pathlengthPair = antiDau1Helix.pathLengths(antiDau2Helix);// return value in centimenters

	    //   StThreeVectorD dcaTwoTracks=dau1Helix.at(pathlengthPair.first)-dau2Helix.at(pathlengthPair.second);

	    //      std::vector<std::pair<StMuTrack*, StMuTrack*> > pairTracksVector;
	    //      std::vector<std::pair<Double_t, Double_t> > pathlengthPairVector;
	    // Apply Long's code: helix1, helix2, magnetic field, pv------event's primary vertex,xv0------v0's position, op1------dau1's momentum at v0 position, op2------......

	    StThreeVectorF xv0, op1, op2;
	    Double_t dcaTwoTracks = closestDistance(antiDau1Helix, antiDau2Helix, mV0Dst.magneticField, primaryVertexPosition, xv0, op1, op2); 

	    StThreeVectorF xV02PV = xv0 - primaryVertexPosition;
	    StThreeVectorD pV0 = op1 + op2;

	    Double_t rDotP = xV02PV.dot(pV0);
	    Double_t dcaV02PV = sqrt(xV02PV.mag2() -(rDotP * rDotP)/(pV0 * pV0));
            StPhysicalHelixD helixv0(pV0, xv0, 0, 0);//TODO: check the initialization.

	    if(!pCutDcaV02PV(dcaV02PV)) continue;
	    if(!pCutDecLength(xV02PV.mag())) continue;
	    if(pCutTwoTracksDca(dcaTwoTracks)) {
#ifdef DEBUG
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
		Double_t antiV0InvariantMass = sqrt((protonE + pionE) * (protonE + pionE) - momentumSumTwoTracks.mag2());
                // Fill the antiparticle invariant mass spectrum
                //hV0Pt -> Fill(pV0.perp());
                hAntiV0InvariantMass -> Fill(antiV0InvariantMass);
		// Assign the values to the mV0Dst, 0 -> AntiLambda
		mV0Dst.v0Type[numberOfV0] = 0;
		mV0Dst.v0Mass[numberOfV0] = antiV0InvariantMass;
                mV0Dst.v0Rapidity[numberOfV0] = TMath::Log((sqrt(antiV0InvariantMass * antiV0InvariantMass + pV0.mag2()))/(antiV0InvariantMass * antiV0InvariantMass + pV0.perp2()));
                mV0Dst.v0Eta[numberOfV0] = 0.5 * TMath::Log((pV0.mag()+pV0.z())/(pV0.mag()-pV0.z()));
		mV0Dst.v0Pt[numberOfV0] = pV0.perp();
		mV0Dst.v0X[numberOfV0] = xv0.x();
		mV0Dst.v0Y[numberOfV0] = xv0.y();
		mV0Dst.v0Z[numberOfV0] = xv0.z();
		mV0Dst.v0Px[numberOfV0] = pV0.x();
		mV0Dst.v0Py[numberOfV0] = pV0.y();
		mV0Dst.v0Pz[numberOfV0] = pV0.z();

		mV0Dst.v0DecLength[numberOfV0] = xV02PV.mag();
		mV0Dst.v0Dca2PV[numberOfV0] = dcaV02PV;
                mV0Dst.v0Dca2PV2d[numberOfV0] = helixv0.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
		mV0Dst.v0Dca2Dau[numberOfV0] = dcaTwoTracks;
               

                mV0Dst.dau1TrackId[numberOfV0] = antiDau1TracksVector[i] -> id(); 
                mV0Dst.dau1Charge[numberOfV0] = antiDau1TracksVector[i] -> charge();
                mV0Dst.dau1Dca2PV[numberOfV0] = dau1dca2PV.mag();
                mV0Dst.dau1Dca2PV2d[numberOfV0] = antiDau1Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mV0Dst.dau1NHits[numberOfV0] = antiDau1TracksVector[i] -> nHits();
                mV0Dst.dau1NSigmaPion[numberOfV0] = antiDau1TracksVector[i] -> nSigmaPion();
                mV0Dst.dau1NSigmaProton[numberOfV0] = antiDau1TracksVector[i] -> nSigmaProton();
                mV0Dst.dau1Pt[numberOfV0] = op1.perp();
                mV0Dst.dau1Px[numberOfV0] = op1.x();
                mV0Dst.dau1Py[numberOfV0] = op1.y();
                mV0Dst.dau1Pz[numberOfV0] = op1.z();
                mV0Dst.dau1NHitsTpc[numberOfV0] = antiDau1TracksVector[i] -> nHitsFit(kTpcId);
                mV0Dst.dau1NHitsSsd[numberOfV0] = antiDau1TracksVector[i] -> nHitsFit(kSsdId);
                mV0Dst.dau1NHitsSvt[numberOfV0] = antiDau1TracksVector[i] -> nHitsFit(kSvtId);
                mV0Dst.dau1TofFlag[numberOfV0] = antiDau1TracksVector[i] -> btofPidTraits().matchFlag();
                mV0Dst.dau1Tof[numberOfV0] = antiDau1TracksVector[i] -> btofPidTraits().timeOfFlight();
              
                mV0Dst.dau2TrackId[numberOfV0] = antiDau2TracksVector[j] -> id(); 
                mV0Dst.dau2Charge[numberOfV0] = antiDau2TracksVector[j] -> charge();
                mV0Dst.dau2Dca2PV[numberOfV0] = dau2dca2PV.mag();
                mV0Dst.dau2Dca2PV2d[numberOfV0] = antiDau2Helix.geometricSignedDistance(primaryVertexPosition.x(), primaryVertexPosition.y());
                mV0Dst.dau2NHits[numberOfV0] = antiDau2TracksVector[j] -> nHits();
                mV0Dst.dau2NSigmaPion[numberOfV0] = antiDau2TracksVector[j] -> nSigmaPion();
                mV0Dst.dau2NSigmaProton[numberOfV0] = antiDau2TracksVector[j] -> nSigmaProton();
                mV0Dst.dau2Pt[numberOfV0] = op2.perp();
                mV0Dst.dau2Px[numberOfV0] = op2.x();
                mV0Dst.dau2Py[numberOfV0] = op2.y();
                mV0Dst.dau2Pz[numberOfV0] = op2.z();
                mV0Dst.dau2NHitsTpc[numberOfV0] = antiDau2TracksVector[j] -> nHitsFit(kTpcId);
                mV0Dst.dau2NHitsSsd[numberOfV0] = antiDau2TracksVector[j] -> nHitsFit(kSsdId);
                mV0Dst.dau2NHitsSvt[numberOfV0] = antiDau2TracksVector[j] -> nHitsFit(kSvtId);
                mV0Dst.dau2TofFlag[numberOfV0] = antiDau2TracksVector[j] -> btofPidTraits().matchFlag();
                mV0Dst.dau2Tof[numberOfV0] = antiDau2TracksVector[j] -> btofPidTraits().timeOfFlight();

		numberOfV0++;
	    }
	}
    }


    //  Int_t sizeOfPairTracksVector = pairTracksVector.size();
    //  Int_t sizeOfAntiPairTracksVector = antiPairTracksVector.size();
#ifdef DEBUG
    //printf("the number of tracks passing dca cut is %d\n", sizeOfPairTracksVector); 
#endif

    /*  for(Int_t i = 0; i < sizeOfPairTracksVector; i++){
    //StThreeVectorF momentumProton=pairTracksVector[i].first->helix().momentumAt(pathlengthPairVector[i].first, magneticField);
    //StThreeVectorF momentumPion=pairTracksVector[i].second->helix().momentumAt(pathlengthPairVector[i].second, magneticField);


    //Double_t protonE=sqrt(mProtonMass*mProtonMass+momentumProton.dot(momentumProton));
    Double_t protonE = sqrt(mProtonMass * mProtonMass + (momentumPairVector[i].first).mag2());
    Double_t pionE = sqrt(mPionMass * mPionMass + (momentumPairVector[i].second).mag2());

    StThreeVectorF momentumSumTwoTracks = momentumPairVector[i].first + momentumPairVector[i].second;// what is the difference between helix's and track's momenta
    //StThreeVectorF momentumSumTwoTracks=momentumProton+momentumPion;// what is the difference between helix's and track's momenta
    Double_t invV0Mass = sqrt((protonE + pionE) * (protonE + pionE) - momentumSumTwoTracks.mag2());
    hInvariantV0Mass -> Fill(invV0Mass);

#ifdef DEBUG
printf("Mass is %f\n", invV0Mass); 
#endif
}*/

    /*  for(Int_t i = 0; i < sizeOfAntiPairTracksVector; i++){
    //StThreeVectorF momentumProton=pairTracksVector[i].first->helix().momentumAt(pathlengthPairVector[i].first, magneticField);
    //StThreeVectorF momentumPion=pairTracksVector[i].second->helix().momentumAt(pathlengthPairVector[i].second, magneticField);


    //Double_t protonE=sqrt(mAntiProtonMass*mAntiProtonMass+momentumProton.dot(momentumProton));
    Double_t antiProtonE=sqrt(mAntiProtonMass*mAntiProtonMass+(antiMomentumPairVector[i].first).mag2());
    Double_t pionPlusE=sqrt(mPionPlusMass*mPionPlusMass+(antiMomentumPairVector[i].second).mag2());

    StThreeVectorF momentumSumTwoTracks=momentumPairVector[i].first+momentumPairVector[i].second;// what is the difference between helix's and track's momenta
    //StThreeVectorF momentumSumTwoTracks=momentumProton+momentumPion;// what is the difference between helix's and track's momenta
    Double_t invV0Mass=sqrt((protonE+pionE)*(protonE+pionE)-momentumSumTwoTracks.mag2());
    hInvariantV0Mass->Fill(invV0Mass);

#ifdef DEBUG
printf("Mass is %f\n", invV0Mass); 
#endif
}*/
    mV0Dst.numberOfV0 = numberOfV0;

    mV0Tree -> Fill();
    return kStOK; 
}

Int_t StV0Maker::Finish(){

    mOutHistoFile -> Write();
#ifdef DEBUG_HISTO
//    gDirectory -> pwd();
//    if(mOutHistoFile -> Write()) printf("Histo write OK!\n");
#endif
#ifdef DEBUG_HISTO
//    if(mV0PicoDstFile -> Write()) printf("picoDst write OK!\n");
#endif
    mV0PicoDstFile -> Write();
    printf("File Closed!\n");
    return kStOK;

}

StThreeVectorF StV0Maker::GetTrackDca2PV(StPhysicalHelixD& helix, StThreeVectorF& primaryVertexPosition){

    Double_t pathLength = helix.pathLength( primaryVertexPosition, true); 
    StThreeVectorF dca = helix.at( pathLength) - primaryVertexPosition;
    return dca;
}

Bool_t StV0Maker::SetOutHistoFileName(TString outHistoFileName){ mOutHistoFileName = outHistoFileName; return 1;}

Bool_t StV0Maker::SetOutDSTFileName(TString outDSTFileName){ mV0PicoDstFileName = outDSTFileName; return 1;}

Bool_t StV0Maker::pCutProtonNSigma(StMuTrack* singleTrack){

    Double_t nSigmaProton = singleTrack -> nSigmaProton();   
    if (fabs(nSigmaProton) < mCutNSigmaProtonLe) return true;//kStOK has been defined in StTypeDefs.h
    else return false;
}

Bool_t StV0Maker::pCutAntiProtonNSigma(StMuTrack* singleTrack){
    Double_t nSigmaProton = singleTrack -> nSigmaProton();
    if(fabs(nSigmaProton) < mCutNSigmaProtonLe) return true;
    else return false;
}

Bool_t StV0Maker::pCutPionMinusNSigma(StMuTrack* singleTrack){
    Double_t nSigmaPion = singleTrack -> nSigmaPion();   
    if (fabs(nSigmaPion) < mCutNSigmaPionLe) return true;//kStOK has been defined in StTypeDefs.h
    else return false;
}

Bool_t StV0Maker::pCutPionPlusNSigma(StMuTrack* singleTrack){
    Double_t nSigmaPion = singleTrack -> nSigmaPion();
    if(fabs(nSigmaPion) < mCutNSigmaPionLe) return true;
    else return false;
}

Bool_t StV0Maker::pCutGoodEvents(StMuEvent* event){
    if(!event) return false;
    else return true;
}

Bool_t StV0Maker::pCutEventsTrigger(StMuEvent* event){
    StTriggerId triggerId =event -> triggerIdCollection().nominal();
    if(!triggerId.isTrigger(350003) && !triggerId.isTrigger(350013) && !triggerId.isTrigger(350023) && !triggerId.isTrigger(350033) && !triggerId.isTrigger(350043)) return false;
    else return true;
}

Bool_t StV0Maker::pCutEventsPVZ(StMuEvent* event){
    if(fabs(event -> primaryVertexPosition().z()) > mCutAbsPrimaryVertexZLeEq ) return false;
    else return true;
}

Bool_t StV0Maker::pCutNHits(Int_t nHits){
    if(nHits < mCutNHitsGrEq) return false;
    else return true;
}

Bool_t StV0Maker::pCutPt(Double_t pt){
    if(pt < mCutPtGrEq) return false;
    else return true;
}
 
Bool_t StV0Maker::pCutDau1Dca2PV(Double_t dau1dca2PV){
    if(dau1dca2PV < mCutDau1Dca2PVGrEq) return false;
    else return true;
}

Bool_t StV0Maker::pCutDau2Dca2PV(Double_t dau2dca2PV){
    if(dau2dca2PV < mCutDau2Dca2PVGrEq) return false;
    else return true;
}

Bool_t StV0Maker::pCutTwoTracksDca(Double_t twoTracksDca){
    if(twoTracksDca > mCutTwoTracksDcaLeEq) return false;
    else return true;
}

Bool_t StV0Maker::pCutDecLength(Double_t decayLength){
    if(decayLength < mCutDecLengthGrEq) return false;
    else return true;
}

Bool_t StV0Maker::pCutDcaV02PV(Double_t v0Dca2PV){
    if(v0Dca2PV > mCutDcaV02PVLe) return false;
    else return true;
}


