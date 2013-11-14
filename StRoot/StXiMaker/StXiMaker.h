#ifndef StXiMaker_def// Generally, class cannot be defined (not the implementation) twice, so this prevent this occasion.
#define StXiMaker_def

#include "StMaker.h"// Why we dont need to include StMuDstMaker.h??-----Because below "class StMuDstMaker" is called a forward declaration, and we only used the pointer or the reference of this type so we don't need to define it.
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"// StV0Maker's parameter 
#include "StMuDSTMaker/COMMON/StMuEvent.h"// StV0Maker's parameter 
#include "StMuDSTMaker/COMMON/StMuTrack.h"// StV0Maker's parameter 
#include "TString.h"
#include "StV0Maker.h"
//#include "StV0Type.h"
//#include "StV0Dst.h"
#include "StXiDst.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include <utility>
#include <vector> 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
class StXiMaker: public StMaker{

    public:

	StXiMaker(StMuDstMaker* muDstMaker, StV0Maker* mV0Maker, TString outHistoFileName, TString outPicoFileName ); 
	virtual ~StXiMaker(){;}

	void SetOutFileName(TString outFileName);
	void SetRefMultiplicity(Int_t, Int_t);     // Actually I have no idea what does the reference multiplicity mean. I think maybe the multiplicity interval's upper and lower bounds.

	Int_t Init();
	Int_t Make();
	Int_t Finish();

    private:
	//-------------------General Variables-----------------------------------
        TFile* mOutHistoFile;
        TFile* mXiPicoDstFile;
	TString mOutHistoFileName;
	TString mXiPicoDstFileName;
	TString mBuffer;    // Temp buffer to store the string
	//TString mV0ParticleName;   // Define the V0 particle name
	Double_t mRotation;    // Background Estimation
	Double_t mSameSign;
	Double_t mProtonMass;   // Define the daughters' masses and V0 mass 
	Double_t mAntiProtonMass;
	Double_t mPionMinusMass; 
	Double_t mPionPlusMass;
	Double_t mXiZeroMass;
	Double_t mXiMinusMass;
	Double_t mXiMassWidth;
	TString mHInvariantMassName;    // Invariant mass spectrum name defined by user code
	//StXiType mXiType;
        //TClonesArray* mGlobalTracks;// To store the global tracks in each event being processed.
        //TClonesArray* mPrimaryTracks;// To store the primary tracks in each event being processed.
	//std::vector<std::pair<StMuTrack*, StMuTrack*> > pairTracksVector;
	//std::vector<std::pair<StMuTrack*, StMuTrack*> > antiPairTracksVector;
	//std::vector<std::pair<StThreeVectorF, StThreeVectorF> > momentumPairVector;
	std::vector<std::pair<StThreeVectorF, StThreeVectorF> > antiMomentumPairVector;
	StMuDstMaker* mMuDstMaker;
        StV0Maker* mV0Maker;
	Int_t mEventFailed;    // Record the events read unsuccessfully and successfully 
	Int_t mEventPassed; 	
	StXiDst mXiDst;	// PicoDst Data
	TTree* mXiTree;  // Data Tree


        //--------------------Cuts------------------------------------------------
	Int_t mCutTriggerIdEq;    // TriggerCuts
	Double_t mCutAbsPrimaryVertexZLeEq;     // Event Primary Vertex Position Cuts
	      //-------------------------TrackCuts
	unsigned short mCutNHitsGrEq;
	Double_t mCutPtGrEq;
	Double_t mCutNSigmaPionLe;
	Double_t mCutNSigmaProtonLe;
	Double_t mCutNSigmaElectronLe;
	Double_t mCutNSigmaKaonLe;
	      //-------------------------- DcaCuts
	Double_t mCutDau1Dca2PVGrEq;   // Single helix dca to primary vertex
	Double_t mCutDau2Dca2PVGrEq;  
	Double_t mCutTwoTracksDcaLeEq;    // Two tracks dca cuts
	Double_t mCutDecLengthGrEq;    // Cuts on decay length    
	Double_t mCutDcaXi2PVLe;      // Cuts on v0 to PV distance

        //---------------------Histograms----------------------------------------
	//TH1F* hNPV;
	TH1F* hVz; // Vertex position
	TH1F* hRefMult; // Referenced multiplicity----between a certain eta range the charged particle number 
	    //TH1F* hSelectNRefMult;
	    //TH1F* hSelectNRefMultM1;
	    //TH1F* hSelectNRefMultM2;
	    //TH1F* hSelectNRefMultFtpcE;
	TH1F* hPtRaw;
	TH1F* hEtaRaw;
	TH1F* hPhiRaw;    
        TH1F* hPt;
	TH1F* hEta;  // Compared to raw histograms, these are cutted after nhits and charge and tofflag
	TH1F* hPhi;
            //TH1F* hPhiLowPt;
	    //TH1F* hPhiHighPt;
	TH2F* hDeDxP;// Raw 
	TH1F* hNHitsFit;
	TH1F* hNHits;// Raw
	TH1F* hNSigmaPion;// Raw
	TH1F* hNSigmaProton;// Raw
	    //TH1F* hNSigmaKaon;
	    //TH1F* hPtDiff;
	    //TH1F* hOrDiff;
            //TH1F* hPDiff;
	    //TH1F* hDcaDiff;
	TH1F* hXiPt;
	    //TH1F* hDau1Dca;
	    //TH1F* hDau2Dca;
            //TH1F* hDau1Dca2;
	    //TH1F* hDau2Dca2;
	//TH1F* hXiZeroInvariantMass;
	TH1F* hXiMinusInvariantMass;

	//---------------------Member functions------------------------ 
	StThreeVectorF GetTrackDca2PV(StPhysicalHelixD&, StThreeVectorF&);    
	Bool_t SetOutHistoFileName(TString);    
	Bool_t SetOutDSTFileName(TString);
        Bool_t pCutLambdaMass();// to be implemented
	Bool_t pCutProtonNSigma(StMuTrack*);
	Bool_t pCutAntiProtonNSigma(StMuTrack*);
	Bool_t pCutPionMinusNSigma(StMuTrack*);
	Bool_t pCutPionPlusNSigma(StMuTrack*);
	    // Bool_t SetRotation(){ mRotation=1;}
	    // Bool_t SetSameSign(){ mSameSign=1;}
	Bool_t pCutGoodEvents(StMuEvent*);
	Bool_t pCutEventsTrigger(StMuEvent*);
	Bool_t pCutEventsPVZ(StMuEvent*);
        Bool_t pCutNHits(Int_t);
        Bool_t pCutPt(Double_t);
        Bool_t pCutDau1Dca2PV(Double_t);
        Bool_t pCutDau2Dca2PV(Double_t);
        Bool_t pCutDcaXi2PV(Double_t);
        Bool_t pCutTwoTracksDca(Double_t);
        Bool_t pCutDecLength(Double_t);
        Bool_t pCutV0Mass(Double_t);
	//------------------Macro for CINT compatability-------------------
	ClassDef(StXiMaker,1) // protected: no use because it has no derived classes.
    
};
#endif
