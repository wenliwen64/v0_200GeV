#ifndef StXiDst_Def
#define StXiDst_Def
//#include "StXiType.h"
//#include "StXiDauType.h"
#define MAX_NUM_Xi 2000
struct StXiDst{

    int runNumber;
    int eventId;
    int triggerMode;
    int nRefMult;
    int nRefMultTof; 
// problematic 
    //float zdcadc0;
    //int bbadcsumeast;
    //double bbccirate;

    float magneticField;

    float PVX;
    float PVY;
    float PVZ;

    int numberOfXi;
// 1-Lambda, 0-AntiLambda
    int xiType[MAX_NUM_Xi];

    float xiMass[MAX_NUM_Xi];
    float xiPt[MAX_NUM_Xi];
    float xiRapidity[MAX_NUM_Xi];
    float xiEta[MAX_NUM_Xi];
    float xiX[MAX_NUM_Xi];
    float xiY[MAX_NUM_Xi];
    float xiZ[MAX_NUM_Xi];
    float xiPx[MAX_NUM_Xi];
    float xiPy[MAX_NUM_Xi];
    float xiPz[MAX_NUM_Xi];
 
    float xiDecLength[MAX_NUM_Xi];
    float xiDca2PV[MAX_NUM_Xi];
    float xiDca2PV2d[MAX_NUM_Xi];
    float xiDca2Dau[MAX_NUM_Xi];
//?? float xiPathLength[MAX_NUM_Xi];

// I dount the tracks Id's use
    //int dau1TrackId[MAX_NUM_Xi];
    //int dau1Charge[MAX_NUM_Xi];
    //float dau1Dca2PV[MAX_NUM_Xi];
    //float dau1Dca2PV2d[MAX_NUM_Xi];
    //int dau1NHits[MAX_NUM_Xi];
//    float dau1DeDx[MAX_NUM_Xi];
    //float dau1NSigmaPion[MAX_NUM_Xi];
    //float dau1NSigmaProton[MAX_NUM_Xi];
//    float dau1Eta[MAX_NUM_Xi];
    float v0Pt[MAX_NUM_Xi];
    float v0Px[MAX_NUM_Xi];
    float v0Py[MAX_NUM_Xi];
    float v0Pz[MAX_NUM_Xi];
    int v0NHitsTpc[MAX_NUM_Xi];
    int v0NHitsSsd[MAX_NUM_Xi];
    int v0NHitsSvt[MAX_NUM_Xi];
    int v0TofFlag[MAX_NUM_Xi];
    float v0Tof[MAX_NUM_Xi];
//    float dau1PathLength[MAX_NUM_Xi];

    int pionMinusTrackId[MAX_NUM_Xi];
    int pionMinusCharge[MAX_NUM_Xi];
    float pionMinusDca2PV[MAX_NUM_Xi];
    float pionMinusDca2PV2d[MAX_NUM_Xi];
    int pionMinusNHits[MAX_NUM_Xi];
//    float dau2DeDx[MAX_NUM_Xi];
    float pionMinusNSigmaPion[MAX_NUM_Xi];
    //float dau2NSigmaProton[MAX_NUM_Xi];
//    float dau2Eta[MAX_NUM_Xi];
    float pionMinusPt[MAX_NUM_Xi];
    float pionMinusPx[MAX_NUM_Xi];
    float pionMinusPy[MAX_NUM_Xi];
    float pionMinusPz[MAX_NUM_Xi];
    int pionMinusNHitsTpc[MAX_NUM_Xi];
    int pionMinusNHitsSsd[MAX_NUM_Xi];
    int pionMinusNHitsSvt[MAX_NUM_Xi];
    int pionMinusTofFlag[MAX_NUM_Xi];
    float pionMinusTof[MAX_NUM_Xi];
//    float dau2PathLength[MAX_NUM_Xi];
//  short type; // all global tracks???
//  short flag;// negative=bad, positive=good, we only need to store the good ones? but there is also a function called bad(), what does it mean?
//  unsigned int flagExtension;// ?? 
//  int bad;// ?
//  int index2Global; // ?? 
//  int index2RichSpectra;// ??
//  int index2BTofHit;// ??
//  int vertexIndex;// ??
//  unsigned short nHits;// Total hits number
//  unsigned short nHitsPoss;// possibel hits number  
//  double nSigmaElectron;
//  double nSigmaPion;
//  double nSigmaProton;
//  double nSigmaKaon;
/*  double dEdx;
  double pt;
  double eta;
  double charge;
  //StThreeVectorF momentum;// this is for local analysis, so there would not be any StThreeVectorF.h files on giant or quarf.
  //StThreeVectorF p;
  StXiType xiType;
  double xiMass;
  double xiPX;
  double xiPY;
  double xiPZ;
  double dcaXiPV;
  

  dauType dau1Type;
  short dau1id;
  unsigned short dau1NHits;
  double dau1Charge;
  double dau1DEDx;
  short dau1id;
  double dau1PX;
  double dau1PY;
  double dau1PZ;
  double dcaDau1PV;


  dauType dau2type;
  short dau2id;
  unsigned short dau2NHits;
  double dau2Charge;
  double dau2DEDx;
  short dau2id;
  double dau2PX;
  double dau2PY;
  double dau2PZ;
  double dcaDau2PV

  double dcaDau1Dau2;
*/
};
#endif
