#ifndef StV0Dst_Def
#define StV0Dst_Def
#include "StV0Type.h"
#include "StV0DauType.h"
#define MAX_NUM_V0 2000
struct StV0Dst{

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

    int numberOfV0;
// 1-Lambda, 0-AntiLambda
    int v0Type[MAX_NUM_V0];

    float v0Mass[MAX_NUM_V0];
    float v0Pt[MAX_NUM_V0];
    float v0Rapidity[MAX_NUM_V0];
    float v0Eta[MAX_NUM_V0];
    float v0X[MAX_NUM_V0];
    float v0Y[MAX_NUM_V0];
    float v0Z[MAX_NUM_V0];
    float v0Px[MAX_NUM_V0];
    float v0Py[MAX_NUM_V0];
    float v0Pz[MAX_NUM_V0];
 
    float v0DecLength[MAX_NUM_V0];
    float v0Dca2PV[MAX_NUM_V0];
    float v0Dca2PV2d[MAX_NUM_V0];
    float v0Dca2Dau[MAX_NUM_V0];
//?? float v0PathLength[MAX_NUM_V0];

// I dount the tracks Id's use
    int dau1TrackId[MAX_NUM_V0];
    int dau1Charge[MAX_NUM_V0];
    float dau1Dca2PV[MAX_NUM_V0];
    float dau1Dca2PV2d[MAX_NUM_V0];
    int dau1NHits[MAX_NUM_V0];
//    float dau1DeDx[MAX_NUM_V0];
    float dau1NSigmaPion[MAX_NUM_V0];
    float dau1NSigmaProton[MAX_NUM_V0];
//    float dau1Eta[MAX_NUM_V0];
    float dau1Pt[MAX_NUM_V0];
    float dau1Px[MAX_NUM_V0];
    float dau1Py[MAX_NUM_V0];
    float dau1Pz[MAX_NUM_V0];
    int dau1NHitsTpc[MAX_NUM_V0];
    int dau1NHitsSsd[MAX_NUM_V0];
    int dau1NHitsSvt[MAX_NUM_V0];
    int dau1TofFlag[MAX_NUM_V0];
    float dau1Tof[MAX_NUM_V0];
//    float dau1PathLength[MAX_NUM_V0];

    int dau2TrackId[MAX_NUM_V0];
    int dau2Charge[MAX_NUM_V0];
    float dau2Dca2PV[MAX_NUM_V0];
    float dau2Dca2PV2d[MAX_NUM_V0];
    int dau2NHits[MAX_NUM_V0];
//    float dau2DeDx[MAX_NUM_V0];
    float dau2NSigmaPion[MAX_NUM_V0];
    float dau2NSigmaProton[MAX_NUM_V0];
//    float dau2Eta[MAX_NUM_V0];
    float dau2Pt[MAX_NUM_V0];
    float dau2Px[MAX_NUM_V0];
    float dau2Py[MAX_NUM_V0];
    float dau2Pz[MAX_NUM_V0];
    int dau2NHitsTpc[MAX_NUM_V0];
    int dau2NHitsSsd[MAX_NUM_V0];
    int dau2NHitsSvt[MAX_NUM_V0];
    int dau2TofFlag[MAX_NUM_V0];
    float dau2Tof[MAX_NUM_V0];
//    float dau2PathLength[MAX_NUM_V0];
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
  StV0Type v0Type;
  double v0Mass;
  double v0PX;
  double v0PY;
  double v0PZ;
  double dcaV0PV;
  

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
