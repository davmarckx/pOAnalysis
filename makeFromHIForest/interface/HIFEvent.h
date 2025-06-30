#ifndef _hifevent_h_
#define _hifevent_h_

#include "TTree.h"

struct HIFEvent_t
{
  HIFEvent_t()
  {
    ntrk=0;
  }

  static const int MAXTRACKS     =  400;
  static const int MAXGENTRACKS  =  4000;  
  static const int MAXPROTONS    =  4;

  Bool_t isData;
  Int_t run,lumi;
  Int_t event;

  Float_t weight;
  
  // Vertex info
  Int_t nvtx;

  Int_t ntrk;

  //track info
  std::vector<int> *trk_nMeasure[MAXTRACKS], *trk_nSaturMeasure[MAXTRACKS], *trk_nMeasureLayer[MAXTRACKS], *trk_numberOfPixelHits[MAXTRACKS], *trk_numberOfHits[MAXTRACKS];
  std::vector<int> *trk_q[MAXTRACKS];
  std::vector<float> *trk_dxy[MAXTRACKS], *trk_dz[MAXTRACKS];
  std::vector<float> *trk_pt_error[MAXTRACKS], *trk_pt[MAXTRACKS], *trk_eta[MAXTRACKS], *trk_phi[MAXTRACKS], *trk_dedx[MAXTRACKS];

  // Gen info
  Int_t typevt;
  
};

void attachToHIFEventTree(TTree *t,HIFEvent_t &ev);

#endif

