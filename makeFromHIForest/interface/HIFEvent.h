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
  std::vector<int> *trk_nMeasure;
  std::vector<int> *trk_nSaturMeasure;
  std::vector<int> *trk_nMeasureLayer;
  std::vector<int> *trk_numberOfPixelHits; 
  std::vector<int> *trk_numberOfHits;

  std::vector<char> *trk_q;


  std::vector<float> *trk_dxy;
  std::vector<float> *trk_dz;
  std::vector<float> *trk_pt_error;
  std::vector<float> *trk_pt;
  std::vector<float> *trk_eta;
  std::vector<float> *trk_phi;
  std::vector<float> *trk_dedx;

  // Gen info
  Int_t typevt;
  
};

void attachToHIFEventTree(TTree *t,HIFEvent_t &ev);

#endif

