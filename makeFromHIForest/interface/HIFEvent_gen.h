#ifndef _hifgenevent_h_
#define _hifgenevent_h_

#include "TTree.h"

struct HIFEventgen_t
{
  HIFEventgen_t()
  {
    gen_ntrk=0;
  }

  static const int MAXTRACKS     =  400;
  static const int MAXGENTRACKS  =  4000;  
  static const int MAXPROTONS    =  4;

  //gen track info
  Float_t gen_ntrk;
  Int_t gen_trk_charge[MAXGENTRACKS];
  std::vector<float> *gen_trk_pt[MAXGENTRACKS], *gen_trk_phi[MAXGENTRACKS], *gen_trk_eta[MAXGENTRACKS];
  Int_t gen_trk_id[MAXGENTRACKS];
  
  
  // Gen info
  //Int_t typevt;
  
};

void attachToHIFEventTree_gen(TTree *t, HIFEventgen_t &ev);

#endif

