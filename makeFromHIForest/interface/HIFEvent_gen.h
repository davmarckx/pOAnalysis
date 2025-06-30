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
  std::vector<char> *gen_trk_charge;
  std::vector<float> *gen_trk_pt, *gen_trk_phi, *gen_trk_eta;
  std::vector<int> *gen_trk_id;
  
  
  // Gen info
  //Int_t typevt;
  
};

void attachToHIFEventTree_gen(TTree *t, HIFEventgen_t &ev);

#endif

