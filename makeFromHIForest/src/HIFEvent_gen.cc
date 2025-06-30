#include "../interface/HIFEvent_gen.h"
#include <iostream>

void attachToHIFEventTree_gen(TTree *t,HIFEventgen_t &ev)
{

  //event header
  //t->SetBranchAddress("isData",    &ev.isData);
  //t->SetBranchAddress("nRun",       &ev.run);
  //t->SetBranchAddress("nEv",     &ev.event);
  //t->SetBranchAddress("nLumi",      &ev.lumi);
  //t->SetBranchAddress("weight",    &ev.weight);
  //t->SetBranchAddress("typevt",    &ev.typevt);
  
  t->SetBranchAddress("npart",        &ev.gen_ntrk);
  t->SetBranchAddress("pt", ev.gen_trk_pt);
  t->SetBranchAddress("chg",     ev.gen_trk_charge);
  t->SetBranchAddress("phi",     ev.gen_trk_phi);
  t->SetBranchAddress("eta",     ev.gen_trk_eta);
  t->SetBranchAddress("pdg", ev.gen_trk_id);
}

