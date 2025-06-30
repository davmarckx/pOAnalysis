#include "../interface/HIFEvent.h"
#include <iostream>

void attachToHIFEventTree(TTree *t,HIFEvent_t &ev)
{

  //event header
  //t->SetBranchAddress("isData",    &ev.isData);
  t->SetBranchAddress("nRun",       &ev.run);
  t->SetBranchAddress("nEv",     &ev.event);
  t->SetBranchAddress("nLumi",      &ev.lumi);

  //t->SetBranchAddress("weight",    &ev.weight);
  //t->SetBranchAddress("typevt",    &ev.typevt);
  
  t->SetBranchAddress("nTrk",        &ev.ntrk);
  t->SetBranchAddress("trkPt",      &ev.trk_pt);
  t->SetBranchAddress("trkPtError",      &ev.trk_pt_error);
  t->SetBranchAddress("trkEta",     &ev.trk_eta);
  t->SetBranchAddress("trkPhi",     &ev.trk_phi);
  t->SetBranchAddress("trkDxyAssociatedVtx",     &ev.trk_dxy);
  t->SetBranchAddress("trkDzAssociatedVtx",      &ev.trk_dz);
  t->SetBranchAddress("trkCharge",       &ev.trk_q);
  t->SetBranchAddress("dedxAllLikelihood",    &ev.trk_dedx);
  t->SetBranchAddress("trkNPixHits",      &ev.trk_numberOfPixelHits);
  t->SetBranchAddress("trkNHits",      &ev.trk_numberOfHits);
}

