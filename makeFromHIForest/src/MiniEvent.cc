#include "../interface/MiniEvent.h"
#include <iostream>

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{

  //event header
  t->Branch("isData",    &ev.isData,   "isData/O");
  t->Branch("run",       &ev.run,      "run/i");
  t->Branch("lumi",      &ev.lumi,     "lumi/i");
  t->Branch("event",     &ev.event,    "event/l");

  t->Branch("weight",     &ev.weight,    "weight/F");
  t->Branch("typevt",        &ev.typevt,          "typevt/I");

  t->Branch("ntrk",        &ev.ntrk,          "ntrk/I");
  t->Branch("trk_p",        ev.trk_p,         "trk_p[ntrk]/F");
  t->Branch("trk_pt",       ev.trk_pt,        "trk_pt[ntrk]/F");
  t->Branch("trk_eta",      ev.trk_eta,       "trk_eta[ntrk]/F");
  t->Branch("trk_phi",      ev.trk_phi,       "trk_phi[ntrk]/F");
  t->Branch("trk_dxy",      ev.trk_dxy,       "trk_dxy[ntrk]/F");
  t->Branch("trk_dz",       ev.trk_dz,        "trk_dz[ntrk]/F");
  t->Branch("trk_q",        ev.trk_q,         "trk_q[ntrk]/I");
  t->Branch("trk_dedx",     ev.trk_dedx,      "trk_dedx[ntrk]/F");
  t->Branch("trk_dedxerr",  ev.trk_dedxerr,   "trk_dedxerr[ntrk]/F");
  t->Branch("trk_hasPF",       ev.trk_hasPF,        "trk_hasPF[ntrk]/I");
  t->Branch("trk_isPi",       ev.trk_isPi,        "trk_isPi[ntrk]/I");
  t->Branch("trk_isK",       ev.trk_isK,        "trk_isK[ntrk]/I");
  t->Branch("trk_isP",       ev.trk_isP,        "trk_isP[ntrk]/I");
  t->Branch("trk_isD",       ev.trk_isD,        "trk_isD[ntrk]/I");

  t->Branch("trk_isPi_loose",       ev.trk_isPi_loose,        "trk_isPi_loose[ntrk]/I");
  t->Branch("trk_isK_loose",       ev.trk_isK_loose,        "trk_isK_loose[ntrk]/I");
  t->Branch("trk_isP_loose",       ev.trk_isP_loose,        "trk_isP_loose[ntrk]/I");
  t->Branch("trk_isD_loose",       ev.trk_isD_loose,        "trk_isD_loose[ntrk]/I");
  t->Branch("trk_isPi_binned",       ev.trk_isPi_binned,        "trk_isPi_binned[ntrk]/I");
  t->Branch("trk_isK_binned",       ev.trk_isK_binned,        "trk_isK_binned[ntrk]/I");
  t->Branch("trk_isP_binned",       ev.trk_isP_binned,        "trk_isP_binned[ntrk]/I");

  t->Branch("trk_genPt",       ev.trk_genPt,        "trk_genPt[ntrk]/F");
  t->Branch("trk_dRmatch",       ev.trk_dRmatch,    "trk_dRmatch[ntrk]/F");
  t->Branch("trk_matchedPdgId",  ev.trk_matchedPdgId,  "trk_matchedPdgId[ntrk]/I");
  t->Branch("trk_numberOfPixelHits",  ev.trk_numberOfPixelHits,  "trk_numberOfPixelHits[ntrk]/I");
  t->Branch("trk_numberOfHits",  ev.trk_numberOfHits,  "trk_numberOfHits[ntrk]/I");

  t->Branch("gen_ntrk",        &ev.gen_ntrk,          "gen_ntrk/I");
  t->Branch("gen_nD",        &ev.gen_nD,          "gen_nD/I");
  t->Branch("gen_trk_pt",        ev.gen_trk_pt,         "gen_trk_pt[gen_ntrk]/F");
  t->Branch("gen_trk_p",        ev.gen_trk_p,         "gen_trk_p[gen_ntrk]/F");
  t->Branch("gen_trk_E",        ev.gen_trk_E,         "gen_trk_E[gen_ntrk]/F");
  t->Branch("gen_trk_id",       ev.gen_trk_id,        "gen_trk_id[gen_ntrk]/I");
  t->Branch("gen_trk_eta",      ev.gen_trk_eta,       "gen_trk_eta[gen_ntrk]/F");
  t->Branch("gen_trk_phi",      ev.gen_trk_phi,       "gen_trk_phi[gen_ntrk]/F");
  t->Branch("gen_trk_hasReco",    ev.gen_trk_hasReco,     "gen_trk_hasReco[gen_ntrk]/I");

}

void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{

  //event header
  t->SetBranchAddress("isData",    &ev.isData);
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  t->SetBranchAddress("weight",    &ev.weight);
  t->SetBranchAddress("typevt",    &ev.typevt);
  
  t->SetBranchAddress("ntrk",        &ev.ntrk);
  t->SetBranchAddress("trk_p",       ev.trk_p);
  t->SetBranchAddress("trk_pt",      ev.trk_pt);
  t->SetBranchAddress("trk_eta",     ev.trk_eta);
  t->SetBranchAddress("trk_phi",     ev.trk_phi);
  t->SetBranchAddress("trk_dxy",     ev.trk_dxy);
  t->SetBranchAddress("trk_dz",      ev.trk_dz);
  t->SetBranchAddress("trk_q",       ev.trk_q);
  t->SetBranchAddress("trk_dedx",    ev.trk_dedx);
  t->SetBranchAddress("trk_dedxerr", ev.trk_dedxerr);
  t->SetBranchAddress("trk_hasPF",       ev.trk_hasPF);
  t->SetBranchAddress("trk_isPi",       ev.trk_isPi);
  t->SetBranchAddress("trk_isK",       ev.trk_isK);
  t->SetBranchAddress("trk_isP",       ev.trk_isP);
  t->SetBranchAddress("trk_isD",       ev.trk_isD);

  t->SetBranchAddress("trk_isPi_loose",       ev.trk_isPi_loose);
  t->SetBranchAddress("trk_isK_loose",       ev.trk_isK_loose);
  t->SetBranchAddress("trk_isP_loose",       ev.trk_isP_loose);
  t->SetBranchAddress("trk_isPi_binned",       ev.trk_isPi_binned);
  t->SetBranchAddress("trk_isK_binned",       ev.trk_isK_binned);
  t->SetBranchAddress("trk_isP_binned",       ev.trk_isP_binned);


  t->SetBranchAddress("trk_genPt",      ev.trk_genPt);
  t->SetBranchAddress("trk_dRmatch",      ev.trk_dRmatch);
  t->SetBranchAddress("trk_matchedPdgId",      ev.trk_matchedPdgId);
  t->SetBranchAddress("trk_numberOfPixelHits",      ev.trk_numberOfPixelHits);
  t->SetBranchAddress("trk_numberOfHits",      ev.trk_numberOfHits);

  t->SetBranchAddress("gen_ntrk",        &ev.gen_ntrk);
  t->SetBranchAddress("gen_nD",        &ev.gen_nD);
  t->SetBranchAddress("gen_trk_pt", ev.gen_trk_pt);
  t->SetBranchAddress("gen_trk_p", ev.gen_trk_p);
  t->SetBranchAddress("gen_trk_E", ev.gen_trk_E);
  t->SetBranchAddress("gen_trk_eta",     ev.gen_trk_eta);
  t->SetBranchAddress("gen_trk_phi",     ev.gen_trk_phi);
  t->SetBranchAddress("gen_trk_hasReco",        &ev.gen_trk_hasReco);
  t->SetBranchAddress("gen_trk_id", ev.gen_trk_id);
}

