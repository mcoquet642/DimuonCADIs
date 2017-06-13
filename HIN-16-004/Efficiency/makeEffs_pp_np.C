#include "oniaEff.C"
#include "oniaEff_TnPToyStudy.C"
#include "oniaEff_pTShapeVary.C"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TDirectory.h>

#include <iostream>

using namespace std;

void makeEffs_pp_np(bool ispbpb=false, bool isprompt=false, bool isacc=false) {
   TChain *tch_npjpsi_pp = new TChain("hionia/myTree");
   if (isacc)
     tch_npjpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
   else
     tch_npjpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_Extended.root");

   // make the efficiency histos
   string dir = "files/eff"; // output files will be stored under this directory
   gSystem->mkdir(Form("./%s",dir.c_str()), kTRUE);

   cout << "Efficiencies for pp non-prompt Jpsi" << endl;
   oniaEff obj_npjpsi_pp(tch_npjpsi_pp);
   obj_npjpsi_pp.Loop(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,obj_npjpsi_pp.trg,isacc);

//   oniaEff_pTShapeVary obj_npjpsi_pp(tch_npjpsi_pp);
//   obj_npjpsi_pp.LoopVary(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,obj_npjpsi_pp.trg_ptWeighting);

//   oniaEff_TnPToyStudy obj_npjpsi_pp(tch_npjpsi_pp);
//   obj_npjpsi_pp.LoopVary(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,obj_npjpsi_pp.tnpTypes::trg_toy);

}
