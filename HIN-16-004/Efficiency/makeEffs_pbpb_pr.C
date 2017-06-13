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

void makeEffs_pbpb_pr(bool ispbpb=true, bool isprompt=true, bool isacc=false) {
   TChain *tch_jpsi_pbpb = new TChain("hionia/myTree");
   if (isacc) {
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
   } else {
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   }

   // make the efficiency histos
   string dir = "files/eff"; // output files will be stored under this directory
   gSystem->mkdir(Form("./%s",dir.c_str()), kTRUE);

   cout << "Efficiencies for pbpb prompt Jpsi" << endl;
   oniaEff obj_jpsi_pbpb(tch_jpsi_pbpb);
   obj_jpsi_pbpb.Loop(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,obj_jpsi_pbpb.trg,isacc);
   
//   oniaEff_pTShapeVary obj_jpsi_pbpb(tch_jpsi_pbpb);
//   obj_jpsi_pbpb.LoopVary(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,obj_jpsi_pbpb.trg_ptWeighting);

//   oniaEff_TnPToyStudy obj_jpsi_pbpb(tch_jpsi_pbpb);
//   obj_jpsi_pbpb.LoopVary(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,obj_jpsi_pbpb.tnpTypes::trg_toy);

}
