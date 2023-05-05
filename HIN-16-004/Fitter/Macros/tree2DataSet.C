// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program creates two RooDataSets (opposite-sign and same-sign dimuons) from an Onia Tree.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 26 19:08 CET 2016
//
//

#include "Utilities/initOniaTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"

#include "TObjArray.h"

map<int, double>   fCentMap; // map for centrality-Ncoll mapping
double             fCentBinning[200];
int                fCentBins;
TObjArray*         fcorrArray = NULL; // Array with the 2D correction for weighting

string  findMyTree(string FileName);
bool    getTChain(TChain* fChain, vector<string> FileNames, string TreeName="O2rtdimuonall");
void    iniBranch(TChain* fChain,bool isMC=false, string TreeName="O2rtdimuonall");
bool    checkDS(RooDataSet* DS, string DSName);
double  deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
bool    isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR=0.03);
double  getCorr(Double_t rapidity, Double_t pt, Double_t mass);
bool    readCorrection(const char* file);


bool tree2DataSet(RooWorkspace& Workspace, vector<string> InputFileNames, string DSName, string OutputFileName, bool UpdateDS=false)
{
  string TreeName("O2rtdimuonall");
//  string TreeName("DimuonsAll");
  RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL; RooDataSet* dataOSNoBkg = NULL;
  
  bool isMC = false;
  if (DSName.find("MC")!=std::string::npos) isMC =true;
  
  int triggerIndex_PP   = PP::HLT_HIL1DoubleMu0_v1;
  int CentFactor = 1;
  
  bool usePeriPD = false;
  
  bool applyWeight = false;
  
  bool isPureSDataset = false;
  if (OutputFileName.find("_PureS")!=std::string::npos) isPureSDataset = true;

  bool applyWeight_Corr = false;
  if ( (OutputFileName.find("_AccEff")!=std::string::npos) || (OutputFileName.find("_lJpsiEff")!=std::string::npos) ) applyWeight_Corr = true;
  if(applyWeight == true) applyWeight_Corr = false;
  
  TString corrName = "";
  TString corrFileName = "";
  if (OutputFileName.find("_AccEff")!=std::string::npos)
  {
    corrFileName = "correction_AccEff.root";
    corrName = "AccEff";
  }
  else if (OutputFileName.find("_lJpsiEff")!=std::string::npos)
  {
    corrFileName = "correction_lJpsiEff.root";
    corrName = "lJpsiEff";
  }
  
  bool createDS = ( gSystem->AccessPathName(OutputFileName.c_str()) || UpdateDS );
  if ( !gSystem->AccessPathName(OutputFileName.c_str()) ) {
    cout << "[INFO] Loading RooDataSet from " << OutputFileName << endl;
    
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"READ");
    if (isMC && isPureSDataset) {
      dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg", DSName.c_str()));
      if (checkDS(dataOSNoBkg, DSName)==false) { createDS = true; }
    } 
    else if (applyWeight_Corr) {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
    }
    else {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s", DSName.c_str()));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
      dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s", DSName.c_str()));
      if (checkDS(dataSS, DSName)==false) { createDS = true; }
    }
    DBFile->Close(); delete DBFile;
  }

  if (createDS) {
    cout << "[INFO] Creating " << (isPureSDataset ? "pure signal " : "") << "RooDataSet for " << DSName << endl;
//    TreeName = findMyTree(InputFileNames[0]); if(TreeName==""){return false;} 
    TChain* theTree = new TChain(TreeName.c_str(),"");
//    string TreeNameExtra = TreeName + "Extra";
//    TChain* theTreeExtra = new TChain(TreeNameExtra.c_str(),"");
    cout << "[INFO] Creating TChain" << endl;
    if(!getTChain(theTree, InputFileNames, TreeName)){ return false; }     // Import files to TChain
//    if(!getTChain(theTreeExtra, InputFileNames, TreeNameExtra)){ return false; }     // Import files to TChain
//    theTree->AddFriend(theTreeExtra);
    initOniaTree(theTree, TreeName);                                       // Initialize the Onia Tree
    iniBranch(theTree,isMC,TreeName);                                     // Initialize the Branches
    
    RooRealVar* mass    = new RooRealVar("invMass","#mu#mu mass", 1.0, 6.0, "GeV/c^{2}");
    RooRealVar* ctau    = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauN    = new RooRealVar("ctauN","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauTrue = new RooRealVar("ctauTrue","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauNRes = new RooRealVar("ctauNRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauRes = new RooRealVar("ctauRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauErr = new RooRealVar("ctauErr","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ptQQ    = new RooRealVar("pt","#mu#mu p_{T}", -1.0, 10000.0, "GeV/c");
    RooRealVar* rapQQ   = new RooRealVar("rap","#mu#mu y", -2.5, 2.5, "");
    RooRealVar* chi21   = new RooRealVar("chi21","#chi^2_{MFT-MCH}", -1000,1000000, "");
    RooRealVar* chi22   = new RooRealVar("chi22","#chi^2_{MFT-MCH}", -1000,1000000, "");
    RooRealVar* cent    = new RooRealVar("cent","centrality", -1.0, 1000.0, "");
    RooRealVar* weight  = new RooRealVar("weight","MC weight", 0.0, 10000000.0, "");
    RooRealVar* weightCorr   = new RooRealVar("weightCorr","Data correction weight", 0.0, 10000000.0, "");
    RooArgSet*  cols    = NULL;
    
    if (applyWeight_Corr)
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weightCorr, *chi21, *chi22);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
      } else {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weightCorr, *chi21, *chi22);
        cols->add(*ctauN);
      }
      if (!readCorrection(Form("%s/Input/%s",gSystem->ExpandPathName(gSystem->pwd()),corrFileName.Data()))){ return false; }
      dataOS = new RooDataSet(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()), "dOS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      //      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      cout<<"--- 1./applyWeight_Corr applied---"<<endl;
    }
    else
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *chi21, *chi22);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
      } else {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *chi21, *chi22);
        cols->add(*ctauN);
      }  
      dataOS = new RooDataSet(Form("dOS_%s", DSName.c_str()), "dOS", *cols, StoreAsymError(*mass));
      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, StoreAsymError(*mass));
      if (isMC && isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg", DSName.c_str()), "dOSNoBkg", *cols, StoreAsymError(*mass));
    }
    
    Long64_t nentries = theTree->GetEntries();
    //nentries = 50000;
    
    float normF = 0.;
    
    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries << endl;
      
      if (theTree->LoadTree(jentry)<0) break;
      if (theTree->GetTreeNumber()!=fCurrent) {
        fCurrent = theTree->GetTreeNumber();
//        cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
//        cout << "[INFO] Processing Root File: " << jentry << endl;
      }
/*      Reco_QQ_4mom->Clear();
      Reco_QQ_mumi_4mom->Clear();
      Reco_QQ_mupl_4mom->Clear();
      if (isMC) {
        Gen_QQ_mumi_4mom->Clear();
        Gen_QQ_mupl_4mom->Clear();
      }*/
      theTree->GetEntry(jentry);
      
      //for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
//        TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
//        mass->setVal(RecoQQ4mom->M());
              ROOT::Math::PtEtaPhiMVector v1_part(fPt1, fEtaMC1, fPhiMC1, 0.105658);
              ROOT::Math::PtEtaPhiMVector v2_part(fPt2, fEtaMC2, fPhiMC2, 0.105658);
              ROOT::Math::PtEtaPhiMVector v12_part = v1_part + v2_part;
//        mass->setVal(v12_part.M());
        mass->setVal(fMass);

/*        if (theTree->GetBranch("Reco_QQ_ctau3D")) { ctau->setVal(Reco_QQ_ctau3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctau")) { ctau->setVal(Reco_QQ_ctau[iQQ]); }
        else { cout << "[ERROR] No ctau information found in the Onia Tree" << endl; }
        if (theTree->GetBranch("Reco_QQ_ctauErr3D")) { ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctauErr")) { ctauErr->setVal(Reco_QQ_ctauErr[iQQ]); }
        else { cout << "[ERROR] No ctauErr information found in the Onia Tree" << endl; }
*/        

        float ct= fTauz*299792458.e-7*10.;
//	cout << "ctau = " << ct << ", tau = " << fTauz << endl;
	float ctErr= fTauzErr*299792458.e-7*10.;
	ctau->setVal(ct);
	ctauErr->setVal(ctErr);
	
//	cout << "[DEBUG] ctau : " << ct << ", ctauErr : " << ctErr << ", ctauN : " << ctau->getVal()/ctauErr->getVal() << endl;
//        ctauN->setVal(ctau->getVal()/ctauErr->getVal());
        ctauN->setVal(ctau->getVal());
	ROOT::Math::PtEtaPhiMVector v12(fPt, fEta, fPhi, fMass);
        
        ptQQ->setVal(fPt);
        chi21->setVal(fChi2MatchMCHMFT1);
        chi22->setVal(fChi2MatchMCHMFT2);
        rapQQ->setVal(v12.Rapidity());
        if (isMC && fMcDecision) {
//          if (theTree->GetBranch("Reco_QQ_ctauTrue3D")) { ctauTrue->setVal(Reco_QQ_ctauTrue3D[iQQ]); }
//          else if (theTree->GetBranch("Reco_QQ_ctauTrue")) { ctauTrue->setVal(Reco_QQ_ctauTrue[iQQ]); }
//          else { cout << "[ERROR] No ctauTrue information found in the Onia Tree" << endl; }

              ROOT::Math::PtEtaPhiMVector v1MC(fPtMC1, fEtaMC1, fPhiMC1, 0.105658);
              ROOT::Math::PtEtaPhiMVector v2MC(fPtMC2, fEtaMC2, fPhiMC2, 0.105658);
              ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
              float Tauz1MC = (fMCPosZ - fVz1) * v12MC.M()*1e1 / TMath::Abs(v12MC.Pz());
              float Tauz2MC = (fMCPosZ - fVz2) * v12MC.M()*1e1 / TMath::Abs(v12MC.Pz());


          ctauTrue->setVal((Tauz1MC+Tauz2MC)/2);
//          cout << "!!!!![DEBUG]!!!!! ctauTrue = " << ctauTrue->getVal() << endl;;
          ctauNRes->setVal( (ctauN->getValV() - ctauTrue->getValV())/(ctauErr->getValV()) );
          ctauRes->setVal( (ctau->getValV() - ctauTrue->getValV()) );
          ptQQ->setVal(v12MC.Pt());
          rapQQ->setVal(v12MC.Rapidity());

        }


        else if (applyWeight_Corr)
        {
          double Corr = getCorr(fEta,fPt,fMass);
          double wCorr = 1/Corr;
          weightCorr->setVal(wCorr);
        }
        float pplus = 0;
        float pminus = 0;
	if (fSign1 > 0 && fSign2 < 0){
		pplus = fPt1*cosh(fEta1);
		pminus = fPt2*cosh(fEta2);
	}else if (fSign1 < 0 && fSign2 > 0){
		pminus = fPt1*cosh(fEta1);
		pplus = fPt2*cosh(fEta2);
	}
	float delta = pplus - pminus;
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 1 && fPt2 > 1 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6){
	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -3.0 && fEta1 > - 3.6 && fEta2 < -3.0 && fEta2 > - 3.6){
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && !fIsAmbig1 && !fIsAmbig2){

//	if (fTauzErr>0 && !isnan(fTauzErr) && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fPt1 > 0.5 && fPt2 > 0.5 && fChi2MatchMCHMFT1 < 0 && fChi2MatchMCHMFT2 < 0){
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fChi2MatchMCHMFT1 > 0 && fChi2MatchMCHMFT2 > 0){
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fMcDecision){
//	if (fTauzErr>0 && !isnan(fTauzErr) &&  fChi2MatchMCHMFT1 < 25 && fChi2MatchMCHMFT2 < 25 && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6){
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fMcDecision && !(fMcMask1 & (0x1 << 7)) && !(fMcMask2 & (0x1 << 7))){
//	if (fTauzErr>0 && !isnan(fTauzErr) && fPt1 > 0.5 && fPt2 > 0.5 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fMcDecision){
//	if (fTauzErr>0 && !isnan(fTauzErr) && !(fMcMask1 & (0x1 << 7)) && !(fMcMask2 & (0x1 << 7)) && fPt1 > 0.5 && fPt2 > 0.5){
//	if (!(fMcMask1 & (0x1 << 7)) && !(fMcMask2 & (0x1 << 7)) && fMcDecision){
//	if (fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 < -2.5 && fEta2 > - 3.6 && fChi2MatchMCHMFT1 < 5 && fChi2MatchMCHMFT2 < 5){
//	if (fMcDecision && fPt1 > 1. && fPt2 > 1.){
//	&& !fIsAmbig1 && !fIsAmbig2
//
          if (fSign==0) { // Opposite-Sign dimuons
//		cout << "delta = " << delta << endl;
            if (isMC && isPureSDataset && isMatchedRecoDiMuon(0)) dataOSNoBkg->add(*cols, (applyWeight ? weight->getVal() : 1.0)); // Signal-only dimuons
            else if (applyWeight_Corr) dataOS->add(*cols,weightCorr->getVal()); //Signal and background dimuons
            else dataOS->add(*cols, ( applyWeight ? weight->getVal() : 1.0)); //Signal and background dimuons            
          }
          else { // Like-Sign dimuons
            if (!isPureSDataset && !applyWeight_Corr ) dataSS->add(*cols, ( applyWeight  ? weight->getVal() : 1.0));
          }
      }
    }
    // Close the TChain and all its pointers
//    delete Reco_QQ_4mom; delete Reco_QQ_mumi_4mom; delete Reco_QQ_mupl_4mom; delete Gen_QQ_mumi_4mom; delete Gen_QQ_mupl_4mom;
    theTree->Reset(); delete theTree;
    
    // Save all the datasets
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"RECREATE");
    DBFile->cd();
    if (isMC && isPureSDataset) dataOSNoBkg->Write(Form("dOS_%s_NoBkg", DSName.c_str()));
    else if (applyWeight_Corr) dataOS->Write(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()));
    else
    {
      dataOS->Write(Form("dOS_%s", DSName.c_str()));
      dataSS->Write(Form("dSS_%s", DSName.c_str()));
    }
    DBFile->Write(); DBFile->Close(); delete DBFile;
  }
  
  // Import datasets to workspace
  if (isMC && isPureSDataset)
  {
    if (!dataOSNoBkg) { cout << "[ERROR] " << DSName << "_NoBkg was not found" << endl; return false; }
    Workspace.import(*dataOSNoBkg);
  }
  else if (applyWeight_Corr)
  {
    if(!dataOS) { cout << "[ERROR] " << DSName << "_" << corrName.Data() << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
  }
  else
  {
    if(!dataOS || !dataSS) { cout << "[ERROR] " << DSName << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
    Workspace.import(*dataSS);
  }
  
  // delete the local datasets
  delete dataSS; delete dataOS; delete dataOSNoBkg;
  
  // delete the correction array
  if (fcorrArray) delete fcorrArray;
  
  return true;
};

string findMyTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  string name = "";
  if(f->GetListOfKeys()->Contains("hionia")){ name = "hionia/myTree"; }
  else if(f->GetListOfKeys()->Contains("O2rtdimuonall")){ name = "O2rtdimuonall"; }
//  else if(f->GetListOfKeys()->Contains("DimuonsAll")){ name = "DimuonsAll"; }
  else { cout << "[ERROR] myTree was not found in: " << FileName << endl; }
  f->Close(); delete f;
  return name;
};

bool getTChain(TChain *fChain, vector<string> FileNames, string TreeName)
{
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;

  for (vector<string>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
  TFile *inputFile = TFile::Open(FileName->c_str());
  cout << "opening file " << FileName->c_str() << endl;
  TIter keyList(inputFile->GetListOfKeys());
	cout << "got list of keys" << endl;
  TKey *key;

  while ((key = (TKey*)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TDirectoryFile")) continue;
    string dir = key->GetName();
	cout << "dir : " << dir.c_str() << endl;
    string nb = dir;
    nb.erase(0,9);
//    int nb_int = std::stoi(nb);
    //if (nb_int%2!=0 && nb_int < 9400000){
    cout << "file name : " << dir.c_str() << endl;
    cout << "[INFO] Adding TFile " << FileName->c_str() << dir.c_str() << endl;
    fChain->Add(Form("%s/%s/%s",FileName->c_str(),dir.c_str(),TreeName.c_str()));
    //}
  }
  }
  if (!fChain) { cout << "[ERROR] fChain was not created, some input files are missing" << endl; return false; }
  return true;
};

void iniBranch(TChain* fChain, bool isMC, string TreeName)
{
  cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
/*  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(false); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(false); }
  if (isMC) {
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->GetBranch("Gen_QQ_mupl_4mom")->SetAutoDelete(false); }
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->GetBranch("Gen_QQ_mumi_4mom")->SetAutoDelete(false); }
  }
  fChain->SetBranchStatus("*",0);
  RecoQQ::iniBranches(fChain);
  if (fChain->GetBranch("Centrality"))        { fChain->SetBranchStatus("Centrality",1);        }
  if (fChain->GetBranch("Reco_QQ_size"))      { fChain->SetBranchStatus("Reco_QQ_size",1);      }
  if (fChain->GetBranch("Reco_QQ_sign"))      { fChain->SetBranchStatus("Reco_QQ_sign",1);      }
  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->SetBranchStatus("Reco_QQ_4mom",1);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_ctau3D"))    { fChain->SetBranchStatus("Reco_QQ_ctau3D",1);    }
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) { fChain->SetBranchStatus("Reco_QQ_ctauErr3D",1); }
  if (fChain->GetBranch("Reco_QQ_ctau"))      { fChain->SetBranchStatus("Reco_QQ_ctau",1);      }
  if (fChain->GetBranch("Reco_QQ_ctauErr"))   { fChain->SetBranchStatus("Reco_QQ_ctauErr",1);   }
  if (isMC)
  {
    if (fChain->GetBranch("Gen_QQ_size"))      { fChain->SetBranchStatus("Gen_QQ_size",1);      }
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1); }
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue3D")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue3D",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue",1); }
  }
*/
 if (fChain->GetBranch("fMass"))   fChain->GetBranch("fMass")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt"))   fChain->GetBranch("fPt")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta"))   fChain->GetBranch("fEta")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi"))   fChain->GetBranch("fPhi")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign"))   fChain->GetBranch("fSign")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauz"))  fChain->GetBranch("fTauz")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauzErr"))   fChain->GetBranch("fTauzErr")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauxy"))  fChain->GetBranch("fTauxy")->SetAutoDelete(false);
 if (fChain->GetBranch("fTauxyErr"))  fChain->GetBranch("fTauxyErr")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosX"))   fChain->GetBranch("fPosX")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosY"))   fChain->GetBranch("fPosY")->SetAutoDelete(false);
 if (fChain->GetBranch("fPosZ"))   fChain->GetBranch("fPosZ")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosX"))   fChain->GetBranch("fMCPosX")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosY"))   fChain->GetBranch("fMCPosY")->SetAutoDelete(false);
 if (fChain->GetBranch("fMCPosZ"))   fChain->GetBranch("fMCPosZ")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt1"))   fChain->GetBranch("fPt1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta1"))   fChain->GetBranch("fEta1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi1"))   fChain->GetBranch("fPhi1")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign1"))   fChain->GetBranch("fSign1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPt2"))   fChain->GetBranch("fPt2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEta2"))   fChain->GetBranch("fEta2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhi2"))   fChain->GetBranch("fPhi2")->SetAutoDelete(false);
 if (fChain->GetBranch("fSign2"))   fChain->GetBranch("fSign2")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcMask1"))   fChain->GetBranch("fMcMask1")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcMask2"))   fChain->GetBranch("fMcMask2")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMID1"))   fChain->GetBranch("fChi2MatchMCHMID1")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMID2"))   fChain->GetBranch("fChi2MatchMCHMID2")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMFT1"))   fChain->GetBranch("fChi2MatchMCHMFT1")->SetAutoDelete(false);
 if (fChain->GetBranch("fChi2MatchMCHMFT2"))   fChain->GetBranch("fChi2MatchMCHMFT2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPtMC1"))   fChain->GetBranch("fPtMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEtaMC1"))   fChain->GetBranch("fEtaMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhiMC1"))   fChain->GetBranch("fPhiMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fEMC1"))   fChain->GetBranch("fEMC1")->SetAutoDelete(false);
 if (fChain->GetBranch("fPtMC2"))   fChain->GetBranch("fPtMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEtaMC2"))   fChain->GetBranch("fEtaMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fPhiMC2"))   fChain->GetBranch("fPhiMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fEMC2"))   fChain->GetBranch("fEMC2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVx1"))   fChain->GetBranch("fVx1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVy1"))   fChain->GetBranch("fVy1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVz1"))   fChain->GetBranch("fVz1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVt1"))   fChain->GetBranch("fVt1")->SetAutoDelete(false);
 if (fChain->GetBranch("fVx2"))   fChain->GetBranch("fVx2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVy2"))   fChain->GetBranch("fVy2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVz2"))   fChain->GetBranch("fVz2")->SetAutoDelete(false);
 if (fChain->GetBranch("fVt2"))   fChain->GetBranch("fVt2")->SetAutoDelete(false);
 if (fChain->GetBranch("fMcDecision"))   fChain->GetBranch("fMcDecision")->SetAutoDelete(false);

 if (fChain->GetBranch("fIsAmbig1"))   fChain->GetBranch("fIsAmbig1")->SetAutoDelete(false);
 if (fChain->GetBranch("fIsAmbig2"))   fChain->GetBranch("fIsAmbig2")->SetAutoDelete(false);


 if (fChain->GetBranch("fMass"))   fChain->SetBranchStatus("fMass", 1);
 if (fChain->GetBranch("fPt"))   fChain->SetBranchStatus("fPt", 1);
 if (fChain->GetBranch("fEta"))   fChain->SetBranchStatus("fEta", 1);
 if (fChain->GetBranch("fPhi"))   fChain->SetBranchStatus("fPhi", 1);
 if (fChain->GetBranch("fSign"))   fChain->SetBranchStatus("fSign", 1);
 if (fChain->GetBranch("fTauz"))  fChain->SetBranchStatus("fTauz", 1);
 if (fChain->GetBranch("fTauzErr"))   fChain->SetBranchStatus("fTauzErr", 1);
 if (fChain->GetBranch("fTauxy"))  fChain->SetBranchStatus("fTauxy", 1);
 if (fChain->GetBranch("fTauxyErr"))  fChain->SetBranchStatus("fTauxyErr", 1);
 if (fChain->GetBranch("fPosX"))   fChain->SetBranchStatus("fPosX", 1);
 if (fChain->GetBranch("fPosY"))   fChain->SetBranchStatus("fPosY", 1);
 if (fChain->GetBranch("fPosZ"))   fChain->SetBranchStatus("fPosZ", 1);
 if (fChain->GetBranch("fMCPosX"))   fChain->SetBranchStatus("fMCPosX", 1);
 if (fChain->GetBranch("fMCPosY"))   fChain->SetBranchStatus("fMCPosY", 1);
 if (fChain->GetBranch("fMCPosZ"))   fChain->SetBranchStatus("fMCPosZ", 1);
 if (fChain->GetBranch("fPt1"))   fChain->SetBranchStatus("fPt1", 1);
 if (fChain->GetBranch("fEta1"))   fChain->SetBranchStatus("fEta1", 1);
 if (fChain->GetBranch("fPhi1"))   fChain->SetBranchStatus("fPhi1", 1);
 if (fChain->GetBranch("fSign1"))   fChain->SetBranchStatus("fSign1", 1);
 if (fChain->GetBranch("fPt2"))   fChain->SetBranchStatus("fPt2", 1);
 if (fChain->GetBranch("fEta2"))   fChain->SetBranchStatus("fEta2", 1);
 if (fChain->GetBranch("fPhi2"))   fChain->SetBranchStatus("fPhi2", 1);
 if (fChain->GetBranch("fSign2"))   fChain->SetBranchStatus("fSign2", 1);
 if (fChain->GetBranch("fMcMask1"))   fChain->SetBranchStatus("fMcMask1", 1);
 if (fChain->GetBranch("fMcMask2"))   fChain->SetBranchStatus("fMcMask2", 1);
 if (fChain->GetBranch("fChi2MatchMCHMID1"))   fChain->SetBranchStatus("fChi2MatchMCHMID1", 1);
 if (fChain->GetBranch("fChi2MatchMCHMID2"))   fChain->SetBranchStatus("fChi2MatchMCHMID2", 1);
 if (fChain->GetBranch("fChi2MatchMCHMFT1"))   fChain->SetBranchStatus("fChi2MatchMCHMFT1", 1);
 if (fChain->GetBranch("fChi2MatchMCHMFT2"))   fChain->SetBranchStatus("fChi2MatchMCHMFT2", 1);
 if (fChain->GetBranch("fPtMC1"))   fChain->SetBranchStatus("fPtMC1", 1);
 if (fChain->GetBranch("fEtaMC1"))   fChain->SetBranchStatus("fEtaMC1", 1);
 if (fChain->GetBranch("fPhiMC1"))   fChain->SetBranchStatus("fPhiMC1", 1);
 if (fChain->GetBranch("fEMC1"))   fChain->SetBranchStatus("fEMC1", 1);
 if (fChain->GetBranch("fPtMC2"))   fChain->SetBranchStatus("fPtMC2", 1);
 if (fChain->GetBranch("fEtaMC2"))   fChain->SetBranchStatus("fEtaMC2", 1);
 if (fChain->GetBranch("fPhiMC2"))   fChain->SetBranchStatus("fPhiMC2", 1);
 if (fChain->GetBranch("fEMC2"))   fChain->SetBranchStatus("fEMC2", 1);
 if (fChain->GetBranch("fVx1"))   fChain->SetBranchStatus("fVx1", 1);
 if (fChain->GetBranch("fVy1"))   fChain->SetBranchStatus("fVy1", 1);
 if (fChain->GetBranch("fVz1"))   fChain->SetBranchStatus("fVz1", 1);
 if (fChain->GetBranch("fVt1"))   fChain->SetBranchStatus("fVt1", 1);
 if (fChain->GetBranch("fVx2"))   fChain->SetBranchStatus("fVx2", 1);
 if (fChain->GetBranch("fVy2"))   fChain->SetBranchStatus("fVy2", 1);
 if (fChain->GetBranch("fVz2"))   fChain->SetBranchStatus("fVz2", 1);
 if (fChain->GetBranch("fVt2"))   fChain->SetBranchStatus("fVt2", 1);
 if (fChain->GetBranch("fMcDecision"))   fChain->SetBranchStatus("fMcDecision", 1);

 if (fChain->GetBranch("fIsAmbig1"))   fChain->SetBranchStatus("fIsAmbig1", 1);
 if (fChain->GetBranch("fIsAmbig2"))   fChain->SetBranchStatus("fIsAmbig2", 1);

};

bool checkDS(RooDataSet* DS, string DSName)
{
  bool incCtauTrue = (DSName.find("MC")!=std::string::npos);
  const RooArgSet* row = DS->get();
  if (
      (row->find("invMass")!=0) &&
      (row->find("pt")!=0)      &&
      (row->find("ctau")!=0)    &&
      (row->find("ctauErr")!=0) &&
      (incCtauTrue ? row->find("ctauTrue")!=0 : true) &&
      (incCtauTrue ? row->find("ctauRes")!=0 : true) &&
      (incCtauTrue ? row->find("ctauNRes")!=0 : true) &&
      (incCtauTrue ? true : row->find("ctauN")!=0)
      ) 
    { return true; }
  else 
    { cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << endl; }

  return false;
};

double deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

bool isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
{
/*  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);
*/  
  bool isMatched(false);
  int iGenMuon(0);
/*  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
  {
    TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
    TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) isMatched = true;
    iGenMuon++;
  }
*/  
  return isMatched;
};

bool readCorrection(const char* file)
{
  TFile *froot = new TFile(file,"READ");
  if (!froot)
  {
    cout << "[ERROR] File "<< file << " for correction of events not found" << endl;
    return false;
  }
  
  TList* lcorr = froot->GetListOfKeys();
  TIter nextCorr(lcorr);
  
  fcorrArray = new TObjArray();
  fcorrArray->SetOwner(kTRUE);
  
  TObjString* fname(0x0);
  while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
  {
    TH2* h = static_cast<TH2*>(froot->FindObjectAny(fname->GetString().Data()));

    TString sName(h->GetName());
    if ( sName.Contains("hcorr") ) fcorrArray->Add(h->Clone());
    else cout << "[WARNING] Correction histo " << sName.Data() << " not according to naming convention. Not included in correction array" << endl;
  }
  
  if (!(fcorrArray->GetSize()>0))
  {
    cout << "[ERROR] Correction array empty: No corrections found." << endl;
    return false;
  }
  delete lcorr;
//  froot->Close(); delete froot;
  
  return true;
};

double getCorr(Double_t rapidity, Double_t pt, Double_t mass)
{
  const char* collName = "PP";
  const char* massName = "Interp";
  if (mass<3.3) massName = "Jpsi";
  
  if (!fcorrArray)
  {
    cout << "[ERROR] No correction array exist" << endl;
    return 0;
  }

  Double_t corr = 1.;
  if (!strcmp(massName,"Interp"))
  {
    TH2* corrHistoJpsi = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s",collName)));
    if (!corrHistoJpsi)
    {
      std::cout << "[Error] No histogram provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
      return 1.;
    }
    
    Int_t binJpsi = corrHistoJpsi->FindBin(fabs(rapidity), pt);
    Double_t corrJpsi = corrHistoJpsi->GetBinContent(binJpsi);
    
  }
  else
  {
    TH2* corrHisto = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_%s_%s",massName,collName)));
    if (!corrHisto)
    {
      std::cout << "[Error] No histogram provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
      return 1.;
    }
    
    Int_t bin = corrHisto->FindBin(fabs(rapidity), pt);
    corr = corrHisto->GetBinContent(bin);
  }
  if(corr<0.00001) corr=1.0;
  
  return corr;
};
