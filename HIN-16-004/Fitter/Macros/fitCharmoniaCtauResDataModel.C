#ifndef fitCharmoniaCtauResDataModel_C
#define fitCharmoniaCtauResDataModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauResModel.C"
#include "fitCharmoniaMassModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "drawCtauResDataPlot.C"


void setCtauResDataCutParameters(struct KinCuts& cut);
void setCtauResDataFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut);
void setCtauResDataGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth);
bool setCtauResDataModel( struct OniaModel& model, map<string, string>&  parIni);
bool createSignalCtauDSUsingSPLOT(RooWorkspace& ws, string dsName, map<string, string>  parIni, struct KinCuts cut, bool incJpsi, bool incBkg, bool useSPlot);

int importDataset_Res(RooWorkspace& myws, const RooWorkspace& inputWS, struct KinCuts cut, string label, bool cutSideBand=false)
{
  cout << "Import data set" << endl;
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  if (cutSideBand) {
    indMuonMass =  indMuonMass + "&&" + "((2.0 < invMass && invMass < 2.8) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
  }
  string indMuonRap     = Form("(%.6f <= abs(rap) && abs(rap) < %.6f)",    cut.dMuon.AbsRap.Min,   cut.dMuon.AbsRap.Max);
  string indMuonChi21     = Form("(%.6f <= chi21 && chi21 < %.6f)",    cut.dMuon.Chi2.Min,   cut.dMuon.Chi2.Max);
  string indMuonChi22     = Form("(%.6f <= chi22 && chi22 < %.6f)",    cut.dMuon.Chi2.Min,   cut.dMuon.Chi2.Max);
  string indMuonPt      = Form("(%.6f <= pt && pt < %.6f)",                cut.dMuon.Pt.Min,       cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.6f < ctau && ctau <= %.6f)",            cut.dMuon.ctau.Min,     cut.dMuon.ctau.Max);
  if(cut.dMuon.ctauCut!=""){ indMuonCtau = cut.dMuon.ctauCut; }
  string indMuonCtauErr = Form("(%.12f < ctauErr && ctauErr < %.12f)",     cut.dMuon.ctauErr.Min,  cut.dMuon.ctauErr.Max);
  string inCentrality   = Form("(%d <= cent && cent < %d)",                cut.Centrality.Start,   cut.Centrality.End);
  string indMuonCtauTrue = Form("(%.12f < ctauTrue && ctauTrue < %.12f)",  cut.dMuon.ctauTrue.Min, cut.dMuon.ctauTrue.Max);
  string indMuonCtauRes = Form("(%.12f < ctauRes && ctauRes < %.12f)",     cut.dMuon.ctauRes.Min,  cut.dMuon.ctauRes.Max);
  string indMuonCtauNRes = Form("(%.12f < ctauNRes && ctauNRes < %.12f)",  cut.dMuon.ctauNRes.Min, cut.dMuon.ctauNRes.Max);
  string indMuonCtauN   = Form("(%.12f < ctauN && ctauN < %.12f)",        cut.dMuon.ctauN.Min, cut.dMuon.ctauN.Max);
  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr + "&&" + indMuonChi21 + "&&" + indMuonChi22;
  if (label.find("MC")!=std::string::npos){ strCut = strCut +"&&"+ indMuonCtauTrue +"&&"+ indMuonCtauNRes +"&&"+ indMuonCtauRes;  }
  else { strCut = strCut +"&&"+ indMuonCtauN;  }

  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return -1;
  }
  strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt + "&&" + indMuonChi21 + "&&" + indMuonChi22;
  cout << "[INFO] Importing local RooDataSet with cuts: " << strCut << endl;
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataOS->sumEntries()==0){
    cout << "[ERROR] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return -1;
  }
  myws.import(*dataOS);
  delete dataOS;

  if (label.find("NoBkg")==std::string::npos && label.find("AccEff")==std::string::npos && label.find("lJpsiEff")==std::string::npos) // Don't try to find SS dataset if label contais NoBkg or correction
  {
    if (!(inputWS.data(Form("dSS_%s", label.c_str())))){
      cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
      return -1;
    }
    RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
    if (dataSS->sumEntries()==0){
      cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    }
    myws.import(*dataSS);
    delete dataSS;
  }

  const RooArgSet* rowOS = myws.data(Form("dOS_%s", label.c_str()))->get();
  ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);
  ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
  ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);
  ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
  ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);
  ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
  if (label.find("MC")!=std::string::npos){
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
    ((RooRealVar*)rowOS->find("ctauRes"))->setMin(cut.dMuon.ctauRes.Min);
    ((RooRealVar*)rowOS->find("ctauRes"))->setMax(cut.dMuon.ctauRes.Max);
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMin(cut.dMuon.ctauNRes.Min);
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMax(cut.dMuon.ctauNRes.Max);
  }
  else {
    ((RooRealVar*)rowOS->find("ctauN"))->setMin(cut.dMuon.ctauN.Min);
    ((RooRealVar*)rowOS->find("ctauN"))->setMax(cut.dMuon.ctauN.Max);
  }

  if (myws.data(Form("dOS_%s_SPLOT_INPUT", label.c_str()))){
    RooDataSet* dataOS = (RooDataSet*)myws.data(Form("dOS_%s_SPLOT_INPUT", label.c_str()))->reduce(strCut.c_str());
    if (dataOS) {
      const RooArgSet* rowOS = dataOS->get();
      ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);
      ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
      ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);
      ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
      ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);
      ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
      ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);
      ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
      if (label.find("MC")!=std::string::npos){
        ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);
        ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
        ((RooRealVar*)rowOS->find("ctauRes"))->setMin(cut.dMuon.ctauRes.Min);
        ((RooRealVar*)rowOS->find("ctauRes"))->setMax(cut.dMuon.ctauRes.Max);
        ((RooRealVar*)rowOS->find("ctauNRes"))->setMin(cut.dMuon.ctauNRes.Min);
        ((RooRealVar*)rowOS->find("ctauNRes"))->setMax(cut.dMuon.ctauNRes.Max);
      }
      else {
        ((RooRealVar*)rowOS->find("ctauN"))->setMin(cut.dMuon.ctauN.Min);
        ((RooRealVar*)rowOS->find("ctauN"))->setMax(cut.dMuon.ctauN.Max);
      }
      if (dataOS->sumEntries()==0){ cout << "[ERROR] No events from dataset " <<  Form("dOS_%s_SPLOT_INPUT", label.c_str()) << " passed the kinematic cuts!" << endl; }
      else if (!isCompatibleDataset(*dataOS, *(RooDataSet*)myws.data(Form("dOS_%s", label.c_str())))){ cout << "[ERROR] sPlot and Original Datasets are inconsistent!" << endl; delete dataOS; return -1; }
      else {
        myws.import(*dataOS, Rename(Form("dOS_%s_SPLOT", label.c_str())));
        if (myws.data(Form("dOS_%s_SPLOT", label.c_str()))) { cout << "[INFO] RooDataset " << Form("dOS_%s_SPLOT", label.c_str()) << " was imported!" << endl; }
        else { cout << "[ERROR] Importing RooDataset " << Form("dOS_%s_SPLOT", label.c_str()) << " failed!" << endl; delete dataOS; return -1; }
        cout << "[INFO] SPlotDS Events: " << dataOS->sumEntries() << " , origDS Events: " << myws.data(Form("dOS_%s", label.c_str()))->sumEntries() << std::endl;
      }
      delete dataOS;
    }
  }

  myws.var("invMass")->setMin(cut.dMuon.M.Min);
  myws.var("invMass")->setMax(cut.dMuon.M.Max);
  myws.var("pt")->setMin(cut.dMuon.Pt.Min);
  myws.var("pt")->setMax(cut.dMuon.Pt.Max);
  myws.var("rap")->setMin(cut.dMuon.AbsRap.Min);
  myws.var("rap")->setMax(cut.dMuon.AbsRap.Max);
  myws.var("chi21")->setMin(cut.dMuon.Chi2.Min);
  myws.var("chi21")->setMax(cut.dMuon.Chi2.Max);
  myws.var("chi22")->setMin(cut.dMuon.Chi2.Min);
  myws.var("chi22")->setMax(cut.dMuon.Chi2.Max);
  myws.var("ctau")->setMin(cut.dMuon.ctau.Min);
  myws.var("ctau")->setMax(cut.dMuon.ctau.Max);
  myws.var("ctauErr")->setMin(cut.dMuon.ctauErr.Min);
  myws.var("ctauErr")->setMax(cut.dMuon.ctauErr.Max);
  if (label.find("MC")!=std::string::npos){
    myws.var("ctauTrue")->setMin(cut.dMuon.ctauTrue.Min);
    myws.var("ctauTrue")->setMax(cut.dMuon.ctauTrue.Max);
    myws.var("ctauRes")->setMin(cut.dMuon.ctauRes.Min);
    myws.var("ctauRes")->setMax(cut.dMuon.ctauRes.Max);
    myws.var("ctauNRes")->setMin(cut.dMuon.ctauNRes.Min);
    myws.var("ctauNRes")->setMax(cut.dMuon.ctauNRes.Max);
}
  else {
    myws.var("ctauN")->setMin(cut.dMuon.ctauN.Min);
    myws.var("ctauN")->setMax(cut.dMuon.ctauN.Max);
  }
  cout << "[INFO] Analyzing bin: " << Form(
                                           "%.3f < pt < %.3f, %.3f < rap < %.3f, %d < cent < %d, %.3f < chi2 < %.3f",
                                           cut.dMuon.Pt.Min,
                                           cut.dMuon.Pt.Max,
                                           cut.dMuon.AbsRap.Min,
                                           cut.dMuon.AbsRap.Max,
                                           cut.Centrality.Start,
                                           cut.Centrality.End,
                                           cut.dMuon.Chi2.Min,
                                           cut.dMuon.Chi2.Max
                                           ) << endl;

  return 1;
};


bool loadCtauErrRange_ctaures(string FileName, struct KinCuts& cut)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[ERROR] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
        cout << "File accessed" << endl;
  TFile *file = new TFile(FileName.c_str());
        cout << "File created " << FileName << endl;
  if (!file) return false;
  RooWorkspace *wksp = (RooWorkspace*) file->Get("workspace");
  cout << "Ws created" << endl;
  if (!wksp) {
    cout << "[ERROR] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (wksp->var("ctauErr")) {
        cout << "var accessed" << endl;
    cut.dMuon.ctauErr.Min = wksp->var("ctauErr")->getMin();
    cut.dMuon.ctauErr.Max = wksp->var("ctauErr")->getMax();
  } else {
    cout << Form("[ERROR] ctauErr was not found!") << endl;
    delete wksp;
    file->Close(); delete file;
    return false;
  }

  delete wksp;
  file->Close(); delete file;
  return true;
};

bool fitCharmoniaCtauResDataModel( RooWorkspace& myws,             // Local Workspace
                               const RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                               struct KinCuts& cut,            // Variable containing all kinematic cuts
                               map<string, string>&  parIni,   // Variable containing all initial parameters
                               struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                               string outputDir,               // Path to output directory
                               // Select the type of datasets to fit
                               string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                               bool importDS      = true,      // Select if the dataset is imported in the local workspace
                               // Select the type of object to fit
                               bool incJpsi       = true,      // Includes Jpsi model
                               bool incBkg        = true,      // Includes Bkg model
                               // Select the fitting options
                               bool useSPlot      = true,      // If yes, then use SPlot technique, if no, use mass range
                               bool useTotctauErrPdf = false,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                               bool doFit         = true,      // Flag to indicate if we want to perform the fit
                               bool loadFitResult = false,     // Load previous fit results
                               map<string, string> inputFitDir={},// User-defined Location of the fit results
                               int  numCores      = 2,         // Number of cores used for fitting
                               // Select the drawing options
                               bool setLogScale   = true,      // Draw plot with log scale
                               bool incSS         = false,     // Include Same Sign data
                               map<string, double> binWidth={} // User-defined Location of the fit results
                               )  
{

//  bool usePerEventError = true;
  bool usePerEventError = false;
  incBkg=true;
//  useSPlot=false;

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  if (DSTAG.find("MC")!=std::string::npos) {
    cout << "[ERROR] We can only fit resolution from Data using Data" << endl; return false;
  }
  
  string COLL = "PP";
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));

  if (importDS) {
	cout << "import DS" << endl;
    setMassCutParameters(cut, incJpsi, false, false);
    setCtauResDataCutParameters(cut);
    if (usePerEventError) {
      // check if we have already done the ctauErr fits. If yes, load their parameters
      string FileName = "";
      string pdfName = "pdfCTAUERR_Tot_PP";
      string plotLabel = "";
      bool incJpsi = true;
      bool fitSideBand = false;
      if (incJpsi)  { plotLabel = plotLabel + "_Jpsi";     }
      plotLabel = plotLabel + "_Bkg";
      setCtauErrFileName(FileName, (inputFitDir["CTAUERR"]=="" ? outputDir : inputFitDir["CTAUERR"]), "DATA", plotLabel, cut, fitSideBand);
      bool foundFit = false;
      if ( loadCtauErrRange_ctaures(FileName, cut) ) { foundFit = true; }
      if (foundFit) { cout << "[INFO] The ctauErr fit was found and I'll load the ctau Error range used." << endl; }
//      else { cout << "[ERROR] The ctauErr fit was not found!" << endl; return false; }
      else { cout << "[ERROR] The ctauErr fit was not found!" << endl; }
      // Importing SPLOT DS from ctauErr results
      string sPlotDSName = Form("dOS_%s_SPLOT", label.c_str());
      if (useSPlot) { loadSPlotDS(myws, FileName, sPlotDSName); }
    }
    setCtauErrCutParameters(cut);
  }

  // Import the local datasets
  double numEntries = 1000000;
  string dsName = Form("dOS_%s", label.c_str());
  if (importDS) {
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset_Res(myws, inputWorkspace, cut, label);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
  }
  else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }

  string dsName2Fit = dsName;
  if (incJpsi) dsName2Fit += "_JPSI";
  else if (incBkg) dsName2Fit += "_BKG";
  string dsName2FitCut = dsName2Fit+"_CTAUNCUT";
  cout << "!!!!![DEBUG]!!!!! dsName2Fit= " << dsName2Fit << "; dsName2FitCut=" << dsName2FitCut << endl;

  // Define pdf and plot names
  string plotLabel = "";
  plotLabel = plotLabel + Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
  string pdfName = "pdfCTAURES_Jpsi_PP";

  string FileName = "";
  setCtauResDataFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    if (loadFitResult) return false;
    setCtauResDataFileName(FileName, outputDir, DSTAG, plotLabel, cut);
  }

  // Check if we have already made the Signal DataSet
  vector<string> dsNames = { dsName2Fit };
  bool createSignalDS = ( !isSPlotDSAlreadyFound(myws, FileName, dsNames, true) );
  bool compDS = loadYields(myws, FileName, dsName, pdfName);
  if (!compDS && (incJpsi || incBkg)) {
    // Setting extra input information needed by each fitter
    string iMassFitDir = inputFitDir["MASS"];
    double ibWidth = binWidth["MASS"];
    bool loadMassFitResult = true;
    bool doMassFit = false;
    bool importDS = true;
    bool getMeanPT = false;
    const char* applyCorr = "";
    bool cutCtau = false;
    bool doConstrFit = false;
    cout << "[DEBUG] In mass fit, include Bkg ? " << incBkg << endl;    
    if ( !fitCharmoniaMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir,
                                 DSTAG, importDS,
                                 true, incBkg,
                                 doMassFit, cutCtau, doConstrFit, false, applyCorr, loadMassFitResult, iMassFitDir, numCores,
                                 false, incSS, ibWidth, getMeanPT
                                 )
         ) { return false; }
    // Let's set all mass parameters to constant except the yields
    if (myws.pdf("pdfMASS_Tot_PP")) {
      cout << "[INFO] Setting mass parameters to constant!" << endl;
      myws.pdf("pdfMASS_Tot_PP")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
    } else { cout << "[ERROR] Mass PDF was not found!" << endl; return false; }
    std::vector< std::string > objs = {"Bkg", "Jpsi"};
    for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), "PP"))){ setConstant( myws, Form("N_%s_%s", obj.c_str(), "PP"), false); cout << Form("N_%s_%s", obj.c_str(), "PP") << " set constant to " << myws.var(Form("N_%s_%s", obj.c_str(), "PP"))->getVal() << endl;}}
    // Let's set the minimum value of the yields to zero
    for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), "PP"))) myws.var(Form("N_%s_%s", obj.c_str(), "PP"))->setMin(0.0); }
  }



  if (createSignalDS) {
    if (importDS && !myws.data(dsName2Fit.c_str())) {
      string pdfMassName = "pdfMASS_Tot_PP";
      if (!createSignalCtauDSUsingSPLOT(myws, dsName, parIni, cut, incJpsi, incBkg, useSPlot)){ cout << "[ERROR] Creating the Signal Ctau Resolution Templates using sPLOT failed" << endl; return false; }
    }
  }
  if (importDS) {
    // Set global parameters
    setCtauResDataGlobalParameterRange(myws, parIni, cut, dsName2Fit, binWidth["CTAURES"]);
    if (!myws.data(dsName2FitCut.c_str())) {
      RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName2Fit.c_str())->reduce(parIni["CtauNRange_Cut"].c_str()))->Clone(dsName2FitCut.c_str());
      myws.import(*dataToFit, Rename(dsName2FitCut.c_str()));
      cout << "!!!!![DEBUG]!!!!! CtauNRange_Cut : " << parIni["CtauNRange_Cut"].c_str() << endl; 
      cout << "!!!!![DEBUG]!!!!! dsName2Fit : " << dsName2Fit <<  ", num entries : " << myws.data(dsName2Fit.c_str())->sumEntries() << endl;
      cout << "!!!!![DEBUG]!!!!! dsName2FitCut : " << dsName2FitCut <<  ", num entries : " << myws.data(dsName2FitCut.c_str())->sumEntries() << endl;
    }
  }
//  myws.var("N_Jpsi_PP")->setVal(100);
//  myws.var("N_Jpsi_PP")->setVal(myws.data(dsName2FitCut.c_str())->sumEntries());
//  myws.var("N_Jpsi_PP")->setMax(myws.data(dsName2FitCut.c_str())->sumEntries()*0.5);

  if (!loadFitResult) {
    // Set models based on initial parameters
    struct OniaModel model;
    if (!setCtauResDataModel(model, parIni)) { return false; }

    //// LOAD CTAU ERROR PDF
    if (usePerEventError) {
      // Setting extra input information needed by each fitter
      bool loadCtauErrFitResult = true;
      bool doCtauErrFit = true;
      bool importDS = true;
      bool incJpsi = true;
      string DSTAG = "DATA_PP";
  
      if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                      DSTAG, importDS, 
                                      incJpsi, incBkg, 
                                      doCtauErrFit, false, loadCtauErrFitResult, inputFitDir, numCores,
                                      true, incSS, binWidth
                                      ) 
           ) { return false; }
    }

    // Build the Fit Model
    if (!buildCharmoniaCtauResModel(myws, model.PP, parIni, dsName2FitCut, "ctau", "pdfCTAURES", usePerEventError, useTotctauErrPdf, numEntries))  { return false; }
    if (!buildCharmoniaCtauResModel(myws, model.PP, parIni, dsName2FitCut, "ctauN", "pdfCTAUNRES", false, useTotctauErrPdf, numEntries))  { return false; }
    cout << "Building res model" << endl;
//    if (!buildCharmoniaCtauResModel(myws, model.PP, parIni, dsName2FitCut, "ctau", "pdfCTAURES", false, false, numEntries))  { return false; }
//    if (!buildCharmoniaCtauResModel(myws, model.PP, parIni, dsName2FitCut, "ctauN", "pdfCTAUNRES", false, false, numEntries))  { return false; }

    // save the initial values of the model we've just created
    RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauN"), *myws.var("ctauErr"), *myws.var("ctau")));
    myws.saveSnapshot(Form("%s_parIni", pdfName.c_str()),*params,kTRUE);
    delete params;
  }

  // check if we have already done this fit. If yes, do nothing and return true.
  bool found =  true; bool skipFit = !doFit;
  RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"), *myws.var("ctauN")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (loadFitResult) {
    if ( loadPreviousFitResult(myws, FileName, DSTAG, false, false) ) { skipFit = true; } else  { skipFit = false; }
    if (skipFit) { cout << "[INFO] This ctau fit was already done, so I'll load the fit results." << endl; }
    myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()),*newpars,kTRUE);
  } else if (found) {
    cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
    return true;
  }
  
  // Fit the Datasets
  if (skipFit==false) {
    bool isWeighted = myws.data(dsName2FitCut.c_str())->isWeighted();
    cout << "!!!!![DEBUG]!!!!! fitting pdf " << pdfName.c_str() << endl;
    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName2FitCut.c_str()), Extended(kTRUE), NumCPU(numCores), SumW2Error(isWeighted), Save());
    fitResult->Print("v");
    myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));

    // Draw the mass plot
    if(!drawCtauResDataPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, dsName2Fit.c_str(), setLogScale, incSS, binWidth["CTAURES"])) { return false; }
    // Save the results
    string FileName = ""; setCtauResDataFileName(FileName, outputDir, DSTAG, plotLabel, cut);
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctauRes/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);
  }
  
  return true;
};


bool setCtauResDataModel( struct OniaModel& model, map<string, string>&  parIni)
{
  map< string , CtauModel > CtauModelDictionary = {
  {"InvalidModel",             CtauModel::InvalidModel},
  {"QuadrupleGaussianResolution", CtauModel::QuadrupleGaussianResolution},
  {"TripleGaussianResolution", CtauModel::TripleGaussianResolution},
  {"DoubleGaussianResolution", CtauModel::DoubleGaussianResolution},
  {"SingleGaussianResolution", CtauModel::SingleGaussianResolution},
  {"TripleDecay",              CtauModel::TripleDecay},
  {"QuadrupleDecay",           CtauModel::QuadrupleDecay},
  {"DoubleSingleSidedDecay",   CtauModel::DoubleSingleSidedDecay},
  {"SingleSidedDecay",         CtauModel::SingleSidedDecay},
  {"Delta",                    CtauModel::Delta},
  {"SymPwrGaussianResolution", CtauModel::SymPwrGaussianResolution},
  {"DoubleGaussianExp",        CtauModel::DoubleGaussianExp}

};
 
    if (parIni.count("Model_CtauRes_PP")>0) {
      model.PP.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PP"]];
      if (model.PP.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PP was not found in the initial parameters!" << endl; return false;
    }

  return true;
};


void setCtauResDataGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string dsName, double binWidth)
{
  Double_t ctauNMax; Double_t ctauNMin;
  myws.data(dsName.c_str())->getRange(*myws.var("ctauN"), ctauNMin, ctauNMax);
  ctauNMin -= 0.00000001; ctauNMax += 0.00000001;
  cout << "[INFO] Range from data: ctauNMin: " << ctauNMin << "  ctauNMax: " << ctauNMax << endl;
  int nBins = min(int( round((ctauNMax - ctauNMin)/binWidth) ), 1000);
  TH1D* hTot = (TH1D*)myws.data(dsName.c_str())->createHistogram("TMP", *myws.var("ctauN"), Binning(nBins, ctauNMin, ctauNMax));
  vector<double> rangeCtauN;
  getRange(hTot, (int)(ceil(5)), rangeCtauN);
  hTot->Delete();
  ctauNMin = rangeCtauN[0];
//___________________________________________
//  if (ctauNMin<-0.5){ ctauNMin = -0.5; }
  if (ctauNMin<-2.){ ctauNMin = -2.; }
  ctauNMin = -2.;
  if (ctauNMax>  0.05){ ctauNMax =   0.05; }
//  if (ctauNMax>  0.2){ ctauNMax =   0.2; }
//________________________________________________
  if (parIni.count("ctauNCut")>0 && parIni["ctauNCut"]!="") {
    parIni["ctauNCut"].erase(parIni["ctauNCut"].find("["), string("[").length());
    parIni["ctauNCut"].erase(parIni["ctauNCut"].find("]"), string("]").length());
    double ctauNCut = 0.0;
    try {
      ctauNCut = std::stod(parIni["ctauNCut"].c_str());
    } catch (const std::invalid_argument&) {
      cout << "[WARNING] ctauNCut is not a number, will ignore it!" << endl;
    }
    ctauNCut = abs(ctauNCut); if(ctauNCut>0.0) { ctauNMin = -1.0*abs(ctauNCut); }
  }
//  ctauNMin = -20.0;
//  ctauNMax =   -0.001;
  cout << "[INFO] Selected range: ctauNMin: " << ctauNMin << "  ctauNMax: " << ctauNMax << endl;
  myws.var("ctauN")->setRange("CtauNWindow", ctauNMin, ctauNMax);
  parIni["CtauNRange_Cut"]   = Form("(%.12f <= ctauN && ctauN <= %.12f)", ctauNMin, ctauNMax);
  cut.dMuon.ctauN.Max = ctauNMax;
  cut.dMuon.ctauN.Min = ctauNMin;
  myws.var("ctauN")->setRange("FullWindow", cut.dMuon.ctauN.Min, cut.dMuon.ctauN.Max);
  
  Double_t ctauMax; Double_t ctauMin;
  myws.data(dsName.c_str())->getRange(*myws.var("ctau"), ctauMin, ctauMax);
  ctauMin -= 0.00000001; ctauMax += 0.00000001;
  if (ctauMax >  0.0){ ctauMax =   0.0; }
  cout << "[INFO] Range from data: ctauMin: " << ctauMin << "  ctauMax: " << ctauMax << endl;
  myws.var("ctau")->setRange("CtauWindow", ctauMin, ctauMax);
  parIni["CtauRange_Cut"]   = Form("(%.12f <= ctau && ctau <= %.12f)", ctauMin, ctauMax);
  cut.dMuon.ctau.Max = ctauMax;
  cut.dMuon.ctau.Min = ctauMin;
  myws.var("ctau")->setRange("FullWindow", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);

  return;
};


void setCtauResDataFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauRes/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), TAG.c_str(), "CTAURES", TAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0));
  
  return;
};
 

void setCtauResDataCutParameters(struct KinCuts& cut)
{
  
  return;
};


bool createSignalCtauDSUsingSPLOT(RooWorkspace& ws, string dsName, map<string, string>  parIni, struct KinCuts cut, bool incJpsi, bool incBkg, bool useSPlot)
{
  
  if (dsName.find("MC")!=std::string::npos)   { return false;  }  // Only accept data
  
  string pdfMassName = "pdfMASS_Tot_PP";
  RooArgList yieldList;
  if (incJpsi)  { cout << "Include Jpsi yield" << endl; yieldList.add( *ws.var("N_Jpsi_PP")) ;  }
  if (incBkg) { yieldList.add( *ws.var("N_Bkg_PP") ); }
  RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->Clone("TMP_DATA");
  cout << "data cloned" << endl;
  RooAbsPdf* pdf = clone(*ws.pdf(pdfMassName.c_str()));
  cout << "cloned mass pdf" << endl;

  if (useSPlot) {
    if (!ws.data((dsName+"_SPLOT").c_str())) {
      RooStats::SPlot sData = RooStats::SPlot("sData","An SPlot", *data, pdf, yieldList);
      ws.import(*data, Rename((dsName+"_SPLOT").c_str()));
      if (incJpsi) {
        cout <<  "[INFO] Jpsi yield -> Mass Fit: " << ws.var("N_Jpsi_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Jpsi_PP") << std::endl;
      }
      if (incBkg) {
        cout <<  "[INFO] Bkg yield -> Mass Fit: " << ws.var("N_Bkg_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Bkg_PP") << std::endl;
      }
//      if (!isCompatibleDataset(*data, *(RooDataSet*)ws.data(dsName.c_str()))){ cout << "[ERROR] sPlot and original Datasets are inconsistent!" << endl; return false; }
    }
    if (incBkg) {
      RooDataSet* dataw_Bkg  = new RooDataSet("TMP_BKG_DATA","TMP_BKG_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("invMass"), *ws.var("ctau"), *ws.var("ctauN"), *ws.var("ctauErr"), *ws.var("N_Bkg_PP_sw")), 0, "N_Bkg_PP_sw");
//      RooDataSet* dataw_Bkg_reduced = (RooDataSet*)dataw_Bkg->reduce("(invMass<3.4)&&(ctau<=0)");
      RooDataSet* dataw_Bkg_reduced = (RooDataSet*)dataw_Bkg->reduce("(invMass<3.4)");
      ws.import(*dataw_Bkg_reduced, Rename((dsName+"_BKG").c_str()));
      cout <<  "[INFO] sPLOT_BKG_DS Events: " << ws.data((dsName+"_BKG").c_str())->sumEntries() << std::endl;
	cout << "New SPLOT Data set : " << (dsName+"_BKG").c_str() << endl;
      delete dataw_Bkg;
    }
    if (incJpsi) {
      RooDataSet* dataw_Jpsi  = new RooDataSet("TMP_JPSI_DATA","TMP_JPSI_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("invMass"), *ws.var("ctau"), *ws.var("ctauN"), *ws.var("ctauErr"), *ws.var("N_Jpsi_PP_sw")), 0, "N_Jpsi_PP_sw");
//      RooDataSet* dataw_Jpsi_reduced = (RooDataSet*)dataw_Jpsi->reduce("(invMass<3.4)&&(ctau<=0)");
      RooDataSet* dataw_Jpsi_reduced = (RooDataSet*)dataw_Jpsi->reduce("(invMass<3.4)");
//!!!!!!!!!!!!!Only signal      
    ws.import(*dataw_Jpsi_reduced, Rename((dsName+"_JPSI").c_str()));
//      ws.import(*data, Rename((dsName+"_JPSI").c_str()));
      cout <<  "[INFO] sPLOT_JPSI_DS Events: " << ws.data((dsName+"_JPSI").c_str())->sumEntries() << std::endl;
	cout << "New SPLOT Data set : " << (dsName+"_JPSI").c_str() << endl;
      delete dataw_Jpsi;
    }
  }
  else {
    if (incJpsi) {
      cout <<  "[INFO] Using mass window: " << (parIni["JpsiMassRange_Cut"]+"&&(ctau<=0)") <<  std::endl;
      cout <<  "[INFO] Number of entries in Jpsi mass window: " << ((RooDataSet*)data->reduce((parIni["JpsiMassRange_Cut"]+"&&(ctau<=0)").c_str()))->numEntries() <<  std::endl;
      RooDataSet* dataw_Jpsi_reduced = (RooDataSet*)data->reduce((parIni["JpsiMassRange_Cut"]+"&&(ctau<=0)").c_str());
      ws.var("N_Jpsi_PP")->setVal(((RooDataSet*)data->reduce((parIni["JpsiMassRange_Cut"]+"&&(ctau<=0)").c_str()))->numEntries());
      ws.import(*dataw_Jpsi_reduced, Rename((dsName+"_JPSI").c_str()));
      cout <<  "[INFO] MASSCUT_JPSI_DS Events: " << ws.data((dsName+"_JPSI").c_str())->sumEntries() << std::endl;
    }
    if (incBkg) {
      cout <<  "[INFO] Using mass window: " << (parIni["BkgMassRange_FULL_Cut"]+"&&(ctau<=0)") <<  std::endl;
      cout <<  "[INFO] Number of entries in sideband mass window: " << ((RooDataSet*)data->reduce((parIni["BkgMassRange_FULL_Cut"]+"&&(ctau<=0)").c_str()))->numEntries() <<  std::endl;
      RooDataSet* dataw_Bkg_reduced = (RooDataSet*)data->reduce((parIni["BkgMassRange_FULL_Cut"]+"&&(ctau<=0)").c_str());
      ws.var("N_Bkg_PP")->setVal(((RooDataSet*)data->reduce((parIni["BkgMassRange_FULL_Cut"]+"&&(ctau<=0)").c_str()))->numEntries());
      ws.import(*dataw_Bkg_reduced, Rename((dsName+"_BKG").c_str()));
      cout <<  "[INFO] MASSCUT_BKG_DS Events: " << ws.data((dsName+"_BKG").c_str())->sumEntries() << std::endl;
    }
  }
  delete data; delete pdf;
  
  return true;
};


#endif // #ifndef fitCharmoniaCtauResDataModel_C
