#ifndef fitCharmoniaCtauModel_test_C
#define fitCharmoniaCtauModel_test_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauModel.C"
#include "fitCharmoniaMassModel.C"
#include "fitCharmoniaCtauResModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "drawCtauPlot.C"

bool setCtauModel2( struct OniaModel& model, map<string, string>&  parIni, bool incJpsi, bool incBkg, bool incPrompt, bool incNonPrompt );
void setCtauFileName2(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool fitSideBand, bool usectauBkgTemplate);
void setCtauGlobalParameterRange2(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool fitCtauRes=false, bool is2DFit=false);
void setCtauCutParameters2(struct KinCuts& cut, bool incNonPrompt);
bool isCtauBkgPdfAlreadyFound(RooWorkspace& myws, string FileName, string pdfName, bool loadCtauBkgPdf=false);
bool createCtauDSUsingSPLOT2(RooWorkspace& ws, string dsName, map<string, string>  parIni, struct KinCuts cut, bool incJpsi, bool incBkg);

bool loadCtauErrRange_ctau2(string FileName, struct KinCuts& cut)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[ERROR] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
        cout << "File accessed" << endl;
  TFile *file = new TFile(FileName.c_str());
        cout << "File created " << FileName << endl;
  if (!file) return false;
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  cout << "Ws created" << endl;
  if (!ws) {
    cout << "[ERROR] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (ws->var("ctauErr")) {
        cout << "var accessed" << endl;
    cut.dMuon.ctauErr.Min = ws->var("ctauErr")->getMin();
    cut.dMuon.ctauErr.Max = ws->var("ctauErr")->getMax();
  } else {
    cout << Form("[ERROR] ctauErr was not found!") << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }

  delete ws;
  file->Close(); delete file;
  return true;
};


bool fitCharmoniaCtauModel_test( RooWorkspace& myws,             // Local Workspace
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
                            bool incBkg        = true,      // Includes Background model
                            bool incPrompt     = true,      // Includes Prompt ctau model
                            bool incNonPrompt  = false,     // Includes NonPrompt ctau model
                            bool useTotctauErrPdf = false,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                            bool usectauBkgTemplate = false,// If yes use a template for Bkg ctau instead of the fitted Pdf
                            // Select the fitting options
                            bool useSPlot      = true,      // If yes, then use SPlot technique, if no, use mass range
                            bool doFit         = true,      // Flag to indicate if we want to perform the fit
                            bool wantPureSMC   = false,     // Flag to indicate if we want to fit pure signal MC
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
  
  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  bool fitSideBand = true;

  string fitType = "CTAUSB";
  string label = Form("%s_PP", DSTAG.c_str());

    setMassCutParameters(cut, false, false, true);
    setCtauCutParameters2(cut, incNonPrompt);
    if (usePerEventError) {
      // check if we have already done the ctauErr fits. If yes, load their parameters
      string FileName = "";
      string pdfName = "pdfCTAUERR_Tot_PP";
      string plotLabel = "";
      bool incJpsi = true;
      bool fitSideBand = false;
      if (incJpsi)  { plotLabel = plotLabel + "_Jpsi";     }
      plotLabel = plotLabel + "_Bkg";
	cout << "plotLabel : " << plotLabel << endl;
      setCtauErrFileName(FileName, (inputFitDir["CTAUERR"]=="" ? outputDir : inputFitDir["CTAUERR"]), "DATA", plotLabel, cut, fitSideBand);
      bool foundFit = false;
      if ( loadCtauErrRange_ctau2(FileName, cut) ) { foundFit = true; }
      if (foundFit) { cout << "[INFO] The ctauErr fit was found and I'll load the ctau Error range used." << endl; }
      else { cout << "[ERROR] The ctauErr fit was not found!" << endl; return false; }
      // Importing SPLOT DS from ctauErr results
      string sPlotDSName = Form("dOS_%s_SPLOT", label.c_str());
      loadSPlotDS(myws, FileName, sPlotDSName);
    }
    setCtauErrCutParameters(cut);
  // Import the local datasets
  double numEntries = 1000000;
  string dsName = Form("dOS_%s", label.c_str());
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset(myws, inputWorkspace, cut, label, false);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }

    // Set global parameters
    cout << "setting global parameters" << endl;
    setCtauErrGlobalParameterRange(myws, parIni, cut, "", binWidth["CTAUERR"], true);
    setCtauGlobalParameterRange2(myws, parIni, cut, label, binWidth[fitType.c_str()], false);

  string dsName2Fit = dsName;
  dsName2Fit += "_BKG";
  string dsName2FitCut = dsName2Fit+"_CTAUCUT";

  string plotLabel = "";
  map<string, bool> plotLabels = {{"BkgNoPR", (incNonPrompt)},
                                  {"CtauRes", (true)}};
  for (map<string, bool>::iterator iter = plotLabels.begin(); iter!=plotLabels.end(); iter++) {
    string obj = iter->first;
    bool cond = iter->second;
    if (cond && parIni.count(Form("Model_%s_PP", obj.c_str()))>0) { 
      plotLabel = plotLabel + Form("_%s_%s", obj.c_str(), parIni[Form("Model_%s_PP", obj.c_str())].c_str()); 
    }
  }
  cout << "plotLabel2 : " << plotLabel  << endl;
  string pdfName = "pdfCTAU_Tot_PP";

  string FileName = "";
  setCtauFileName2(FileName, (inputFitDir[fitType.c_str()]=="" ? outputDir : inputFitDir[fitType.c_str()]), DSTAG, plotLabel, cut, fitSideBand, false);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir[fitType.c_str()]!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    setCtauFileName2(FileName, outputDir, DSTAG, plotLabel, cut, fitSideBand, false);
  }

  // check if we have already done this fit. If yes, do nothing and return true.
  bool found =  true; bool skipFit = !doFit;
    vector<string> pdfNames; pdfNames.push_back("pdfCTAUCOND_Bkg_PP");
    found = found && isPdfAlreadyFound(myws, FileName, pdfNames, false);
    if (found) {
        cout << "[INFO] This ctau Bkg Pdf was already made, so I'll just go to the next one." << endl;
      return true;
    }

  if (true) {
    // Check if we have already made the Background DataSet
    vector<string> dsNames = { dsName2Fit };
    bool createSPLOTDS = ( !isSPlotDSAlreadyFound(myws, FileName, dsNames, true) );
    bool compDS = loadYields(myws, FileName, dsName, pdfName);
    if (!compDS) {
	cout << "Yields not loaded" << endl;
      // Setting extra input information needed by each fitter
      string iMassFitDir = inputFitDir["MASS"];
      double ibWidth = binWidth["MASS"];
      bool loadMassFitResult = true;
      bool doMassFit = false;
      bool impDS = false;
      bool getMeanPT = false;
      const char* applyCorr = "";
      bool cutCtau = false;
      bool doConstrFit = false;
      
      if ( !fitCharmoniaMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir,
                                   DSTAG, impDS,
                                   true, true,
                                   doMassFit, cutCtau, doConstrFit, false, applyCorr, loadMassFitResult, iMassFitDir, numCores,
                                   setLogScale, incSS, ibWidth, getMeanPT
                                   )
           ) { return false; }
      // Let's set all mass parameters to constant except the yields
      if (myws.pdf("pdfMASS_Tot_PP")) {
        cout << "[INFO] Setting mass parameters to constant!" << endl;
        myws.pdf("pdfMASS_Tot_PP")->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
      } else { cout << "[ERROR] Mass PDF was not found!" << endl; return false; }
      std::vector< std::string > objs = {"Bkg", "Jpsi", "Psi2S"};
      for (auto obj : objs) { if (myws.var(Form("N_%s_PP", obj.c_str())))  setConstant( myws, Form("N_%s_PP", obj.c_str()), false); }
      // Let's set the minimum value of the yields to zero
      for (auto obj : objs) { if (myws.var(Form("N_%s_PP", obj.c_str()))) myws.var(Form("N_%s_PP", obj.c_str()))->setMin(0.0); }
    }else{

	cout << "Yields loaded, no mass fit" << endl;
    }
    if (createSPLOTDS) {
	cout << "Creating SPlot for Ctau" << endl;
      if (!myws.data(dsName2Fit.c_str())) {
        if (!createCtauDSUsingSPLOT2(myws, dsName, parIni, cut, false, true)){ cout << "[ERROR] Creating the Ctau Templates using sPLOT failed" << endl; return false; }
      }
    }
  }
  if (!myws.data(dsName2FitCut.c_str())) {
    // Cut the RooDataSet
    RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName2Fit.c_str())->reduce(parIni["CtauRange_Cut"].c_str()))->Clone(dsName2FitCut.c_str());
    myws.import(*dataToFit, Rename(dsName2FitCut.c_str()));
    double lostEvents = (myws.data(dsName2Fit.c_str())->numEntries()-myws.data(dsName2FitCut.c_str())->numEntries());
    double lostPerc = lostEvents*100./(myws.data(dsName2Fit.c_str())->numEntries());
    cout << "[INFO] After applying the cut: " << parIni["CtauRange_Cut"] << " , we lost " << lostEvents << " events ( " << lostPerc << " % ) " << endl;
  }
  
  if (true) {
	cout << "No load fit result" << endl;
    // Set models based on initial parameters
    struct OniaModel model;
	cout << "No use bkg template" << endl;
      if (!setCtauModel2(model, parIni, false, true, incPrompt, incNonPrompt)) { return false; }
	cout << "CTAU model set" << endl;

    //// LOAD CTAU ERROR PDF
    if (usePerEventError) {
	cout << "Use per event error" << endl;
      // Setting extra input information needed by each fitter
      bool loadCtauErrFitResult = true;
      bool doCtauErrFit = true;
      bool importDS = false;
      bool incJpsi = true;
      string DSTAG = "DATA_PP";
        
      if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir,
                                      DSTAG, importDS,
                                      incJpsi, incBkg,
                                      doCtauErrFit, false, loadCtauErrFitResult, inputFitDir, numCores,
                                      setLogScale, incSS, binWidth
                                      )
           ) { return false; }
	cout << "CTAUERR model fitted" << endl;
    }

    // Build the Fit Model
    if (!buildCharmoniaCtauModel(myws, model.PP, parIni, dsName2FitCut, cut, true, false, incPrompt, incNonPrompt, useTotctauErrPdf, false, binWidth["CTAUSB"], numEntries))  { return false; }
cout << "CTAU model built" << endl;

    //// LOAD CTAU RESOLUTION PDF
    if (fitSideBand) {
      // check if we have already done the resolution fits. If yes, load their results
      string FileName = "";
      string plotLabel = Form("_CtauRes_%s", parIni["Model_CtauRes_PP"].c_str());
	cout << "plotLabel3 : " << plotLabel << endl;
      string DSTAG = "DATA_PP";
      setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut);
      bool found = false;
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
        plotLabel = Form("_CtauRes_%s", parIni["Model_CtauRes_PP"].c_str());
        setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut);
      } else if (inputFitDir["CTAURES"]!="") { found = true; }
      if ( !loadPreviousFitResult(myws, FileName, DSTAG, false, false) ) {
        cout << "[ERROR] The ctau resolution fit results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau resolution fits were found, so I'll load the fit results." << endl; 
      }
      std::vector< std::string > objs = {"Bkg", "Jpsi"}; bool hasPromptPdf = false;
      for (auto obj : objs) {
        if (myws.pdf(Form("pdfCTAU_%sPR_PP", obj.c_str()))) {
          cout << "[INFO] Setting Prompt " << obj << " parameters to constant!" << endl; hasPromptPdf = true;
          myws.pdf(Form("pdfCTAU_%sPR_PP", obj.c_str()))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
        }
      }
//      if (!hasPromptPdf) { cout << "[ERROR] Prompt Ctau PDF was not found!" << endl; return false; }
      if (!setConstant(myws, "s1_CtauRes_PP", true)) { return false; }
    }

    // save the initial values of the model we've just created
    RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("invMass"), *myws.var("ctauErr")));
    myws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
    delete params;
  }

  RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")));
  
    found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
    if (found) {
      cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
      return true;
    }

  // Fit the Datasets
  if (skipFit==false) {
      bool isWeighted = myws.data(dsName2FitCut.c_str())->isWeighted();
      RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName2FitCut.c_str()), Extended(kTRUE), NumCPU(numCores), Optimize(kFALSE), SumW2Error(isWeighted), Save(), PrintLevel(1));
      fitResult->Print("v");
      myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
    
    // Draw the mass plot
    drawCtauPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, false, true, incPrompt, incNonPrompt, true, false, setLogScale, incSS, binWidth[fitType.c_str()]);
    // Save the results
    string FileName = ""; setCtauFileName2(FileName, outputDir, DSTAG, plotLabel, cut, fitSideBand, false);
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctauSB/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);
  }

  return true;
};


bool setCtauModel2( struct OniaModel& model, map<string, string>&  parIni, bool incJpsi, bool incBkg, bool incPrompt, bool incNonPrompt )
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
  {"Delta",                    CtauModel::Delta}
};

    if (parIni.count("Model_CtauRes_PP")>0) {
      model.PP.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PP"]];
      if (model.PP.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PP was not found in the initial parameters!" << endl; return false;
    }
  if (incBkg && incNonPrompt) {
    if (parIni.count("Model_BkgNoPR_PP")>0) {
      model.PP.Bkg.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_BkgNoPR_PP"]];
      if (model.PP.Bkg.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The background non-prompt ctau model: " << parIni["Model_BkgNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
//      cout << "[ERROR] Background non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (incBkg && incPrompt) {
    if (parIni.count("Model_BkgPR_PP")>0) {
      model.PP.Bkg.Ctau.Prompt = CtauModelDictionary[parIni["Model_BkgPR_PP"]];
      if (model.PP.Bkg.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The background prompt ctau model: " << parIni["Model_BkgPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_BkgPR_PP"] = "Delta";
      model.PP.Bkg.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (incJpsi && incNonPrompt) {
    if (parIni.count("Model_JpsiNoPR_PP")>0) {
      model.PP.Jpsi.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_JpsiNoPR_PP"]];
      if (model.PP.Jpsi.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi non-prompt ctau model: " << parIni["Model_JpsiNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (incJpsi && incPrompt) {
    if (parIni.count("Model_JpsiPR_PP")>0) {
      model.PP.Jpsi.Ctau.Prompt = CtauModelDictionary[parIni["Model_JpsiPR_PP"]];
      if (model.PP.Jpsi.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi prompt ctau model: " << parIni["Model_JpsiPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_JpsiPR_PP"] = "Delta";
      model.PP.Jpsi.Ctau.Prompt=CtauModel::Delta;
    }
  }

  return true;
};


void setCtauGlobalParameterRange2(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool optimizeRange, bool is2DFit)
{
  Double_t ctauMax; Double_t ctauMin;
  myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctau"), ctauMin, ctauMax);
  ctauMin -= 0.00001;  ctauMax += 0.00001;
cout << "ctaumin = " << ctauMin << ", ctaumax = " << ctauMax << endl;
  int nBins = min(int( round((ctauMax - ctauMin)/binWidth) ), 1000);
  if (optimizeRange) {
    TH1D* hTot = (TH1D*)((RooDataSet*)myws.data(Form("dOS_%s", label.c_str()))->Clone("dTMP"))->createHistogram("hTMP", *myws.var("ctau"), Binning(nBins, ctauMin, ctauMax));
    vector<double> rangeCtau; 
    getRange(hTot, 4, rangeCtau);
    hTot->Delete();
    ctauMin = rangeCtau[0];
    if (ctauMin<cut.dMuon.ctau.Min) { ctauMin = cut.dMuon.ctau.Min; }
    if (ctauMax>cut.dMuon.ctau.Max) { ctauMax = cut.dMuon.ctau.Max; }
    if (ctauMin < -4.0) { ctauMin = -4.0; }
    cout << " optimized range : ctaumin = " << ctauMin << ", ctaumax = " << ctauMax << endl;
  }
  else if (!is2DFit && parIni.count("ctauCut_PP")>0 && parIni["ctauCut_PP"]!="") {
    parIni["ctauCut_PP"].erase(parIni["ctauCut_PP"].find("["), string("[").length());
    parIni["ctauCut_PP"].erase(parIni["ctauCut_PP"].find("]"), string("]").length());
    TString sctauCut(parIni["ctauCut_PP"].c_str());
    TObjArray* actauCut = sctauCut.Tokenize(",");    
    ctauMin = (static_cast<TObjString*>(actauCut->At(0)))->GetString().Atof();
    ctauMax = (static_cast<TObjString*>(actauCut->At(1)))->GetString().Atof();
    delete actauCut;
  }

    ctauMin =  -1.5;
    ctauMax =  1.5;
  cout << "[INFO] Range from data: ctauMin: " << ctauMin << "  ctauMax: " << ctauMax << endl;
  myws.var("ctau")->setRange("CtauWindow", ctauMin, ctauMax);
  parIni["CtauRange_Cut"]   = Form("(%.12f <= ctau && ctau < %.12f)", ctauMin, ctauMax);
  cut.dMuon.ctau.Max = ctauMax;
  cut.dMuon.ctau.Min = ctauMin;
  myws.var("ctau")->setRange("CtauFullWindow", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("FullWindow", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandTOP_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max); 
  myws.var("ctau")->setRange("SideBandMID_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandBOT_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max); 
  myws.var("ctau")->setRange("SideBandMID_JPSI", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandBOT_JPSI", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);

  return;
};


void setCtauFileName2(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool fitSideBand, bool usectauBkgTemplate)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  if (!usectauBkgTemplate) {FileName = Form("%sctau%s/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), (fitSideBand?"SB":""), TAG.c_str(), "CTAU", TAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0));}
  else {FileName = Form("%sctau%sTemp/%s/result/FIT_%s_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), (fitSideBand?"SB":""), TAG.c_str(), "CTAU", TAG.c_str(), "PP",
                        "Bkg", (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0));}

  return;
};
 

void setCtauCutParameters2(struct KinCuts& cut, bool incNonPrompt)
{
  // Define the ctau range
  if (cut.dMuon.ctau.Min==-1000. && cut.dMuon.ctau.Max==1000.) { 
    // Default ctau values, means that the user did not specify a ctau range
    if (incNonPrompt) {
      cut.dMuon.ctau.Min = -30.0;
      cut.dMuon.ctau.Max = 100.0;
//      cut.dMuon.ctau.Min = -3.0;
//      cut.dMuon.ctau.Max = 3.0;
    } else {
      cut.dMuon.ctau.Min = -30.0;
      cut.dMuon.ctau.Max = 100.0;
//      cut.dMuon.ctau.Min = -3.0;
//      cut.dMuon.ctau.Max = 3.0;
    }
  }
        
  cout << "[INFO] Setting ctau range to min: " << cut.dMuon.ctau.Min << " and max " << cut.dMuon.ctau.Max << endl;
  
  return;
};


bool createCtauDSUsingSPLOT2(RooWorkspace& ws, string dsName, map<string, string>  parIni, struct KinCuts cut, bool incJpsi, bool incBkg)
{  
    string pdfMassName = "pdfMASS_Tot_PP";
    RooArgList yieldList;
    yieldList.add( *ws.var("N_Jpsi_PP") );
    yieldList.add( *ws.var("N_Bkg_PP") ); // Always add background
    RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->Clone("TMP_DATA");
    RooAbsPdf* pdf = clone(*ws.pdf(pdfMassName.c_str()));
  if (!ws.data((dsName+"_SPLOT").c_str())) {
    RooStats::SPlot sData = RooStats::SPlot("sData","An SPlot", *data, pdf, yieldList);
    ws.import(*data, Rename((dsName+"_SPLOT").c_str()));
    cout <<  "[INFO] Bkg yield -> Mass Fit: " << ws.var("N_Bkg_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Bkg_PP") << std::endl;
    cout <<  "[INFO] Jpsi yield -> Mass Fit: " << ws.var("N_Jpsi_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Jpsi_PP") << std::endl;
//    if (!isCompatibleDataset(*data, *(RooDataSet*)ws.data(dsName.c_str()))){ cout << "[ERROR] sPlot and original Datasets are inconsistent!" << endl; return false; }
    delete data; delete pdf;
  }
  
    RooDataSet* dataw_Bkg  = new RooDataSet("TMP_BKG_DATA","TMP_BKG_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("invMass"), *ws.var("ctau"), *ws.var("ctauN"), *ws.var("ctauErr"), *ws.var("N_Bkg_PP_sw")), 0, "N_Bkg_PP_sw");
    ws.import(*dataw_Bkg, Rename((dsName+"_BKG").c_str()));
    cout <<  "[INFO] sPLOT_BKG_DS Events: " << ws.data((dsName+"_BKG").c_str())->sumEntries() << std::endl;
    delete dataw_Bkg;
  
  return true;
};


#endif // #ifndef fitCharmoniaCtauModel_test_C
