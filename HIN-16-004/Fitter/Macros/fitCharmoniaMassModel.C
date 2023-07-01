#ifndef fitCharmoniaMassModel_C
#define fitCharmoniaMassModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaMassModel.C"
#include "drawMassPlot.C"

void setCtauCuts(struct KinCuts& cut);
bool setMassModel( struct OniaModel& model, map<string, string>&  parIni, bool incJpsi, bool incBkg );
void setMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool cutSideBand=false) ;
void setMassGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, bool incJpsi, bool incBkg, bool wantPureSMC=false);
void setMassCutParameters(struct KinCuts& cut, bool incJpsi, bool isMC=false, bool useForCtauFits=false);


bool fitCharmoniaMassModel( RooWorkspace& myws,            // Local Workspace
                            const RooWorkspace& inputWorkspace,  // Workspace with all the input RooDatasets
                            struct KinCuts& cut,           // Variable containing all kinematic cuts
                            map<string, string>&  parIni,  // Variable containing all initial parameters
                            struct InputOpt& opt,          // Variable with run information (kept for legacy purpose)
                            string outputDir,              // Path to output directory
                            // Select the type of datasets to fit
                            string DSTAG,                  // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                            bool importDS    = true,       // Select if the dataset is imported in the local workspace
                            // Select the type of object to fit
                            bool incJpsi     = true,       // Includes Jpsi model
                            bool incBkg      = true,       // Includes Background model
                            // Select the fitting options
                            bool doFit       = true,       // Flag to indicate if we want to perform the fit
                            bool cutCtau     = false,      // Apply prompt ctau cuts
                            bool doConstrFit   = false,    // Do constrained fit
                            bool wantPureSMC = false,      // Flag to indicate if we want to fit pure signal MC
                            const char* applyCorr ="",     // Flag to indicate if we want corrected dataset and which correction
                            uint loadFitResult = false,    // Load previous fit results
                            string inputFitDir = "",       // Location of the fit results
                            int  numCores    = 2,          // Number of cores used for fitting
                            // Select the drawing options
                            bool setLogScale = true,       // Draw plot with log scale
                            bool incSS       = false,      // Include Same Sign data
                            double  binWidth = 0.05,       // Bin width used for plotting
                            bool getMeanPT   = false       // Compute the mean PT (NEED TO FIX)
                            )  
{

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    isMC = true;
  }
  wantPureSMC = (isMC && wantPureSMC);
  bool cutSideBand = (incBkg &&  !incJpsi);
  bool applyWeight_Corr = ( strcmp(applyCorr,"") );
  
  // Define the mass range
  setMassCutParameters(cut, incJpsi, isMC);
  parIni["invMassNorm"] = Form("RooFormulaVar::%s('( -1.0 + 2.0*( @0 - @1 )/( @2 - @1) )', {%s, mMin[%.6f], mMax[%.6f]})", "invMassNorm", "invMass", cut.dMuon.M.Min, cut.dMuon.M.Max );
  cout << "[INFO] RooFormula defined" << endl;
  // Apply the ctau cuts to reject non-prompt charmonia
//  if (cutCtau) { setCtauCuts(cut); }
  
  string COLL = "PP";
  string plotLabelPP;

    // Set models based on initial parameters
    cout << "[INFO] fit pp" << endl;
    struct OniaModel model;
    if (!setMassModel(model, parIni, incJpsi, incBkg)) { 
	cout << "Setting mass model failed" << endl;
	return false; }

    cout << "[INFO] Mass model set" << endl;
    // Import the local datasets
    double numEntries = 1000000;
    string label = ((DSTAG.find("PP")!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), "PP"));
    if (wantPureSMC) label += "_NoBkg";
    if (applyWeight_Corr) label += Form("_%s",applyCorr);
    string dsName = Form("dOS_%s", label.c_str());

    cout << "[INFO] Label : " << label << endl;
    cout << "[INFO] DataSet name : " << dsName.c_str() << endl;
    if (importDS) {
      if ( !(myws.data(dsName.c_str())) ) {
        int importID = importDataset(myws, inputWorkspace, cut, label, cutSideBand);
        if (importID<0) { return false; }
        else if (importID==0) { doFit = false; }
      }
      numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
    }
    else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }
    if (myws.data(dsName.c_str())) numEntries = myws.data(dsName.c_str())->sumEntries();

    // Set global parameters
    setMassGlobalParameterRange(myws, parIni, cut, incJpsi, incBkg, wantPureSMC);

    // Build the Fit Model
    if (!buildCharmoniaMassModel(myws, model.PP, parIni, doConstrFit, incBkg, incJpsi, numEntries))  {

	cout << "Building mass model failed" << endl;
	 return false; }

    // Define plot names
    if (incJpsi)  { plotLabelPP += Form("_Jpsi_%s", parIni["Model_Jpsi_PP"].c_str());   } 
    if (incBkg)   { plotLabelPP += Form("_Bkg_%s", parIni["Model_Bkg_PP"].c_str());     }
    if (wantPureSMC) plotLabelPP +="_NoBkg";
    if (applyWeight_Corr) plotLabelPP +=Form("_%s",applyCorr);


    // Define pdf and plot names
    string pdfName = Form("pdfMASS_Tot_%s", COLL.c_str());
    string plotLabel = plotLabelPP;

    if (wantPureSMC) label += "_NoBkg";
    if (applyWeight_Corr) label += Form("_%s",applyCorr);
      
    // check if we have already done this fit. If yes, do nothing and return true.
    string FileName = "";
    setMassFileName(FileName, (inputFitDir=="" ? outputDir : inputFitDir), DSTAG, plotLabel, cut, cutSideBand);
    if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir!="") {
      cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
      if (loadFitResult) return false;
      setMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, cutSideBand);
    }
    bool found =  true; bool skipFit = !doFit;
    RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("invMass")));
	cout << "parameters got" << endl;
    found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
	cout << "fit found" << endl;
    if (loadFitResult) {
      if ( loadPreviousFitResult(myws, FileName, DSTAG, (!isMC && !cutSideBand && loadFitResult==1), loadFitResult==1) ) { skipFit = true; } else { skipFit = false; } 
      if (skipFit) { cout << "[INFO] This mass fit was already done, so I'll load the fit results." << endl; }
      myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()), *newpars, kTRUE);
    } else if (found) {
      cout << "[INFO] This mass fit was already done, so I'll just go to the next one." << endl;
      return true;
    }

    // Fit the Datasets
    if (skipFit==false) {
      bool isWeighted = myws.data(dsName.c_str())->isWeighted();
      RooFitResult* fitResult(0x0);
      if (doConstrFit)
      {
        cout << "[INFO] Performing constrained fit" << endl;
          cout << "[INFO] Constrained variables: alpha, n, ratio of sigmas" << endl;
          fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(isWeighted), Range(cutSideBand ? parIni["BkgMassRange_FULL_Label"].c_str() : "MassWindow"), ExternalConstraints(RooArgSet(*(myws.pdf("sigmaAlphaConstr")),*(myws.pdf("sigmaNConstr")))), NumCPU(numCores), Save());
      }
      else
      {
	cout << "No constr fit" << endl;
       fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(isWeighted), Range(cutSideBand ? parIni["BkgMassRange_FULL_Label"].c_str() : "MassWindow"), NumCPU(numCores), Save());
      }
      fitResult->Print("v"); 
      myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str())); 
      // Create the output files
      int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
      drawMassPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, incJpsi, incBkg, cutCtau, wantPureSMC, setLogScale, incSS, nBins, getMeanPT);
      // Save the results
      string FileName = ""; setMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, cutSideBand);
      myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()), *newpars, kTRUE);
      saveWorkSpace(myws, Form("%smass%s/%s/result", outputDir.c_str(), (cutSideBand?"SB":""), DSTAG.c_str()), FileName);
    }

  return true;
};


void setCtauCuts(struct KinCuts& cut) 
{
  if (cut.dMuon.AbsRap.Max<=1.6) {
    cut.dMuon.ctauCut = "( ctau < (0.012 + (0.23/pt)) )";
  }
  if (cut.dMuon.AbsRap.Min>=1.6) {
    cut.dMuon.ctauCut = "( ctau < (0.014 + (0.28/pt)) )";
  }

  return;
};


bool setMassModel( struct OniaModel& model, map<string, string>&  parIni, bool incJpsi, bool incBkg)
{

map< string , MassModel > MassModelDictionary = {
  {"InvalidModel",            MassModel::InvalidModel},
  {"SingleGaussian",          MassModel::SingleGaussian},
  {"DoubleGaussian",          MassModel::DoubleGaussian},
  {"SingleCrystalBall",       MassModel::SingleCrystalBall},
  {"ExtendedCrystalBall",     MassModel::ExtendedCrystalBall},
  {"DoubleCrystalBall",       MassModel::DoubleCrystalBall},
  {"GaussianAndCrystalBall",  MassModel::GaussianAndCrystalBall},
  {"Uniform",                 MassModel::Uniform},
  {"Chebychev1",              MassModel::Chebychev1},
  {"Chebychev2",              MassModel::Chebychev2},
  {"Chebychev3",              MassModel::Chebychev3},
  {"Chebychev4",              MassModel::Chebychev4},
  {"Chebychev5",              MassModel::Chebychev5},
  {"Chebychev6",              MassModel::Chebychev6},
  {"ExpChebychev1",           MassModel::ExpChebychev1},
  {"ExpChebychev2",           MassModel::ExpChebychev2},
  {"ExpChebychev3",           MassModel::ExpChebychev3},
  {"ExpChebychev4",           MassModel::ExpChebychev4},
  {"ExpChebychev5",           MassModel::ExpChebychev5},
  {"ExpChebychev6",           MassModel::ExpChebychev6},
  {"Exponential",             MassModel::Exponential},
  {"NA60",             MassModel::NA60},
  {"VWG",             MassModel::VWG}
};


  if (incBkg) {
    if (parIni.count("Model_Bkg_PP")>0) {
      model.PP.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PP"]];
      if (model.PP.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (incJpsi) {
    if (parIni.count("Model_Jpsi_PP")>0) {
      model.PP.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PP"]];
	cout << "Mass model set to " << parIni["Model_Jpsi_PP"] << " " << endl;
      if (model.PP.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};


void setMassGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, bool incJpsi, bool incBkg,  bool wantPureSMC)
{
  if (wantPureSMC)
    {
      if (incJpsi)
        {
          if (cut.dMuon.AbsRap.Min >= 1.6) {
            myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.32);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.32);
          }
          else {
            myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.26);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.26);
          }
        }
    }
  else {
    myws.var("invMass")->setRange("InclusiveMassRegion", 2.6, cut.dMuon.M.Max);
    myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
    myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
    parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, cut.dMuon.M.Max);
    if (incJpsi) {
      myws.var("invMass")->setRange("JpsiWindow", 2.9, 3.3);
      parIni["JpsiMassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 2.9, 3.3);
    }
    if (incBkg) {
      myws.var("invMass")->setRange("SideBandMID_FULL",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.5)?3.5:cut.dMuon.M.Max));
      myws.var("invMass")->setRange("SideBandMID_JPSI",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.5)?3.5:cut.dMuon.M.Max));
      parIni["BkgMassRange_FULL_Label"]  = "SideBandMID_FULL";
      parIni["BkgMassRange_JPSI_Label"]  = "SideBandMID_JPSI";
      if (cut.dMuon.M.Min < 2.8) {
        myws.var("invMass")->setRange("SideBandBOT_FULL", cut.dMuon.M.Min, 2.8);
        myws.var("invMass")->setRange("SideBandBOT_JPSI", ((cut.dMuon.M.Min<2.0)?2.0:cut.dMuon.M.Min), 2.8);
        parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandBOT_FULL";
        parIni["BkgMassRange_JPSI_Label"] = parIni["BkgMassRange_JPSI_Label"] + "," + "SideBandBOT_JPSI";
      }
      if (cut.dMuon.M.Max > 3.9) {
        myws.var("invMass")->setRange("SideBandTOP_FULL", 3.9, cut.dMuon.M.Max);
        myws.var("invMass")->setRange("SideBandTOP_PSI2S", 3.9, ((cut.dMuon.M.Max>4.2)?4.2:cut.dMuon.M.Max));
        parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandTOP_FULL";
      }
      parIni["BkgMassRange_FULL_Cut"]  = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
      parIni["BkgMassRange_FULL_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "((2.0 < invMass && invMass < 2.8) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
      parIni["BkgMassRange_JPSI_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "(2.0 < invMass && invMass < 3.5)";
      parIni["BkgMassRange_FULL_Cut"]  = "("+parIni["BkgMassRange_FULL_Cut"]+")";
      parIni["BkgMassRange_JPSI_Cut"]  = "("+parIni["BkgMassRange_JPSI_Cut"]+")";
    }
  }

  return;
};

void setMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool cutSideBand) 
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%smass%s/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), (cutSideBand?"SB":""), TAG.c_str(), "MASS", TAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0));
};


void setMassCutParameters(struct KinCuts& cut, bool incJpsi, bool isMC, bool useForCtauFits)
{
  // Define the mass range
  if (cut.dMuon.M.Max==5 && cut.dMuon.M.Min==2) { 
    // Default mass values, means that the user did not specify a mass range
    if ( incJpsi) {
      if (isMC && !useForCtauFits){
        cut.dMuon.M.Min = 2.6;
        cut.dMuon.M.Max = 3.5;
      } else {
        cut.dMuon.M.Min = 2.6;
        cut.dMuon.M.Max = 3.5;
      }
    }
    else {
      cut.dMuon.M.Min = 2.6;
      cut.dMuon.M.Max = 3.5;
    }
  }
  cout << "[INFO] Setting mass range to min: " << cut.dMuon.M.Min << " and max " << cut.dMuon.M.Max << endl;

  return;
};


#endif // #ifndef fitCharmoniaMassModel_C
