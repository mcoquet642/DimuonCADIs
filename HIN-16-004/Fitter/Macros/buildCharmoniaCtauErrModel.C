#ifndef buildCharmoniaCtauErrModel_C
#define buildCharmoniaCtauErrModel_C

#include "Utilities/initClasses.h"
#include "RooStats/SPlot.h"

using namespace RooStats;

void setCtauErrDefaultParameters(map<string, string> &parIni, double numEntries);
bool histToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, vector<double> rangeErr);
bool makeCtauErrPdf(RooWorkspace& ws, vector<TH1D*>& outputHist, vector<string> objectColl, vector<string> rangeColl, vector<string> cutColl, string dsName, struct KinCuts cut, double binWidth);
bool createCtauErrTemplateUsingMatrix(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut,  map<string,string> parIni, bool incJpsi, double binWidth);
bool createCtauErrTemplateUsingSPLOT(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, double binWidth);
bool addCtauErrModel(RooWorkspace& ws, string object, string pdfType, map<string,string> parIni);
TH1* rebinhist(TH1 *hist, double xmin=1e99, double xmax=-1e99);
void MatrixInverse2x2(TMatrixD& Inverse);


bool buildCharmoniaCtauErrModel(RooWorkspace& ws, map<string, string>  parIni, 
                                struct KinCuts cut,          // Variable containing all kinematic cuts
                                string dsName,               // Name of current input dataset
                                bool incJpsi,                // Include Jpsi model
                                double  binWidth,            // Bin width
                                double  numEntries = 300000. // Number of entries in the dataset
                                )
{


  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) { isMC = true; }

  // If the initial parameters are empty, set defaul parameter values
  setCtauErrDefaultParameters(parIni, numEntries);

  // C r e a t e   m o d e l 

  string pdfType = "pdfCTAUERR";
  if(!createCtauErrTemplateUsingSPLOT(ws, dsName, pdfType, cut, incJpsi, binWidth)) { cout << "[ERROR] Creating the Ctau Error Templates using sPLOT failed" << endl; return false; }
  //if(!createCtauErrTemplateUsingMatrix(ws, dsName, pdfType, cut, parIni, incJpsi, binWidth)) { cout << "[ERROR] Creating the Ctau Error Templates failed" << endl; return false; }
  cout << "Template created w/ sPlot" << endl;
  if (incJpsi)  {  if(!addCtauErrModel(ws, "Jpsi", pdfType, parIni))  { return false; }  }
  if (!isMC)    {  if(!addCtauErrModel(ws, "Bkg", pdfType, parIni))   { return false; }  }

cout << "Models added " << endl;
  // Total PDF
  
  RooArgList pdfList;
  if (incJpsi)  { pdfList.add( *ws.pdf(Form("%sTot_Jpsi_PP", pdfType.c_str())) );  }
  if (!isMC)    { pdfList.add( *ws.pdf(Form("%sTot_Bkg_PP", pdfType.c_str())) );   }
  if (!incJpsi && isMC) {
    cout << "[ERROR] User did not include any model for MC, please fix your input settings!" << endl; return false;
  }
cout << "Pdf added " << endl;
  string pdfName = Form("%s_Tot_PP", pdfType.c_str());
  RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList);
  ws.import(*themodel);
cout << "Model imported " << endl;
  ws.pdf(pdfName.c_str())->setNormRange("CtauErrWindow");
  
  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("ctauErr")));
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);

  return true;
};


bool addCtauErrModel(RooWorkspace& ws, string object, string pdfType, map<string,string> parIni)
{
  if (ws.pdf(Form("%sTot_%s_PP", pdfType.c_str(), object.c_str()))) { 
    cout << Form("[ERROR] The %s Extended Ctau Error Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  if (!ws.var(Form("N_%s_PP", object.c_str()))){ ws.factory( parIni[Form("N_%s_PP", object.c_str())].c_str() ); } 
  // create the Extended PDF
  ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_%s_PP", pdfType.c_str(), object.c_str()),
                  Form("%s_%s_PP", pdfType.c_str(), object.c_str()),
                  Form("N_%s_PP", object.c_str())
                  ));

  return true;
};


bool createCtauErrTemplateUsingSPLOT(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");

  if (dsName.find("MC")!=std::string::npos)   { cout << "[ERROR] sPlot can only be run on data!" << endl; return false;  }  // Only accept data

  string pdfMassName = "pdfMASS_Tot_PP";
  RooArgList yieldList;
  if (incJpsi)  { yieldList.add( *ws.var("N_Jpsi_PP") );  }
  yieldList.add( *ws.var("N_Bkg_PP") ); // Always add background
  RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->Clone("TMP_DATA");
  RooAbsPdf* pdf = clone(*ws.pdf(pdfMassName.c_str()));

  RooStats::SPlot sData = RooStats::SPlot("sData","An SPlot", *data, pdf, yieldList);
  ws.import(*data, Rename((dsName+"_SPLOT").c_str()));
  if (incJpsi) {
    cout <<  "[INFO] Jpsi yield -> Mass Fit: " << ws.var("N_Jpsi_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Jpsi_PP") << std::endl;
  }
  cout <<  "[INFO] Bkg yield -> Mass Fit: " << ws.var("N_Bkg_PP")->getVal() << " , sWeights: " << sData.GetYieldFromSWeight("N_Bkg_PP") << std::endl;
  
  // create weighted data sets
  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  vector<double> rangeErr; rangeErr.push_back(cut.dMuon.ctauErr.Min); rangeErr.push_back(cut.dMuon.ctauErr.Max);
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);

  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_Tot_PP", hType.c_str()), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  if ( !histToPdf(ws, hTot, Form("%sTot_Tot_PP", pdfType.c_str()), rangeErr)) { return false; }
  hTot->Delete();

  RooDataSet* dataw_Bkg  = new RooDataSet("TMP_BKG_DATA","TMP_BKG_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("ctauErr"), *ws.var("N_Bkg_PP_sw")), 0, "N_Bkg_PP_sw");
  TH1D* hBkg = (TH1D*)dataw_Bkg->createHistogram(Form("%s_Bkg_PP_sWEIGHT", hType.c_str()), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  if ( !histToPdf(ws, hBkg, Form("%s_Bkg_PP", pdfType.c_str()), rangeErr)) { return false; }
  delete hBkg; delete dataw_Bkg;
  
  if (incJpsi) {
    RooDataSet* dataw_Jpsi  = new RooDataSet("TMP_JPSI_DATA","TMP_JPSI_DATA", (RooDataSet*)ws.data((dsName+"_SPLOT").c_str()), RooArgSet(*ws.var("ctauErr"), *ws.var("N_Jpsi_PP_sw")), 0, "N_Jpsi_PP_sw");
    TH1D* hJpsi = (TH1D*)dataw_Jpsi->createHistogram(Form("%s_Jpsi_PP_sWEIGHT", hType.c_str()), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
    if ( !histToPdf(ws, hJpsi, Form("%s_Jpsi_PP", pdfType.c_str()), rangeErr)) { return false; }
    delete hJpsi; delete dataw_Jpsi;
  }

  return true;
};
  

bool createCtauErrTemplateUsingMatrix(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut,  map<string,string> parIni, bool incJpsi, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");

  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) { isMC = true; }

  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);
  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_Tot_PP", hType.c_str()), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  vector<double> rangeErr; rangeErr.push_back(cut.dMuon.ctauErr.Min); rangeErr.push_back(cut.dMuon.ctauErr.Max);

  // Set the range of the ctau Error
  if (isMC) {
    TH1D* hSig = (TH1D*)hTot->Clone("SIGTMP");
    delete hSig;
  }
  else {
    TH1D* hBkg = (TH1D*)hTot->Clone("BKGTMP");
    if (incJpsi) {
      vector<string> rangeColl; vector<string> cutColl; vector<string> objectColl;
      rangeColl.push_back("JpsiWindow"); cutColl.push_back(parIni["JpsiMassRange_Cut"]);  objectColl.push_back("Jpsi");
      rangeColl.push_back(parIni["BkgMassRange_JPSI_Label"]); cutColl.push_back(parIni["BkgMassRange_JPSI_Cut"]);   objectColl.push_back("Bkg");
      vector<TH1D*> outputHist_JPSI;
      if ( !makeCtauErrPdf(ws, outputHist_JPSI, objectColl, rangeColl, cutColl, dsName, cut, binWidth) ) { return false; }
      hBkg->Add(outputHist_JPSI.at(0), -1.0);
      if ( !histToPdf(ws, outputHist_JPSI.at(0), Form("%s_Jpsi_PP", pdfType.c_str()), rangeErr)) { return false; }
      if ( !histToPdf(ws, hBkg/*outputHist_JPSI.at(1)*/, Form("%s_Bkg_PP", pdfType.c_str()), rangeErr)) { return false; }
      delete outputHist_JPSI.at(0); delete outputHist_JPSI.at(1);
    }
    delete hBkg;
  }
  if ( !histToPdf(ws, hTot, Form("%sTot_Tot_PP", pdfType.c_str()), rangeErr)) { return false; }
  hTot->Delete();

  return true;
};


bool makeCtauErrPdf(RooWorkspace& ws, vector<TH1D*>& outputHist, vector<string> objectColl, vector<string> rangeColl, vector<string> cutColl, string dsName, struct KinCuts cut, double binWidth)
{
  if (objectColl.size()!=rangeColl.size() || objectColl.size()<1) {
    cout << "[ERROR] Object and range collections are wrong, CtauErrPdf maker failed!" << endl; return false;
  }

  double ctauErrMax = cut.dMuon.ctauErr.Max, ctauErrMin = cut.dMuon.ctauErr.Min;
  int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);


  // Extract the PDFs and the Number of events for each physics object
  vector<Double_t>    N;
  vector<RooAbsPdf*>  pdfMass;
  for (unsigned int i=0; i<objectColl.size(); i++) {
    if (!ws.pdf(Form("pdfMASS_%s_PP", objectColl.at(i).c_str()))) {
      cout << Form("[ERROR] Mass PDF for %s is missing!", objectColl.at(i).c_str()) << endl; return false;
    } else {
      pdfMass.push_back(ws.pdf(Form("pdfMASS_%s_PP", objectColl.at(i).c_str())));
    }
    if (ws.var(Form("N_%s_PP", objectColl.at(i).c_str()))) { 
      N.push_back(ws.var(Form("N_%s_PP", objectColl.at(i).c_str()))->getValV());
    } else if (ws.function(Form("N_%s_PP", objectColl.at(i).c_str()))) {
      N.push_back(ws.function(Form("N_%s_PP", objectColl.at(i).c_str()))->getValV());
    } else {
     cout << Form("[ERROR] N parameter for %s is missing!", objectColl.at(i).c_str()) << endl; return false;
    } 
  }

  // Extract the ctau error histograms for each mass range 
  vector<TH1*> histCtauErr;
  for (unsigned int j=0; j<rangeColl.size(); j++) {
    RooDataSet* data = (RooDataSet*)ws.data(dsName.c_str())->reduce(*ws.var("ctauErr"), cutColl.at(j).c_str());
    histCtauErr.push_back( (TH1*)data->createHistogram(Form("hCtauErr_Range_%s_PP", rangeColl.at(j).c_str()), *ws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax)) );
    histCtauErr.at(j)->Scale(1/histCtauErr.at(j)->GetEntries());
    delete data;
  }

  // Get the mass yield coefficients by integrating the mass PDFs
  TMatrixD CoefficientMatrix = TMatrixD(objectColl.size(), rangeColl.size());
  for (unsigned int j=0; j<rangeColl.size(); j++) {
    Double_t Norm = 0;
    for (unsigned int i=0; i<objectColl.size(); i++) {
      CoefficientMatrix(j, i) = pdfMass.at(i)->createIntegral(*ws.var("invMass"), NormSet(*ws.var("invMass")), Range(rangeColl.at(j).c_str()))->getValV();
      CoefficientMatrix(j, i) = N.at(i) * CoefficientMatrix(j, i);
      Norm = Norm + CoefficientMatrix(j, i);
    }
    if (Norm>0.0) {
      double testNorm = 0.0;
      for (unsigned int i=0; i<objectColl.size(); i++) {
        CoefficientMatrix(j, i) = CoefficientMatrix(j, i) / Norm;
        testNorm += CoefficientMatrix(j, i);
      }
      if (abs(testNorm-1.0)>0.000000000001) { cout << "[ERROR] Sum of alphas is: " << Form("%.10f",testNorm) << " !" << endl; return false; }
    } else {
      cout << "[ERROR] Normalization factor for Range: " << rangeColl.at(j) << " is invalid: " << Form("%.4f",Norm) << endl; return false;
    }
  }

  // Get the inverse of the coefficients
  CoefficientMatrix.Print();
  Double_t Det = CoefficientMatrix.Determinant();
  TMatrixD InvMatrix = CoefficientMatrix;
  if (Det>0) {
    InvMatrix.Invert();
    TMatrixD TestMatrix = InvMatrix * CoefficientMatrix;
    double testDet = TestMatrix.Determinant();
    if (fabs(testDet-1.0)>0.00000000001) { cout << Form("[ERROR] Determinant of invMatrix * CoefficienctMatrix is: %.10f", testDet) << endl; return false; }
  } else {
    cout << "[ERROR] Determinant of Coefficient Matrix is invalid: " << Form("%.4f",Det) << endl; return false;
  }

  // Get the ctau error histrograms for each physics object
  for (unsigned int j=0; j<objectColl.size(); j++) {
    outputHist.push_back(new TH1D(Form("h_%s_PP", objectColl.at(j).c_str()), "", nBins, ctauErrMin, ctauErrMax));
    for (unsigned int i=0; i<rangeColl.size(); i++) {
      outputHist.at(j)->Add(histCtauErr.at(i), InvMatrix(j, i));
    }
    outputHist.at(j)->Scale(N.at(j)*pdfMass.at(j)->createIntegral(*ws.var("invMass"), NormSet(*ws.var("invMass")), Range("FullWindow"))->getValV());
  }

  for (unsigned int i=0; i<rangeColl.size(); i++) {
    histCtauErr.at(i)->Delete();
  }
   
  return true;
  
};


bool histToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, vector<double> rangeErr)
{
  if (ws.pdf(pdfName.c_str())) { 
    cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << endl;
    return true; 
  }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0) { hist->SetBinContent(i, 0.0000000001); } }
  // 2) Reduce the range of histogram and rebin it
  TH1* hClean = rebinhist(hist, rangeErr[0], rangeErr[1]);

  cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << endl;

  string dataName = pdfName;
  dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
  RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *ws.var("ctauErr"), hClean);
  if (dataHist==NULL) { cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << endl; return false; } 
  if (dataHist->sumEntries()==0) { cout << "[ERROR] DataHist used to create " << pdfName << " is empty!" << endl; return false; } 
  if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { cout << "[ERROR] DataHist used to create " << pdfName << "  " << " is invalid!  " << endl; return false; } 
  ws.import(*dataHist);

  RooHistPdf* pdf = new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *((RooDataHist*)ws.data(dataName.c_str())));
  //RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *((RooDataSet*)ws.data(dataName.c_str())),RooKeysPdf::NoMirror, isPbPb?0.4:0.4);
  if (pdf==NULL) { cout << "[ERROR] RooKeysPDF " << pdfName << " is NULL!" << endl; return false; } 
  ws.import(*pdf);

  delete pdf;
  delete dataHist;
 
  return true;
};


TH1* rebinhist(TH1 *hist, double xmin, double xmax)
{
  TH1 *hcopy = (TH1*) hist->Clone("hcopy");

  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();

  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.1) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.00000001) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          hcopy->SetBinContent(i,hcopy->GetBinContent(i)/nrebin);
          break;
        }
        nrebin++;
      }
    }
  }

  if (xmin < newbins[1]) newbins[0] = xmin;
  if (xmax > newbins[newbins.size()-2]) newbins[newbins.size()-1] = xmax;

  TH1 *ans = hcopy->Rebin(newbins.size()-1,"hnew",newbins.data());

  delete hcopy;
  return ans;
};


void setCtauErrDefaultParameters(map<string, string> &parIni, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count("N_Jpsi_PP")==0 || parIni["N_Jpsi_PP"]=="") { 
    parIni["N_Jpsi_PP"]  = Form("%s[%.12f,%.12f,%.12f]", "N_Jpsi_PP", numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count("N_Bkg_PP")==0 || parIni["N_Bkg_PP"]=="") { 
    parIni["N_Bkg_PP"]  = Form("%s[%.12f,%.12f,%.12f]", "N_Bkg_PP", numEntries, 0.0, numEntries*2.0);
  }
cout << "[build ctauerr] default param set" << endl;

  return;
};


void MatrixInverse2x2(TMatrixD& Inverse)
{
  Double_t A = Inverse(0,0)  ,  B = Inverse(1,0);
  Double_t C = Inverse(0,1)  ,  D = Inverse(1,1);

  Double_t det = ( ( A*D ) - ( B*C ) );

  Double_t Ainv = ( D/det )       ,  Binv = -1.0*( B/det );
  Double_t Cinv = -1.0*( C/det )  ,  Dinv = ( A/det );

  Inverse(0,0) = Ainv;   Inverse(1,0) = Binv;
  Inverse(0,1) = Cinv;   Inverse(1,1) = Dinv;

  return;
};


#endif // #ifndef buildCharmoniaCtauErrModel_C
