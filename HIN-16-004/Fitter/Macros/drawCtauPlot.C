#ifndef drawCtauPlot_C
#define drawCtauPlot_C

#include "Utilities/initClasses.h"

void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, vector<double> rangeErr, double excEvts=0.0);
void printCtauParameters(RooWorkspace myws, TPad* Pad, string pdfName, bool isWeighted);

void drawCtauPlot(RooWorkspace& myws,   // Local workspace
                  string outputDir,     // Output directory
		  struct InputOpt opt,  // Variable with run information (kept for legacy purpose)
                  struct KinCuts cut,   // Variable with current kinematic cuts
                  map<string, string>  parIni,   // Variable containing all initial parameters
                  string plotLabel,     // The label used to define the output file name
                  // Select the type of datasets to fit
                  string DSTAG,         // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                  // Select the type of object to fit
                  bool incJpsi,         // Includes Jpsi model
                  bool incBkg,          // Includes Background model        
                  bool incPrompt,       // Includes Prompt ctau model       
                  bool incNonPrompt,    // Includes Non-Prompt ctau model
                  // Select the fitting options
                  bool useSPlot,        // If yes, then use SPlot ds, if no, use mass range one
                  bool plotPureSMC,     // Flag to indicate if we want to fit pure signal MC
                  // Select the drawing options
                  bool setLogScale,     // Draw plot with log scale
                  bool incSS,           // Include Same Sign data
                  double binWidth       // Bin width
                  ) 
{

//  bool usectauBkgTemplate = incBkg && !(myws.pdf(Form("pdfCTAU_BkgPR_%s", "PP"))) && !(myws.pdf(Form("pdfCTAU_BkgNoPR_%s", "PP")));
  bool usectauBkgTemplate = false;
bool usePerEventError = false;
  
  RooMsgService::instance().getStream(0).removeTopic(Caching);  
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  string pdfTotName  = Form("pdfCTAU_Tot_%s", "PP");
  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), "PP");
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), "PP");
  string dsOSNameCut = dsOSName+"_CTAUCUT";
  string dsOSName2Fit = dsOSName;
  if (useSPlot && incBkg) dsOSName2Fit += "_BKG";
  else if (useSPlot && incJpsi) dsOSName2Fit += "_JPSI";
  string dsOSName2FitCut = dsOSName2Fit+"_CTAUCUT";
  string hOSName = Form("dhCTAUERRTot_Tot_%s", "PP");
  string hOSNameBkg  = Form("dhCTAUERR_Bkg_%s", "PP");
  string hOSNameJpsi = Form("dhCTAUERR_Jpsi_%s", "PP");
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), "PP");

  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();
  bool isMC = (DSTAG.find("MC")!=std::string::npos);
  vector<double> range; range.push_back(cut.dMuon.ctau.Min); range.push_back(cut.dMuon.ctau.Max);

  double minRange = (double)(floor(range[0]*10.)/10.)+0.1;
  double maxRange = (double)(ceil(range[1]*10.)/10.)+0.1;
  if (!incNonPrompt) {
    if (abs(maxRange)>abs(minRange)) { minRange = -1.0*abs(maxRange); } else { maxRange = abs(minRange); }
  } else {
    minRange = -4.0;
    maxRange = 7.0;
  }
  Double_t numTot = myws.data(dsOSName2Fit.c_str())->sumEntries();
  Double_t outTot = myws.data(dsOSName2Fit.c_str())->numEntries();
  Double_t outErr = myws.data(dsOSName2Fit.c_str())->reduce(Form("(ctau>%.6f || ctau<%.6f)", range[1], range[0]))->numEntries();
  int nBins = min(int( round((maxRange - minRange)/binWidth) ), 1000);

  double normDSTot   = 1.0;  if (myws.data(dsOSName2FitCut.c_str()))  { normDSTot   = myws.data(dsOSName2Fit.c_str())->sumEntries()/myws.data(dsOSName2FitCut.c_str())->sumEntries();  }
  double normJpsi  = 1.0;  if (myws.data(hOSNameJpsi.c_str()))  { normJpsi  = myws.data(dsOSName2Fit.c_str())->sumEntries()*normDSTot/myws.data(hOSNameJpsi.c_str())->sumEntries();  }
  double normBkg   = 1.0;  if (myws.data(hOSNameBkg.c_str()))   { normBkg   = myws.data(dsOSName2Fit.c_str())->sumEntries()*normDSTot/myws.data(hOSNameBkg.c_str())->sumEntries();   }
  double normTot   = 1.0;  if (myws.data(hOSName.c_str()))  { normTot   = myws.data(dsOSName2Fit.c_str())->sumEntries()*normDSTot/myws.data(hOSName.c_str())->sumEntries();  }
cout << "normBkg : " << normBkg << ", normTot : " << normTot << ", numTot : " <<  numTot << endl;

  // Create the main plot of the fit
  RooPlot*   frame = myws.var("ctau")->frame(Bins(nBins), Range(minRange, maxRange));
// // FUntionality not working yet
//  RooPlot*   frame(0x0);
//  if (!usectauBkgTemplate) frame = myws.var("ctau")->frame(Bins(nBins), Range(minRange, maxRange));
//  else frame = myws.var("ctau")->frame(Binning("TemplateBinning"),Range(minRange, maxRange));
  
  frame->updateNormVars(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr"))) ;
  myws.data(dsOSName2Fit.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  if (incBkg && !incJpsi) {
    if (!usectauBkgTemplate)
    {
	if (usePerEventError){
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(32),
//                                           ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
//                                           FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), Precision(1e-4)
                                           LineColor(kAzure-9), Precision(1e-4)
                                           );
      myws.data(dsOSName2Fit.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
      if (myws.pdf(Form("pdfCTAU_BkgPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_BkgPR_%s", "PP")))),
//                                                                                                          ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                                                                                          Normalization(normBkg, RooAbsReal::NumEvent),
                                                                                                          LineColor(kRed+2), Precision(1e-4), NumCPU(32)
                                                                                                          );
      }else{ cout << "pdfCTAUCOND_BkgPR nt found" << endl;}
      if (myws.pdf(Form("pdfCTAU_BkgNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_BkgNoPR_%s", "PP")))),
//                                                                                                            ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                                                                                            Normalization(normBkg, RooAbsReal::NumEvent),
                                                                                                            LineColor(kOrange+10), Precision(1e-4), NumCPU(32)
                                                                                                            );
      }else{ cout << "pdfCTAUCOND_BkgNoPR nt found" << endl;}
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(32),
 //                                          ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                           LineColor(kBlack), Precision(1e-4)
                                           );
	}else{

      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),  Normalization(numTot, RooAbsReal::NumEvent), NumCPU(32),
//                                           FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), Precision(1e-4)
                                           LineColor(kAzure-9), Precision(1e-4)
                                           );
      myws.data(dsOSName2Fit.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
      if (myws.pdf(Form("pdfCTAUCOND_BkgPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUCOND_BkgPR_%s", "PP")))),
                                                                                                          Normalization(numTot, RooAbsReal::NumEvent),
                                                                                                          LineColor(kRed+2), Precision(1e-4), NumCPU(32)
                                                                                                          );
      }else{ cout << "pdfCTAUCOND_BkgPR nt found" << endl;}
/*      if (myws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", "PP")))),
                                                                                                            Normalization(numTot, RooAbsReal::NumEvent),
                                                                                                            LineColor(kOrange+10), Precision(1e-4), NumCPU(32)
                                                                                                            );
      }else{ cout << "pdfCTAUCOND_BkgNoPR nt found" << endl;}*/
      if (myws.pdf(Form("pdfCTAUDSS_BkgNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPRDSS"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUDSS_BkgNoPR_%s", "PP")))),
                                                                                                            Normalization(numTot, RooAbsReal::NumEvent),
                                                                                                            LineColor(kYellow+10), Precision(1e-4), NumCPU(32)
                                                                                                            );
      }else{ cout << "pdfCTAUCOND_BkgNoPR nt found" << endl;}
      if (myws.pdf(Form("pdfCTAUDF_BkgNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPRDF"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUDF_BkgNoPR_%s", "PP")))),
                                                                                                            Normalization(numTot, RooAbsReal::NumEvent),
                                                                                                            LineColor(kGreen+10), Precision(1e-4), NumCPU(32)
                                                                                                            );
      }else{ cout << "pdfCTAUCOND_BkgNoPR nt found" << endl;}
      if (myws.pdf(Form("pdfCTAUDDS_BkgNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPRDDS"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUDDS_BkgNoPR_%s", "PP")))),
                                                                                                            Normalization(numTot, RooAbsReal::NumEvent),
                                                                                                            LineColor(kGreen+10), Precision(1e-4), NumCPU(32)
                                                                                                            );
      }else{ cout << "pdfCTAUCOND_BkgNoPR nt found" << endl;}
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(32),
                                           LineColor(kBlack), Precision(1e-4)
                                           );
	}
    }
    else
    {
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(numTot, RooAbsReal::NumEvent), NumCPU(32),
                                           LineColor(kBlack), Precision(1e-4)
                                           );
    }
  }
  if (useSPlot) {
    if (incJpsi) {
	cout << "hOSNameJpsi.c_str() : " << hOSNameJpsi.c_str() << endl;
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSI"),  Normalization(normJpsi, RooAbsReal::NumEvent), NumCPU(32),
//                                           ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameJpsi.c_str()), kTRUE),
                                           FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), Precision(1e-4)
                                           );
      myws.data(dsOSName2Fit.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
      if (myws.pdf(Form("pdfCTAU_JpsiPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSIPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_JpsiPR_%s", "PP")))),
//                                                                                                           ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameJpsi.c_str()), kTRUE),
                                                                                                           Normalization(normJpsi, RooAbsReal::NumEvent),
                                                                                                           LineColor(kRed+2), Precision(1e-4), NumCPU(32)
                                                                                                           );
      }
      if (myws.pdf(Form("pdfCTAU_JpsiNoPR_%s", "PP"))) {myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_JpsiNoPR_%s", "PP")))),
//                                                                                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameJpsi.c_str()), kTRUE),
                                                                                                             Normalization(normJpsi, RooAbsReal::NumEvent),
                                                                                                             LineColor(kOrange+10), Precision(1e-4), NumCPU(32)
                                                                                                             );
      }
      myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(normJpsi, RooAbsReal::NumEvent), NumCPU(32),
//                                           ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameJpsi.c_str()), kTRUE),
                                           LineColor(kBlack), Precision(1e-4)
                                           );
    }
  }
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName2Fit.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  // set the CMS style
  setTDRStyle();

  // Create the pull distribution of the fit 
  RooHist *hpull = frame->pullHist(0, "PDF", true);
  hpull->SetName("hpull");
  RooPlot* frame2 = myws.var("ctau")->frame(Title("Pull Distribution"), Bins(nBins), Range(minRange, maxRange));
// // FUntionality not working yet
//  RooPlot* frame2(0x0);
//  if (!usectauBkgTemplate) frame2 = myws.var("ctau")->frame(Title("Pull Distribution"), Bins(nBins), Range(minRange, maxRange));
//  else frame2 = myws.var("ctau")->frame(Title("Pull Distribution"),Binning("TemplateBinning"),Range(minRange, maxRange));
  frame2->addPlotable(hpull, "PX");
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cCtauFig_%s", "PP"), "cCtauFig",800,800);
  TPad    *pad1  = new TPad(Form("pad1_%s", "PP"),"",0,0.23,1,1);
  TPad    *pad2  = new TPad(Form("pad2_%s", "PP"),"",0,0,1,.228);
  TLine   *pline = new TLine(minRange, 0.0, maxRange, 0.0);
  
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetTitleOffset(3);
  frame->GetXaxis()->SetLabelOffset(3);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.7);
  frame->GetYaxis()->SetTitleFont(42);
  setCtauRange(myws, frame, dsOSName2FitCut, setLogScale, range, outErr);
 
  cFig->cd();
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000); 
  pad1->SetBottomMargin(0.015); 
  //plot fit
  pad1->Draw();
  pad1->cd(); 
  frame->Draw();

  if (!usectauBkgTemplate) printCtauParameters(myws, pad1, pdfTotName, isWeighted);
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f < #chi^{2}_{MFT-MCH} < %.1f",cut.dMuon.Chi2.Min,cut.dMuon.Chi2.Max)); dy+=0.045;
  if (outErr>0.0) {
      t->DrawLatex(0.21, 0.86-dy, Form("Excl: (%.4f%%) %.0f evts", (outErr*100.0/outTot), outErr)); dy+=1.5*0.045;
  }

  // Drawing the Legend
  double ymin = 0.7802;
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if(frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),(usectauBkgTemplate&&!incJpsi)?"Bkg template":"Total fit","l"); }
  if((incBkg && incJpsi) && frame->findObject("BKG")) { leg->AddEntry(frame->findObject("BKG"),"Background","fl");  }
  if(incBkg && incJpsi && frame->findObject("JPSI")) { leg->AddEntry(frame->findObject("JPSI"),"J/#psi PDF","l"); }
  leg->Draw("same");

  //Drawing the title
  TString label;
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  
  gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  //---plot pull
  pad2->Draw();
  pad2->cd();
    
  frame2->SetTitle("");
  frame2->GetYaxis()->CenterTitle(kTRUE);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.1);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetXaxis()->CenterTitle(kTRUE);
  frame2->GetXaxis()->SetTitleOffset(1);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetXaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

  frame2->Draw(); 
  
  // *** Print chi2/ndof 
  if (!usectauBkgTemplate) printChi2(myws, pad2, frame, "ctau", dsOSName2Fit.c_str(), pdfTotName.c_str(), nBins);
  
  pline->Draw("same");
  pad2->Update();
  
  bool SB = (incBkg&&!incJpsi);
  // Save the plot in different formats
  if (!usectauBkgTemplate)
  {
    gSystem->mkdir(Form("%sctau%s/%s/plot/root/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%s/%s/plot/root/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
    gSystem->mkdir(Form("%sctau%s/%s/plot/png/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%s/%s/plot/png/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.png", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
    gSystem->mkdir(Form("%sctau%s/%s/plot/pdf/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%s/%s/plot/pdf/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.pdf", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
  }
  else
  {
    gSystem->mkdir(Form("%sctau%sTemp/%s/plot/root/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%sTemp/%s/plot/root/PLOT_%s_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.root", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", "Bkg", (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
    gSystem->mkdir(Form("%sctau%sTemp/%s/plot/png/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%sTemp/%s/plot/png/PLOT_%s_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.png", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", "Bkg", (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
    gSystem->mkdir(Form("%sctau%sTemp/%s/plot/pdf/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
    cFig->SaveAs(Form("%sctau%sTemp/%s/plot/pdf/PLOT_%s_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_chi2%.0f%.0f.pdf", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), "PP", "Bkg", (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, (cut.dMuon.Chi2.Min*10.0), (cut.dMuon.Chi2.Max*10.0)));
  }
  cFig->Clear();
  cFig->Close();
  
}

void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, vector<double> rangeErr, double excEvts)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("ctau"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));

  Double_t Yup(0.),Ydown(0.);
  
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) isMC = true;

  if(setLogScale)
  {
    if (isMC) {
      Yup = YMax*TMath::Power((YMax/YMin), (0.5/(1.0-0.5-0.2)));
      Ydown = YMin/(TMath::Power((YMax/YMin), (0.2/(1.0-0.5-0.2))));
    } 
    else {
      Yup = YMax*TMath::Power((YMax/0.01), 0.5);
      Ydown = 0.01;
    }
  }
  else
  {
    Yup = YMax+(YMax-0.0)*0.5;
    Ydown = 0.0;
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;

  if (excEvts>0.0) {
    TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.4)):(Ydown + (Yup-Ydown)*0.4)));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    frame->addObject(minline);
    TLine   *maxline = new TLine(rangeErr[1], 0.0, rangeErr[1], (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.4)):(Ydown + (Yup-Ydown)*0.4)));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    frame->addObject(maxline);
  }
};


void printCtauParameters(RooWorkspace myws, TPad* Pad, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.045; 
  RooArgSet* Parameters = (RooArgSet*)myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->selectByAttrib("Constant",kFALSE);
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="invMass"){continue;} else if(s1=="ctau"){continue;} else if(s1=="MassRatio"){continue;}   
    else if(s1=="One"){continue;} else if(s1=="mMin"){continue;} else if(s1=="mMax"){continue;}
    else if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){
      s1=Form("#%s",s1.c_str());
    }

    if(s2=="CtauRes")  { s2="Res";   } 
    else if(s2=="JpsiNoPR")  { s2="J/#psi[NoPR]";   } 
    else if(s2=="JpsiPR")  { s2="J/#psi[PR]";   } 
    else if(s2=="Jpsi" && (s1=="N" || s1=="b"))  { s2="J/#psi";   } 
    else if(s2=="BkgNoPR")   { s2="bkg[NoPR]";      } 
    else if(s2=="BkgPR")   { s2="bkg[PR]";      } 
    else if(s2=="Bkg" && (s1=="N" || s1=="b"))   { s2="bkg";      }
    else {continue;}
    if(s3!=""){
      label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());
    }
    // Print the parameter's results
    if(s1=="N"){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.0f#pm%.0f ", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1=="b"){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f ", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("sigma")!=std::string::npos) { 
      if (s2.find("Res")!=std::string::npos) {
        t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f", label.c_str(), it->getValV(), it->getError())); dy+=0.045;
      } else {
        t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f mm", label.c_str(), it->getValV(), it->getError())); dy+=0.045;
      }
    }
    else if(s1.find("lambda")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("ctau")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f mm", (label.insert(1, string("#"))).c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else { 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
  }
};


#endif // #ifndef drawCtauPlot_C
