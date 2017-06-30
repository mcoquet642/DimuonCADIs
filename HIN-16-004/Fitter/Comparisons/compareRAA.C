#include "../Macros/CMS/CMS_lumi.C"
#include "../Macros/CMS/tdrstyle.C"
#include "datasheet.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>



using namespace std;


void compareRAA(
    const bool doprompt=true,
    const int compare=0, //0: all, 1: D meson, 2: B meson, 3: charged particles
    const int drawopt=2, //0: 4 rap bins, 1: 3 cent bins, 2: all integrated
    const bool logx=false
    ) {

  setTDRStyle();

  // Other results
  TGraphAsymmErrors *DmesonRAA_cent_0_100_1 = 
  DmesonRAA_cent_0_100_1 = new TGraphAsymmErrors(5,DmesonRAA_fx3001,DmesonRAA_fy3001,DmesonRAA_felx3001,DmesonRAA_fehx3001,DmesonRAA_fely3001,DmesonRAA_fehy3001);
  DmesonRAA_cent_0_100_1->SetName("DmesonRAA_cent_0_100_1");
  DmesonRAA_cent_0_100_1->SetTitle("DmesonRAA_cent_0_100");
  
  TGraphAsymmErrors *DmesonRAA_cent_0_100_2 = 
  DmesonRAA_cent_0_100_2 = new TGraphAsymmErrors(9,DmesonRAA_fx3002,DmesonRAA_fy3002,DmesonRAA_felx3002,DmesonRAA_fehx3002,DmesonRAA_fely3002,DmesonRAA_fehy3002);
  DmesonRAA_cent_0_100_2->SetName("DmesonRAA_cent_0_100_2");
  DmesonRAA_cent_0_100_2->SetTitle("DmesonRAA_cent_0_100");
  
  DmesonRAA_cent_0_100_1->SetMarkerStyle(kFullCircle);
  DmesonRAA_cent_0_100_1->SetMarkerColor(kOrange+2);
  DmesonRAA_cent_0_100_1->SetLineColor(kOrange+2);
  DmesonRAA_cent_0_100_1->SetFillColorAlpha(kOrange+2, 0.25);
  DmesonRAA_cent_0_100_2->SetMarkerStyle(kFullCircle);
  DmesonRAA_cent_0_100_2->SetMarkerColor(kOrange+2);
  DmesonRAA_cent_0_100_2->SetLineColor(kOrange+2);
  DmesonRAA_cent_0_100_2->SetFillColorAlpha(kOrange+2, 0.25);

  TGraphAsymmErrors *BmesonRAA_cent_0_100 = 
  BmesonRAA_cent_0_100 = new TGraphAsymmErrors(5,BmesonRAA_fx3002,BmesonRAA_fy3002,BmesonRAA_felx3002,BmesonRAA_fehx3002,BmesonRAA_fely3002,BmesonRAA_fehy3002);
  BmesonRAA_cent_0_100->SetName("BmesonRAA_cent_0_100");
  BmesonRAA_cent_0_100->SetTitle("BmesonRAA_cent_0_100");

  BmesonRAA_cent_0_100->SetMarkerStyle(kFullCircle);
  BmesonRAA_cent_0_100->SetMarkerColor(kGreen+2);
  BmesonRAA_cent_0_100->SetLineColor(kGreen+2);
  BmesonRAA_cent_0_100->SetFillColorAlpha(kGreen+2, 0.25);

  TH1D *hChgRAA_cent_0_100 = new TH1D("h_ChgRAA_cent_0_100",";p_{T} [GeV/c];R_{AA}",37,chgRAA_xAxis);
  for (int i=0; i<37; i++) {
    hChgRAA_cent_0_100->SetBinContent(i+1,chgRAA_yAxis[i]);
    hChgRAA_cent_0_100->SetBinError(i+1,chgRAA_yError[i]);
  }
  TGraphErrors *ChgRAA_cent_0_100 = new TGraphErrors(hChgRAA_cent_0_100);
  ChgRAA_cent_0_100->SetName("ChgRAA_cent_0_100");
  ChgRAA_cent_0_100->SetTitle("ChgRAA_cent_0_100");

  ChgRAA_cent_0_100->SetMarkerStyle(kFullCircle);
  ChgRAA_cent_0_100->SetMarkerColor(kViolet+9);
  ChgRAA_cent_0_100->SetLineColor(kViolet+9);
  ChgRAA_cent_0_100->SetFillColorAlpha(kViolet+9, 0.25);

  // 16025 results
  string pr_fname[3] = {
    "../resultPlots/result_JPsi_RAA_pt_prompt_JPsi_0_linearX_accEffCorr.root",
    "../resultPlots/result_JPsi_RAA_pt_prompt_JPsi_1_linearX_accEffCorr.root",
    "../resultPlots/result_JPsi_RAA_pt_prompt_JPsi_2_linearX_accEffCorr.root"
  };
  string np_fname[3] = {
    "../resultPlots/result_JPsi_RAA_pt_nonprompt_0_linearX_accEffCorr.root",
    "../resultPlots/result_JPsi_RAA_pt_nonprompt_1_linearX_accEffCorr.root",
    "../resultPlots/result_JPsi_RAA_pt_nonprompt_2_linearX_accEffCorr.root"
  };
  
  TFile inputf( doprompt ? pr_fname[drawopt].c_str() : np_fname[drawopt].c_str() );
  TFile inputf0( doprompt ? pr_fname[0].c_str() : np_fname[0].c_str() ); // to grab the lowest 3 pt bins when drawing integrated RAA vs. pt

  // Load histograms
  TCanvas *c1 = (TCanvas*)inputf.Get("c1");
  TGraphAsymmErrors *bin_0 = (TGraphAsymmErrors *)c1->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst = (TGraphAsymmErrors *)c1->FindObject("bin_0_syst");
  bin_0->SetName("bin_0");
  bin_0_syst->SetName("bin_0_syst");
  TH1F *haxes = (TH1F*)c1->FindObject("haxes");

  TCanvas *c0 = (TCanvas*)inputf0.Get("c1");
  TGraphAsymmErrors *bin_lowpt_ = (TGraphAsymmErrors *)c0->FindObject("bin_3");
  TGraphAsymmErrors *bin_lowpt_syst_ = (TGraphAsymmErrors *)c0->FindObject("bin_3_syst");

  int npoint_lowpt = bin_lowpt_->GetN();
  Double_t *x_lowpt = bin_lowpt_->GetX();
  Double_t *y_lowpt = bin_lowpt_->GetY();
  Double_t *exlow_lowpt = bin_lowpt_->GetEXlow();
  Double_t *eylow_lowpt = bin_lowpt_->GetEYlow();
  Double_t *exhigh_lowpt = bin_lowpt_->GetEXhigh();
  Double_t *eyhigh_lowpt = bin_lowpt_->GetEYhigh();
  Double_t *eylow_syst_lowpt = bin_lowpt_syst_->GetEYlow();
  Double_t *eyhigh_syst_lowpt = bin_lowpt_syst_->GetEYhigh();

  Double_t g_x_lowpt[3], g_y_lowpt[3], g_exl_lowpt[3], g_exh_lowpt[3], g_eyl_lowpt[3], g_eyh_lowpt[3], g_eyl_syst_lowpt[3], g_eyh_syst_lowpt[3];
  for (int i=0; i<3; i++) {
    g_x_lowpt[i] = x_lowpt[i];
    g_y_lowpt[i] = y_lowpt[i];
    g_exl_lowpt[i] = exlow_lowpt[i];
    g_exh_lowpt[i] = exhigh_lowpt[i];
    g_eyl_lowpt[i] = eylow_lowpt[i];
    g_eyh_lowpt[i] = eyhigh_lowpt[i];
    g_eyl_syst_lowpt[i] = eylow_syst_lowpt[i];
    g_eyh_syst_lowpt[i] = eyhigh_syst_lowpt[i];
  }
  
  TGraphAsymmErrors *bin_lowpt = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_lowpt,g_eyh_lowpt);
  TGraphAsymmErrors *bin_lowpt_syst = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_syst_lowpt,g_eyh_syst_lowpt);
  bin_lowpt->SetName("bin_lowpt");
  bin_lowpt_syst->SetName("bin_lowpt_syst");
  bin_lowpt->SetMarkerSize(1.5);
  bin_lowpt->SetMarkerStyle(kFullCross);
  bin_lowpt->SetMarkerColor(kCyan+2);
  bin_lowpt->SetLineColor(kCyan+2);
  bin_lowpt_syst->SetLineColor(kCyan+2);
  bin_lowpt_syst->SetFillColorAlpha(kCyan,0.35);


  double xpos1 = haxes->GetXaxis()->GetBinLowEdge(1);
  int lastbin = haxes->GetNbinsX();
  double xpos2 = 60;
//  double xpos2 = haxes->GetXaxis()->GetBinLowEdge(lastbin) + haxes->GetXaxis()->GetBinWidth(lastbin);
  double xaxis[2] = {0};
  xaxis[0] = logx ? 1 : xpos1;
  xaxis[1] = 60;

  TH1F *gaxes = new TH1F("gaxes","",1,xaxis);
  gaxes->GetXaxis()->CenterTitle();
  gaxes->GetXaxis()->SetTitle(haxes->GetXaxis()->GetTitle());
  if (logx) {
    gaxes->GetXaxis()->SetTitleSize(0.05);
    gaxes->GetXaxis()->SetTitleOffset(1.15);
    gaxes->GetXaxis()->SetLabelSize(0.04);
    gaxes->GetXaxis()->SetLabelOffset(0);
  }
  gaxes->GetYaxis()->SetTitle(haxes->GetYaxis()->GetTitle());
  gaxes->GetYaxis()->SetRangeUser(0, 1.6);

  TLegend leg(0.18,0.65,0.7,0.90,NULL,"brNDC");
  leg.SetMargin(0.1);
  leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.03);
  leg.SetLineColor(1);
  leg.SetLineStyle(1);
  leg.SetLineWidth(1);
  leg.SetFillColor(10);
  leg.SetFillStyle(0);

  //0: 4 rap bins, 1: 3 cent bins, 2: all integrated
  TGraphAsymmErrors *bin_1, *bin_1_syst, *bin_2, *bin_2_syst, *bin_3, *bin_3_syst;
  if (drawopt<2) {
    bin_1 = (TGraphAsymmErrors *)c1->FindObject("bin_1");
    bin_1->SetName("bin_1");
    bin_1_syst = (TGraphAsymmErrors *)c1->FindObject("bin_1_syst");
    bin_2 = (TGraphAsymmErrors *)c1->FindObject("bin_2");
    bin_2->SetName("bin_2");
    bin_2_syst = (TGraphAsymmErrors *)c1->FindObject("bin_2_syst");
    if (drawopt<1) {
      bin_3 = (TGraphAsymmErrors *)c1->FindObject("bin_3");
      bin_3->SetName("bin_3");
      bin_3_syst = (TGraphAsymmErrors *)c1->FindObject("bin_3_syst");
    }
  }

  // Draw graphs into 1 canvas
  gStyle->SetEndErrorSize(5);
  TCanvas canv("canv","canv",600,600);
  canv.SetLogx(logx);
  gaxes->Draw("");
  
  bin_0_syst->Draw("5 same");
  bin_0->Draw("p same");
  bin_lowpt_syst->Draw("5 same");
  bin_lowpt->Draw("p same");
  if (drawopt<2) {
    bin_1_syst->Draw("5 same");
    bin_1->Draw("p same");
    bin_2_syst->Draw("5 same");
    bin_2->Draw("p same");
    if (drawopt<1) {
      bin_3_syst->Draw("5 same");
      bin_3->Draw("p same");
    }
  }

  TLegendEntry *entry;

  if (drawopt==0) {
    entry = leg.AddEntry("bin_0","|y| < 0.6","p");
    entry = leg.AddEntry("bin_1","0.6 < |y| < 1.2","p");
    entry = leg.AddEntry("bin_2","1.2 < |y| < 1.8","p");
    entry = leg.AddEntry("bin_3","1.8 < |y| < 2.4","p");
  }
  else if (drawopt==1) {
    entry = leg.AddEntry("bin_0","Cent. 0-10%","p");
    entry = leg.AddEntry("bin_1","Cent. 10-30%","p");
    entry = leg.AddEntry("bin_2","Cent. 30-100%","p");
  }
  else if (drawopt==2) {
    string jpsi = doprompt? "Prompt J/#psi" : "Nonprompt J/#psi";
    entry = leg.AddEntry("bin_lowpt",Form("%s (1.8 < |y| < 2.4)",jpsi.c_str()),"p");
    entry = leg.AddEntry("bin_0",Form("%s (|y| < 2.4)",jpsi.c_str()),"p");
  }

  if (compare==1) { 
    DmesonRAA_cent_0_100_1->Draw("p 5 same");
    DmesonRAA_cent_0_100_2->Draw("p 5 same");
    entry = leg.AddEntry("DmesonRAA_cent_0_100_1","#splitline{D^{0} (|y| < 1)}{#scale[0.7]{CMS PAS HIN-16-001}}","p");
  }
  else if (compare==2) {
    BmesonRAA_cent_0_100->Draw("p 5 same");
    entry = leg.AddEntry("BmesonRAA_cent_0_100","#splitline{B^{+} (|y| < 2.4)}{#scale[0.7]{arXiv:1705.04727}}","p");
  }
  else if (compare==3) {
    ChgRAA_cent_0_100->Draw("p 5 same");
    entry = leg.AddEntry("ChgRAA_cent_0_100","#splitline{Charged hadrons (|#eta| < 1)}{#scale[0.7]{JHEP 04 (2017) 039}}","p");
  }
  else {
    DmesonRAA_cent_0_100_1->Draw("p 5 same");
    DmesonRAA_cent_0_100_2->Draw("p 5 same");
    BmesonRAA_cent_0_100->Draw("p 5 same");
    ChgRAA_cent_0_100->Draw("p 5 same");
    entry = leg.AddEntry("DmesonRAA_cent_0_100_1","D^{0} (|y| < 1) #scale[0.67]{CMS PAS HIN-16-001}","p");
    entry = leg.AddEntry("BmesonRAA_cent_0_100","B^{+} (|y| < 2.4) #scale[0.67]{arXiv:1705.04727}","p");
    entry = leg.AddEntry("ChgRAA_cent_0_100","Charged hadrons (|#eta| < 1) #scale[0.67]{JHEP 04 (2017) 039}","p");
  }

  // Draw again the main points to be on top of every comparison
  bin_0_syst->Draw("5 same");
  bin_0->Draw("p same");
  bin_lowpt_syst->Draw("5 same");
  bin_lowpt->Draw("p same");
  if (drawopt<2) {
    bin_1_syst->Draw("5 same");
    bin_1->Draw("p same");
    bin_2_syst->Draw("5 same");
    bin_2->Draw("p same");
    if (drawopt<1) {
      bin_3_syst->Draw("5 same");
      bin_3->Draw("p same");
    }
  }

  // TEXTS
  TLatex tl;
  tl.SetNDC(); tl.SetTextAlign(32);

  if (drawopt!=2) {
    if (doprompt) tl.DrawLatexNDC(0.92,0.69,"Prompt J/#psi");
    else tl.DrawLatexNDC(0.92,0.69,"Nonprompt J/#psi");
  }
  CMS_lumi( (TPad*) gPad, 106, 33, "" ); 

  // draw a line at RAA=1
  TLine line; line.SetLineColor(kBlack); line.SetNDC();
  line.DrawLine(xpos1,1,xpos2,1);
  canv.Update();
  leg.Draw("same");

  // Global uncertainty boxes
  TBox *box0, *box1, *box2, *box3, *box_lowpt;
  if (drawopt==0) {
    if (doprompt) {
      box0 = new TBox(0, 0.951740, 2.5, 1.044238);
      box1 = new TBox(2.5, 0.951740, 5, 1.044238);
      box2 = new TBox(5, 0.951740, 7.5, 1.044238);
      box3 = new TBox(7.5, 0.951740, 10, 1.044238);
    } else {
      box0 = new TBox(0, 1, 2.5, 1);
      box1 = new TBox(2.5, 1, 5, 1);
      box2 = new TBox(5, 1, 7.5, 1);
      box3 = new TBox(7.5, 1, 10, 1);
    }
    box0->SetLineColor(kRed);
    box0->SetFillColorAlpha(kRed,0.35);
    box0->Draw("lf");
    box1->SetLineColor(kBlue+3);
    box1->SetFillColorAlpha(kBlue+3,0.35);
    box1->Draw("lf");
    box2->SetLineColor(kGreen);
    box2->SetFillColorAlpha(kGreen,0.35);
    box2->Draw("lf");
    box3->SetLineColor(kCyan);
    box3->SetFillColorAlpha(kCyan,0.35);
    box3->Draw("lf");
  } else if (drawopt==1) {
    if (doprompt) {
      box0 = new TBox(0, 0.951740, 2.5, 1.044238);
      box1 = new TBox(2.5, 0.951740, 5, 1.044238);
      box2 = new TBox(5, 0.951740, 7.5, 1.044238);
    } else {
      box0 = new TBox(0, 1, 2.5, 1);
      box1 = new TBox(2.5, 1, 5, 1);
      box2 = new TBox(5, 1, 7.5, 1);
    }
    box0->SetLineColor(kRed);
    box0->SetFillColorAlpha(kRed,0.35);
    box0->Draw("lf");
    box1->SetLineColor(kBlue+3);
    box1->SetFillColorAlpha(kBlue+3,0.35);
    box1->Draw("lf");
    box2->SetLineColor(kGreen);
    box2->SetFillColorAlpha(kGreen,0.35);
    box2->Draw("lf");
  } else if (drawopt==2) {
    if (doprompt) {
      box_lowpt = new TBox(0, 0.951740, 2.5, 1.044238);
      box0 = new TBox(2.5, 0.951740, 5, 1.044238);
    } else {
      box_lowpt = new TBox(0, 1, 2.5, 1);
      box0 = new TBox(2.5, 1, 5, 1);
    }
    box_lowpt->SetLineColor(kCyan);
    box_lowpt->SetFillColorAlpha(kCyan,0.35);
    box_lowpt->Draw("lf");
    box0->SetLineColor(kRed+2);
    box0->SetFillColorAlpha(kRed,0.35);
    box0->Draw("lf");
  }
  canv.RedrawAxis();

  // Save canvas
  string outputname = "RAA_pt_";
  outputname += Form("%s",doprompt?"pr_":"np_");
  if (drawopt==0) outputname += "4raps_";
  else if (drawopt==1) outputname += "3cents_";
  else if (drawopt==2) outputname += "integrated_";
  if (compare==0) outputname += "all";
  else if (compare==1) outputname += "D";
  else if (compare==2) outputname += "B";
  else if (compare==3) outputname += "Chg";
  if (logx) outputname += "_log.pdf";
  else outputname += "_lin.pdf";

  cout << outputname << endl;
  canv.SaveAs(outputname.c_str());


  return ;
}
