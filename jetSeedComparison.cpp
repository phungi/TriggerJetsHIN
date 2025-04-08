#if !defined(__CINT__) || defined(__MAKECINT__)                                                       
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#endif

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06)
{
  TLatex* latex = new TLatex(x, y, text);                                                                                                
  latex->SetNDC();                                                                                                                       
  latex->SetTextSize(textSize);                                                                                                          
  latex->SetTextColor(color);                                                                                                            
  latex->SetTextFont(42);                                                                                                                
  latex->Draw();                                                                                                                         
}                                                                                                  


void jetSeedComparison(){
  Int_t ci = 2555;
  // colours and names of the sampleskaaskopPink
  new TColor(++ci, 54/256., 161/256., 24/256.); // dark green
  new TColor(++ci, 252/256., 161/256., 3/256.); // dark orange
  new TColor(++ci, 173/256., 16/256., 204/256.); // purple
  new TColor(++ci, 132/256., 240/256., 209/256.); // cyan
  new TColor(++ci, 36/256., 116/256., 156/256.); // deep blue
  new TColor(++ci, 242/256., 235/256., 36/256.); // yellow  
  new TColor(++ci, 200/256., 135/256., 136/256.); //   
  new TColor(++ci, 136/256., 200/256., 36/256.); //   
  ci = 2555;
  std::vector < Int_t > colours = { ci+1, ci+2, ci+3, ci+4, ci+5, ci+6, ci+7 };
  // TString fname1="data_23_ppring_newjec.root";
  // TString fname2="data_23_ppring_nonewjec.root";
  TH1::SetDefaultSumw2();
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TH1F *dummy_histo = new TH1F("","",100,0,200);
  dummy_histo->GetYaxis()->SetTitleOffset(0.8);
  dummy_histo->GetXaxis()->SetTitleOffset(0.95);
  dummy_histo->GetYaxis()->SetTitle("Trigger efficiency");
  dummy_histo->GetXaxis()->SetTitle("corrected p_{T,PF}^{lead} [GeV]");
  // dummy_histo->GetXaxis()->SetTitle("corrected p_{T,PF}^{subl.} [GeV]");
  dummy_histo->GetXaxis()->SetLabelFont(42);
  dummy_histo->GetYaxis()->SetLabelFont(42);
  dummy_histo->GetXaxis()->SetLabelSize(0.04);
  dummy_histo->GetYaxis()->SetLabelSize(0.04);
  dummy_histo->GetXaxis()->SetTitleFont(42);
  dummy_histo->GetYaxis()->SetTitleFont(42);
  dummy_histo->GetXaxis()->SetTitleSize(0.045);
  dummy_histo->GetYaxis()->SetTitleSize(0.045);
  dummy_histo->GetYaxis()->SetRangeUser(0,1.25);
  dummy_histo->GetXaxis()->SetRangeUser(0,80);

  TString out_file = "SingleJet_efficiencies_PbPb.pdf";
  std::vector<Double_t> thresholds = {0.5,1,2,4,6,8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  std::vector<TString> triggers = {};
  for(size_t i{0}; i< thresholds.size(); ++i){
    TString a = Form("L1_SingleJet%i_BptxAND",thresholds.at(i));
    triggers.push_back(a);
  }
  TCanvas *c11 = new TCanvas();
  c11->Print( out_file + "(","pdf");
  // TString fname1="ChunkyRerun.root";
  // TString fname2="2023_pp_ZeroBias_new_.root";
  std::cout << "before opening file" << std::endl;
  // TFile *input1 = new TFile(fname1, "READ");
  TString in_dir = "/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/";
  std::vector<TString> names = {
    // "jS=0.5 w/o JEC @ L1",
    // "jS=1.0 w/o JEC @ L1",
    // "jS=1.5 w/o JEC @ L1",
    // "jS=2.0 w/o JEC @ L1",
    // "jS=2.5 w/o JEC @ L1",
    // "jS=1.5 w/ JEC",
    // "jS=2.0 w/ JEC",
    // "jS=2.5 w/ JEC @ L1",
    // "jS=2.5 w/ #phi-ring",
    // "both_L1_2p5_noPhi_all.root"
    // "jS=2.5 w/o JEC",
    // "jS=2.5 w/ JEC no #phi",
    // "jS=2.5 w/ JEC with #phi",
    "h"
  };
  // TString ntuple_dir = "/eos/cms/store/group/phys_heavyions/vavladim/L1_ntuples_jS/";
  TString ntuple_dir = "/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/";

  std::vector<TString> in_files = {
    //double jet triggers, no JEC on L1, fewer bins
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_0p5jS_noJEC.root_doubleL1triggers_histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_1p0jS_noJEC.root_doubleL1triggers_histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_1p5jS_noJEC.root_doubleL1triggers_histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_2p0jS_noJEC.root_doubleL1triggers_histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_2p5jS_noJEC.root_doubleL1triggers_histos.root",
    // ntuple_dir + "L1Ntuple_2p5_no_Phi_allEvents.root_doubleL1triggers_histos.root"

    //single jet triggers, no JEC on L1, fewer bins
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_0p5jS_noJEC.root.histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_1p0jS_noJEC.root.histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_1p5jS_noJEC.root.histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_2p0jS_noJEC.root.histos.root",
    // ntuple_dir + "L1Ntuple_gg2uu_NoPhi_2p5jS_noJEC.root.histos.root",
    // ntuple_dir + "L1Ntuple_2p5_no_Phi_allEvents.root.histos2.root",

    "PbPb_387879.root"

    //single jet triggers below
    //some of those are pp MC, with/without phi, JEC on L1
    // ntuple_dir + "L1Ntuple_2p5_incl_Phi_allEvents.root.histos.root",
    // "SJ_1p5_noPhi.root",
    // "SJ_1p5_Phi.root",
    // "SJ_2p0_noPhi.root",
    // "SJ_2p0_Phi.root",
    // "SJ_2p5_noPhi.root",
    // "SJ_2p5_Phi.root"
    // "1p5_noPhi_SJet_gg2uu.root",
    // "1p5_withPhi_SJet_gg2uu.root",
    // "2p0_withPhi_SJet_gg2uu.root",
    // "2p0_noPhi_SJet_gg2uu.root",
    // "2p5_noPhi_SJet_gg2uu.root",
    // "2p5_withPhi_SJet_gg2uu.root"
    // "both_L1_1p5_noPhi_all.root",
    // "both_L1_1p5_withPhi_all.root",
    // "both_L1_2p0_noPhi_all.root",
    // "both_L1_2p0_withPhi_all.root",
    // "both_L1_2p5_noPhi_all.root",
    // "both_L1_2p5_withPhi_all.root"

    // "1p5_noPhi_all.root",
    // "1p5_withPhi_all.root",
    // "2p0_noPhi_all.root",
    // "2p0_withPhi_all.root",
    // "two_jet_trigger_gg2uu/2p5_noPhi_all.root",
    // "2p5_withPhi_all.root",
    // "both_L1_2p5_noPhi_all.root"
    // "Seed2_withPhi.root",
    // "2p5_noPhi_all.root",
    // "2p5_withPhi_all.root"
  };
  std::vector<TFile *> inputs = {};
  for(size_t i{0}; i < in_files.size(); ++i){
    // in_files.at(i) = ntuple_dir + in_files.at(i);
    TFile *iF = new TFile(in_files.at(i), "READ");
    inputs.push_back(iF);
  }
  std::cout << "out of file loop " << std::endl;
  TCanvas *can = new TCanvas();
  can->cd();
  can->SetTicks();
  
  for(Int_t t{0}; t < triggers.size(); ++t){
    std::vector<TH1F*> histo = {};
    std::vector<TH1F*> hjet_all = {};
    std::vector<TEfficiency *> eff = {};
    TLegend *lego;
    lego = new TLegend(0.65, 0.45, 0.75, 0.55);
    // lego->SetHeader("     no phiRing PUS","L");
    lego->SetBorderSize(0);
    lego->SetTextSize(0.03);
    lego->SetTextFont(42);
    lego->SetBorderSize(0);
    lego->SetFillStyle(0);
    lego->SetFillColor(0);
    can->cd();
    dummy_histo->Draw();
    TLine *Line=new TLine(0, 1, 80, 1);
    Line->SetLineColorAlpha(kRed,0.2);
    Line->SetLineStyle(9);
    Line->SetLineWidth(2);
    Line->Draw("SAME");
    for(size_t i{0}; i < inputs.size(); ++i){
      std::cout << "here?" << std::endl;
      TH1F *denom = (TH1F*)inputs.at(i)->Get("jetpt_all_pf");
      denom->GetXaxis()->SetRangeUser(1, 20);
      hjet_all.push_back(denom);
      TH1F *num = (TH1F*)inputs.at(i)->Get(Form("jetpt_trig_calo%i",t));
      num->GetXaxis()->SetRangeUser(1, 20);
      histo.push_back(num);
      if(TEfficiency::CheckConsistency(*(TH1*)num,*(TH1*)denom)){
        TEfficiency *effi = new TEfficiency(*(TH1*)num,*(TH1*)denom);
        eff.push_back(effi);
        eff.at(i)->SetMarkerSize(0.5);
        eff.at(i)->SetMarkerStyle(20);
        eff.at(i)->SetMarkerColor(colours.at(i));
        eff.at(i)->SetLineColor(colours.at(i));
        lego->AddEntry(eff.at(i), names.at(i), "LP" );
        can->cd();
        eff.at(i)->Draw("SAME");
      }
    }
    lego->Draw();
    DrawLatex(0.15,0.82,1,Form("PbPb, PF#rightarrowL1 matching"),0.05);
    // DrawLatex(0.15,0.77,1,Form("double jet trigger, p_{T}>%0.1f GeV", thresholds.at(t)),0.05);
    DrawLatex(0.15,0.77,1,Form("single jet trigger, p_{T}>%0.1f GeV", thresholds.at(t)),0.05);
    
    can->Print( out_file,"pdf");
    can->Clear();
  }
  can->Print( out_file + ")","pdf");
  exit(0);
}
