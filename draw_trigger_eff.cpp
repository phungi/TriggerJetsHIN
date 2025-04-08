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


void draw_trigger_eff(){
  Int_t ci = 2555;
  // colours and names of the sampleskaaskopPink
  new TColor(++ci, 54/256., 161/256., 24/256.); // dark green
  new TColor(++ci, 252/256., 161/256., 3/256.); // dark orange
  new TColor(++ci, 173/256., 16/256., 204/256.); // purple
  new TColor(++ci, 132/256., 240/256., 209/256.); // cyan
  new TColor(++ci, 36/256., 116/256., 156/256.); // deep blue
  new TColor(++ci, 242/256., 235/256., 36/256.); // yellow  
  ci = 2555;
  std::vector < Int_t > colours = { ci+1, ci+2, ci+3, ci+4, ci+5, ci+6 };
  // TString fname1="data_23_ppring_newjec.root";
  // TString fname2="data_23_ppring_nonewjec.root";
  TH1::SetDefaultSumw2();
  gROOT->SetBatch(kTRUE);  
  std::string out_file = "2024_PbPb_L1_turnons_387973_highest_full_lead_only_CHF.pdf";
  std::vector<Int_t> thresholds = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  std::vector<TString> triggers = {};
  for(size_t i{0}; i< thresholds.size(); ++i){
    TString a = Form("L1_SingleJet%i_BptxAND",thresholds.at(i));
    triggers.push_back(a);
  }
  TCanvas *c11 = new TCanvas();
  c11->Print( (out_file + "(").c_str(),"pdf");
  // TString fname1="ChunkyRerun.root";
  // TString fname2="2023_pp_ZeroBias_new_.root";
  std::cout << "before opening file" << std::endl;
  // TFile *input1 = new TFile(fname1, "READ");
  TFile *input2 = new TFile("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/CHF_cut_full_run.root","READ");
  // TFile *input3;
  std::cout << "list?" << std::endl;
  // int triggers.size()=14;
  TH1F *dummy_histo = new TH1F("","",100,0,200);
  TH1F *histo[1][22];
  TH1F *hjet_all[1]; 
  TH1F *histo2[1][22];
  TEfficiency *eff[1][22];
  // TEfficiency *eff[22];
  gStyle->SetOptTitle(0);//dsa
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(0);
  // hjet_all[0]=(TH1F*)input1->Get("jetpt_all_calo");
  hjet_all[0]=(TH1F*)input2->Get("jetpt_all_pf");
  hjet_all[0]->GetXaxis()->SetRangeUser(1, 20);
  for(int i=0;i<triggers.size();i++){
    // histo[0][i]=(TH1F*)input1->Get(Form("jetpt_trig_calo%i",i));
    histo[0][i]=(TH1F*)input2->Get(Form("jetpt_trig_calo%i",i));
    histo[0][i]->GetXaxis()->SetRangeUser(1, 20);
  }

  for(int j=0;j<1;j++){
    for(int i=0;i<triggers.size();i++){
      histo2[j][i]=(TH1F*)histo[j][i]->Clone(Form("histo2[%i][%i]",j,i));
      histo2[j][i]->SetDirectory(0);
      // TEfficiency *eff = 0;
      // eff[j][i]->SetDirectory(0);
      if(TEfficiency::CheckConsistency(*(TH1*)histo2[j][i],*(TH1*)hjet_all[j])){
        eff[j][i] = new TEfficiency(*(TH1*)histo2[j][i],*(TH1*)hjet_all[j]);
      }
      histo2[j][i]->Divide(histo[j][i],hjet_all[j],1,1,"b");
      // histo2[j][i]->GetXaxis()->SetRangeUser(1, 20);
      // eff[i]= new TEfficiency(*histo[j][i],*hjet_all[j]);
      // histo2[j][i]->Divide(histo[j][i],hjet_all[j],1,1,"pois");
      // std::cout << "errors " << 
      // histo2[j][i]->SetLineWidth(3);
      // cout<<"what "<<histo2[j][i]<<endl;
    }
  }
  TCanvas *can = new TCanvas();
  can->cd();
  // TPad *pad = new TPad("pad0","This is pad0",0.,0.,1,1.);                                                                                     
  TLegend *lego;
  lego = new TLegend(0.65, 0.74, 0.75, 0.9);
  lego->SetBorderSize(0);
  lego->SetTextSize(0.03);
  lego->SetTextFont(42);
  lego->SetBorderSize(0);
  lego->SetFillStyle(0);
  can->SetTicks();
  dummy_histo->GetYaxis()->SetTitleOffset(0.8);
  dummy_histo->GetXaxis()->SetTitleOffset(0.95);
  dummy_histo->GetYaxis()->SetTitle("Trigger efficiency");
  dummy_histo->GetXaxis()->SetTitle("uncorr. lead jet p_{T,PF   }[GeV]");
  dummy_histo->GetXaxis()->SetLabelFont(42);
  dummy_histo->GetYaxis()->SetLabelFont(42);
  dummy_histo->GetXaxis()->SetLabelSize(0.04);
  dummy_histo->GetYaxis()->SetLabelSize(0.04);
  dummy_histo->GetXaxis()->SetTitleFont(42);
  dummy_histo->GetYaxis()->SetTitleFont(42);
  dummy_histo->GetXaxis()->SetTitleSize(0.045);
  dummy_histo->GetYaxis()->SetTitleSize(0.045);
  dummy_histo->GetYaxis()->SetRangeUser(0,1.25);
  dummy_histo->GetXaxis()->SetRangeUser(0,150);
  for(int k=0;k<triggers.size();k++){
    can->cd();
    // histo2[0][k]->GetYaxis()->SetTitleOffset(0.8);
    // histo2[0][k]->GetXaxis()->SetTitleOffset(0.95);
    // histo2[0][k]->GetYaxis()->SetTitle("Trigger efficiency");
    // histo2[0][k]->GetXaxis()->SetTitle("p_{T,PF   }[GeV]");
    // histo2[0][k]->GetXaxis()->SetLabelFont(42);
    // histo2[0][k]->GetYaxis()->SetLabelFont(42);
    // histo2[0][k]->GetXaxis()->SetLabelSize(0.04);
    // histo2[0][k]->GetYaxis()->SetLabelSize(0.04);
    // histo2[0][k]->GetXaxis()->SetTitleFont(42);
    // histo2[0][k]->GetYaxis()->SetTitleFont(42);
    // histo2[0][k]->GetXaxis()->SetTitleSize(0.045);
    // histo2[0][k]->GetYaxis()->SetTitleSize(0.045);
    // histo2[0][k]->GetYaxis()->SetRangeUser(0,1.25);
    // histo2[0][k]->GetXaxis()->SetRangeUser(0,150);
    
    eff[0][k]->SetMarkerSize(0.5);
    eff[0][k]->SetMarkerStyle(20);
    histo2[0][k]->SetMarkerSize(0.5);
    histo2[0][k]->SetMarkerStyle(20);
    if(k==0){
      eff[0][k]->SetMarkerColor(kBlack);
      eff[0][k]->SetLineColor(kBlack);
      histo2[0][k]->SetMarkerColor(kBlack);
      histo2[0][k]->SetLineColor(kBlack);
    }
    else{
      eff[0][k]->SetMarkerColor(colours.at(k%4));
      eff[0][k]->SetLineColor(colours.at(k%4));
      histo2[0][k]->SetMarkerColor(colours.at(k%4));
      histo2[0][k]->SetLineColor(colours.at(k%4));
    }
    if((k%4==1 && k!=1 ) || k==0 ){
      dummy_histo->Draw("E1");
      // histo2[0][k]->Draw("E1,SAME");
      eff[0][k]->Draw("E1,SAME");
      std::cout << "first drawn curve for " << k << std::endl;
    }
    else{
      eff[0][k]->Draw("E1,SAME");
      // histo2[0][k]->Draw("E1,SAME");
    }
    TLine *l1=new TLine(0, 1, 150, 1);
    l1->SetLineColor(kRed);
    l1->Draw("SAME");
    // histo2[1][k]->SetLineColor(2);
    // histo2[1][k]->SetMarkerColor(2);
    // histo2[1][k]->SetMarkerSize(0.5);
    // histo2[1][k]->SetMarkerStyle(20);
    // histo2[1][k]->Draw("E1,same");
    // eff[k]->Draw("same");
    DrawLatex(0.15,0.82,1,"2024 PbPb 387973",0.05);
    lego->AddEntry(histo2[0][k],triggers.at(k), "LPE");
    // lego->AddEntry(histo2[1][k],"2023 Zero Bias", "L");
    // lego->AddEntry(eff[k],"TEfficiency", "L");
    lego->Draw("");
    lego->SetFillColor(0);
    can->cd();
    // can->SaveAs(Form("PP_2023ZeroBias/fig_trig%i_pp_newJEC.pdf",k));
    if(k!=0 && k%4==0){
      can->Print( (out_file).c_str(),"pdf");
      // pad->Clear();
      lego->Clear();
      can->Clear();
    }
    // can->Clear();
    // pad->Clear();
  }
  c11->Print( (out_file + ")").c_str(),"pdf");
}
