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


void plotResults(){
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
  TH1::SetDefaultSumw2( );
  gROOT->SetBatch(kTRUE);  
  std::string out_file = "pp_2023_ZB_cuts_JEC.pdf";
  std::vector<TString> triggers = {"L1_SingleJet8_BptxAND","L1_SingleJet16_BptxAND","L1_SingleJet20_BptxAND","L1_SingleJet24_BptxAND","L1_SingleJet28_BptxAND","L1_SingleJet32_BptxAND","L1_SingleJet35_BptxAND","L1_SingleJet40_BptxAND","L1_SingleJet44_BptxAND","L1_SingleJet48_BptxAND","L1_SingleJet50_BptxAND","L1_SingleJet56_BptxAND","L1_SingleJet60_BptxAND","L1_SingleJet80_BptxAND","L1_SingleJet90_BptxAND","L1_SingleJet120_BptxAND","L1_SingleJet140_BptxAND","L1_SingleJet150_BptxAND","L1_SingleJet160_BptxAND","L1_SingleJet170_BptxAND","L1_SingleJet180_BptxAND","L1_SingleJet200_BptxAND"};
  TCanvas *c = new TCanvas();
  c->Print( (out_file + "(").c_str(),"pdf");
  std::vector<TString> filenames ={"pp2023ZB_cuts_JEC.root"};
  std::vector<TFile*> files = {};
  for(size_t i{0}; i < filenames.size(); ++i){
    TFile *a = new TFile(filenames.at(i),"READ");
    files.push_back(a);
  }
// TFile *input1 = new TFile(fname1, "READ");
  // TFile *input2 = new TFile("2023_pp_ZeroBias.root","READ");
// TFile *input3;
  std::vector<TH1F*> histo;
  std::vector<TH1F*> histo2;
  TH1F *hjet_all;
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
// hjet_all[0]=(TH1F*)input1->Get("jetpt_all_calo");
  hjet_all=(TH1F*)input2->Get("jetpt_all_calo");
  hjet_all[0]->GetXaxis()->SetRangeUser(1, 20);
  for(int i=0;i<nbins;i++){
// histo[0][i]=(TH1F*)input1->Get(Form("jetpt_trig_calo%i",i));
    histo[0][i]=(TH1F*)input2->Get(Form("jetpt_trig_calo%i",i));
    histo[0][i]->GetXaxis()->SetRangeUser(1, 20);
  }

  for(int j=0;j<1;j++){
    for(int i=0;i<nbins;i++){
      
      
      histo2[j][i]=(TH1F*)histo[j][i]->Clone(Form("histo2[%i][%i]",j,i));
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
  lego = new TLegend(0.65, 0.25, 0.75, 0.55);
  lego->SetBorderSize(0);
  lego->SetTextSize(0.03);
  lego->SetTextFont(42);
  lego->SetBorderSize(0);
  lego->SetFillStyle(0);
  can->SetTicks();
// can[k]->cd();
  // pad->SetFillColor(0);
  // pad->SetMargin(0.15,0.9,0.25,0.9);
  // pad->Draw();
  // pad->SetTicks(1,1);
  // pad->cd();
  for(int k=0;k<nbins;k++){
  // can->Clear();
  // TCanvas *can[k];                                                                                   
  // TPad *pad;                                                                                      
  // TLegend *lego;
  // can[k]= new TCanvas(Form("canvas1%d",k),Form("canvas1%d",k) ,1100,1100);
  // can[k]->SetTicks();
  // can[k]->cd();
  // pad = new TPad("pad0","This is pad0",0.,0.,1,1.);
  // pad->SetFillColor(0);
  // pad->SetMargin(0.15,0.9,0.25,0.9);
  // pad->Draw();
  // pad->SetTicks(1,1);
    can->cd();
    histo2[0][k]->GetYaxis()->SetTitleOffset(0.8);
    histo2[0][k]->GetXaxis()->SetTitleOffset(0.95);
    histo2[0][k]->GetYaxis()->SetTitle("Trigger efficiency");
    histo2[0][k]->GetXaxis()->SetTitle("p_{T,pf}[GeV]");
    histo2[0][k]->GetXaxis()->SetLabelFont(42);
    histo2[0][k]->GetYaxis()->SetLabelFont(42);
    histo2[0][k]->GetXaxis()->SetLabelSize(0.04);
    histo2[0][k]->GetYaxis()->SetLabelSize(0.04);

    histo2[0][k]->GetXaxis()->SetTitleFont(42);
    histo2[0][k]->GetYaxis()->SetTitleFont(42);
    histo2[0][k]->GetXaxis()->SetTitleSize(0.045);
    histo2[0][k]->GetYaxis()->SetTitleSize(0.045);
    histo2[0][k]->GetYaxis()->SetRangeUser(0,1.25);
    histo2[0][k]->GetXaxis()->SetRangeUser(0,100);
    histo2[0][k]->SetMarkerSize(0.5);
    histo2[0][k]->SetMarkerStyle(20);
    if(k==0){
      histo2[0][k]->SetMarkerColor(colours.at(0));
      histo2[0][k]->SetLineColor(colours.at(0));
    }
    else{
      histo2[0][k]->SetMarkerColor(colours.at(k%4));
      histo2[0][k]->SetLineColor(colours.at(k%4));
    }
    if(k%4==1 || k==0 )
      histo2[0][k]->Draw("E1");
    else
      histo2[0][k]->Draw("E1,SAME");
    TLine *l1=new TLine(0, 1, 100, 1);
    l1->SetLineColor(kRed);
    l1->Draw("SAME");
// histo2[1][k]->SetLineColor(2);
// histo2[1][k]->SetMarkerColor(2);
// histo2[1][k]->SetMarkerSize(0.5);
// histo2[1][k]->SetMarkerStyle(20);
// histo2[1][k]->Draw("E1,same");
// eff[k]->Draw("same");
    DrawLatex(0.15,0.82,1,"2023 pp ZeroBias",0.05);



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
  c->Print( (out_file + ")").c_str(),"pdf");
}
