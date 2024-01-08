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


void plotResults_pp_summaryJec(){                                                                                          
                                                                               
 
  // TString fname1="data_23_ppring_newjec.root";
  // TString fname2="data_23_ppring_nonewjec.root";

  TString fname1="datapp_22_newjec_pf.root";
  TString fname2="datapp_22_nonewjec_pf.root";  
   TFile *input1;
   input1=TFile::Open(fname1);
   TFile *input2;                                                                                                      input2=TFile::Open(fname2);
   TFile *input3;                                          
   

   int nbins=22;
   TH1F *histo[2][22];
   TH1F *hjet_all[2]; 
   TH1F *histo2[2][22];
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    hjet_all[0]=(TH1F*)input1->Get("jetpt_all_calo");
    hjet_all[1]=(TH1F*)input2->Get("jetpt_all_calo");

    
    
    for(int i=0;i<nbins;i++){
    histo[0][i]=(TH1F*)input1->Get(Form("jetpt_trig_calo%i",i));
    histo[1][i]=(TH1F*)input2->Get(Form("jetpt_trig_calo%i",i));
     }

    for(int j=0;j<2;j++){
      for(int i=0;i<nbins;i++){
     histo2[j][i]=(TH1F*)histo[j][i]->Clone(Form("histo2[%j][%i]",j,i));
     histo2[j][i]->Divide(histo[j][i],hjet_all[j],1,1,"b");
     histo2[j][i]->SetLineWidth(3);
     cout<<"what "<<histo2[j][i]<<endl;}}

   
    for(int k=0;k<nbins;k++){
    TCanvas *can[k];                                                                                   
    TPad *pad;                                                                                      
    TLegend *lego;
    can[k]= new TCanvas(Form("canvas1%d",k),Form("canvas1%d",k) ,1100,1100);
   can[k]->SetTicks();
   can[k]->cd();
   pad = new TPad("pad0","This is pad0",0.,0.,1,1.);
   pad->SetFillColor(0);
  pad->SetMargin(0.15,0.9,0.25,0.9);
   pad->Draw();
   pad->SetTicks(1,1);
   pad->cd();
   histo2[0][k]->GetYaxis()->SetTitleOffset(0.9);
  histo2[0][k]->GetXaxis()->SetTitleOffset(0.9);
  histo2[0][k]->GetYaxis()->SetTitle("Trigger efficiency");
   histo2[0][k]->GetXaxis()->SetTitle("p_{T,pf}");
   histo2[0][k]->GetXaxis()->SetLabelFont(42);
   histo2[0][k]->GetYaxis()->SetLabelFont(42);
   histo2[0][k]->GetXaxis()->SetLabelSize(0.04);
   histo2[0][k]->GetYaxis()->SetLabelSize(0.04);

   histo2[0][k]->GetXaxis()->SetTitleFont(42);
   histo2[0][k]->GetYaxis()->SetTitleFont(42);
   histo2[0][k]->GetXaxis()->SetTitleSize(0.065);
   histo2[0][k]->GetYaxis()->SetTitleSize(0.065);
   histo2[0][k]->SetMarkerSize(1);
   histo2[0][k]->SetMarkerSize(1);
   histo2[0][k]->GetYaxis()->SetRangeUser(0,1.25);
   histo2[0][k]->GetXaxis()->SetRangeUser(0,400);
   
 
   histo2[0][k]->SetLineColor(1);
   histo2[0][k]->Draw("");
   histo2[1][k]->SetLineColor(2);
   histo2[1][k]->Draw("same");
   



   if(k==0) DrawLatex(0.2,0.85,1,"L1_SinglJet8_BptxAND",0.05);
   if(k==1) DrawLatex(0.2,0.85,1,"L1_SinglJet16_BptxAND",0.05);
   if(k==2) DrawLatex(0.2,0.85,1,"L1_SinglJet20_BptxAND",0.05);
   if(k==3) DrawLatex(0.2,0.85,1,"L1_SinglJet24_BptxAND",0.05);
   if(k==4) DrawLatex(0.2,0.85,1,"L1_SinglJet28_BptxAND",0.05);
   if(k==5) DrawLatex(0.2,0.85,1,"L1_SinglJet32_BptxAND",0.05);
   if(k==6) DrawLatex(0.2,0.85,1,"L1_SinglJet35_BptxAND",0.05);
   if(k==7) DrawLatex(0.2,0.85,1,"L1_SinglJet40_BptxAND",0.05);
   if(k==8) DrawLatex(0.2,0.85,1,"L1_SinglJet44_BptxAND",0.05);
   if(k==9) DrawLatex(0.2,0.85,1,"L1_SinglJet48_BptxAND",0.05);
   if(k==10) DrawLatex(0.2,0.85,1,"L1_SinglJet50_BptxAND",0.05);
   if(k==11) DrawLatex(0.2,0.85,1,"L1_SinglJet56_BptxAND",0.05);
   if(k==12) DrawLatex(0.2,0.85,1,"L1_SinglJet60_BptxAND",0.05);
   if(k==13) DrawLatex(0.2,0.85,1,"L1_SinglJet80_BptxAND",0.05);
   if(k==14) DrawLatex(0.2,0.85,1,"L1_SinglJet90_BptxAND",0.05);
   if(k==15) DrawLatex(0.2,0.85,1,"L1_SinglJet120_BptxAND",0.05);
   if(k==16) DrawLatex(0.2,0.85,1,"L1_SinglJet140_BptxAND",0.05);
   if(k==17) DrawLatex(0.2,0.85,1,"L1_SinglJet150_BptxAND",0.05);
   if(k==18) DrawLatex(0.2,0.85,1,"L1_SinglJet160_BptxAND",0.05);
   if(k==19) DrawLatex(0.2,0.85,1,"L1_SinglJet170_BptxAND",0.05);
   if(k==20) DrawLatex(0.2,0.85,1,"L1_SinglJet180_BptxAND",0.05);
   if(k==21) DrawLatex(0.2,0.85,1,"L1_SinglJet200_BptxAND",0.05);
   


   
   
  lego = new TLegend(0.75, 0.25, 0.85, 0.45);
  lego->SetBorderSize(0);
  lego->SetTextSize(0.015);
  lego->SetTextFont(42);


   lego->AddEntry(histo2[0][k],"corrected", "L");
   lego->AddEntry(histo2[1][k],"uncorrected", "L");
   

	lego->Draw("");
  lego->SetFillColor(0);
  can[k]->SaveAs(Form("fig_trig%i_pp_newJEC.pdf",k));
  
    }
}
