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




void plotResults_PbPb_data(int flag=0){                                                                                          
                                                                               
 
   TString fname="data_mb_corr.root";
   if(flag==1) fname="data_mb_uncorr.root";
     
         if(flag==2) fname="data_PbPb_pf_unmatched_corr.root";
	  if(flag==3) fname="data_2018_corr.root";  
	  if(flag==4) fname="data_2018_uncorr.root";
	  if(flag==5) fname="mc_18.root";
	  if(flag==6) fname="data_22_newjec.root";
            	   if(flag==7) fname="data_22_nonewjec.root";

         	   if(flag==8) fname="data_22_newjec_030.root";
		   	   if(flag==9) fname="data_22_newjec_3070.root";
			   	   if(flag==10) fname="data_22_newjec_3090.root";
                                     if(flag==11) fname="data_22_newjec_pf.root";
				      if(flag==12) fname="data_22_nonewjec_pf.root";
				        if(flag==13) fname="data_22_newjec_pf_matchedcalo.root";
					if(flag==14) fname="dataPbPb_22_newjec_pf_nonfhcut.root";
					if(flag==15) fname="dataPbPb_22_newjec_pf_deltaRmatch_periph.root";
					if(flag==16)fname="dataPbPb_22_newjec_pf_matchcalo_periph.root";

	   TFile *input;
   input=TFile::Open(fname);
   int nbins=14;
   TH1F *histo[15];
   TH1F *hjet_all;
    TH1F *hjet_closest;
     TH1F *hjet_hardest;   
     gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    hjet_all=(TH1F*)input->Get("jetpt_all_calo");
    //   histo[0]=(TH1F*)input->Get(Form("jetpt_trig_calo%i",0));
   for(int i=0;i<nbins;i++){
     
    histo[i]=(TH1F*)input->Get(Form("jetpt_trig_calo%i",i));
   
   
    histo[i]->SetLineColor(i+1);
   
    histo[i]->Divide(histo[i],hjet_all,1,1,"b");
    if(i==8) histo[i]->SetLineColor(kGray+1);
    if(i==9) histo[i]->SetLineColor(kRed-2);
    if(i==10) histo[i]->SetLineColor(kGreen-2);
     if(i==11) histo[i]->SetLineColor(kBlue-2);
     if(i==12) histo[i]->SetLineColor(kMagenta-2);
      if(i==13) histo[i]->SetLineColor(49);
     histo[i]->SetLineWidth(3);
     cout<<"what "<<histo[i]<<endl;
   }



   // if(flag<2)hjet_closest->Divide(hjet_hardest);



   
    TCanvas *can;                                                                                   
   TPad *pad;                                                                                      
   TLegend *lego;
  can= new TCanvas(Form("canvas1"),Form("canvas1") ,1100,1100);
  can->SetTicks();
  can->cd();
  pad = new TPad("pad0","This is pad0",0.,0.,1,1.);
  pad->SetFillColor(0);
  pad->SetMargin(0.15,0.9,0.25,0.9);
  pad->Draw();
  pad->SetTicks(1,1);
  pad->cd();
  histo[0]->GetYaxis()->SetTitleOffset(0.9);
  histo[0]->GetXaxis()->SetTitleOffset(0.9);
 
  histo[0]->GetXaxis()->SetLabelFont(42);
  histo[0]->GetYaxis()->SetLabelFont(42);
  histo[0]->GetXaxis()->SetLabelSize(0.04);
  histo[0]->GetYaxis()->SetLabelSize(0.04);

  histo[0]->GetXaxis()->SetTitleFont(42);
  histo[0]->GetYaxis()->SetTitleFont(42);
  histo[0]->GetXaxis()->SetTitleSize(0.065);
  histo[0]->GetYaxis()->SetTitleSize(0.065);
  histo[0]->SetMarkerSize(1);
  histo[0]->SetMarkerSize(1);
  histo[0]->GetYaxis()->SetRangeUser(0.1,1.1);
  histo[0]->GetXaxis()->SetRangeUser(0,140);
  if(flag==0) histo[0]->GetXaxis()->SetTitle("p_{T,jet}^{calo,corr}");
   if(flag==1 || flag==4 || flag==5 || flag==7) histo[0]->GetXaxis()->SetTitle("p_{T,jet}^{calo,uncorr}");
    if(flag==2 || flag==11 || flag==13 || flag==14 || flag==15) histo[0]->GetXaxis()->SetTitle("p_{T,jet}^{pf,corr}");
     if(flag==12) histo[0]->GetXaxis()->SetTitle("p_{T,jet}^{pf,uncorr}");
    if(flag==3 || flag==6 ||flag==8||flag==9||flag==10) histo[0]->GetXaxis()->SetTitle("p_{T,jet}^{calo,corr}");  
    histo[0]->GetYaxis()->SetTitle("Trigger efficiency");
 
  histo[0]->Draw("");

  for(int n=1;n<nbins;n++){
   
    histo[n]->Draw("same");
  }

  lego = new TLegend(0.5, 0.3, 0.7, 0.5);
  lego->SetBorderSize(0);
  lego->SetTextSize(0.015);
  lego->SetTextFont(42);



  lego->AddEntry(histo[0],"L1_SinglJet8_BptxAND", "L");
   lego->AddEntry(histo[1],"L1_SinglJet16_BptxAND", "L");
    lego->AddEntry(histo[2],"L1_SinglJet24_BptxAND", "L");
     lego->AddEntry(histo[3],"L1_SinglJet28_BptxAND", "L");
      lego->AddEntry(histo[4],"L1_SinglJet32_BptxAND", "L");
       lego->AddEntry(histo[5],"L1_SinglJet36_BptxAND", "L");
        lego->AddEntry(histo[6],"L1_SinglJet40_BptxAND", "L");
	if(nbins>7){
	lego->AddEntry(histo[7],"L1_SinglJet44_BptxAND", "L");
	  lego->AddEntry(histo[8],"L1_SinglJet48_BptxAND", "L");
	   lego->AddEntry(histo[9],"L1_SinglJet56_BptxAND", "L");
	    lego->AddEntry(histo[10],"L1_SinglJet60_BptxAND", "L");
	     lego->AddEntry(histo[11],"L1_SinglJet64_BptxAND", "L");
                lego->AddEntry(histo[12],"L1_SinglJet72_BptxAND", "L");
		lego->AddEntry(histo[13],"L1_SinglJet80_BptxAND", "L");}
  lego->Draw("");
  lego->SetFillColor(0);

  if(flag==1) can->SaveAs(Form("L1TriggerEfficiency_data_calo_old.pdf"));
  if(flag==0) can->SaveAs(Form("L1TriggerEfficiency_data_calo_test.pdf"));    








 //  if(flag<2){




 //  TCanvas *canb;                                                                                                                                                              
 //   TPad *padb;                                                                                                                                                                  
 //   TLegend *legob;                                                                                                                                                              
 //  canb= new TCanvas(Form("canvasb"),Form("canvasb") ,1100,1100);                                                                                                                
 //  canb->SetTicks();                                                                                                                                                             
 //  canb->cd();                                                                                                                                                                   
 //  padb = new TPad("padb","This is padb",0.,0.,1,1.);                                                                                                                            
 //  padb->SetFillColor(0);                                                                                                                                                        
 //  padb->SetMargin(0.15,0.9,0.25,0.9);                                                                                                                                           
 //  padb->Draw();                                                                                                                                                                 
 //  padb->SetTicks(1,1);                                                                                                                                                          
 //  padb->cd();                                                                                                                                          hjet_closest->GetYaxis()->SetTitle("Efficiency hardest PF");
 //  hjet_closest->GetXaxis()->SetTitle("p_{T,jet}^{PF}");
 // hjet_closest->GetYaxis()->SetTitleOffset(0.9);                                                                                                                                   
 // hjet_closest->GetXaxis()->SetTitleOffset(0.9);                                                                                                                                   
 // hjet_closest->Draw("");



 //  }


  
}
