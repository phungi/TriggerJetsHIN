


#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TDirectory.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>                                                                                                                                                                                                         
#include <cstdio>   // needed for io       
using namespace std;
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
#include <TEfficiency.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TMath.h>
#include <vector>
#include <TGraphAsymmErrors.h>

#endif


Float_t getDPHI(Float_t phi1, Float_t phi2){
  Float_t dphi = phi1 - phi2;

  if(dphi > TMath::Pi())
    dphi = dphi - 2. * TMath::Pi();
  if(dphi <= -TMath::Pi())
    dphi = dphi + 2. * TMath::Pi();

  if(TMath::Abs(dphi) > TMath::Pi()) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
  }

  return dphi;
}


Float_t getDR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2){
  Float_t theDphi = getDPHI(phi1, phi2);
  Float_t theDeta = eta1 - eta2;
  return TMath::Sqrt(theDphi*theDphi + theDeta*theDeta);
}

int Plot_Final_datacalo_JEC(std::string filename){


  //int bitmin=274;
  // int bitmax=281;
  // int binnum=7;
  int jec_correct=0;
  int bitmin=259;
  int bitmax=273;
  int binnum=14;

  double centmin=0;
  double centmax=60;
  const Double_t limits[14] = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  TH1F* hjet_all_calo;
  TH1F* hjet_all_calo_matched;
  TH1F* hjet_fired_calo[15];
  TH1F* hjet_fired_calo_ratio[15];
  TH1F* hcentrality;
  TH1F* hsumenergy;
  hcentrality=new TH1F("centrality","centrality",200,0,200);
  hsumenergy=new TH1F("sumenergy","sumenergy",500,0,20000);
  hjet_all_calo=new TH1F(Form("jetpt_all_calo"),Form("jetpt_all_calo"),150,0,300);
  hjet_all_calo_matched=new TH1F(Form("jetpt_all_calo_matched"),Form("jetpt_all_calo_matched"),150,0,300);
  hjet_all_calo->Sumw2();
  hjet_all_calo_matched->Sumw2();  
  for(int j=0;j<binnum;j++){

   
    hjet_fired_calo[j]=new TH1F(Form("jetpt_trig_calo%i",j),Form("jetpt_trig_calo%i",j),150,0,300);
    hjet_fired_calo[j]->Sumw2();}

  TFile* fout = new TFile("out.root","recreate"); 
 

  gSystem->Load("libFWCoreFWLite.so");
  JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/l/lcunquei/TestTrigger23/ParallelMC_L2Relative_AK4Calo_offline.txt");  
  vector<JetCorrectorParameters> vPar;
   vPar.push_back(*L2JetPar);

  FWLiteEnabler::enable(); 
  FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);
  cout<<JEC<<endl;

   
  TFile *input = TFile::Open(filename.data());
 
    TDirectory* emuDir = input->GetDirectory("l1UpgradeEmuTree");
    TDirectoryFile *hiEvtAnalyzer=(TDirectoryFile*)input->Get("hiEvtAnalyzer");
    TDirectoryFile *jetbranch=(TDirectoryFile*)input->Get("akCs4PFJetAnalyzer");



    TTree *jets=(TTree*)jetbranch->Get("t");
    TTree *ev=(TTree*)hiEvtAnalyzer->Get("HiTree");
    TTree *emuTree=(TTree*)emuDir->Get("L1UpgradeTree");


    emuTree->AddFriend(jets);
    emuTree->AddFriend(ev);


    TTreeReader emuReader(emuTree);

    TTreeReaderValue<Int_t> mycent(emuReader,"hiBin");
    TTreeReaderValue<Float_t> myhf(emuReader,"hiHF");  
    TTreeReaderValue<Int_t> nref(emuReader,"nref");
    TTreeReaderValue<Int_t> ncalo(emuReader,"ncalo");
    TTreeReaderArray<Float_t> jtpt(emuReader,"jtpt");
    TTreeReaderArray<Float_t> jteta(emuReader,"jteta");
    TTreeReaderArray<Float_t> jtphi(emuReader,"jtphi");
    TTreeReaderArray<Float_t> calopt(emuReader,"calopt");

    TTreeReaderArray<Float_t> caloeta(emuReader,"caloeta");
    TTreeReaderArray<Float_t> calophi(emuReader,"calophi");
    

    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
    TTreeReaderValue<vector<float>> emuJetEta(emuReader, "jetEta");
    TTreeReaderValue<vector<float>> emuJetPhi(emuReader, "jetPhi");

       
    if(!emuTree)return 1;
    if(!jets) return 1;



    while (emuReader.Next()) {

       hcentrality->Fill(*mycent);
      
      if(*mycent<centmin || *mycent>centmax) continue;
      hsumenergy->Fill(*myhf);
      //if(*myhf > 5000) continue;
      double calomax=-999;
      double etacalomax=0;
      double phicalomax=0;
     
      for(int l=0;l<*ncalo;l++){

	
	JEC->setJetPt(calopt[l]);
	JEC->setJetEta(caloeta[l]);
	JEC->setJetPhi(calophi[l]);

	double Correction = JEC->getCorrection();
      
	double CorrectedPT=calopt[l];
        if(jec_correct)CorrectedPT = calopt[l]*JEC->getCorrection(); 
	if(caloeta[l]>2 || caloeta[l]<-2) continue;
        if(CorrectedPT>calomax){

        calomax=CorrectedPT;
         etacalomax=caloeta[l];
         phicalomax=calophi[l];}
      }
     
      
        if(calomax==-999) continue;
      
     

      hjet_all_calo->Fill(calomax);
  
      double l1pt=0;
      double l1eta=0;
      double l1phi=0;
      double drmin=10;
      
     
     
      for(long unsigned n=0;n<(*emuJetPt).size();n++){
      double dr=getDR((*emuJetEta)[n],(*emuJetPhi)[n],etacalomax,phicalomax); 

        if(dr<drmin){
          drmin=dr;
          l1pt=(*emuJetPt)[n];
          l1eta=(*emuJetEta)[n];
          l1phi=(*emuJetPhi)[n];
        }
      }
    
     

      
      if(drmin>0.4) continue;
      if(l1pt==0) continue;
     
      hjet_all_calo_matched->Fill(calomax);


   
     
      for(int l=0;l<binnum;l++){

	if(l1pt>limits[l]){ 
	  hjet_fired_calo[l]->Fill(calomax);}
	  

      }
   
     
                                                               
   



             
                

        
    }



    
   
    delete jets;
    delete ev;
    delete emuTree;
    input->Close();
    delete input;

   
 
    fout->cd();
    hcentrality->Write();
    hsumenergy->Write();
  
    hjet_all_calo->Write();
    hjet_all_calo_matched->Divide(hjet_all_calo);
    hjet_all_calo_matched->Write();
    for(int k=bitmin;k<bitmax;k++){
    
      
      hjet_fired_calo[k-bitmin]->Write();
      hjet_fired_calo_ratio[k-bitmin]=(TH1F*)hjet_fired_calo[k-bitmin]->Clone(Form("hjet_fired_calo_ratio[%i]",k-bitmin));
      hjet_fired_calo_ratio[k-bitmin]->Write();
    }
   
       fout->Close();


    return 1;

}



