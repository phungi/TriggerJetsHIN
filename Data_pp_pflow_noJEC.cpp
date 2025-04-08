


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

int Data_pp_pflow_noJEC(std::string filename){


//int bitmin=274;
//int bitmax=281;
//int binnum=7;
  int jec_correct=0;
  int bitmin=259;
  int bitmax=273;
  int binnum=14;

  double centmin=60;
  double centmax=180;
  const Double_t limits[14] = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  TH1F* hjet_all_calo;
  TH1F* hjet_all_calo_matched;
  TH1F* hjet_fired_calo[15];
  TH1F* hjet_fired_calo_ratio[15];
  TH1F* hcentrality;
  TH1F* hjetcalopf;

  hcentrality=new TH1F("centrality","centrality",200,0,200);
  hjet_all_calo=new TH1F(Form("jetpt_all_calo"),Form("jetpt_all_calo"),150,0,300);
  hjet_all_calo_matched=new TH1F(Form("jetpt_all_calo_matched"),Form("jetpt_all_calo_matched"),150,0,300);
  hjetcalopf=new TH1F(Form("jetcalopf"),Form("jetcalopf"),150,0,3.16);
  hjet_all_calo->Sumw2();
  hjet_all_calo_matched->Sumw2();  
  hjetcalopf->Sumw2();
  for(int j=0;j<binnum;j++){


    hjet_fired_calo[j]=new TH1F(Form("jetpt_trig_calo%i",j),Form("jetpt_trig_calo%i",j),150,0,300);
    hjet_fired_calo[j]->Sumw2();}

    TFile* fout = new TFile("out.root","recreate"); 


    gSystem->Load("libFWCoreFWLite.so");
    JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
    vector<JetCorrectorParameters> vPar;
    vPar.push_back(*L2JetPar);

    FWLiteEnabler::enable(); 
    FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);
    cout<<JEC<<endl;
//cout<<indexdir<<" "<<indexfile<<indexrun1<<" "<<indexrun2<<endl;  
// TFile *input = TFile::Open(Form("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/mitaylor/L1MenuStudies/TestRun2022/MinimumBias/HITestRaw%d_MinBias_HIRun2022A_1260pre1_v2/HITestRaw%d/MinBias_HIRun2022A_1260pre1_v2/230%d_%d/0000/L1Ntuple_%d.root",indexdir,indexdir,indexrun1,indexrun2,indexfile));
//TFile *input = TFile::Open("L1Ntuple_86.root");   
//if(!input || input->IsZombie()) { cout <<"The file could not be opened!"; return 1;}

    TFile *input = TFile::Open(filename.data());
    // std::cout << "Done with the file" << std::endl;
    TDirectory* emuDir = input->GetDirectory("l1UpgradeTree");
    TDirectoryFile *hiEvtAnalyzer=(TDirectoryFile*)input->Get("hiEvtAnalyzer");
// TDirectoryFile *jetbranch=(TDirectoryFile*)input->Get("akCs4PFJetAnalyzer");
    TDirectoryFile *jetbranch=(TDirectoryFile*)input->Get("ak4PFJetAnalyzer");
// TDirectory* emuDirdec = input->GetDirectory("l1uGTEmuTree");


    TTree *jets=(TTree*)jetbranch->Get("t");
    TTree *ev=(TTree*)hiEvtAnalyzer->Get("HiTree");
    TTree *emuTree=(TTree*)emuDir->Get("L1UpgradeTree");
//TTree *emdec=(TTree*)emuDirdec->Get("L1uGTTree"); 

    
    emuTree->AddFriend(jets);
    emuTree->AddFriend(ev);
//emuTree->AddFriend(emdec);
    // std::cout << "Done with the trees" << std::endl;
    TTreeReader emuReader(emuTree);
// TTreeReaderValue<std::vector<bool>> emuDecision(emuReader, "m_algoDecisionInitial");
// TTreeReaderValue<Int_t> mycent(emuReader,"hiBin");
    TTreeReaderArray<Float_t> NHF(emuReader, "jtPfNHF");
    TTreeReaderArray<Float_t> CEMF(emuReader, "jtPfCEF");
    TTreeReaderArray<Float_t> CHF(emuReader, "jtPfCHF");
    TTreeReaderArray<Float_t> NEMF(emuReader, "jtPfNEF");
    TTreeReaderArray<Float_t> MUF(emuReader, "jtPfMUF");

    TTreeReaderArray<Int_t>   CHM(emuReader, "jtPfCHM");
    TTreeReaderArray<Int_t>   NHM(emuReader, "jtPfNHM");
    TTreeReaderArray<Int_t>   CEM(emuReader, "jtPfCEM");
    TTreeReaderArray<Int_t>   NEM(emuReader, "jtPfNEM");
    TTreeReaderArray<Int_t>   MUM(emuReader, "jtPfMUM");

    TTreeReaderValue<Int_t> nref(emuReader,"nref");
    TTreeReaderArray<Float_t> jtpt(emuReader,"jtpt");
    TTreeReaderArray<Float_t> jteta(emuReader,"jteta");
    TTreeReaderArray<Float_t> jtphi(emuReader,"jtphi");

    TTreeReaderValue<Int_t> ncalo(emuReader,"ncalo");
    TTreeReaderArray<Float_t> calopt(emuReader,"calopt");
    TTreeReaderArray<Float_t> caloeta(emuReader,"caloeta");
    TTreeReaderArray<Float_t> calophi(emuReader,"calophi");


    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
    TTreeReaderValue<vector<float>> emuJetEta(emuReader, "jetEta");
    TTreeReaderValue<vector<float>> emuJetPhi(emuReader, "jetPhi");

    // std::cout << "Done with the arrays" << std::endl;
    if(!emuTree)return 1;
    if(!jets) return 1;

//vector<bool>* emuDecisionVec;

//cout<<"hello"<<endl;   

    while (emuReader.Next()) {
// emuDecisionVec = emuDecision.Get();      
// hcentrality->Fill(*mycent);
// if(*mycent<centmin || *mycent>centmax) continue;
      double calomax=-999;
      double etacalomax=0;
      double phicalomax=0;
      double calomaxb=-999;
      double etacalomaxb=0;
      double phicalomaxb=0;  

      for(int l=0;l<*nref;l++){
//cout<<jtPfMUF[l]<<endl;
//if(jtPfNHF[l]==0) continue;
        // std::cout << "Setting JEC vars " << jtpt[l] << " " << jteta[l] << " " << jtphi[l] << std::endl;
        //Below cuts are taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
        Bool_t passSelections = false;
        Int_t totalParticles = CHM[l] + NHM[l] + CEM[l] + NEM[l] + MUM[l];
        Int_t NumNeutralParticles = NHM[l] + NEM[l];
        if( abs(jteta[l])<=2.6 && CEMF[l]<0.8 && CHM[l]>0 && CHF[l]>0.01 && NEMF[l]<0.9 && MUF[l] <0.8 && NHF[l] < 0.9 && totalParticles > 1 ) passSelections = true;
        // if( abs(jteta[l])>2.6 && abs(jteta[l])<=2.7 && CEMF[l]<0.8 && CHM[l]>0 && NEMF[l]<0.99 && MUF[l] <0.8 && NHF[l] < 0.9 ) passSelections = true;
        // if( abs(jteta[l])>2.7 && abs(jteta[l])<=3.0 && NEMF[l]<0.99 && NumNeutralParticles > 1 ) passSelections = true;
        // if( abs(jteta[l])>3.0 && NEMF[l]<0.90 && NumNeutralParticles > 10 ) passSelections = true;
        //if jet doesn't satisfy above conditions, skip
        if(!passSelections) continue;
        JEC->setJetPt(jtpt[l]);
        JEC->setJetEta(jteta[l]);
        JEC->setJetPhi(jtphi[l]);

        double Correction = JEC->getCorrection();
        // std::cout << Correction << " jet correction" << std::endl;
        double CorrectedPT = jtpt[l];
        if(jec_correct) CorrectedPT = jtpt[l]*JEC->getCorrection();
        if(jteta[l]>2 || jteta[l]<-2) continue;
        if(CorrectedPT>calomax){
          calomax=CorrectedPT;
          etacalomax=jteta[l];
          phicalomax=jtphi[l];}
        }


        if(calomax==-999) continue;



        for(int m=0;m<*ncalo;m++){
          if(caloeta[m]>2 || caloeta[m]<-2) continue;
          double CorrectedCaloPT=calopt[m];
          if(CorrectedCaloPT>calomaxb){

            calomaxb=CorrectedCaloPT;
            etacalomaxb=caloeta[m];
            phicalomaxb=calophi[m];}

          }

          double drjets=getDR(etacalomaxb,phicalomaxb,etacalomax,phicalomax);
          hjetcalopf->Fill(drjets);

          if(drjets>0.4) continue; 

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
//      cout<<drmin<<" "<<l1pt<<endl;



//      if(drmin>0.4) continue;
          if(l1pt==0) continue;

          hjet_all_calo_matched->Fill(calomax);



//cout<<binnum<<endl;
          for(int l=0;l<binnum;l++){

            if(l1pt>limits[l]){ 
              hjet_fired_calo[l]->Fill(calomax);}


            }

//cout<<"emuloop"<<endl;









          }




//cout<<"here"<<endl;


//    delete emuDir;
//delete hiEvtAnalyzer;
// delete jetbranch;
          delete jets;
          delete ev;
          delete emuTree;
          input->Close();
          delete input;

//cout<<"there"<<endl;

          fout->cd();
//cout<<"and here"<<endl;
          hjet_all_calo->Write();
          hjet_all_calo_matched->Divide(hjet_all_calo);
          hjet_all_calo_matched->Write();
          hjetcalopf->Write();
          for(int k=bitmin;k<bitmax;k++){


            hjet_fired_calo[k-bitmin]->Write();
            hjet_fired_calo_ratio[k-bitmin]=(TH1F*)hjet_fired_calo[k-bitmin]->Clone(Form("hjet_fired_calo_ratio[%i]",k-bitmin));
            hjet_fired_calo_ratio[k-bitmin]->Write();
          }

          fout->Close();


          return 0;

        }



