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

int Data_pp_pflow_JEC_pflow_first(std::string filename){
  TH1::SetDefaultSumw2();
  //int bitmin=274;
  //int bitmax=281;
  //int binnum=7;
  int jec_correct=1;
  int bitmin=259;
  int bitmax=273;
  // int binnum=14;

  double centmin=60;
  double centmax=180;
  std::vector<Double_t> limits = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  TH1F* hjet_all_pf;
  TH1F* hjet_all_pf_matched;
  TH1F* hjet_fired_calo[15];
  TH1F* hjet_fired_calo_ratio[15];
  TH1F* hcentrality;
  TH1F* hjetl1pfdR;
  TH2F *dR_vs_PFpt = new TH2F("dR_vs_PFpt", "#DeltaR vs PFlow pt",50,0,5, 100, 0, 300);
  hcentrality=new TH1F("centrality","centrality",200,0,200);
  hjet_all_pf=new TH1F(Form("jetpt_all_pf"),Form("jetpt_all_pf"),150,0,300);
  hjet_all_pf_matched=new TH1F(Form("jetpt_all_pf_matched"),Form("jetpt_all_pf_matched"),150,0,300);
  hjetl1pfdR=new TH1F(Form("dR_PF_to_L1"),Form("#DeltaR between hardest PF and closest L1"),150,0,3.16);
  for(int j=0;j<limits.size();j++){
    hjet_fired_calo[j]=new TH1F(Form("jetpt_trig_calo%i",j),Form("jetpt_trig_calo%i",j),150,0,300);
    // hjet_fired_calo[j]->Sumw2();
  }
  TFile* fout = new TFile("out.root","recreate"); 
  gSystem->Load("libFWCoreFWLite.so");
  JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L2JetPar);

  FWLiteEnabler::enable(); 
  FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);
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
  if(!emuTree) return 1;
  if(!jets) return 1;

  Int_t term_count=0;
  while (emuReader.Next()) {
    // std::cout << "pf jets: " << *nref << " l1 jets: " << (*emuJetPt).size() << " calo jets: " << *ncalo << std::endl;
      // term_count++;
      // if(term_count==2000) break;
    // emuDecisionVec = emuDecision.Get();      
    // hcentrality->Fill(*mycent);
    // if(*mycent<centmin || *mycent>centmax) continue;
    double L1_match_pt= -std::numeric_limits<double>::max();
    double L1_match_eta=-std::numeric_limits<double>::max();
    double L1_match_phi=-std::numeric_limits<double>::max();
    double calomaxb=-std::numeric_limits<double>::max();
    double etacalomaxb=-std::numeric_limits<double>::max();
    double phicalomaxb=-std::numeric_limits<double>::max();  
    double PFpt=-std::numeric_limits<double>::max();
    double PFeta=-std::numeric_limits<double>::max();
    double PFphi=-std::numeric_limits<double>::max();
    double drmin=std::numeric_limits<double>::max();
    for(int l=0;l<*nref;l++){
      //Below cuts are taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
      if(jteta[l]>2 || jteta[l]<-2) continue;
      Bool_t passSelections = false;
      Int_t totalParticles = CHM[l] + NHM[l] + CEM[l] + NEM[l] + MUM[l];
      Int_t NumNeutralParticles = NHM[l] + NEM[l];
      if( abs(jteta[l])<=2.6 && CEMF[l]<0.8 && CHM[l]>0 && CHF[l]>0.01 && NEMF[l]<0.9 && MUF[l] <0.8 && NHF[l] < 0.9 && totalParticles > 1 ) passSelections = true;
      //below cuts are removed since the selection is for <2 eta anyway
      // if( abs(jteta[l])>2.6 && abs(jteta[l])<=2.7 && CEMF[l]<0.8 && CHM[l]>0 && NEMF[l]<0.99 && MUF[l] <0.8 && NHF[l] < 0.9 ) passSelections = true;
      // if( abs(jteta[l])>2.7 && abs(jteta[l])<=3.0 && NEMF[l]<0.99 && NumNeutralParticles > 1 ) passSelections = true;
      // if( abs(jteta[l])>3.0 && NEMF[l]<0.90 && NumNeutralParticles > 10 ) passSelections = true;
      //if jet doesn't satisfy above conditions, skip
      if(!passSelections) continue;
      JEC->setJetPt(jtpt[l]);
      JEC->setJetEta(jteta[l]);
      double Correction = JEC->getCorrection();
      double CorrectedPT = jtpt[l];
      if(jec_correct) CorrectedPT = jtpt[l]*JEC->getCorrection();
      
      if(CorrectedPT > PFpt){
        PFpt = CorrectedPT;
        PFeta = jteta[l];
        PFphi = jtphi[l];
      }
    }
    if(PFpt < 0) continue;
    for(long unsigned n=0;n<(*emuJetPt).size();n++){
      double dr=getDR(PFeta,PFphi,emuJetEta->at(n),emuJetPhi->at(n)); 
      if(dr < drmin){
        L1_match_pt=emuJetPt->at(n);
        L1_match_eta=emuJetEta->at(n);
        L1_match_phi=emuJetPhi->at(n);
        drmin = dr;
      }
    }
    hjetl1pfdR->Fill(drmin);
    hjet_all_pf->Fill(PFpt);
    dR_vs_PFpt->Fill(drmin, PFpt);
    if(L1_match_pt==-std::numeric_limits<double>::max()) continue;
        // double drjets=getDR(l1eta,l1phi,PF_match_eta,PF_match_phi);  
    if(drmin > 0.4) continue;
    hjet_all_pf_matched->Fill(PFpt);
    for(int l=0;l<limits.size();l++){
      if(L1_match_pt>=limits.at(l)){ 
        hjet_fired_calo[l]->Fill(PFpt);
      }
    }
  }
  delete jets;
  delete ev;
  delete emuTree;
  input->Close();
  delete input;
  fout->cd();
  hjet_all_pf->Write();
  hjet_all_pf_matched->Divide(hjet_all_pf);
  hjet_all_pf_matched->Write();
  hjetl1pfdR->Write();
  dR_vs_PFpt->Write();
  for(int k=bitmin;k<bitmax;k++){
    hjet_fired_calo[k-bitmin]->Write();
    hjet_fired_calo_ratio[k-bitmin]=(TH1F*)hjet_fired_calo[k-bitmin]->Clone(Form("hjet_fired_calo_ratio[%i]",k-bitmin));
    hjet_fired_calo_ratio[k-bitmin]->Write();
  }
  fout->Close();
  return 0;
}



