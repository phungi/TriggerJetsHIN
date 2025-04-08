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

void gammagamma_jets_pf(TString input){
  TH1::SetDefaultSumw2();
  //int bitmin=274;
  //int bitmax=281;
  //int binnum=7;
  int jec_correct=1;
  // int bitmin=259;
  // int bitmax=273;
  // int binnum=14;

  double centmin=60;
  double centmax=180;
  std::vector<Double_t> threshold = {0.5,1,2,4,6,8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  std::vector<TH1F> emuMatchedHistos;
  std::vector<TH1F> emuHistos;

  int nbins = 75;
  float min = 0;
  float max = 150;

  TH1F* hjet_all_pf;
  TH1F* hjet_all_pf_matched;
  TH1F* hjet_fired_calo[15];
  TH1F* hjet_fired_calo_ratio[15];
  TH1F* hcentrality;
  TH1F* hjetl1pfdR;
  TH2F *dR_vs_PFpt = new TH2F("dR_vs_PFpt", "#DeltaR vs PFlow pt",50,0,5, 100, 0, 300);
  hcentrality=new TH1F("centrality","centrality",200,0,200);
  hjet_all_pf=new TH1F(Form("jetpt_all_pf"),Form("jetpt_all_pf"), nbins, min, max);
  hjet_all_pf_matched=new TH1F(Form("jetpt_all_pf_matched"),Form("jetpt_all_pf_matched"), nbins, min, max);
  hjetl1pfdR=new TH1F(Form("dR_PF_to_L1"),Form("#DeltaR between hardest PF and closest L1"),150,0,3.16);
  for(Int_t i{0}; i<threshold.size(); ++i){
        TH1F a(Form("jetpt_trig_calo_probably_dont_use%i",i), "", nbins, min, max);
        TH1F b(Form("jetpt_trig_calo%i",i), "", nbins, min, max);
        emuMatchedHistos.push_back(b);
        emuHistos.push_back(a);

    }
  // for(int j=0;j<threshold.size();j++){
  //   hjet_fired_calo[j]=new TH1F(Form("jetpt_trig_calo%i",j),Form("jetpt_trig_calo%i",j),150,0,300);
  //   // hjet_fired_calo[j]->Sumw2();
  // }
  TFile* fout = new TFile(input + "_doubleL1triggers_histos.root","recreate"); 
  gSystem->Load("libFWCoreFWLite.so");
  JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L2JetPar);

  FWLiteEnabler::enable(); 
  FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);
  TFile *in_f = new TFile(input, "READ");
  // std::cout << "Done with the file" << std::endl;
  TDirectory* emuDir = in_f->GetDirectory("l1UpgradeEmuTree");
  TDirectoryFile *hiEvtAnalyzer=(TDirectoryFile*)in_f->Get("hiEvtAnalyzer");
  // TDirectoryFile *jetbranch=(TDirectoryFile*)in_f->Get("akCs4PFJetAnalyzer");
  TDirectoryFile *jetbranch=(TDirectoryFile*)in_f->Get("akFlowPuCs4PFJetAnalyzer");
  // TDirectoryFile *jetbranch=(TDirectoryFile*)in_f->Get("ak4PFJetAnalyzer");
  // TDirectory* emuDirdec = in_f->GetDirectory("l1uGTEmuTree");

  TTree *jets=(TTree*)jetbranch->Get("t");
  TTree *ev=(TTree*)hiEvtAnalyzer->Get("HiTree");
  TTree *emuTree=(TTree*)emuDir->Get("L1UpgradeTree");
  //TTree *emdec=(TTree*)emuDirdec->Get("L1uGTTree"); 


  emuTree->AddFriend(jets);
  emuTree->AddFriend(ev);
  //emuTree->AddFriend(emdec);
  // std::cout << "Done with the trees" << std::endl;
  TTreeReader emuReader(emuTree);
  TTreeReaderValue<Float_t> weight(emuReader, "weight");
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
  if(!emuTree) return;
  if(!jets) return;

  Int_t term_count=0;
  while (emuReader.Next()) {
    // std::cout << "pf jets: " << *nref << " l1 jets: " << (*emuJetPt).size() << " calo jets: " << *ncalo << std::endl;
      // term_count++;
      // if(term_count==2000) break;
    // emuDecisionVec = emuDecision.Get();      
    // hcentrality->Fill(*mycent);
    // if(*mycent<centmin || *mycent>centmax) continue;
    vector<Double_t> L1_match_pt(2,-std::numeric_limits<double>::max());
    vector<Double_t> L1_match_eta(2,-std::numeric_limits<double>::max());
    vector<Double_t> L1_match_phi(2,-std::numeric_limits<double>::max());
    // double L1_match_eta=-std::numeric_limits<double>::max();
    // double L1_match_phi=-std::numeric_limits<double>::max();
    // double L1_match_pt_subl= -std::numeric_limits<double>::max();
    // double L1_match_eta_subl=-std::numeric_limits<double>::max();
    // double L1_match_phi_subl=-std::numeric_limits<double>::max();

    vector<Double_t> PFpt(2,-std::numeric_limits<double>::max());
    vector<Double_t> PFeta(2,-std::numeric_limits<double>::max());
    vector<Double_t> PFphi(2,-std::numeric_limits<double>::max());
    // double PFeta=-std::numeric_limits<double>::max();
    // double PFphi=-std::numeric_limits<double>::max();
    // double PFpt_subl=-std::numeric_limits<double>::max();
    // double PFeta_subl=-std::numeric_limits<double>::max();
    // double PFphi_subl=-std::numeric_limits<double>::max();
    vector<Double_t> drmin(2,std::numeric_limits<double>::max());
    // double drmin=std::numeric_limits<double>::max();
    // double drmin_subl=std::numeric_limits<double>::max();
    for(int l=0;l<*nref;l++){
      //Below cuts are taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
      if( abs(jteta[l])>2.5 ) continue;
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
      //if greater than leading, set subleading to leading, take new value for leading
      // std::cout << "jet number " << l << " with pt " << CorrectedPT << " and eta " << jteta[l] << std::endl;
      if(CorrectedPT > PFpt.at(0)){
        PFpt.at(1) = PFpt.at(0);
        PFeta.at(1) = PFeta.at(0);
        PFphi.at(1) = PFphi.at(0);
        PFpt.at(0) = CorrectedPT;
        PFeta.at(0) = jteta[l];
        PFphi.at(0) = jtphi[l];
      }
      //else if greater than subleading, replace it
      else if(CorrectedPT > PFpt.at(1)){
        PFpt.at(1) = CorrectedPT;
        PFeta.at(1) = jteta[l];
        PFphi.at(1) = jtphi[l];
      }
    }
    // std::cout << "Leading jet: " << PFpt.at(0) << " subleading jet: " << PFpt.at(1) << std::endl;
    TLorentzVector leading_PF;
    leading_PF.SetPtEtaPhiM(PFpt.at(0), PFeta.at(0), PFphi.at(0), 0);
    TLorentzVector leading_L1;
    leading_L1.SetPtEtaPhiM(0, 0, 0, 0);
    TLorentzVector subleading_PF;
    subleading_PF.SetPtEtaPhiM(PFpt.at(1), PFeta.at(1), PFphi.at(1), 0);
    TLorentzVector subleading_L1;
    subleading_L1.SetPtEtaPhiM(0, 0, 0, 0);
    vector<TLorentzVector> L1s = {leading_L1, subleading_L1};
    // std::cout << "angle between jets " << abs(leading_PF.DeltaPhi(subleading_PF)) << std::endl;
    if( PFpt.at(1) < 0 || abs(leading_PF.DeltaPhi(subleading_PF)) < 2. ) continue;
    for(long unsigned n=0; n < (*emuJetPt).size(); n++){
      if(n!=0){
        TLorentzVector one, two;
        one.SetPtEtaPhiM(emuJetPt->at(n), emuJetEta->at(n), emuJetPhi->at(n), 0);
        two.SetPtEtaPhiM(emuJetPt->at(n-1), emuJetEta->at(n-1), emuJetPhi->at(n-1), 0);
        if(one.DeltaR(two) < 0.4){
          std::cout << one.DeltaR(two) << " distance between two l1 objects" << std::endl;
        }
      }
      // std::cout << emuJetPt->at(n) << " L1 jet pt for " << n << std::endl;
      std::vector<double> dr_PF_L1 = {getDR(PFeta.at(0), PFphi.at(0), emuJetEta->at(n), emuJetPhi->at(n)), getDR(PFeta.at(1), PFphi.at(1), emuJetEta->at(n), emuJetPhi->at(n))}; 
      // double dr_sublead = getDR(PFeta.at(1), PFphi.at(1), emuJetEta->at(n), emuJetPhi->at(n));
      for(size_t r{0}; r < dr_PF_L1.size(); ++r){
        if(dr_PF_L1.at(r) < drmin.at(r)){
          L1s.at(r).SetPtEtaPhiM(emuJetPt->at(n), emuJetEta->at(n), emuJetPhi->at(n), 0);
          drmin.at(r) = dr_PF_L1.at(r);
        }
      }
      // if(dr_lead < drmin){
      //   L1_match_pt=emuJetPt->at(n);
      //   L1_match_eta=emuJetEta->at(n);
      //   L1_match_phi=emuJetPhi->at(n);
      //   drmin = dr;
      // }
    }
    hjetl1pfdR->Fill(drmin.at(1), *weight);
    hjet_all_pf->Fill(PFpt.at(1), *weight);
    dR_vs_PFpt->Fill(drmin.at(1), PFpt.at(1), *weight);
    // std::cout << " L1 leading jet: " << L1s.at(0).Pt() << " subleading jet: " << L1s.at(1).Pt() << std::endl;
    if(L1s.at(1).Pt() == 0) continue;
        // double drjets=getDR(l1eta,l1phi,PF_match_eta,PF_match_phi);  
    if(drmin.at(0) > 0.4 or drmin.at(1) > 0.4) continue;
    hjet_all_pf_matched->Fill(PFpt.at(1), *weight);
    for(Int_t l{0}; l < threshold.size(); l++){
      if( L1s.at(0).Pt() >= threshold.at(l) && L1s.at(1).Pt() >= threshold.at(l) && abs(L1s.at(1).DeltaPhi(L1s.at(0))) > 2.*TMath::Pi()/3.){ 
        emuMatchedHistos.at(l).Fill(PFpt.at(1), *weight);
        // std::cout << "Matched jets have pts " << L1s.at(0).Pt() << " and " << L1s.at(1).Pt() << " L1s and for PF " << leading_PF.Pt() << " and " << subleading_PF.Pt() << std::endl;
        // std::cout << "Matched jets have etas " << L1s.at(0).Eta() << " and " << L1s.at(1).Eta() << " L1s and for PF " << leading_PF.Eta() << " and " << subleading_PF.Eta() << std::endl;
        // std::cout << "Matched jets have phis " << L1s.at(0).Phi() << " and " << L1s.at(1).Phi() << " L1s and for PF " << leading_PF.Phi() << " and " << subleading_PF.Phi() << std::endl;
      }
    }
  }
  delete jets;
  delete ev;
  delete emuTree;
  in_f->Close();
  delete in_f;
  fout->cd();
  hjet_all_pf->Write();
  hjet_all_pf_matched->Divide(hjet_all_pf);
  hjet_all_pf_matched->Write();
  hjetl1pfdR->Write();
  dR_vs_PFpt->Write();
  for(size_t i{0}; i<threshold.size(); ++i){
    // emuHistos.at(i).Write();
    emuMatchedHistos.at(i).Write();
  }
  // for(int k=bitmin;k<bitmax;k++){
  //   hjet_fired_calo[k-bitmin]->Write();
  //   hjet_fired_calo_ratio[k-bitmin]=(TH1F*)hjet_fired_calo[k-bitmin]->Clone(Form("hjet_fired_calo_ratio[%i]",k-bitmin));
  //   hjet_fired_calo_ratio[k-bitmin]->Write();
  // }
  fout->Close();
  exit(0);
}



