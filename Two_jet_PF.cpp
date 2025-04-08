/*
Input: Folder of L1Ntuples
Output: A plot of the jet turn-ons with and with out L1 dR matching vs PF jet pT
*/

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TChain.h"

#include "TMath.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

double dr(float eta1, float phi1, float eta2, float phi2) {
    float deta = TMath::Abs(eta1 - eta2);
    float dphi = TMath::Abs(phi1 - phi2);
    if (dphi > TMath::Pi()) dphi = TMath::Abs(dphi - 2*TMath::Pi());

    return TMath::Sqrt(dphi*dphi + deta*deta);
}

// void GetFiles(char const* input, vector<string>& files) {
//     TSystemDirectory dir(input, input);
//     TList *list = dir.GetListOfFiles();

//     if (list) {
//         TSystemFile *file;
//         string fname;
//         TIter next(list);
//         while ((file = (TSystemFile*) next())) {
//             fname = file->GetName();

//             if (file->IsDirectory() && (fname.find(".") == string::npos)) {
//                 string newDir = string(input) + fname + "/";
//                 GetFiles(newDir.c_str(), files);
//             }
//             else if ((fname.find(".root") != string::npos)) {
//                 files.push_back(string(input) + fname);
//                 cout << files.back() << endl;
//             }
//         }
//     }

//     return;
// }

// void FillChain(TChain& chain, vector<string>& files) {
//     for (auto file : files) {
//         chain.Add(file.c_str());
//     }
// }

int Efficiency(TString input) {
    /* read in all files in the input folder */
    // vector<string> files;
    // GetFiles(input, files);
    // TFile *input = TFile::Open(filename.data());
    /* read in reco jet information */
    TChain offChain("akFlowPuCs4PFJetAnalyzer/t");
    // FillChain(offChain, files);
    offChain.Add(input);
    TTreeReader offReader(&offChain);
    TTreeReaderValue<int>   jetN(offReader, "nref");
    TTreeReaderArray<float> jetPt(offReader, "jtpt");
    TTreeReaderArray<float> jetEta(offReader, "jteta");
    TTreeReaderArray<float> jetPhi(offReader, "jtphi");

    TTreeReaderArray<Float_t> NHF(offReader, "jtPfNHF");
    TTreeReaderArray<Float_t> CEMF(offReader, "jtPfCEF");
    TTreeReaderArray<Float_t> CHF(offReader, "jtPfCHF");
    TTreeReaderArray<Float_t> NEMF(offReader, "jtPfNEF");
    TTreeReaderArray<Float_t> MUF(offReader, "jtPfMUF");

    TTreeReaderArray<Int_t>   CHM(offReader, "jtPfCHM");
    TTreeReaderArray<Int_t>   NHM(offReader, "jtPfNHM");
    TTreeReaderArray<Int_t>   CEM(offReader, "jtPfCEM");
    TTreeReaderArray<Int_t>   NEM(offReader, "jtPfNEM");
    TTreeReaderArray<Int_t>   MUM(offReader, "jtPfMUM");
    /* read in emulated jet information */
    TChain emuChain("l1UpgradeTree/L1UpgradeTree");
    // FillChain(emuChain, files);
    emuChain.Add(input);
    TTreeReader emuReader(&emuChain);
    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
    TTreeReaderValue<vector<float>> emuJetEta(emuReader, "jetEta");
    TTreeReaderValue<vector<float>> emuJetPhi(emuReader, "jetPhi");

    string seed = "L1_SingleJet64_BptxAND";
    std::vector<float> threshold = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
    std::vector<TH1F> emuMatchedHistos;
    std::vector<TH1F> emuHistos;

    gSystem->Load("libFWCoreFWLite.so");
    JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
    vector<JetCorrectorParameters> vPar;
    vPar.push_back(*L2JetPar);

    FWLiteEnabler::enable(); 
    FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);

    /* create histograms for efficiency plots */
    int nbins = 150;
    float min = 0;
    float max = 300;

    TH1F emuHist("emuHist", "", nbins, min, max);
    TH1F emuMatchedHist("emuMatchedHist", "", nbins, min, max);
    TH1F recoHist("recoHist", "", nbins, min, max);
    TH2F *dR_vs_PFpt = new TH2F("dR_vs_PFpt", "#DeltaR vs PFlow pt",50,0,5, 100, 0, 300);

    for(Int_t i{0}; i<threshold.size(); ++i){
        TH1F a(Form("jetpt_trig_calo_%i",i), "", nbins, min, max);
        TH1F b(Form("jetpt_trig_calo_matched%i",i), "", nbins, min, max);
        emuMatchedHistos.push_back(b);
        emuHistos.push_back(a);

    }
    Long64_t totalEvents = emuReader.GetEntries(true);

    /* read in information from TTrees */
    for (Long64_t i = 0; i < totalEvents; i++) {
        emuReader.Next(); offReader.Next();

        if (i % 20000 == 0) { 
            cout << "Entry: " << i << " / " <<  totalEvents << endl; 
        }

        float maxJetPt = -999;
        float maxJetPhi = -999;
        float maxJetEta = -999;
        float sublJetPt = -999;
        float sublJetPhi = -999;


        float emuMaxJetPt = -999;
        float emuMatchedJetPt = -999;
        float sublemuMatchedJetPt = -999;
        float emuMatchedJetPhi = -999;
        float sublemuMatchedJetPhi = -999;
        float minDR = 10;
        float subminDR = 10;
        TLorentzVector PF_jet(0,0,0,0);
        TLorentzVector PF_subljet(0,0,0,0);
        TLorentzVector emu_jet(0,0,0,0);
        TLorentzVector emu_subljet(0,0,0,0);
        /* iterate through jets and find the jet with max pT */
        for (int i = 0; i < *jetN; ++i) {
            if (TMath::Abs(jetEta[i]) > 2.5) { continue; }
            Bool_t passSelections = false;
            Int_t totalParticles = CHM[i] + NHM[i] + CEM[i] + NEM[i] + MUM[i];
            Int_t NumNeutralParticles = NHM[i] + NEM[i];
            if( abs(jetEta[i])<=2.6 && CEMF[i]<0.8 && CHM[i]>0 && CHF[i]>0.01 && NEMF[i]<0.9 && MUF[i] <0.8 && NHF[i] < 0.9 && totalParticles > 1 ) passSelections = true;
            // if( abs(jteta[l])>2.6 && abs(jteta[l])<=2.7 && CEMF[l]<0.8 && CHM[l]>0 && NEMF[l]<0.99 && MUF[l] <0.8 && NHF[l] < 0.9 ) passSelections = true;
            // if( abs(jteta[l])>2.7 && abs(jteta[l])<=3.0 && NEMF[l]<0.99 && NumNeutralParticles > 1 ) passSelections = true;
            // if( abs(jteta[l])>3.0 && NEMF[l]<0.90 && NumNeutralParticles > 10 ) passSelections = true;
            //if jet doesn't satisfy above conditions, skip
            if(!passSelections) continue;
            JEC->setJetPt(jetPt[i]);
            JEC->setJetEta(jetEta[i]);
            double Correction = JEC->getCorrection();
            std::cout << "correction " << Correction << std::endl;
            jetPt[i] = jetPt[i]*JEC->getCorrection();            
            if (jetPt[i] > sublJetPt) {
                if(jetPt[i] > maxJetPt){
                    PF_subljet = PF_jet;
                    PF_jet.SetPtEtaPhiM(jetPt[i], jetEta[i], jetPhi[i], 0);
                }
                else{
                    PF_subljet.SetPtEtaPhiM(jetPt[i], jetEta[i], jetPhi[i], 0);
                }
                // maxJetPt = jetPt[i];
                // maxJetEta = jetEta[i];
                // maxJetPhi = jetPhi[i];
            }
        }
        if ( PF_subljet.Pt() > 0 && PF_subljet.DeltaPhi(PF_jet) < 2.*TMath::Pi()/3 ) {
            recoHist.Fill(PF_subljet.Pt());
            // PF_jet.SetPtEtaPhiM(maxJetPt, maxJetEta,maxJetPhi, 0);
            /* iterate through emu jets and find matched and unmatched jets with max pT */
            for (size_t i = 0; i < (*emuJetPt).size(); ++i) {
                if ((*emuJetPt)[i] > emuMaxJetPt) {
                    emuMaxJetPt = (*emuJetPt)[i];
                }
                TLorentzVector L1obj;
                L1obj.SetPtEtaPhiM((*emuJetPt)[i],(*emuJetEta)[i], (*emuJetPhi)[i], 0);
                // std::cout << L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) << "difference" << std::endl;
                if(L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], PF_subljet.Eta(), PF_subljet.Phi()) > 0.01){
                    std::cout << L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) << " difference!" << std::endl;
                }
                
                if (dr((*emuJetEta)[i], (*emuJetPhi)[i], PF_subljet.Eta(), PF_subljet.Phi()) < subminDR) {
                    subminDR = dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi);
                    emu_subljet.SetPtEtaPhiM((*emuJetPt)[i], (*emuJetEta)[i],(*emuJetPhi)[i], 0 );
                    // sublemuMatchedJetPt = (*emuJetPt)[i];
                    // sublemuMatchedJetPhi = (*emuJetPhi)[i];
                }

                if (dr((*emuJetEta)[i], (*emuJetPhi)[i], PF_jet.Eta(), PF_jet.Phi()) < minDR) {
                    minDR = dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi);
                    emu_jet.SetPtEtaPhiM((*emuJetPt)[i], (*emuJetEta)[i],(*emuJetPhi)[i], 0 );
                    // emuMatchedJetPt = (*emuJetPt)[i];
                    // emuMatchedJetPhi = (*emuJetPhi)[i];
                }
            }
            dR_vs_PFpt->Fill(minDR, maxJetPt);
            for(size_t i{0}; i<threshold.size(); i++){
                if (emu_subljet.Pt() >= threshold.at(i)) {
                    emuHistos.at(i).Fill(maxJetPt);
                }
                if (emu_subljet.Pt() >= threshold.at(i) && subminDR < 0.4 && minDR < 0.4 && emu_subljet.DeltaPhi(emu_jet) < 2 ) {
                    emuMatchedHistos.at(i).Fill(PF_subljet.Pt());
                }
            }
        }
    }

    TFile* fout = new TFile("out.root", "recreate");
    for(size_t i{0}; i<threshold.size(); ++i){
        emuHistos.at(i).Write();
        emuMatchedHistos.at(i).Write();
    }
    // emuHist.Write();
    // emuMatchedHist.Write();
    dR_vs_PFpt->Write();
    recoHist.Write();
    fout->Close();
   
    return 0;
}

void Two_jet_PF(TString input) {
    Efficiency(input);
    exit(0);
}