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
    TChain offChain("ak4PFJetAnalyzer/t");
    // FillChain(offChain, files);
    offChain.Add(input);
    TTreeReader offReader(&offChain);
    // TTreeReaderValue<int>   jetN(offReader, "nref");
    // TTreeReaderArray<float> jetPt(offReader, "jtpt");
    // TTreeReaderArray<float> jetEta(offReader, "jteta");
    // TTreeReaderArray<float> jetPhi(offReader, "jtphi");

    TTreeReaderValue<Int_t> jetN(offReader,"ncalo");
    TTreeReaderArray<Float_t> jetPt(offReader,"calopt");
    TTreeReaderArray<Float_t> jetEta(offReader,"caloeta");
    TTreeReaderArray<Float_t> jetPhi(offReader,"calophi");

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
    // JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
    JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4Calo_pp_Reco_v0_12-21-2023.txt");  
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
    TH2F *dR_vs_PFpt = new TH2F("dR_vs_PFpt", "#DeltaR vs PFlow pt",50, 0, 5, 100, 0, 300);
    TH2F *matched_calo_vs_l1_pt = new TH2F("matched_calo_vs_L1_pt", "matched calo vs L1 pt; L1 p_{T}; calo p_{T}", 50, 0, 300, 50, 0, 300);
    TH2F *calo_vs_l1_pt = new TH2F("calo_vs_L1_pt", "calo vs L1 pt; L1 p_{T}; calo p_{T}", 50, 0, 300, 50, 0, 300);
    TH2F *matched_calo_vs_l1_pt_zoom = new TH2F("matched_calo_vs_L1_pt_zoom", "matched calo vs L1 pt; L1 p_{T}; calo p_{T}", 50, 0, 50, 50, 0, 50);
    TH2F *calo_vs_l1_pt_zoom = new TH2F("calo_vs_L1_pt_zoom", "calo vs L1 pt; L1 p_{T}; calo p_{T}", 50, 0, 50, 50, 0, 50);
    TH1F *calo_l1_pt_diff = new TH1F("calo_l1_pt_diff", "calo l1 #Deltap_{T};#Deltap_{T};",50,-25,25);
    TH1F *uncorr_calo_l1_pt_diff = new TH1F("uncorr_calo_l1_pt_diff", "uncorr. calo l1 #Deltap_{T};#Deltap_{T};",50,-25,25);
    for(Int_t i{0}; i<threshold.size(); ++i){
        TH1F a(Form("jetpt_trig_calo_maybe_dont_use_%i",i), "", nbins, min, max);
        TH1F b(Form("jetpt_trig_calo%i",i), "", nbins, min, max);
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
        float maxUncorrCaloPt = -999;

        float emuMaxJetPt = -999;
        float emuMatchedJetPt = -999;
        float minDR = 10;

        float uncorrected_calo = -999;
        TLorentzVector PF_jet;
        /* iterate through jets and find the jet with max pT */
        for (int i = 0; i < *jetN; ++i) {
            if (TMath::Abs(jetEta[i]) > 2) { continue; }
            Bool_t passSelections = true;
            Int_t totalParticles = CHM[i] + NHM[i] + CEM[i] + NEM[i] + MUM[i];
            Int_t NumNeutralParticles = NHM[i] + NEM[i];
            // if( abs(jetEta[i])<=2.6 && CEMF[i]<0.8 && CHM[i]>0 && CHF[i]>0.01 && NEMF[i]<0.9 && MUF[i] <0.8 && NHF[i] < 0.9 && totalParticles > 1 ) passSelections = true;
            // if( abs(jteta[l])>2.6 && abs(jteta[l])<=2.7 && CEMF[l]<0.8 && CHM[l]>0 && NEMF[l]<0.99 && MUF[l] <0.8 && NHF[l] < 0.9 ) passSelections = true;
            // if( abs(jteta[l])>2.7 && abs(jteta[l])<=3.0 && NEMF[l]<0.99 && NumNeutralParticles > 1 ) passSelections = true;
            // if( abs(jteta[l])>3.0 && NEMF[l]<0.90 && NumNeutralParticles > 10 ) passSelections = true;
            //if jet doesn't satisfy above conditions, skip
            if(!passSelections) continue;
            JEC->setJetPt(jetPt[i]);
            JEC->setJetEta(jetEta[i]);
            double Correction = JEC->getCorrection();
            std::cout << "correction " << Correction << std::endl;
            uncorrected_calo = jetPt[i];
            jetPt[i] = jetPt[i]*JEC->getCorrection();            
            if (jetPt[i] > maxJetPt) {
                maxJetPt = jetPt[i];
                maxJetEta = jetEta[i];
                maxJetPhi = jetPhi[i];
                maxUncorrCaloPt = uncorrected_calo;
                
            }
        }
        if (maxJetPt > 0) {
            recoHist.Fill(maxJetPt);
            PF_jet.SetPtEtaPhiM(maxJetPt, maxJetEta,maxJetPhi, 0);
            /* iterate through emu jets and find matched and unmatched jets with max pT */
            for (size_t i = 0; i < (*emuJetPt).size(); ++i) {
                if ((*emuJetPt)[i] > emuMaxJetPt) {
                    emuMaxJetPt = (*emuJetPt)[i];
                }
                // TLorentzVector L1obj;
                // L1obj.SetPtEtaPhiM((*emuJetPt)[i],(*emuJetEta)[i], (*emuJetPhi)[i], 0);
                // std::cout << L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) << "difference" << std::endl;
                // if(L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) > 0.01){
                //     std::cout << L1obj.DeltaR(PF_jet) - dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) << " difference!" << std::endl;
                // }
                
                if (dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) < minDR) {
                    minDR = dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi);
                    emuMatchedJetPt = (*emuJetPt)[i];
                }
            }
            dR_vs_PFpt->Fill(minDR, maxJetPt);
            calo_vs_l1_pt->Fill(emuMaxJetPt, maxJetPt);
            calo_vs_l1_pt_zoom->Fill(emuMaxJetPt, maxJetPt);
            if(minDR < 0.4){
                matched_calo_vs_l1_pt->Fill(emuMaxJetPt, maxJetPt);
                matched_calo_vs_l1_pt_zoom->Fill(emuMaxJetPt, maxJetPt);
                calo_l1_pt_diff->Fill(emuMaxJetPt-maxJetPt);
                uncorr_calo_l1_pt_diff->Fill(emuMaxJetPt - maxUncorrCaloPt );
            }
            for(size_t i{0}; i<threshold.size(); i++){
                if (emuMaxJetPt >= threshold.at(i)) {
                    emuHistos.at(i).Fill(maxJetPt);
                }
                if (emuMatchedJetPt >= threshold.at(i) && minDR < 0.4) {
                    emuMatchedHistos.at(i).Fill(maxJetPt);
                }
            }
        }
    }

    // TGraphAsymmErrors emuRecoEff(&emuHist, &recoHist, "cl=0.683 b(1,1) mode");
    // TGraphAsymmErrors emuRecoMatchedEff(&emuMatchedHist, &recoHist, "cl=0.683 b(1,1) mode");

    /* plot the turn ons vs reco jet pt */
    // TCanvas recoCanvas("recoCanvas", "", 0, 0, 500, 500);
    // recoCanvas.cd();

    // emuRecoMatchedEff.GetXaxis()->SetTitle("Reco Jet pT (GeV)");
    // emuRecoMatchedEff.GetXaxis()->CenterTitle(true);
    // emuRecoMatchedEff.GetYaxis()->SetTitle("Efficiency");
    // emuRecoMatchedEff.GetYaxis()->CenterTitle(true);

    // emuRecoMatchedEff.SetMarkerColor(46);
    // emuRecoMatchedEff.SetLineColor(46);
    // emuRecoMatchedEff.SetMarkerSize(0.5);
    // emuRecoMatchedEff.SetMarkerStyle(20);
    // emuRecoMatchedEff.Draw();

    // emuRecoEff.SetMarkerColor(30);
    // emuRecoEff.SetLineColor(30);
    // emuRecoEff.SetMarkerSize(0.5);
    // emuRecoEff.SetMarkerStyle(20);
    // emuRecoEff.Draw("LP SAME");

    // TLegend recoLegend(0.53, 0.12 ,0.88, 0.3);
    // recoLegend.SetTextSize(0.03);
    // recoLegend.SetHeader(seed.c_str());
    // recoLegend.AddEntry(&emuRecoEff, "Not #DeltaR Matched", "lep");
    // recoLegend.AddEntry(&emuRecoMatchedEff, "#DeltaR Matched", "lep");
    // recoLegend.Draw();

    // recoCanvas.SaveAs("PfJetEfficiency.pdf");

    /* save histograms to file so I can look at them */
    TFile* fout = new TFile("out.root", "recreate");
    for(size_t i{0}; i<threshold.size(); ++i){
        emuHistos.at(i).Write();
        emuMatchedHistos.at(i).Write();
    }
    // emuHist.Write();
    // emuMatchedHist.Write();
    matched_calo_vs_l1_pt->Write();
    calo_vs_l1_pt->Write();
    matched_calo_vs_l1_pt_zoom->Write();
    calo_vs_l1_pt_zoom->Write();
    calo_l1_pt_diff->Write();
    uncorr_calo_l1_pt_diff->Write();
    dR_vs_PFpt->Write();
    recoHist.Write();
    fout->Close();
    return 0;
}

int Molly_calo(TString input) {
    return Efficiency(input);
}