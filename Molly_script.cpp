/*
Input: Folder of L1Ntuples
Output: A plot of the jet turn-ons with and with out L1 dR matching vs PF jet pT
*/

// LOOK AT THE ITERATOR NAMES FOR EVENTS AND JETS
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

int Efficiency(TString input) {
    TH1::SetDefaultSumw2();
    TFile *for_ls = new TFile(input, "READ");
    TString AK4PF_string = "ak4PFJetAnalyzer";
    TString AKCS4PF_string = "akCs4PFJetAnalyzer";
    TDirectory *dir_ak4 = for_ls->GetDirectory(AK4PF_string);
    TDirectory *dir_PuCs4 = for_ls->GetDirectory(AKCS4PF_string);
    TString chain_name = "";
    if(dir_PuCs4){
        chain_name = AKCS4PF_string + "/t";
    }
    else if(dir_ak4){
        chain_name = AK4PF_string + "/t";
    }
    else{
        std::cout << "Don't have " << AK4PF_string << " or " << AKCS4PF_string << " directories, exiting!" << std::endl;
        return 1;
    }
    // this is probably wrong, look at the ttreereader
    TChain offChain(chain_name);
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
    TDirectory *dir_L1_emulation = for_ls->GetDirectory("l1UpgradeEmuTree");
    TDirectory *dir_L1objects = for_ls->GetDirectory("l1object");
    TString L1_chain_name = "";

    if(dir_L1objects){
        std::cout << "l1object" << std::endl;
        L1_chain_name = "l1object/L1UpgradeFlatTree";
    }
    else if(dir_L1_emulation){
        std::cout << "l1UpgradeEmuTree" << std::endl;
        L1_chain_name = "l1UpgradeEmuTree/L1UpgradeTree";
    }
    
    TChain emuChain(L1_chain_name);
    // FillChain(emuChain, files);
    emuChain.Add(input);
    TTreeReader emuReader(&emuChain);
    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
    TTreeReaderValue<vector<float>> emuJetEta(emuReader, "jetEta");
    TTreeReaderValue<vector<float>> emuJetPhi(emuReader, "jetPhi");
    TChain weightChain("hiEvtAnalyzer/HiTree");
    weightChain.Add(input);
    TTreeReader weightReader(&weightChain);
    Bool_t isMC = false;
    Float_t weight = 1.;
    if(isMC){
        weightChain.SetBranchAddress("weight",&weight);
    }
    // thresholds for the L1 trigger
    std::vector<float> threshold = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
    std::vector<TH1F> emuMatchedHistos;
    std::vector<TH1F> emuHistos;

    gSystem->Load("libFWCoreFWLite.so");

    //add an if to load the appropriate corrections for jets depending on dir 
    JetCorrectorParameters *L2JetPar=new JetCorrectorParameters("/afs/cern.ch/user/v/vavladim/public/L1trigger/TriggerJetsHIN/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");  
    vector<JetCorrectorParameters> vPar;
    vPar.push_back(*L2JetPar);

    FWLiteEnabler::enable();    
    FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vPar);

    /* create histograms for efficiency plots */
    int nbins = 75;
    float min = 0;
    float max = 150;

    TH1F emuHist("emuHist", "", nbins, min, max);
    TH1F emuMatchedHist("emuMatchedHist", "", nbins, min, max);
    TH1F recoHist("jetpt_all_pf", "", nbins, min, max);
    TH2F *dR_vs_PFpt = new TH2F("dR_vs_PFpt", "#DeltaR vs PFlow pt",50,0,5, 100, 0, 300);

    for(Int_t i{0}; i<threshold.size(); ++i){
        TH1F a(Form("jetpt_trig_calo_probably_dont_use%i",i), "", nbins, min, max);
        TH1F b(Form("jetpt_trig_calo%i",i), "", nbins, min, max);
        emuMatchedHistos.push_back(b);
        emuHistos.push_back(a);

    }
    Long64_t totalEvents = emuReader.GetEntries(true);
    std::cout << emuReader.GetEntries(true) << " " << offReader.GetEntries(true) << " " << weightReader.GetEntries(true) << " entries " << std::endl;
    /* read in information from TTrees */
    for (Long64_t e{0}; e < totalEvents; e++) {
        emuReader.Next(); offReader.Next(); weightReader.Next();
        if (e % 20000 == 0) { 
            cout << "Entry: " << e << " / " <<  totalEvents << endl; 
        }
        // std::cout << "New event" << std::endl;
        float maxJetPt = -999;
        float maxJetPhi = -999;
        float maxJetEta = -999;
        float emuMaxJetPt = -999;
        float emuMatchedJetPt = -999;
        float minDR = 10;
        float l1_within_dR = 0.;
        TLorentzVector PF_jet;
        /* iterate through jets and find the jet with max pT */
        for (int i = 0; i < *jetN; ++i) {
            if(i!=0) continue; //check, but they should be ordered in pT, so no need to loop over all of them // depends on what you want to look
            if( TMath::Abs(jetEta[i]) > 2 ) continue; //look only at jets contained within 2.4>|eta|
            Bool_t passSelections = false;
            Int_t totalParticles = CHM[i] + NHM[i] + CEM[i] + NEM[i] + MUM[i];
            Int_t NumNeutralParticles = NHM[i] + NEM[i];
            // if( abs(jetEta[i])<=2.6 && CEMF[i]<0.8 && CHM[i]>0 && CHF[i]>0.01 && NEMF[i]<0.9 && MUF[i] <0.8 && NHF[i] < 0.9 && totalParticles > 1 ) passSelections = true;
            // original line above, added CHF<0.9 below
            if( abs(jetEta[i])<=2.6 && CEMF[i]<0.8 && CHM[i]>0 && CHF[i]>0.01 && CHF[i]<0.9 && NEMF[i]<0.9 && MUF[i] <0.8 && NHF[i] < 0.9 && totalParticles > 1 ) passSelections = true;
            //don't care about rest as we have an eta cut, but adding them below anyway for completeness
            // if( abs(jetEta[l])>2.6 && abs(jteta[l])<=2.7 && CEMF[l]<0.8 && CHM[l]>0 && NEMF[l]<0.99 && MUF[l] <0.8 && NHF[l] < 0.9 ) passSelections = true;
            // if( abs(jetEta[l])>2.7 && abs(jteta[l])<=3.0 && NEMF[l]<0.99 && NumNeutralParticles > 1 ) passSelections = true;
            // if( abs(jteta[l])>3.0 && NEMF[l]<0.90 && NumNeutralParticles > 10 ) passSelections = true;
            
            if(!passSelections) continue; //if jet doesn't satisfy above conditions, skip
            // JEC->setJetPt(jetPt[i]);
            // JEC->setJetEta(jetEta[i]);
            // double Correction = JEC->getCorrection();
            // std::cout << "correction " << Correction << std::endl;
            // jetPt[i] = jetPt[i]*JEC->getCorrection();
            // std::cout << jetPt[i] << " jet pt" << std::endl;
            //aren't these already ordered in pt? anyway, loop and pick the highest
            if (jetPt[i] > maxJetPt) {
                maxJetPt = jetPt[i];
                maxJetEta = jetEta[i];
                maxJetPhi = jetPhi[i];
                
            }
        }
        // std::cout << "max jet pt " << maxJetPt << std::endl;
        if (maxJetPt > 0) {
            recoHist.Fill(maxJetPt, weight);
            //mass doesn't matter, just declaring a lorentz vector to get dR function
            PF_jet.SetPtEtaPhiM(maxJetPt, maxJetEta,maxJetPhi, 0);
            /* iterate through emu jets and find matched and unmatched jets with max pT */
            // std::cout << (*emuJetPt).size() << " size of L1 objects" << std::endl;
            for (size_t i = 0; i < (*emuJetPt).size(); ++i) {
                // std::cout << (*emuJetPt)[i] << " l1 calo jet pt" << std::endl;
                if ((*emuJetPt)[i] > emuMaxJetPt) {
                    emuMaxJetPt = (*emuJetPt)[i];
                }
                TLorentzVector L1obj;
                L1obj.SetPtEtaPhiM((*emuJetPt)[i],(*emuJetEta)[i], (*emuJetPhi)[i], 0);
                //until now, taking the closest l1 object to the jet axis, change to highest energy l1 object inside R=0.4 below
                // ----------
                // if (dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi) < minDR) {
                //     minDR = dr((*emuJetEta)[i], (*emuJetPhi)[i], maxJetEta, maxJetPhi);
                //     emuMatchedJetPt = (*emuJetPt)[i];
                // }
                // ----------
                if (PF_jet.DeltaR(L1obj) < 0.4 && emuMatchedJetPt < (*emuJetPt)[i]) {
                    //it's no longer a minDR in this case, but doesn't matter
                    minDR = PF_jet.DeltaR(L1obj);
                    emuMatchedJetPt = (*emuJetPt)[i];
                }
                // ----------
            }
            dR_vs_PFpt->Fill(minDR, maxJetPt, weight);
            for(size_t i{0}; i<threshold.size(); i++){
                if (emuMaxJetPt >= threshold.at(i)) {
                    emuHistos.at(i).Fill(maxJetPt, weight);
                }
                if (emuMatchedJetPt >= threshold.at(i) && minDR < 0.4) {
                    emuMatchedHistos.at(i).Fill(maxJetPt, weight);
                }
            }
        }
    }
    TFile* fout = new TFile("out.root", "recreate");
    for(size_t i{0}; i<threshold.size(); ++i){
        emuHistos.at(i).Write();
        emuMatchedHistos.at(i).Write();
    }
    dR_vs_PFpt->Write();
    recoHist.Write();
    fout->Close();
    return 0;
}

int Molly_script(TString input) {
    return Efficiency(input);
}