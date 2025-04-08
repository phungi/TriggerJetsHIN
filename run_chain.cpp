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

void GetFiles(TString input, vector<string>& files) {
    TSystemDirectory dir(input, input);
    TList *list = dir.GetListOfFiles();

    if (list) {
        TSystemFile *file;
        string fname;
        TIter next(list);
        while ((file = (TSystemFile*) next())) {
            fname = file->GetName();

            if (file->IsDirectory() && (fname.find(".") == string::npos)) {
                string newDir = string(input) + fname + "/";
                GetFiles(newDir.c_str(), files);
            }
            else if ((fname.find(".root") != string::npos)) {
                files.push_back(string(input) + fname);
                cout << files.back() << endl;
            }
        }
    }

    return;
}

void FillChain(TChain& chain, vector<string>& files) {
    for (auto file : files) {
        chain.Add(file.c_str());
        std::cout << "adding " << file << " to chain" << std::endl;
        std::cout << chain.GetEntries() << " total number of entries" << std::endl;
    }

}

int Efficiency(TString input){
    vector<string> files;
    GetFiles(input, files);
    TChain offChain("hiEvtAnalyzer/HiTree");
    FillChain(offChain, files);
    // TFile *in = new TFile(input, "READ");
    // TTree *old_tree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    offChain.SetBranchStatus("*",0);
    for(auto activeBranchName: {"run"}){
        offChain.SetBranchStatus(activeBranchName,1);
    }
    TFile *outf = new TFile("out.root", "recreate");
    // auto newtree = offChain.GetTree()->CloneTree(0);
    offChain.CloneTree(-1,"fast");
    // newtree->GetBranch("run")->SetFile("out.root");
    // newtree->CopyEntries(offChain.GetTree());
    // cout << newtree->GetEntries() << " new tree entries" << std::endl;
    outf->Write();
    return 0;
    // offChain.Add(input);
   
    return 0;
}

int run_chain(TString input) {
    return Efficiency(input);
}