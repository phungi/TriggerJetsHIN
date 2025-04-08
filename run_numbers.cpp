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
    // ROOT::RDataFrame df("hiEvtAnalyzer/HiTree", input);
    // df.Alias("run","run").Snapshot("tree","out.root",{"run"});
    TFile *in = new TFile(input, "READ");
    TTree *old_tree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    old_tree->SetBranchStatus("*",0);
    for(auto activeBranchName: {"run"}){
        old_tree->SetBranchStatus(activeBranchName,1);
    }
    TFile *outf = new TFile("out10.root", "recreate");
    auto newtree = old_tree->CloneTree(0);
    newtree->GetBranch("run")->SetFile("out10.root");
    newtree->CopyEntries(old_tree);
    outf->Write();
    return 0;
}

int run_numbers(TString input) {
    return Efficiency(input);
}