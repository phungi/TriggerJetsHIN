void simplePlots(){
	TString input = "/eos/cms/store/group/phys_heavyions/cbaldene/HIZeroBias/HIZeroBias2/2023PbPb_jetSeed1p5_v1_0_1_Oct4th_midFill/241004_235101/0000/JetZDC_ZeroBias_22.root";
	TH1::SetDefaultSumw2();
	TFile *for_ls = new TFile(input, "READ");
	TChain emuChain("l1UpgradeEmuTree/L1UpgradeTree");
	emuChain.Add(input);
    TTreeReader emuReader(&emuChain);
    TTreeReaderArray<float> emuJetPt(emuReader, "L1Upgrade.jetEt");
    TTreeReaderArray<float> emuJetEta(emuReader, "L1Upgrade.jetEta");
    TTreeReaderArray<float> emuJetPhi(emuReader, "L1Upgrade.jetPhi");
    TTreeReaderValue<UShort_t>   jetN(emuReader, "L1Upgrade.nJets");
    Long64_t totalEvents = emuReader.GetEntries(true);
    TH2D *Et_vs_phi = new TH2D("Et vs phi, Et>8 GeV", "L1 jet E_{T} vs #phi, E_{T}>8 GeV;#phi;E_{T}", 50, -3.14, 3.14, 50, 0, 300);
    TH2D *Et_vs_eta = new TH2D("Et vs eta, Et>8 GeV", "L1 jet E_{T} vs #eta, E_{T}>8 GeV;#eta;E_{T}", 50, -5, 5, 50, 0, 300);
    TH2D *phi_vs_eta = new TH2D("phi vs eta, Et>8 GeV", "L1 jet #phi vs #eta, E_{T}>8 GeV;#eta;#phi", 50, -5, 5, 50, -3.14, 3.14);
    for (Long64_t i = 0; i < totalEvents; i++){
    	if(i%1000==0) std::cout << i << "/" << totalEvents << std::endl;
    	emuReader.Next();
    	for (size_t j = 0; j < *jetN; ++j){
    		if(emuJetPt[j] < 8) continue;
    		Et_vs_phi->Fill(emuJetPhi[j], emuJetPt[j]);
    		Et_vs_eta->Fill(emuJetEta[j], emuJetPt[j]);
    		phi_vs_eta->Fill(emuJetEta[j], emuJetPhi[j]);
    	}
    }
    TCanvas *c = new TCanvas();
    TString pdffile = "2D_plots.pdf";
    c->cd();
    Et_vs_eta->Draw("colz");
    c->Print( pdffile + "(","pdf");
    c->Clear();
    Et_vs_phi->Draw("colz");
    c->Print( pdffile,"pdf");
    c->Clear();
    phi_vs_eta->Draw("colz");
    c->Print( pdffile + ")","pdf");
    exit(0);
}