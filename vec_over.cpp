int vec_over(TString input) {
    TFile *a = new TFile(input, "READ");
    TTree *t = (TTree*)a->Get("HiTree");
    UInt_t run;
    t->SetBranchAddress("run", &run);
    std::vector<UInt_t> vec = {};
    for(size_t i{0}; i < t->GetEntries(); ++i){
        // std::cout << i << std::endl;
        std::cout << i << std::endl;
        t->GetEntry(i);
        // if(i==1000)
        //     break;
        if ( std::find(vec.begin(), vec.end(), run) == vec.end() )
            vec.push_back(run);
    }
    for(size_t v{0}; v<vec.size();++v){
        std::cout << vec.at(v) << std::endl;
    }
    return 0;
}