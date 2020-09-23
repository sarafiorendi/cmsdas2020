#include "common.h"

void task5_1(int cate, double bdt_min)
{
//     int cate = 0;
//     double bdt_min = 0.2;
    
    RooRealVar bdt("bdt","",0.1,1.0);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds_sig = new RooDataSet("rds_sig","",RooArgSet(bdt,wgt),"wgt");
    RooDataSet *rds_bkg = new RooDataSet("rds_bkg","",RooArgSet(bdt,wgt),"wgt");
    
    TFile *fin = new TFile("/afs/cern.ch/work/k/kfjack/public/cmsdas/bmm/bsmmMc.root");
    TTree *tin = (TTree*)fin->Get("bsmmMc");
    
    int cate_t;
    float m_t, bdt_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("bdt",&bdt_t);
    
    double weight = effyield[cate][N_bsmm]/(double)tin->GetEntries(Form("cate==%d",cate));
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        bdt.setVal(bdt_t);
        wgt.setVal(weight);
        rds_sig->add(RooArgSet(bdt,wgt),weight);
    }
    delete fin;
    
    fin = new TFile("/afs/cern.ch/work/k/kfjack/public/cmsdas/bmm/bmmData-blind.root");
    tin = (TTree*)fin->Get("bmmData");
    
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("bdt",&bdt_t);
    
    weight = 0.25/0.45; // scale from 5.45-5.90 to 5.20-5.45, assuming a uniform mass dist.
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        if (m_t<5.45) continue;
        bdt.setVal(bdt_t);
        wgt.setVal(weight);
        rds_bkg->add(RooArgSet(bdt,wgt),weight);
    }
    delete fin;
    
    cout << "Category: " << cate << endl;
    double ns = rds_sig->sumEntries(Form("bdt>%g",bdt_min));
    double nb = rds_bkg->sumEntries(Form("bdt>%g",bdt_min));
    if (ns>0. && nb>0.) {
        double Z = sqrt(2.*((ns+nb)*log(1.+ns/nb)-ns));
        cout << "bdt_min: " << bdt_min << ", ns: " << ns << ", nb: " << nb << ", Z: " << Z << endl;
    }
}
