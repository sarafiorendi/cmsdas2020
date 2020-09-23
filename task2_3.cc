#include "common.h"

void task2_3()
{
    int cate_t = 0;
    
    RooRealVar m("m","",5.0,5.8);
    RooRealVar cate("cate","",-1,20);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    
    TFile *fin = new TFile("/afs/cern.ch/work/k/kfjack/public/cmsdas/bmm/bupsikData.root");
    TTree *tin = (TTree*)fin->Get("bupsikData");
    
    RooDataSet* rds_data   = new RooDataSet("rds_mc", "rds_mc", tin,  RooArgSet(cate, m, wgt));
    rds_data = (RooDataSet*)rds_data->reduce(RooArgSet(wgt, m), Form("m >= 5.0 && m <= 5.8 && cate==%d", cate_t));
    
/* From task2_1.cc
 *  1  sig_frac     9.65521e-01   1.45707e-03   9.11950e-05   1.03937e+00
 *   2  sig_mean1    5.27839e+00   1.01238e-04   6.31736e-05  -2.17841e-01
 *    3  sig_mean2    5.25550e+00   2.39898e-03   2.99203e-04  -4.61218e-01
 *     4  sig_sigma1   2.99468e-02   9.04925e-05   3.42182e-05  -9.29776e-02
 *      5  sig_sigma2   1.12346e-01   2.50130e-03   2.66008e-04  -9.58226e-02
 *       */
    RooRealVar sig_shift("sig_shift","",0.,-0.02,0.02);
    RooRealVar sig_scale("sig_scale","",1.,0.8,1.2);
    
    RooRealVar sigmc_mean1("sigmc_mean1","",5.27839e+00);
    RooRealVar sigmc_mean2("sigmc_mean2","",5.25550e+00);
    RooRealVar sigmc_sigma1("sigmc_sigma1","",2.99468e-02);
    RooRealVar sigmc_sigma2("sigmc_sigma2","",1.12346e-01);
    
    RooAddition sig_mean1("sig_mean1","",RooArgList(sigmc_mean1,sig_shift));
    RooAddition sig_mean2("sig_mean2","",RooArgList(sigmc_mean2,sig_shift));
    RooProduct sig_sigma1("sig_sigma1","",RooArgList(sigmc_sigma1,sig_scale));
    RooProduct sig_sigma2("sig_sigma2","",RooArgList(sigmc_sigma2,sig_scale));
    RooRealVar sig_frac("sig_frac","",9.65521e-01);
    RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
    RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
    RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));
        
    RooRealVar comb_coeff("comb_coeff","",-1.2,-10.,10.);
    RooExponential pdf_comb("pdf_comb","",m,comb_coeff);
    
    RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.08);
    RooRealVar jpsix_shift("jpsix_shift","",5.13,5.12,5.16);
    RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(m,jpsix_shift,jpsix_scale));
    
    RooRealVar n_sig("n_sig","",100000,0.,1E8);
    RooRealVar n_comb("n_comb","",80000,0.,1E6);
    RooRealVar n_jpsix("n_jpsix","",20000,0.,1E5);
    RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb,pdf_jpsix),RooArgList(n_sig,n_comb,n_jpsix));
    
    model.fitTo(*rds_data, Extended(true), SumW2Error(true));
    
    RooPlot *frame = m.frame(Title(" "),Bins(80));
    rds_data->plotOn(frame, Name("t_rds_data"));
    model.plotOn(frame, Name("t_model"), LineWidth(3));
    model.plotOn(frame, Name("t_pdf_comb"), Components("pdf_comb"), LineWidth(3), LineStyle(2), LineColor(kGray+1));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.15,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();
    
    TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetHeader(Form("Category %d",cate_t));
    leg->AddEntry(frame->findObject("t_rds_data"),"Data","EP");
    leg->AddEntry(frame->findObject("t_model"),"Fit","L");
    leg->AddEntry(frame->findObject("t_pdf_comb"),"Combinatorial bkg.","L");
    leg->Draw();
    
    canvas->Print("task2_3.pdf");
}
