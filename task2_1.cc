#include "common.h"

void task2_1()
{
    int cate = 0;
    
    RooRealVar m("m","",5.0,5.8);
    RooDataSet *rds_mc = new RooDataSet("rds_mc","",RooArgSet(m));
    
    TFile *fin = new TFile("bupsikMc.root");
    TTree *tin = (TTree*)fin->Get("bupsikMc");
    
    int cate_t;
    float m_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);

    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        if (m_t<5.0 || m_t>=5.8) continue;
        m.setVal(m_t);
        rds_mc->add(RooArgSet(m));
    }
    delete fin;
    
    RooRealVar sig_mean1("sig_mean1","",5.28,5.2,5.4);
    RooRealVar sig_mean2("sig_mean2","",5.28,5.2,5.4);
    RooRealVar sig_sigma1("sig_sigma1","",0.030,0.005,0.060);
    RooRealVar sig_sigma2("sig_sigma2","",0.080,0.040,0.200);
    RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
    RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
    RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
    RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));
    
    pdf_sig.fitTo(*rds_mc);
    
    RooPlot *frame = m.frame(Title(" "),Bins(80));
    rds_mc->plotOn(frame, Name("t_rds_mc"));
    pdf_sig.plotOn(frame, Name("t_pdf_sig"), LineWidth(3));
    
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
    leg->SetHeader(Form("Category %d",cate));
    leg->AddEntry(frame->findObject("t_rds_mc"),"Simluation","EP");
    leg->AddEntry(frame->findObject("t_pdf_sig"),"Fit","L");
    leg->Draw();
    
    canvas->Print("task2_1.pdf");
}