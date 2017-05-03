

void nxtal()
{


	TFile fin_with("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/Pi0Tuplizer_QCD_Pt-15to20_EMEnriched_withNxtalCut/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/pi0Ntuple_withNxtalCut_29Apr2017.root");
	TFile fin_without("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/Pi0Tuplizer_QCD_Pt-15to20_EMEnriched_withoutNxtalCut/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/pi0Ntuple_withoutNxtalCut_29Apr2017.root");
	TTree *tree_with = (TTree*)fin_with.Get("ntuples/Pi0Events");
	TTree *tree_without = (TTree*)fin_without.Get("ntuples/Pi0Events");


	double MaxY = 0.0;

        gStyle->SetOptStat(0);
        gStyle->SetPalette(107) ;
	TCanvas *myC = new TCanvas("myC","myC",100,100,800,700);

	//get histograms from tree
	TH1F *h_nxtal_EB_Pi0_G1 = new TH1F("h_nxtal_EB_Pi0_G1","h_nxtal_EB_Pi0_G1",10,0,10);
	tree_without->Draw("nxtalG1_rec>>h_nxtal_EB_Pi0_G1","fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_nxtal_EB_Pi0_G2 = new TH1F("h_nxtal_EB_Pi0_G2","h_nxtal_EB_Pi0_G2",10,0,10);
	tree_without->Draw("nxtalG2_rec>>h_nxtal_EB_Pi0_G2","fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_nxtal_EE_Pi0_G1 = new TH1F("h_nxtal_EE_Pi0_G1","h_nxtal_EE_Pi0_G1",10,0,10);
	tree_without->Draw("nxtalG1_rec>>h_nxtal_EE_Pi0_G1","fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_nxtal_EE_Pi0_G2 = new TH1F("h_nxtal_EE_Pi0_G2","h_nxtal_EE_Pi0_G2",10,0,10);
	tree_without->Draw("nxtalG2_rec>>h_nxtal_EE_Pi0_G2","fromPi0 && abs(etaPi0_rec)>1.5");

	TH1F *h_nxtal_EB_Eta_G1 = new TH1F("h_nxtal_EB_Eta_G1","h_nxtal_EB_Eta_G1",10,0,10);
	tree_without->Draw("nxtalG1_rec>>h_nxtal_EB_Eta_G1","!fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_nxtal_EB_Eta_G2 = new TH1F("h_nxtal_EB_Eta_G2","h_nxtal_EB_Eta_G2",10,0,10);
	tree_without->Draw("nxtalG2_rec>>h_nxtal_EB_Eta_G2","!fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_nxtal_EE_Eta_G1 = new TH1F("h_nxtal_EE_Eta_G1","h_nxtal_EE_Eta_G1",10,0,10);
	tree_without->Draw("nxtalG1_rec>>h_nxtal_EE_Eta_G1","!fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_nxtal_EE_Eta_G2 = new TH1F("h_nxtal_EE_Eta_G2","h_nxtal_EE_Eta_G2",10,0,10);
	tree_without->Draw("nxtalG2_rec>>h_nxtal_EE_Eta_G2","!fromPi0 && abs(etaPi0_rec)>1.5");

	TH1F *h_mPi0_rec_EB_withCut = new TH1F("h_mPi0_rec_EB_withCut","h_mPi0_rec_EB_withCut",100,0.05,0.25);
	tree_with->Draw("mPi0_rec>>h_mPi0_rec_EB_withCut","fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_mPi0_rec_EB_withoutCut = new TH1F("h_mPi0_rec_EB_withoutCut","h_mPi0_rec_EB_withoutCut",100,0.05,0.25);
	tree_without->Draw("mPi0_rec>>h_mPi0_rec_EB_withoutCut","fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_mPi0_rec_EE_withCut = new TH1F("h_mPi0_rec_EE_withCut","h_mPi0_rec_EE_withCut",100,0.05,0.25);
	tree_with->Draw("mPi0_rec>>h_mPi0_rec_EE_withCut","fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_mPi0_rec_EE_withoutCut = new TH1F("h_mPi0_rec_EE_withoutCut","h_mPi0_rec_EE_withoutCut",100,0.05,0.25);
	tree_without->Draw("mPi0_rec>>h_mPi0_rec_EE_withoutCut","fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_mEta_rec_EB_withCut = new TH1F("h_mEta_rec_EB_withCut","h_mEta_rec_EB_withCut",100,0.25,0.8);
	tree_with->Draw("mPi0_rec>>h_mEta_rec_EB_withCut","!fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_mEta_rec_EB_withoutCut = new TH1F("h_mEta_rec_EB_withoutCut","h_mEta_rec_EB_withoutCut",100,0.25,0.8);
	tree_without->Draw("mPi0_rec>>h_mEta_rec_EB_withoutCut","!fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_mEta_rec_EE_withCut = new TH1F("h_mEta_rec_EE_withCut","h_mEta_rec_EE_withCut",100,0.25,0.8);
	tree_with->Draw("mPi0_rec>>h_mEta_rec_EE_withCut","!fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_mEta_rec_EE_withoutCut = new TH1F("h_mEta_rec_EE_withoutCut","h_mEta_rec_EE_withoutCut",100,0.25,0.8);
	tree_without->Draw("mPi0_rec>>h_mEta_rec_EE_withoutCut","!fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_ptPi0_rec_EB_withCut = new TH1F("h_ptPi0_rec_EB_withCut","h_ptPi0_rec_EB_withCut",100,0.0,10.);
	tree_with->Draw("ptPi0_rec>>h_ptPi0_rec_EB_withCut","fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_ptPi0_rec_EB_withoutCut = new TH1F("h_ptPi0_rec_EB_withoutCut","h_ptPi0_rec_EB_withoutCut",100,0.0,10.);
	tree_without->Draw("ptPi0_rec>>h_ptPi0_rec_EB_withoutCut","fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_ptPi0_rec_EE_withCut = new TH1F("h_ptPi0_rec_EE_withCut","h_ptPi0_rec_EE_withCut",100,0.0,10.);
	tree_with->Draw("ptPi0_rec>>h_ptPi0_rec_EE_withCut","fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_ptPi0_rec_EE_withoutCut = new TH1F("h_ptPi0_rec_EE_withoutCut","h_ptPi0_rec_EE_withoutCut",100,0.0,10.);
	tree_without->Draw("ptPi0_rec>>h_ptPi0_rec_EE_withoutCut","fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_ptEta_rec_EB_withCut = new TH1F("h_ptEta_rec_EB_withCut","h_ptEta_rec_EB_withCut",100,0.0,10.);
	tree_with->Draw("ptPi0_rec>>h_ptEta_rec_EB_withCut","!fromPi0 && abs(etaPi0_rec)<1.5");
	
	TH1F *h_ptEta_rec_EB_withoutCut = new TH1F("h_ptEta_rec_EB_withoutCut","h_ptEta_rec_EB_withoutCut",100,0.0,10.);
	tree_without->Draw("ptPi0_rec>>h_ptEta_rec_EB_withoutCut","!fromPi0 && abs(etaPi0_rec)<1.5");

	TH1F *h_ptEta_rec_EE_withCut = new TH1F("h_ptEta_rec_EE_withCut","h_ptEta_rec_EE_withCut",100,0.0,10.0);
	tree_with->Draw("ptPi0_rec>>h_ptEta_rec_EE_withCut","!fromPi0 && abs(etaPi0_rec)>1.5");
	
	TH1F *h_ptEta_rec_EE_withoutCut = new TH1F("h_ptEta_rec_EE_withoutCut","h_ptEta_rec_EE_withoutCut",100,0.0,10.0);
	tree_without->Draw("ptPi0_rec>>h_ptEta_rec_EE_withoutCut","!fromPi0 && abs(etaPi0_rec)>1.5");
	

	//draw to plots
	//myC->SetLogy(1);

	MaxY = 0.0;
	if(h_nxtal_EB_Pi0_G1->GetMaximum() > MaxY) MaxY = h_nxtal_EB_Pi0_G1->GetMaximum();
	if(h_nxtal_EB_Pi0_G2->GetMaximum() > MaxY) MaxY = h_nxtal_EB_Pi0_G2->GetMaximum();
	
	h_nxtal_EB_Pi0_G1->SetLineColor(2);
	h_nxtal_EB_Pi0_G1->SetLineWidth(2);
	
	h_nxtal_EB_Pi0_G2->SetLineColor(4);
	h_nxtal_EB_Pi0_G2->SetLineWidth(2);
	

	h_nxtal_EB_Pi0_G1->Draw();
	h_nxtal_EB_Pi0_G1->GetYaxis()->SetRangeUser(1.0,1.4*MaxY);
	h_nxtal_EB_Pi0_G1->SetTitle("");
	h_nxtal_EB_Pi0_G1->GetYaxis()->SetTitle("Events");
	h_nxtal_EB_Pi0_G1->GetXaxis()->SetTitle("Nxtal");

	h_nxtal_EB_Pi0_G2->Draw("same");

	TLegend *leg_nxtal_EB_Pi0 = new TLegend(0.2,0.7,0.8,0.89);	
	leg_nxtal_EB_Pi0->SetBorderSize(0);
        leg_nxtal_EB_Pi0->SetTextSize(0.05);
        leg_nxtal_EB_Pi0->SetLineColor(1);
        leg_nxtal_EB_Pi0->SetLineStyle(1);
        leg_nxtal_EB_Pi0->SetLineWidth(1);
        leg_nxtal_EB_Pi0->SetFillColor(0);
        leg_nxtal_EB_Pi0->SetFillStyle(1001);
	leg_nxtal_EB_Pi0->AddEntry(h_nxtal_EB_Pi0_G1,"leading photon from EB-#pi^{0}","l");
	leg_nxtal_EB_Pi0->AddEntry(h_nxtal_EB_Pi0_G2,"sub-leading photon from EB-#pi^{0}","l");
	leg_nxtal_EB_Pi0->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Pi0.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Pi0.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Pi0.C");

//
	MaxY = 0.0;
	if(h_nxtal_EE_Pi0_G1->GetMaximum() > MaxY) MaxY = h_nxtal_EE_Pi0_G1->GetMaximum();
	if(h_nxtal_EE_Pi0_G2->GetMaximum() > MaxY) MaxY = h_nxtal_EE_Pi0_G2->GetMaximum();
	
	h_nxtal_EE_Pi0_G1->SetLineColor(2);
	h_nxtal_EE_Pi0_G1->SetLineWidth(2);
	
	h_nxtal_EE_Pi0_G2->SetLineColor(4);
	h_nxtal_EE_Pi0_G2->SetLineWidth(2);
	

	h_nxtal_EE_Pi0_G1->Draw();
	h_nxtal_EE_Pi0_G1->GetYaxis()->SetRangeUser(1.0,1.4*MaxY);
	h_nxtal_EE_Pi0_G1->SetTitle("");
	h_nxtal_EE_Pi0_G1->GetYaxis()->SetTitle("Events");
	h_nxtal_EE_Pi0_G1->GetXaxis()->SetTitle("Nxtal");

	h_nxtal_EE_Pi0_G2->Draw("same");

	TLegend *leg_nxtal_EE_Pi0 = new TLegend(0.2,0.7,0.8,0.89);	
	leg_nxtal_EE_Pi0->SetBorderSize(0);
        leg_nxtal_EE_Pi0->SetTextSize(0.05);
        leg_nxtal_EE_Pi0->SetLineColor(1);
        leg_nxtal_EE_Pi0->SetLineStyle(1);
        leg_nxtal_EE_Pi0->SetLineWidth(1);
        leg_nxtal_EE_Pi0->SetFillColor(0);
        leg_nxtal_EE_Pi0->SetFillStyle(1001);
	leg_nxtal_EE_Pi0->AddEntry(h_nxtal_EE_Pi0_G1,"leading photon from EE-#pi^{0}","l");
	leg_nxtal_EE_Pi0->AddEntry(h_nxtal_EE_Pi0_G2,"sub-leading photon from EE-#pi^{0}","l");
	leg_nxtal_EE_Pi0->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Pi0.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Pi0.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Pi0.C");
//
	MaxY = 0.0;
	if(h_nxtal_EB_Eta_G1->GetMaximum() > MaxY) MaxY = h_nxtal_EB_Eta_G1->GetMaximum();
	if(h_nxtal_EB_Eta_G2->GetMaximum() > MaxY) MaxY = h_nxtal_EB_Eta_G2->GetMaximum();
	
	h_nxtal_EB_Eta_G1->SetLineColor(2);
	h_nxtal_EB_Eta_G1->SetLineWidth(2);
	
	h_nxtal_EB_Eta_G2->SetLineColor(4);
	h_nxtal_EB_Eta_G2->SetLineWidth(2);
	

	h_nxtal_EB_Eta_G1->Draw();
	h_nxtal_EB_Eta_G1->GetYaxis()->SetRangeUser(1.0,1.4*MaxY);
	h_nxtal_EB_Eta_G1->SetTitle("");
	h_nxtal_EB_Eta_G1->GetYaxis()->SetTitle("Events");
	h_nxtal_EB_Eta_G1->GetXaxis()->SetTitle("Nxtal");

	h_nxtal_EB_Eta_G2->Draw("same");

	TLegend *leg_nxtal_EB_Eta = new TLegend(0.2,0.7,0.8,0.89);	
	leg_nxtal_EB_Eta->SetBorderSize(0);
        leg_nxtal_EB_Eta->SetTextSize(0.05);
        leg_nxtal_EB_Eta->SetLineColor(1);
        leg_nxtal_EB_Eta->SetLineStyle(1);
        leg_nxtal_EB_Eta->SetLineWidth(1);
        leg_nxtal_EB_Eta->SetFillColor(0);
        leg_nxtal_EB_Eta->SetFillStyle(1001);
	leg_nxtal_EB_Eta->AddEntry(h_nxtal_EB_Eta_G1,"leading photon from EB-#eta","l");
	leg_nxtal_EB_Eta->AddEntry(h_nxtal_EB_Eta_G2,"sub-leading photon from EB-#eta","l");
	leg_nxtal_EB_Eta->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Eta.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Eta.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EB_Eta.C");

//
	MaxY = 0.0;
	if(h_nxtal_EE_Eta_G1->GetMaximum() > MaxY) MaxY = h_nxtal_EE_Eta_G1->GetMaximum();
	if(h_nxtal_EE_Eta_G2->GetMaximum() > MaxY) MaxY = h_nxtal_EE_Eta_G2->GetMaximum();
	
	h_nxtal_EE_Eta_G1->SetLineColor(2);
	h_nxtal_EE_Eta_G1->SetLineWidth(2);
	
	h_nxtal_EE_Eta_G2->SetLineColor(4);
	h_nxtal_EE_Eta_G2->SetLineWidth(2);
	

	h_nxtal_EE_Eta_G1->Draw();
	h_nxtal_EE_Eta_G1->GetYaxis()->SetRangeUser(1.0,1.4*MaxY);
	h_nxtal_EE_Eta_G1->SetTitle("");
	h_nxtal_EE_Eta_G1->GetYaxis()->SetTitle("Events");
	h_nxtal_EE_Eta_G1->GetXaxis()->SetTitle("Nxtal");

	h_nxtal_EE_Eta_G2->Draw("same");

	TLegend *leg_nxtal_EE_Eta = new TLegend(0.2,0.7,0.8,0.89);	
	leg_nxtal_EE_Eta->SetBorderSize(0);
        leg_nxtal_EE_Eta->SetTextSize(0.05);
        leg_nxtal_EE_Eta->SetLineColor(1);
        leg_nxtal_EE_Eta->SetLineStyle(1);
        leg_nxtal_EE_Eta->SetLineWidth(1);
        leg_nxtal_EE_Eta->SetFillColor(0);
        leg_nxtal_EE_Eta->SetFillStyle(1001);
	leg_nxtal_EE_Eta->AddEntry(h_nxtal_EE_Eta_G1,"leading photon from EE-#eta","l");
	leg_nxtal_EE_Eta->AddEntry(h_nxtal_EE_Eta_G2,"sub-leading photon from EE-#eta","l");
	leg_nxtal_EE_Eta->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Eta.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Eta.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/nxtal_EE_Eta.C");

//
	MaxY = 0.0;
	h_mPi0_rec_EB_withCut->Scale(1.0/h_mPi0_rec_EB_withCut->GetEntries());
	h_mPi0_rec_EB_withoutCut->Scale(1.0/h_mPi0_rec_EB_withoutCut->GetEntries());

	if(h_mPi0_rec_EB_withCut->GetMaximum() > MaxY) MaxY = h_mPi0_rec_EB_withCut->GetMaximum();
	if(h_mPi0_rec_EB_withoutCut->GetMaximum() > MaxY) MaxY = h_mPi0_rec_EB_withoutCut->GetMaximum();
	
	h_mPi0_rec_EB_withCut->SetLineColor(2);
	h_mPi0_rec_EB_withCut->SetLineWidth(2);
	
	h_mPi0_rec_EB_withoutCut->SetLineColor(4);
	h_mPi0_rec_EB_withoutCut->SetLineWidth(2);
	

	h_mPi0_rec_EB_withCut->Draw();
	h_mPi0_rec_EB_withCut->GetYaxis()->SetRangeUser(1.0,0.7*MaxY);
	h_mPi0_rec_EB_withCut->SetTitle("");
	h_mPi0_rec_EB_withCut->GetYaxis()->SetTitle("Events");
	h_mPi0_rec_EB_withCut->GetXaxis()->SetTitle("m_{#gamma#gamma} / GeV");

	h_mPi0_rec_EB_withoutCut->Draw("same");

	TLegend *leg_mPi0_rec_EB = new TLegend(0.2,0.7,0.8,0.89);	
	leg_mPi0_rec_EB->SetBorderSize(0);
        leg_mPi0_rec_EB->SetTextSize(0.05);
        leg_mPi0_rec_EB->SetLineColor(1);
        leg_mPi0_rec_EB->SetLineStyle(1);
        leg_mPi0_rec_EB->SetLineWidth(1);
        leg_mPi0_rec_EB->SetFillColor(0);
        leg_mPi0_rec_EB->SetFillStyle(1001);
	leg_mPi0_rec_EB->AddEntry(h_mPi0_rec_EB_withCut,"with Nxtal>6 cut (EB-#pi^{0})","l");
	leg_mPi0_rec_EB->AddEntry(h_mPi0_rec_EB_withoutCut,"without Nxtal cut (EB-#pi^{0})","l");
	leg_mPi0_rec_EB->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EB.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EB.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EB.C");

//
	MaxY = 0.0;
	h_mPi0_rec_EE_withCut->Scale(1.0/h_mPi0_rec_EE_withCut->GetEntries());
	h_mPi0_rec_EE_withoutCut->Scale(1.0/h_mPi0_rec_EE_withoutCut->GetEntries());

	if(h_mPi0_rec_EE_withCut->GetMaximum() > MaxY) MaxY = h_mPi0_rec_EE_withCut->GetMaximum();
	if(h_mPi0_rec_EE_withoutCut->GetMaximum() > MaxY) MaxY = h_mPi0_rec_EE_withoutCut->GetMaximum();
	
	h_mPi0_rec_EE_withCut->SetLineColor(2);
	h_mPi0_rec_EE_withCut->SetLineWidth(2);
	
	h_mPi0_rec_EE_withoutCut->SetLineColor(4);
	h_mPi0_rec_EE_withoutCut->SetLineWidth(2);
	

	h_mPi0_rec_EE_withCut->Draw();
	h_mPi0_rec_EE_withCut->GetYaxis()->SetRangeUser(1.0,0.7*MaxY);
	h_mPi0_rec_EE_withCut->SetTitle("");
	h_mPi0_rec_EE_withCut->GetYaxis()->SetTitle("Events");
	h_mPi0_rec_EE_withCut->GetXaxis()->SetTitle("m_{#gamma#gamma} / GeV");

	h_mPi0_rec_EE_withoutCut->Draw("same");

	TLegend *leg_mPi0_rec_EE = new TLegend(0.2,0.7,0.8,0.89);	
	leg_mPi0_rec_EE->SetBorderSize(0);
        leg_mPi0_rec_EE->SetTextSize(0.05);
        leg_mPi0_rec_EE->SetLineColor(1);
        leg_mPi0_rec_EE->SetLineStyle(1);
        leg_mPi0_rec_EE->SetLineWidth(1);
        leg_mPi0_rec_EE->SetFillColor(0);
        leg_mPi0_rec_EE->SetFillStyle(1001);
	leg_mPi0_rec_EE->AddEntry(h_mPi0_rec_EE_withCut,"with Nxtal>6 cut (EE-#pi^{0})","l");
	leg_mPi0_rec_EE->AddEntry(h_mPi0_rec_EE_withoutCut,"without Nxtal cut (EE-#pi^{0})","l");
	leg_mPi0_rec_EE->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EE.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EE.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mPi0_rec_EE.C");

//
	MaxY = 0.0;
	h_mEta_rec_EB_withCut->Scale(1.0/h_mEta_rec_EB_withCut->GetEntries());
	h_mEta_rec_EB_withoutCut->Scale(1.0/h_mEta_rec_EB_withoutCut->GetEntries());

	if(h_mEta_rec_EB_withCut->GetMaximum() > MaxY) MaxY = h_mEta_rec_EB_withCut->GetMaximum();
	if(h_mEta_rec_EB_withoutCut->GetMaximum() > MaxY) MaxY = h_mEta_rec_EB_withoutCut->GetMaximum();
	
	h_mEta_rec_EB_withCut->SetLineColor(2);
	h_mEta_rec_EB_withCut->SetLineWidth(2);
	
	h_mEta_rec_EB_withoutCut->SetLineColor(4);
	h_mEta_rec_EB_withoutCut->SetLineWidth(2);
	

	h_mEta_rec_EB_withCut->Draw();
	h_mEta_rec_EB_withCut->GetYaxis()->SetRangeUser(1.0,0.7*MaxY);
	h_mEta_rec_EB_withCut->SetTitle("");
	h_mEta_rec_EB_withCut->GetYaxis()->SetTitle("Events");
	h_mEta_rec_EB_withCut->GetXaxis()->SetTitle("m_{#gamma#gamma} / GeV");

	h_mEta_rec_EB_withoutCut->Draw("same");

	TLegend *leg_mEta_rec_EB = new TLegend(0.2,0.7,0.8,0.89);	
	leg_mEta_rec_EB->SetBorderSize(0);
        leg_mEta_rec_EB->SetTextSize(0.05);
        leg_mEta_rec_EB->SetLineColor(1);
        leg_mEta_rec_EB->SetLineStyle(1);
        leg_mEta_rec_EB->SetLineWidth(1);
        leg_mEta_rec_EB->SetFillColor(0);
        leg_mEta_rec_EB->SetFillStyle(1001);
	leg_mEta_rec_EB->AddEntry(h_mEta_rec_EB_withCut,"with Nxtal>6 cut (EB-#eta)","l");
	leg_mEta_rec_EB->AddEntry(h_mEta_rec_EB_withoutCut,"without Nxtal cut (EB-#eta)","l");
	leg_mEta_rec_EB->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EB.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EB.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EB.C");

//
	MaxY = 0.0;
	h_mEta_rec_EE_withCut->Scale(1.0/h_mEta_rec_EE_withCut->GetEntries());
	h_mEta_rec_EE_withoutCut->Scale(1.0/h_mEta_rec_EE_withoutCut->GetEntries());

	if(h_mEta_rec_EE_withCut->GetMaximum() > MaxY) MaxY = h_mEta_rec_EE_withCut->GetMaximum();
	if(h_mEta_rec_EE_withoutCut->GetMaximum() > MaxY) MaxY = h_mEta_rec_EE_withoutCut->GetMaximum();
	
	h_mEta_rec_EE_withCut->SetLineColor(2);
	h_mEta_rec_EE_withCut->SetLineWidth(2);
	
	h_mEta_rec_EE_withoutCut->SetLineColor(4);
	h_mEta_rec_EE_withoutCut->SetLineWidth(2);
	

	h_mEta_rec_EE_withCut->Draw();
	h_mEta_rec_EE_withCut->GetYaxis()->SetRangeUser(1.0,0.7*MaxY);
	h_mEta_rec_EE_withCut->SetTitle("");
	h_mEta_rec_EE_withCut->GetYaxis()->SetTitle("Events");
	h_mEta_rec_EE_withCut->GetXaxis()->SetTitle("m_{#gamma#gamma} / GeV");

	h_mEta_rec_EE_withoutCut->Draw("same");

	TLegend *leg_mEta_rec_EE = new TLegend(0.2,0.7,0.8,0.89);	
	leg_mEta_rec_EE->SetBorderSize(0);
        leg_mEta_rec_EE->SetTextSize(0.05);
        leg_mEta_rec_EE->SetLineColor(1);
        leg_mEta_rec_EE->SetLineStyle(1);
        leg_mEta_rec_EE->SetLineWidth(1);
        leg_mEta_rec_EE->SetFillColor(0);
        leg_mEta_rec_EE->SetFillStyle(1001);
	leg_mEta_rec_EE->AddEntry(h_mEta_rec_EE_withCut,"with Nxtal>6 cut (EE-#eta)","l");
	leg_mEta_rec_EE->AddEntry(h_mEta_rec_EE_withoutCut,"without Nxtal cut (EE-#eta)","l");
	leg_mEta_rec_EE->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EE.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EE.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/mEta_rec_EE.C");

//
	myC->SetLogy(1);
	MaxY = 0.0;
	//h_ptPi0_rec_EB_withCut->Scale(1.0/h_ptPi0_rec_EB_withCut->GetEntries());
	//h_ptPi0_rec_EB_withoutCut->Scale(1.0/h_ptPi0_rec_EB_withoutCut->GetEntries());

	if(h_ptPi0_rec_EB_withCut->GetMaximum() > MaxY) MaxY = h_ptPi0_rec_EB_withCut->GetMaximum();
	if(h_ptPi0_rec_EB_withoutCut->GetMaximum() > MaxY) MaxY = h_ptPi0_rec_EB_withoutCut->GetMaximum();
	
	h_ptPi0_rec_EB_withCut->SetLineColor(2);
	h_ptPi0_rec_EB_withCut->SetLineWidth(2);
	
	h_ptPi0_rec_EB_withoutCut->SetLineColor(4);
	h_ptPi0_rec_EB_withoutCut->SetLineWidth(2);
	

	h_ptPi0_rec_EB_withCut->Draw();
	h_ptPi0_rec_EB_withCut->GetYaxis()->SetRangeUser(1.0,200.*MaxY);
	h_ptPi0_rec_EB_withCut->SetTitle("");
	h_ptPi0_rec_EB_withCut->GetYaxis()->SetTitle("Events");
	h_ptPi0_rec_EB_withCut->GetXaxis()->SetTitle("pT_{#gamma#gamma} / GeV");

	h_ptPi0_rec_EB_withoutCut->Draw("same");

	TLegend *leg_ptPi0_rec_EB = new TLegend(0.2,0.7,0.8,0.89);	
	leg_ptPi0_rec_EB->SetBorderSize(0);
        leg_ptPi0_rec_EB->SetTextSize(0.05);
        leg_ptPi0_rec_EB->SetLineColor(1);
        leg_ptPi0_rec_EB->SetLineStyle(1);
        leg_ptPi0_rec_EB->SetLineWidth(1);
        leg_ptPi0_rec_EB->SetFillColor(0);
        leg_ptPi0_rec_EB->SetFillStyle(1001);
	leg_ptPi0_rec_EB->AddEntry(h_ptPi0_rec_EB_withCut,"with Nxtal>6 cut (EB-#pi^{0})","l");
	leg_ptPi0_rec_EB->AddEntry(h_ptPi0_rec_EB_withoutCut,"without Nxtal cut (EB-#pi^{0})","l");
	leg_ptPi0_rec_EB->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EB.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EB.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EB.C");

//
	MaxY = 0.0;
	//h_ptPi0_rec_EE_withCut->Scale(1.0/h_ptPi0_rec_EE_withCut->GetEntries());
	//h_ptPi0_rec_EE_withoutCut->Scale(1.0/h_ptPi0_rec_EE_withoutCut->GetEntries());

	if(h_ptPi0_rec_EE_withCut->GetMaximum() > MaxY) MaxY = h_ptPi0_rec_EE_withCut->GetMaximum();
	if(h_ptPi0_rec_EE_withoutCut->GetMaximum() > MaxY) MaxY = h_ptPi0_rec_EE_withoutCut->GetMaximum();
	
	h_ptPi0_rec_EE_withCut->SetLineColor(2);
	h_ptPi0_rec_EE_withCut->SetLineWidth(2);
	
	h_ptPi0_rec_EE_withoutCut->SetLineColor(4);
	h_ptPi0_rec_EE_withoutCut->SetLineWidth(2);
	

	h_ptPi0_rec_EE_withCut->Draw();
	h_ptPi0_rec_EE_withCut->GetYaxis()->SetRangeUser(1.0,200.*MaxY);
	h_ptPi0_rec_EE_withCut->SetTitle("");
	h_ptPi0_rec_EE_withCut->GetYaxis()->SetTitle("Events");
	h_ptPi0_rec_EE_withCut->GetXaxis()->SetTitle("pT_{#gamma#gamma} / GeV");

	h_ptPi0_rec_EE_withoutCut->Draw("same");

	TLegend *leg_ptPi0_rec_EE = new TLegend(0.2,0.7,0.8,0.89);	
	leg_ptPi0_rec_EE->SetBorderSize(0);
        leg_ptPi0_rec_EE->SetTextSize(0.05);
        leg_ptPi0_rec_EE->SetLineColor(1);
        leg_ptPi0_rec_EE->SetLineStyle(1);
        leg_ptPi0_rec_EE->SetLineWidth(1);
        leg_ptPi0_rec_EE->SetFillColor(0);
        leg_ptPi0_rec_EE->SetFillStyle(1001);
	leg_ptPi0_rec_EE->AddEntry(h_ptPi0_rec_EE_withCut,"with Nxtal>6 cut (EE-#pi^{0})","l");
	leg_ptPi0_rec_EE->AddEntry(h_ptPi0_rec_EE_withoutCut,"without Nxtal cut (EE-#pi^{0})","l");
	leg_ptPi0_rec_EE->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EE.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EE.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptPi0_rec_EE.C");

//
	MaxY = 0.0;
	//h_ptEta_rec_EB_withCut->Scale(1.0/h_ptEta_rec_EB_withCut->GetEntries());
	//h_ptEta_rec_EB_withoutCut->Scale(1.0/h_ptEta_rec_EB_withoutCut->GetEntries());

	if(h_ptEta_rec_EB_withCut->GetMaximum() > MaxY) MaxY = h_ptEta_rec_EB_withCut->GetMaximum();
	if(h_ptEta_rec_EB_withoutCut->GetMaximum() > MaxY) MaxY = h_ptEta_rec_EB_withoutCut->GetMaximum();
	
	h_ptEta_rec_EB_withCut->SetLineColor(2);
	h_ptEta_rec_EB_withCut->SetLineWidth(2);
	
	h_ptEta_rec_EB_withoutCut->SetLineColor(4);
	h_ptEta_rec_EB_withoutCut->SetLineWidth(2);
	

	h_ptEta_rec_EB_withCut->Draw();
	h_ptEta_rec_EB_withCut->GetYaxis()->SetRangeUser(1.0,200.*MaxY);
	h_ptEta_rec_EB_withCut->SetTitle("");
	h_ptEta_rec_EB_withCut->GetYaxis()->SetTitle("Events");
	h_ptEta_rec_EB_withCut->GetXaxis()->SetTitle("pT_{#gamma#gamma} / GeV");

	h_ptEta_rec_EB_withoutCut->Draw("same");

	TLegend *leg_ptEta_rec_EB = new TLegend(0.2,0.7,0.8,0.89);	
	leg_ptEta_rec_EB->SetBorderSize(0);
        leg_ptEta_rec_EB->SetTextSize(0.05);
        leg_ptEta_rec_EB->SetLineColor(1);
        leg_ptEta_rec_EB->SetLineStyle(1);
        leg_ptEta_rec_EB->SetLineWidth(1);
        leg_ptEta_rec_EB->SetFillColor(0);
        leg_ptEta_rec_EB->SetFillStyle(1001);
	leg_ptEta_rec_EB->AddEntry(h_ptEta_rec_EB_withCut,"with Nxtal>6 cut (EB-#eta)","l");
	leg_ptEta_rec_EB->AddEntry(h_ptEta_rec_EB_withoutCut,"without Nxtal cut (EB-#eta)","l");
	leg_ptEta_rec_EB->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EB.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EB.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EB.C");

//
	MaxY = 0.0;
	//h_ptEta_rec_EE_withCut->Scale(1.0/h_ptEta_rec_EE_withCut->GetEntries());
	//h_ptEta_rec_EE_withoutCut->Scale(1.0/h_ptEta_rec_EE_withoutCut->GetEntries());

	if(h_ptEta_rec_EE_withCut->GetMaximum() > MaxY) MaxY = h_ptEta_rec_EE_withCut->GetMaximum();
	if(h_ptEta_rec_EE_withoutCut->GetMaximum() > MaxY) MaxY = h_ptEta_rec_EE_withoutCut->GetMaximum();
	
	h_ptEta_rec_EE_withCut->SetLineColor(2);
	h_ptEta_rec_EE_withCut->SetLineWidth(2);
	
	h_ptEta_rec_EE_withoutCut->SetLineColor(4);
	h_ptEta_rec_EE_withoutCut->SetLineWidth(2);
	

	h_ptEta_rec_EE_withCut->Draw();
	h_ptEta_rec_EE_withCut->GetYaxis()->SetRangeUser(1.0,200.*MaxY);
	h_ptEta_rec_EE_withCut->SetTitle("");
	h_ptEta_rec_EE_withCut->GetYaxis()->SetTitle("Events");
	h_ptEta_rec_EE_withCut->GetXaxis()->SetTitle("pT_{#gamma#gamma} / GeV");

	h_ptEta_rec_EE_withoutCut->Draw("same");

	TLegend *leg_ptEta_rec_EE = new TLegend(0.2,0.7,0.8,0.89);	
	leg_ptEta_rec_EE->SetBorderSize(0);
        leg_ptEta_rec_EE->SetTextSize(0.05);
        leg_ptEta_rec_EE->SetLineColor(1);
        leg_ptEta_rec_EE->SetLineStyle(1);
        leg_ptEta_rec_EE->SetLineWidth(1);
        leg_ptEta_rec_EE->SetFillColor(0);
        leg_ptEta_rec_EE->SetFillStyle(1001);
	leg_ptEta_rec_EE->AddEntry(h_ptEta_rec_EE_withCut,"with Nxtal>6 cut (EE-#eta)","l");
	leg_ptEta_rec_EE->AddEntry(h_ptEta_rec_EE_withoutCut,"without Nxtal cut (EE-#eta)","l");
	leg_ptEta_rec_EE->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EE.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EE.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/nxtalCut_QCD/ptEta_rec_EE.C");


}
