

void recHits()
{


	TFile fin("../../../python/pi0Ntuple_hlt.root");
	TTree *tree = (TTree*)fin.Get("ntuples/Pi0Events");

	cout<<tree->GetEntries()<<endl;

	double MaxY = 0.0;

        gStyle->SetOptStat(0);
        gStyle->SetPalette(107) ;
	TCanvas *myC = new TCanvas("myC","myC",100,100,800,700);

	//get histograms from tree
	TH1F *h_N_ebRecHit = new TH1F("h_N_ebRecHit","h_N_ebRecHit",40,0,40);
	tree->Draw("N_ebRecHit>>h_N_ebRecHit");
	
	TH1F *h_N_ebRecHit_Pi0_ = new TH1F("h_N_ebRecHit_Pi0_","h_N_ebRecHit_Pi0_",40,0,40);
	tree->Draw("N_ebRecHit_Pi0_>>h_N_ebRecHit_Pi0_");

	TH1F *h_N_ebRecHit_Eta_ = new TH1F("h_N_ebRecHit_Eta_","h_N_ebRecHit_Eta_",40,0,40);
	tree->Draw("N_ebRecHit_Eta_>>h_N_ebRecHit_Eta_");
	
	TH1F *h_N_eeRecHit = new TH1F("h_N_eeRecHit","h_N_eeRecHit",40,0,40);
	tree->Draw("N_eeRecHit>>h_N_eeRecHit");
	
	TH1F *h_N_eeRecHit_Pi0_ = new TH1F("h_N_eeRecHit_Pi0_","h_N_eeRecHit_Pi0_",40,0,40);
	tree->Draw("N_eeRecHit_Pi0_>>h_N_eeRecHit_Pi0_");

	TH1F *h_N_eeRecHit_Eta_ = new TH1F("h_N_eeRecHit_Eta_","h_N_eeRecHit_Eta_",40,0,40);
	tree->Draw("N_eeRecHit_Eta_>>h_N_eeRecHit_Eta_");
	
	TH1F *h_N_ebeeRecHit = new TH1F("h_N_ebeeRecHit","h_N_ebeeRecHit",40,0,40);
	tree->Draw("N_ebRecHit+N_eeRecHit>>h_N_ebeeRecHit");
	
	TH1F *h_N_ebeeRecHit_Pi0_ = new TH1F("h_N_ebeeRecHit_Pi0_","h_N_ebeeRecHit_Pi0_",40,0,40);
	tree->Draw("N_ebRecHit_Pi0_+N_eeRecHit_Pi0_>>h_N_ebeeRecHit_Pi0_");

	TH1F *h_N_ebeeRecHit_Eta_ = new TH1F("h_N_ebeeRecHit_Eta_","h_N_ebeeRecHit_Eta_",40,0,40);
	tree->Draw("N_ebRecHit_Eta_+N_eeRecHit_Eta_>>h_N_ebeeRecHit_Eta_");


	TH1F *h_N_ebPho_rec = new TH1F("h_N_ebPho_rec","h_N_ebPho_rec",40,0,40);
	tree->Draw("N_ebPho_rec>>h_N_ebPho_rec");
	
	TH1F *h_N_ebPho_rec_Pi0_ = new TH1F("h_N_ebPho_rec_Pi0_","h_N_ebPho_rec_Pi0_",40,0,40);
	tree->Draw("N_ebPho_rec_Pi0_>>h_N_ebPho_rec_Pi0_");

	TH1F *h_N_ebPho_rec_Eta_ = new TH1F("h_N_ebPho_rec_Eta_","h_N_ebPho_rec_Eta_",40,0,40);
	tree->Draw("N_ebPho_rec_Eta_>>h_N_ebPho_rec_Eta_");
	
	TH1F *h_N_eePho_rec = new TH1F("h_N_eePho_rec","h_N_eePho_rec",40,0,40);
	tree->Draw("N_eePho_rec>>h_N_eePho_rec");
	
	TH1F *h_N_eePho_rec_Pi0_ = new TH1F("h_N_eePho_rec_Pi0_","h_N_eePho_rec_Pi0_",40,0,40);
	tree->Draw("N_eePho_rec_Pi0_>>h_N_eePho_rec_Pi0_");

	TH1F *h_N_eePho_rec_Eta_ = new TH1F("h_N_eePho_rec_Eta_","h_N_eePho_rec_Eta_",40,0,40);
	tree->Draw("N_eePho_rec_Eta_>>h_N_eePho_rec_Eta_");
	
	TH1F *h_N_ebeePho_rec = new TH1F("h_N_ebeePho_rec","h_N_ebeePho_rec",40,0,40);
	tree->Draw("N_Pho_rec>>h_N_ebeePho_rec");
	
	TH1F *h_N_ebeePho_rec_Pi0_ = new TH1F("h_N_ebeePho_rec_Pi0_","h_N_ebeePho_rec_Pi0_",40,0,40);
	tree->Draw("N_Pho_rec_Pi0_>>h_N_ebeePho_rec_Pi0_");

	TH1F *h_N_ebeePho_rec_Eta_ = new TH1F("h_N_ebeePho_rec_Eta_","h_N_ebeePho_rec_Eta_",40,0,40);
	tree->Draw("N_Pho_rec_Eta_>>h_N_ebeePho_rec_Eta_");
	
	
	TH1F *h_N_ebeePair_rec = new TH1F("h_N_ebeePair_rec","h_N_ebeePair_rec",10,0,10);
	tree->Draw("N_Pair_rec>>h_N_ebeePair_rec");

	TH1F *h_N_ebeePi0_rec = new TH1F("h_N_ebeePi0_rec","h_N_ebeePi0_rec",10,0,10);
	tree->Draw("N_Pi0_rec>>h_N_ebeePi0_rec");

	TH1F *h_N_ebeeEta_rec = new TH1F("h_N_ebeeEta_rec","h_N_ebeeEta_rec",10,0,10);
	tree->Draw("N_Eta_rec>>h_N_ebeeEta_rec");

	TH1F *h_N_ebPair_rec = new TH1F("h_N_ebPair_rec","h_N_ebPair_rec",10,0,10);
	tree->Draw("N_ebPair_rec>>h_N_ebPair_rec");

	TH1F *h_N_ebPi0_rec = new TH1F("h_N_ebPi0_rec","h_N_ebPi0_rec",10,0,10);
	tree->Draw("N_ebPi0_rec>>h_N_ebPi0_rec");

	TH1F *h_N_ebEta_rec = new TH1F("h_N_ebEta_rec","h_N_ebEta_rec",10,0,10);
	tree->Draw("N_ebEta_rec>>h_N_ebEta_rec");

	TH1F *h_N_eePair_rec = new TH1F("h_N_eePair_rec","h_N_eePair_rec",10,0,10);
	tree->Draw("N_eePair_rec>>h_N_eePair_rec");

	TH1F *h_N_eePi0_rec = new TH1F("h_N_eePi0_rec","h_N_eePi0_rec",10,0,10);
	tree->Draw("N_eePi0_rec>>h_N_eePi0_rec");

	TH1F *h_N_eeEta_rec = new TH1F("h_N_eeEta_rec","h_N_eeEta_rec",10,0,10);
	tree->Draw("N_eeEta_rec>>h_N_eeEta_rec");


	//draw to plots
	//myC->SetLogy(1);

	MaxY = 0.0;
	if(h_N_ebRecHit->GetMaximum() > MaxY) MaxY = h_N_ebRecHit->GetMaximum();
	if(h_N_eeRecHit->GetMaximum() > MaxY) MaxY = h_N_eeRecHit->GetMaximum();
	if(h_N_ebeeRecHit->GetMaximum() > MaxY) MaxY = h_N_ebeeRecHit->GetMaximum();
	
	h_N_ebRecHit->SetLineColor(2);
	h_N_eeRecHit->SetLineColor(3);
	h_N_ebeeRecHit->SetLineColor(4);
	
	h_N_ebRecHit->SetLineWidth(2);
	h_N_eeRecHit->SetLineWidth(2);
	h_N_ebeeRecHit->SetLineWidth(2);
	

	h_N_ebRecHit->Draw();
	h_N_ebRecHit->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebRecHit->SetTitle("");
	h_N_ebRecHit->GetYaxis()->SetTitle("Events");
	h_N_ebRecHit->GetXaxis()->SetTitle("#recHits per event");

	h_N_eeRecHit->Draw("same");
	h_N_ebeeRecHit->Draw("same");

	TLegend *leg_RecHit = new TLegend(0.4,0.7,0.8,0.89);	
	leg_RecHit->SetBorderSize(0);
        leg_RecHit->SetTextSize(0.04);
        leg_RecHit->SetLineColor(1);
        leg_RecHit->SetLineStyle(1);
        leg_RecHit->SetLineWidth(1);
        leg_RecHit->SetFillColor(0);
        leg_RecHit->SetFillStyle(1001);
	leg_RecHit->AddEntry(h_N_ebRecHit,"EB RecHit","l");
	leg_RecHit->AddEntry(h_N_eeRecHit,"EE RecHit","l");
	leg_RecHit->AddEntry(h_N_ebeeRecHit,"EB + EE RecHit","l");
	leg_RecHit->Draw();
	
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_all.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_all.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_all.C");
	//recHit Pi0
	MaxY = 0.0;
	if(h_N_ebRecHit_Pi0_->GetMaximum() > MaxY) MaxY = h_N_ebRecHit_Pi0_->GetMaximum();
	if(h_N_eeRecHit_Pi0_->GetMaximum() > MaxY) MaxY = h_N_eeRecHit_Pi0_->GetMaximum();
	if(h_N_ebeeRecHit_Pi0_->GetMaximum() > MaxY) MaxY = h_N_ebeeRecHit_Pi0_->GetMaximum();
	
	h_N_ebRecHit_Pi0_->SetLineColor(2);
	h_N_eeRecHit_Pi0_->SetLineColor(3);
	h_N_ebeeRecHit_Pi0_->SetLineColor(4);
	
	h_N_ebRecHit_Pi0_->SetLineWidth(2);
	h_N_eeRecHit_Pi0_->SetLineWidth(2);
	h_N_ebeeRecHit_Pi0_->SetLineWidth(2);
	

	h_N_ebRecHit_Pi0_->Draw();
	h_N_ebRecHit_Pi0_->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebRecHit_Pi0_->SetTitle("");
	h_N_ebRecHit_Pi0_->GetYaxis()->SetTitle("Events");
	h_N_ebRecHit_Pi0_->GetXaxis()->SetTitle("#recHits per event");

	h_N_eeRecHit_Pi0_->Draw("same");
	h_N_ebeeRecHit_Pi0_->Draw("same");

	TLegend *leg_RecHit_Pi0_ = new TLegend(0.4,0.7,0.8,0.89);	
	leg_RecHit_Pi0_->SetBorderSize(0);
        leg_RecHit_Pi0_->SetTextSize(0.04);
        leg_RecHit_Pi0_->SetLineColor(1);
        leg_RecHit_Pi0_->SetLineStyle(1);
        leg_RecHit_Pi0_->SetLineWidth(1);
        leg_RecHit_Pi0_->SetFillColor(0);
        leg_RecHit_Pi0_->SetFillStyle(1001);
	leg_RecHit_Pi0_->AddEntry(h_N_ebRecHit_Pi0_,"EB RecHit (Pi0)","l");
	leg_RecHit_Pi0_->AddEntry(h_N_eeRecHit_Pi0_,"EE RecHit (Pi0)","l");
	leg_RecHit_Pi0_->AddEntry(h_N_ebeeRecHit_Pi0_,"EB + EE RecHit (Pi0)","l");
	leg_RecHit_Pi0_->Draw();
	
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Pi0_.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Pi0_.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Pi0_.C");

	//recHit Eta
	MaxY = 0.0;
	if(h_N_ebRecHit_Eta_->GetMaximum() > MaxY) MaxY = h_N_ebRecHit_Eta_->GetMaximum();
	if(h_N_eeRecHit_Eta_->GetMaximum() > MaxY) MaxY = h_N_eeRecHit_Eta_->GetMaximum();
	if(h_N_ebeeRecHit_Eta_->GetMaximum() > MaxY) MaxY = h_N_ebeeRecHit_Eta_->GetMaximum();
	
	h_N_ebRecHit_Eta_->SetLineColor(2);
	h_N_eeRecHit_Eta_->SetLineColor(3);
	h_N_ebeeRecHit_Eta_->SetLineColor(4);
	
	h_N_ebRecHit_Eta_->SetLineWidth(2);
	h_N_eeRecHit_Eta_->SetLineWidth(2);
	h_N_ebeeRecHit_Eta_->SetLineWidth(2);
	

	h_N_ebRecHit_Eta_->Draw();
	h_N_ebRecHit_Eta_->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebRecHit_Eta_->SetTitle("");
	h_N_ebRecHit_Eta_->GetYaxis()->SetTitle("Events");
	h_N_ebRecHit_Eta_->GetXaxis()->SetTitle("#recHits per event");

	h_N_eeRecHit_Eta_->Draw("same");
	h_N_ebeeRecHit_Eta_->Draw("same");

	TLegend *leg_RecHit_Eta_ = new TLegend(0.4,0.7,0.8,0.89);	
	leg_RecHit_Eta_->SetBorderSize(0);
        leg_RecHit_Eta_->SetTextSize(0.04);
        leg_RecHit_Eta_->SetLineColor(1);
        leg_RecHit_Eta_->SetLineStyle(1);
        leg_RecHit_Eta_->SetLineWidth(1);
        leg_RecHit_Eta_->SetFillColor(0);
        leg_RecHit_Eta_->SetFillStyle(1001);
	leg_RecHit_Eta_->AddEntry(h_N_ebRecHit_Eta_,"EB RecHit (Eta)","l");
	leg_RecHit_Eta_->AddEntry(h_N_eeRecHit_Eta_,"EE RecHit (Eta)","l");
	leg_RecHit_Eta_->AddEntry(h_N_ebeeRecHit_Eta_,"EB + EE RecHit (Eta)","l");
	leg_RecHit_Eta_->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Eta_.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Eta_.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hRecHit_Eta_.C");



	MaxY = 0.0;
	if(h_N_ebPho_rec->GetMaximum() > MaxY) MaxY = h_N_ebPho_rec->GetMaximum();
	if(h_N_eePho_rec->GetMaximum() > MaxY) MaxY = h_N_eePho_rec->GetMaximum();
	if(h_N_ebeePho_rec->GetMaximum() > MaxY) MaxY = h_N_ebeePho_rec->GetMaximum();
	
	h_N_ebPho_rec->SetLineColor(2);
	h_N_eePho_rec->SetLineColor(3);
	h_N_ebeePho_rec->SetLineColor(4);
	
	h_N_ebPho_rec->SetLineWidth(2);
	h_N_eePho_rec->SetLineWidth(2);
	h_N_ebeePho_rec->SetLineWidth(2);
	

	h_N_ebPho_rec->Draw();
	h_N_ebPho_rec->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebPho_rec->SetTitle("");
	h_N_ebPho_rec->GetYaxis()->SetTitle("Events");
	h_N_ebPho_rec->GetXaxis()->SetTitle("#Photons per event");

	h_N_eePho_rec->Draw("same");
	h_N_ebeePho_rec->Draw("same");

	TLegend *leg_Pho_rec = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Pho_rec->SetBorderSize(0);
        leg_Pho_rec->SetTextSize(0.04);
        leg_Pho_rec->SetLineColor(1);
        leg_Pho_rec->SetLineStyle(1);
        leg_Pho_rec->SetLineWidth(1);
        leg_Pho_rec->SetFillColor(0);
        leg_Pho_rec->SetFillStyle(1001);
	leg_Pho_rec->AddEntry(h_N_ebPho_rec,"EB Pho_rec","l");
	leg_Pho_rec->AddEntry(h_N_eePho_rec,"EE Pho_rec","l");
	leg_Pho_rec->AddEntry(h_N_ebeePho_rec,"EB + EE Pho_rec","l");
	leg_Pho_rec->Draw();
	
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_all.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_all.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_all.C");
	//recHit Pi0
	MaxY = 0.0;
	if(h_N_ebPho_rec_Pi0_->GetMaximum() > MaxY) MaxY = h_N_ebPho_rec_Pi0_->GetMaximum();
	if(h_N_eePho_rec_Pi0_->GetMaximum() > MaxY) MaxY = h_N_eePho_rec_Pi0_->GetMaximum();
	if(h_N_ebeePho_rec_Pi0_->GetMaximum() > MaxY) MaxY = h_N_ebeePho_rec_Pi0_->GetMaximum();
	
	h_N_ebPho_rec_Pi0_->SetLineColor(2);
	h_N_eePho_rec_Pi0_->SetLineColor(3);
	h_N_ebeePho_rec_Pi0_->SetLineColor(4);
	
	h_N_ebPho_rec_Pi0_->SetLineWidth(2);
	h_N_eePho_rec_Pi0_->SetLineWidth(2);
	h_N_ebeePho_rec_Pi0_->SetLineWidth(2);
	

	h_N_ebPho_rec_Pi0_->Draw();
	h_N_ebPho_rec_Pi0_->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebPho_rec_Pi0_->SetTitle("");
	h_N_ebPho_rec_Pi0_->GetYaxis()->SetTitle("Events");
	h_N_ebPho_rec_Pi0_->GetXaxis()->SetTitle("#Photons per event");

	h_N_eePho_rec_Pi0_->Draw("same");
	h_N_ebeePho_rec_Pi0_->Draw("same");

	TLegend *leg_Pho_rec_Pi0_ = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Pho_rec_Pi0_->SetBorderSize(0);
        leg_Pho_rec_Pi0_->SetTextSize(0.04);
        leg_Pho_rec_Pi0_->SetLineColor(1);
        leg_Pho_rec_Pi0_->SetLineStyle(1);
        leg_Pho_rec_Pi0_->SetLineWidth(1);
        leg_Pho_rec_Pi0_->SetFillColor(0);
        leg_Pho_rec_Pi0_->SetFillStyle(1001);
	leg_Pho_rec_Pi0_->AddEntry(h_N_ebPho_rec_Pi0_,"EB Pho_rec (Pi0)","l");
	leg_Pho_rec_Pi0_->AddEntry(h_N_eePho_rec_Pi0_,"EE Pho_rec (Pi0)","l");
	leg_Pho_rec_Pi0_->AddEntry(h_N_ebeePho_rec_Pi0_,"EB + EE Pho_rec (Pi0)","l");
	leg_Pho_rec_Pi0_->Draw();
	
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Pi0_.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Pi0_.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Pi0_.C");

	//recHit Eta
	MaxY = 0.0;
	if(h_N_ebPho_rec_Eta_->GetMaximum() > MaxY) MaxY = h_N_ebPho_rec_Eta_->GetMaximum();
	if(h_N_eePho_rec_Eta_->GetMaximum() > MaxY) MaxY = h_N_eePho_rec_Eta_->GetMaximum();
	if(h_N_ebeePho_rec_Eta_->GetMaximum() > MaxY) MaxY = h_N_ebeePho_rec_Eta_->GetMaximum();
	
	h_N_ebPho_rec_Eta_->SetLineColor(2);
	h_N_eePho_rec_Eta_->SetLineColor(3);
	h_N_ebeePho_rec_Eta_->SetLineColor(4);
	
	h_N_ebPho_rec_Eta_->SetLineWidth(2);
	h_N_eePho_rec_Eta_->SetLineWidth(2);
	h_N_ebeePho_rec_Eta_->SetLineWidth(2);
	

	h_N_ebPho_rec_Eta_->Draw();
	h_N_ebPho_rec_Eta_->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebPho_rec_Eta_->SetTitle("");
	h_N_ebPho_rec_Eta_->GetYaxis()->SetTitle("Events");
	h_N_ebPho_rec_Eta_->GetXaxis()->SetTitle("#Photons per event");

	h_N_eePho_rec_Eta_->Draw("same");
	h_N_ebeePho_rec_Eta_->Draw("same");

	TLegend *leg_Pho_rec_Eta_ = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Pho_rec_Eta_->SetBorderSize(0);
        leg_Pho_rec_Eta_->SetTextSize(0.04);
        leg_Pho_rec_Eta_->SetLineColor(1);
        leg_Pho_rec_Eta_->SetLineStyle(1);
        leg_Pho_rec_Eta_->SetLineWidth(1);
        leg_Pho_rec_Eta_->SetFillColor(0);
        leg_Pho_rec_Eta_->SetFillStyle(1001);
	leg_Pho_rec_Eta_->AddEntry(h_N_ebPho_rec_Eta_,"EB Pho_rec (Eta)","l");
	leg_Pho_rec_Eta_->AddEntry(h_N_eePho_rec_Eta_,"EE Pho_rec (Eta)","l");
	leg_Pho_rec_Eta_->AddEntry(h_N_ebeePho_rec_Eta_,"EB + EE Pho_rec (Eta)","l");
	leg_Pho_rec_Eta_->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Eta_.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Eta_.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hPho_rec_Eta_.C");






	//Pi0 num
	MaxY = 0.0;
	if(h_N_ebPi0_rec->GetMaximum() > MaxY) MaxY = h_N_ebPi0_rec->GetMaximum();
	if(h_N_eePi0_rec->GetMaximum() > MaxY) MaxY = h_N_eePi0_rec->GetMaximum();
	if(h_N_ebeePi0_rec->GetMaximum() > MaxY) MaxY = h_N_ebeePi0_rec->GetMaximum();
	
	h_N_ebPi0_rec->SetLineColor(2);
	h_N_eePi0_rec->SetLineColor(3);
	h_N_ebeePi0_rec->SetLineColor(4);
	
	h_N_ebPi0_rec->SetLineWidth(2);
	h_N_eePi0_rec->SetLineWidth(2);
	h_N_ebeePi0_rec->SetLineWidth(2);
	

	h_N_ebPi0_rec->Draw();
	h_N_ebPi0_rec->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebPi0_rec->SetTitle("");
	h_N_ebPi0_rec->GetYaxis()->SetTitle("Events");
	h_N_ebPi0_rec->GetXaxis()->SetTitle("#of diphoton pairs");

	h_N_eePi0_rec->Draw("same");
	h_N_ebeePi0_rec->Draw("same");

	TLegend *leg_Pi0_num = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Pi0_num->SetBorderSize(0);
        leg_Pi0_num->SetTextSize(0.04);
        leg_Pi0_num->SetLineColor(1);
        leg_Pi0_num->SetLineStyle(1);
        leg_Pi0_num->SetLineWidth(1);
        leg_Pi0_num->SetFillColor(0);
        leg_Pi0_num->SetFillStyle(1001);
	leg_Pi0_num->AddEntry(h_N_ebPi0_rec,"EB Pi0","l");
	leg_Pi0_num->AddEntry(h_N_eePi0_rec,"EE Pi0","l");
	leg_Pi0_num->AddEntry(h_N_ebeePi0_rec,"EB + EE Pi0","l");
	leg_Pi0_num->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pi0_rec.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pi0_rec.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pi0_rec.C");

	//Eta num
	MaxY = 0.0;
	if(h_N_ebEta_rec->GetMaximum() > MaxY) MaxY = h_N_ebEta_rec->GetMaximum();
	if(h_N_eeEta_rec->GetMaximum() > MaxY) MaxY = h_N_eeEta_rec->GetMaximum();
	if(h_N_ebeeEta_rec->GetMaximum() > MaxY) MaxY = h_N_ebeeEta_rec->GetMaximum();
	
	h_N_ebEta_rec->SetLineColor(2);
	h_N_eeEta_rec->SetLineColor(3);
	h_N_ebeeEta_rec->SetLineColor(4);
	
	h_N_ebEta_rec->SetLineWidth(2);
	h_N_eeEta_rec->SetLineWidth(2);
	h_N_ebeeEta_rec->SetLineWidth(2);
	

	h_N_ebEta_rec->Draw();
	h_N_ebEta_rec->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebEta_rec->SetTitle("");
	h_N_ebEta_rec->GetYaxis()->SetTitle("Events");
	h_N_ebEta_rec->GetXaxis()->SetTitle("#of diphoton pairs");

	h_N_eeEta_rec->Draw("same");
	h_N_ebeeEta_rec->Draw("same");

	TLegend *leg_Eta_num = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Eta_num->SetBorderSize(0);
        leg_Eta_num->SetTextSize(0.04);
        leg_Eta_num->SetLineColor(1);
        leg_Eta_num->SetLineStyle(1);
        leg_Eta_num->SetLineWidth(1);
        leg_Eta_num->SetFillColor(0);
        leg_Eta_num->SetFillStyle(1001);
	leg_Eta_num->AddEntry(h_N_ebEta_rec,"EB Eta","l");
	leg_Eta_num->AddEntry(h_N_eeEta_rec,"EE Eta","l");
	leg_Eta_num->AddEntry(h_N_ebeeEta_rec,"EB + EE Eta","l");
	leg_Eta_num->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Eta_rec.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Eta_rec.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Eta_rec.C");

	//Pair num
	MaxY = 0.0;
	if(h_N_ebPair_rec->GetMaximum() > MaxY) MaxY = h_N_ebPair_rec->GetMaximum();
	if(h_N_eePair_rec->GetMaximum() > MaxY) MaxY = h_N_eePair_rec->GetMaximum();
	if(h_N_ebeePair_rec->GetMaximum() > MaxY) MaxY = h_N_ebeePair_rec->GetMaximum();
	
	h_N_ebPair_rec->SetLineColor(2);
	h_N_eePair_rec->SetLineColor(3);
	h_N_ebeePair_rec->SetLineColor(4);
	
	h_N_ebPair_rec->SetLineWidth(2);
	h_N_eePair_rec->SetLineWidth(2);
	h_N_ebeePair_rec->SetLineWidth(2);
	

	h_N_ebPair_rec->Draw();
	h_N_ebPair_rec->GetYaxis()->SetRangeUser(1.0,1.1*MaxY);
	h_N_ebPair_rec->SetTitle("");
	h_N_ebPair_rec->GetYaxis()->SetTitle("Events");
	h_N_ebPair_rec->GetXaxis()->SetTitle("#of diphoton pairs");

	h_N_eePair_rec->Draw("same");
	h_N_ebeePair_rec->Draw("same");

	TLegend *leg_Pair_num = new TLegend(0.4,0.7,0.8,0.89);	
	leg_Pair_num->SetBorderSize(0);
        leg_Pair_num->SetTextSize(0.04);
        leg_Pair_num->SetLineColor(1);
        leg_Pair_num->SetLineStyle(1);
        leg_Pair_num->SetLineWidth(1);
        leg_Pair_num->SetFillColor(0);
        leg_Pair_num->SetFillStyle(1001);
	leg_Pair_num->AddEntry(h_N_ebPair_rec,"EB Pi0+Eta","l");
	leg_Pair_num->AddEntry(h_N_eePair_rec,"EE Pi0+Eta","l");
	leg_Pair_num->AddEntry(h_N_ebeePair_rec,"EB + EE Pi0+Eta","l");
	leg_Pair_num->Draw();
	
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pair_rec.pdf");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pair_rec.png");
	myC->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/recHits_online/hN_Pair_rec.C");


}
