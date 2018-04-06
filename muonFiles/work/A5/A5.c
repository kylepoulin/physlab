

double the1[5] = {5.2,6.3,6.6,9.0,12.2};
double the2[5] = {6.6,6.3,6.2,5.9,5.5};
double data[5] = {8.0,7.0,4.0,10.0,11.0};

double chi1=0.0;
double chi2=0.0;

void A5(){

	gRandom->SetSeed(12432514);

	for(int i=0; i<5; i++){
		chi1 += pow((data[i]-the1[i]),2)/the1[i];
		chi2 += pow((data[i]-the2[i]),2)/the2[i];
	}
	cout << chi1 << " and " << chi2 << endl;

	double pval1 = TMath::Prob(chi1,5);
	double pval2 = TMath::Prob(chi2,5);
	
	cout << pval1 << " and " << pval2 << endl;
	
	
	TH1D *hchi1 = new TH1D("Croissant","Chi squared WRT Theory1", 30, -1, 30);
	TH1D *hchi2 = new TH1D("Baguette","Chi squared WRT Theory2", 30, -1, 60);
	TH1D *hpval1 = new TH1D("Danish","P value WRT Theory1", 30, -1, 2);
	TH1D *hpval2 = new TH1D("Croissant2","P value WRT Theory2", 30, -1, 2);

	chi1 = 0.0;
	chi2 = 0.0;
	for(int trials=0; trials<10000; trials++){

		chi1 = 0.0;
		chi2 = 0.0;

		double pseudo[5];
		for(int j=0; j<5; j++){
			pseudo[j] = gRandom->Poisson(the1[j]);
		}

		for(int i=0; i<5; i++){
			chi1 += pow((pseudo[i]-the1[i]),2)/the1[i];
			chi2 += pow((pseudo[i]-the2[i]),2)/the2[i];
		}
		
		hchi1->Fill(chi1);
		hchi2->Fill(chi2);	

		pval1 = TMath::Prob(chi1,5);	
		pval2 = TMath::Prob(chi2,5);

		hpval1->Fill(pval1);
		hpval2->Fill(pval2);
	}


	TCanvas *can1 = new TCanvas("C1");
	can1->Divide(2,2);
	can1->cd(1);
	hchi1->Draw();
	can1->cd(2);
	hchi2->Draw();
	can1->cd(3);
	hpval1->Draw();
	can1->cd(4);
	hpval2->Draw();

}
//
