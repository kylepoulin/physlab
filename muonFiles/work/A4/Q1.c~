#include <stdio.h>

TF1 *drawGauss(double mu, double sigma, int col=kRed) {
  TF1 *gaus = new TF1("","TMath::Gaus(x,[0],[1],1)",-100,1000);
  gaus->SetLineColor(col);
  gaus->SetParameters(mu,sigma);
  gaus->Draw("same");
  return gaus;
}

void Q1(){
	//Setup
	double norm = 1;
	double var;
	double sum = 0.0;
	gRandom->SetSeed(12345678);

	
	//Standard Uniform Setup
	TH1D *h_a = new TH1D("SU1","Sum(1): Standard Uniform(0,1)",50,-1,2);
	TH1D *h_a2 = new TH1D("SU2","Sum(2): Standard Uniform(0,1)",50,-1,3);
	TH1D *h_a3 = new TH1D("SU3","Sum(5): Standard Uniform(0,1)",75,-1,6);
	TH1D *h_a4 = new TH1D("SU4","Sum(100): Standard Uniform(0,1)",100,40,60);
	//Poisson Setup
	TH1D *h_a5 = new TH1D("S1","Sum(1): Poisson(1)",50,-2,10);
	TH1D *h_a6 = new TH1D("S2","Sum(2): Poisson(1)",50,-4,16);
	TH1D *h_a7 = new TH1D("S3","Sum(5): Poisson(1)",75,0,20);
	TH1D *h_a8 = new TH1D("S4","Sum(100): Poisson(1)",100,40,160);
	//Exponential Setup
	TH1D *h_a9 = new TH1D("E1","Sum(1): Exponential(1)",50,-4,8);
	TH1D *h_a10 = new TH1D("E2","Sum(2): Exponential(1)",50,-6,15);
	TH1D *h_a11 = new TH1D("E3","Sum(5): Exponential(1)",75,-2,20);
	TH1D *h_a12 = new TH1D("E4","Sum(100): Exponential(1)",100,60,150);
	//Standard Decimal Setup
	TH1D *h_a13 = new TH1D("U1","Sum(1): Uniform(0.4,0.6)",50,0.3,0.7);
	TH1D *h_a14 = new TH1D("U2","Sum(2): Uniform(0.4,0.6)",50,0.5,1.5);
	TH1D *h_a15 = new TH1D("U3","Sum(5): Uniform(0.4,0.6)",75,1.5,3.1);
	TH1D *h_a16 = new TH1D("U4","Sum(100): Uniform(0.4,0.6)",100,46,54);	

	//The meat of it

//Standard Uniform
	TCanvas *can1 = new TCanvas("C1");
	can1->Divide(2,2);
	can1->cd(1);
	for(int i=0; i<10000; i++){
		var = gRandom->Uniform(0,1);
		h_a->Fill(var);
		sum += pow((var-0.5),2);
	}
	//cout << sum << endl;
	sum = sqrt(sum*(1.0/9999));
	h_a->Scale(norm/h_a->Integral("width")); //important
	h_a->SetMaximum(1.6);
	h_a->GetXaxis()->SetTitle("Random Variable X");
	h_a->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a->Draw();
	//cout << sum << endl;
	drawGauss(0.5,sum);

	can1->cd(2);
	sum = 0;
	for(int i=0; i<10000; i++){
		var = gRandom->Uniform(0,1)+gRandom->Uniform(0,1);
		h_a2->Fill(var);
		sum += pow((var-1),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a2->Scale(norm/h_a2->Integral("width"));
	h_a2->SetMaximum(1.2);
	h_a2->GetXaxis()->SetTitle("Random Variable X");
	h_a2->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a2->Draw();
	drawGauss(1.0,sum);

	can1->cd(3);
	sum = 0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<5; j++)
			d+=gRandom->Uniform(0,1);
		h_a3->Fill(d);
		sum += pow((d-2.5),2);
	}
	sum = sqrt(sum*(1.0/9999));
	//cout << sum << endl;
	h_a3->Scale(norm/h_a3->Integral("width"));
	h_a3->SetMaximum(1.0);
	h_a3->GetXaxis()->SetTitle("Random Variable X");
	h_a3->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a3->Draw();
	drawGauss(2.5,sum);

	can1->cd(4);
	sum = 0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<100; j++)
			d+=gRandom->Uniform(0,1);
		h_a4->Fill(d);
		sum += pow((d-50),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a4->Scale(norm/h_a4->Integral("width"));
	h_a4->GetXaxis()->SetTitle("Random Variable X");
	h_a4->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a4->Draw();
	drawGauss(50,sum);

//Poisson

	TCanvas *can2 = new TCanvas("C2");
	sum =0;
	can2->Divide(2,2);
	can2->cd(1);
	for(int i=0; i<10000; i++){
		var = gRandom->Poisson(1);
		h_a5->Fill(var);
		sum+=pow((var-1),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a5->Scale(norm/h_a5->Integral("width"));
	h_a5->GetXaxis()->SetTitle("Random Variable X");
	h_a5->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a5->Draw();
	drawGauss(1,sum);
	

	can2->cd(2);
	sum=0;
	for(int i=0; i<10000; i++){
		var = gRandom->Poisson(1)+gRandom->Poisson(1);
		h_a6->Fill(var);
		sum+=pow((var-2),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a6->Scale(norm/h_a6->Integral("width"));
	h_a6->GetXaxis()->SetTitle("Random Variable X");
	h_a6->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a6->Draw();
	drawGauss(2,sum);

	can2->cd(3);
	sum=0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<5; j++)
			d+=gRandom->Poisson(1);
		h_a7->Fill(d);
		sum+=pow((d-5),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a7->Scale(norm/h_a7->Integral("width"));
	h_a7->GetXaxis()->SetTitle("Random Variable X");
	h_a7->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a7->Draw();
	drawGauss(5,sum);

	can2->cd(4);
	sum=0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<100; j++)
			d+=gRandom->Poisson(1);
		h_a8->Fill(d);
		sum+=pow((d-100),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a8->Scale(norm/h_a8->Integral("width"));
	h_a8->GetXaxis()->SetTitle("Random Variable X");
	h_a8->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a8->Draw();
	drawGauss(100,sum);

//Exponential
	TCanvas *can3 = new TCanvas("C3");
	can3->Divide(2,2);
	can3->cd(1);
	sum=0;
	for(int i=0; i<10000; i++){
		var = gRandom->Exp(1);
		h_a9->Fill(var);
		sum+=pow((var-1),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a9->Scale(norm/h_a9->Integral("width"));
	h_a9->GetXaxis()->SetTitle("Random Variable X");
	h_a9->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a9->Draw();
	drawGauss(1,sum);

	can3->cd(2);
	sum=0;
	for(int i=0; i<10000; i++){
		var=gRandom->Exp(1)+gRandom->Exp(1);
		h_a10->Fill(var);
		sum+=pow((var-2),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a10->Scale(norm/h_a10->Integral("width"));
	h_a10->GetXaxis()->SetTitle("Random Variable X");
	h_a10->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a10->Draw();
	drawGauss(2,sum);

	can3->cd(3);
	sum=0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<5; j++)
			d+=gRandom->Exp(1);
		h_a11->Fill(d);
		sum+=pow((d-5),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a11->Scale(norm/h_a11->Integral("width"));
	h_a11->GetXaxis()->SetTitle("Random Variable X");
	h_a11->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a11->Draw();
	drawGauss(5,sum);

	can3->cd(4);
	sum=0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<100; j++)
			d+=gRandom->Exp(1);
		h_a12->Fill(d);
		sum+=pow((d-100),2);
	}
	sum = sqrt(sum*(1.0/9999));
	h_a12->Scale(norm/h_a12->Integral("width"));
	h_a12->GetXaxis()->SetTitle("Random Variable X");
	h_a12->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a12->Draw();
	drawGauss(100,sum);

//Uniform Decimal
	TCanvas *can4 = new TCanvas("C4");
	can4->Divide(2,2);
	can4->cd(1);
	sum = 0;
	for(int i=0; i<10000; i++){
		var = gRandom->Uniform(0.4,0.6);
		h_a13->Fill(var);
		sum += pow(var-0.5,2);
	}
	h_a13->Scale(norm/h_a13->Integral("width"));
	sum = sqrt(sum*(1.0/9999));
	h_a13->SetMaximum(8);
	h_a13->GetXaxis()->SetTitle("Random Variable X");
	h_a13->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a13->Draw();
	drawGauss(0.5,sum);

	can4->cd(2);
	sum=0;
	for(int i=0; i<10000; i++){
		var = gRandom->Uniform(0.4,0.6)+gRandom->Uniform(0.4,0.6);
		h_a14->Fill(var);
		sum += pow(var-1,2);
	}
	h_a14->Scale(norm/h_a14->Integral("width"));
	sum = sqrt(sum*(1.0/9999));
	h_a14->GetXaxis()->SetTitle("Random Variable X");
	h_a14->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a14->Draw();
	drawGauss(1,sum);

	can4->cd(3);
	sum = 0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<5; j++)
			d+=gRandom->Uniform(0.4,0.6);
		h_a15->Fill(d);
		sum += pow(d-2.5,2);
	}
	h_a15->Scale(norm/h_a15->Integral("width"));
	sum = sqrt(sum*(1.0/9999));
	h_a15->GetXaxis()->SetTitle("Random Variable X");
	h_a15->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a15->Draw();
	drawGauss(2.5,sum);

	can4->cd(4);
	sum=0;
	for(int i=0; i<10000; i++){
		double d = 0;
		for(int j=0; j<100; j++)
			d+=gRandom->Uniform(0.4,0.6);
		h_a16->Fill(d);
		sum += pow(d-50,2);
	}
	h_a16->Scale(norm/h_a16->Integral("width"));
	sum = sqrt(sum*(1.0/9999));
	h_a16->GetXaxis()->SetTitle("Random Variable X");
	h_a16->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_a16->Draw();
	drawGauss(50,sum);
}



