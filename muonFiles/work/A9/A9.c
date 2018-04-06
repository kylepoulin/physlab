#include <stdio.h>

const double e = 2.71828182845904523536;

void A9(){
	gRandom->SetSeed(12345);
	TH1D *NormHist = new TH1D("","", 50, -2, 10);
	TH1D *NormHist2 = new TH1D("","", 50, -2, 10);
	for(int trials=0; trials < 10000; trials++){
		double xi=0;
		for(int i = 0; i<5; i++){
			xi+= gRandom->Exp(1.0);
		}
		double tau = xi/5;
		
		NormHist->Fill(tau);
		///////DO THE SAME FOR TAU2
		xi=0;
		for(int i = 0; i<5; i++){
			xi+= gRandom->Exp(2.5);
		}
		tau = xi/5;
		NormHist2->Fill(tau);
	}
	NormHist->Scale(1/NormHist->Integral("width"));
	NormHist2->Scale(1/NormHist2->Integral("width"));
	
	TCanvas *can1 = new TCanvas("C1");
	NormHist->GetYaxis()->SetTitle("Normalized number of bins filled");
	NormHist->GetXaxis()->SetTitle("Value of Tau-hat with Tau=1.0");
	NormHist->GetYaxis()->CenterTitle();
	NormHist->GetXaxis()->CenterTitle();
	NormHist->Draw();

	

	TF1 *fa1 = new TF1("fa1","pow([p0],[p0])/TMath::Factorial([p0]-1)*pow(x,[p0]-1)/pow([p1],[p0])*pow([p2],-[p0]*x/[p1])",0,20);
	fa1->SetParameters(5.0,1.0,e);
	fa1->Draw("same");
	
	TCanvas *can2 = new TCanvas("C2");
	NormHist2->GetYaxis()->SetTitle("Normalized number of bins filled");
	NormHist2->GetXaxis()->SetTitle("Value of Tau-hat with Tau=2.5");
	NormHist2->GetYaxis()->CenterTitle();
	NormHist2->GetXaxis()->CenterTitle();
	NormHist2->Draw();
	TF1 *fa2 = new TF1("fa2","pow([p0],[p0])/TMath::Factorial([p0]-1)*pow(x,[p0]-1)/pow([p1],[p0])*pow([p2],-[p0]*x/[p1])",0,20);
	fa2->SetParameters(5.0,2.5,e);
	fa2->Draw("same");



}
