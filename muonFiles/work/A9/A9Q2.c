#include <stdio.h>
TF1 *drawGauss(double mu, double sigma, int col=kRed) {
  TF1 *gaus = new TF1("","TMath::Gaus(x,[0],[1],1)",-100,1000);
  gaus->SetLineColor(col);
  gaus->SetParameters(mu,sigma);
  gaus->Draw("same");
  return gaus;
}

void A9Q2(){
	TH1D *NormHist1 = new TH1D("Sample Variance","Sample Variance", 50, 0, 100);
	TH1D *NormHist2 = new TH1D("Variance/Real-Variance*(n-1)","Variance/Real-Variance*(n-1)", 50, 0, 50);
	TH1D *NormHist3 = new TH1D("t-value","t-value", 50, -4, 4);
	

	double realMean=5;
	double realVariance=9;
	double sampMean;
	double sampVariance;
	int n = 5;
	for(int trials = 0; trials <10000; trials++){
		sampMean=0;
		sampVariance=0;
		double obs;
		for(int i=0; i<n; i++){
			obs = gRandom->Gaus(realMean,sqrt(realVariance));
			sampMean+= obs;
			sampVariance+= pow(obs-realMean,2);
		}
		sampMean=sampMean/n;
		sampVariance = sampVariance/(n-1);
		NormHist1->Fill(sampVariance);
		NormHist2->Fill(sampVariance/(realVariance)*(n-1));
		NormHist3->Fill((sampMean-realMean)/sqrt(sampVariance/n));
	}
	
	NormHist1->Scale(1/NormHist1->Integral("width"));
	NormHist2->Scale(1/NormHist2->Integral("width"));
	NormHist3->Scale(1/NormHist3->Integral("width"));

	TCanvas *can1 = new TCanvas("C1");
	NormHist1->GetYaxis()->SetTitle("Normalized number of bins filled");
	NormHist1->GetXaxis()->SetTitle("Sample Variance");
	NormHist1->GetYaxis()->CenterTitle();
	NormHist1->GetXaxis()->CenterTitle();
	NormHist1->Draw();
	//TF1 *chi2_1 = new TF1("chi2_1dofF","ROOT::Math::chisquared_pdf(x,10,0)",0,50);
	//chi2_1->SetLineColor(kRed);
	//chi2_1->Draw("same");
	TF1 *fits = drawGauss(9,3);

	TCanvas *can2 = new TCanvas("C2");
	NormHist2->GetYaxis()->SetTitle("Normalized number of bins filled");
	NormHist2->GetXaxis()->SetTitle("z value");
	NormHist2->GetYaxis()->CenterTitle();
	NormHist2->GetXaxis()->CenterTitle();
	NormHist2->Draw();
	TF1 *chi2_2 = new TF1("chi2_4dofF","ROOT::Math::chisquared_pdf(x,4,0)",0,25);
	chi2_2->SetLineColor(kRed);
	chi2_2->Draw("same");

	TCanvas *can3 = new TCanvas("C3");
	NormHist3->GetYaxis()->SetTitle("Normalized number of bins filled");
	NormHist3->GetXaxis()->SetTitle("t-value");
	NormHist3->GetYaxis()->CenterTitle();
	NormHist3->GetXaxis()->CenterTitle();
	NormHist3->Draw();
	TF1 *fit1 = drawGauss(0,1);
	cout << "When changing the mean, we don't see any difference in any of the calculated parameters. When we change the real variance, we see a difference only in the sample variance." << endl;
}
