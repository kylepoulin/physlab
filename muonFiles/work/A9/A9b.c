#include <stdio.h>

void A9b(){
	cout << TMath::ChisquareQuantile(1-0.1585,10);
	cout << TMath::ChisquareQuantile(0.1585,10);
	TF1 *fa1 = new TF1("fa1","x*TMath::ChisquareQuantile(1-[0],10)/10",0,5);
	TF1 *fa2 = new TF1("fa2","x*TMath::ChisquareQuantile([0],10)/10",0,5);
	//choosing a central confidence interval alpha = beta = gamma/2
	//where 1-gamma = 0.683 for this first case
	fa1->SetParameter(0,0.1585);
	fa2->SetParameter(0,0.1585);
	TCanvas *can1 = new TCanvas("C1");
	fa1->Draw();
	fa2->Draw("same");
	
	TF1 *fa3 = new TF1("fa3","x*TMath::ChisquareQuantile(1-[0],10)/10",0,5);
	TF1 *fa4 = new TF1("fa4","x*TMath::ChisquareQuantile([0],10)/10",0,5);
	fa3->SetParameter(0,0.025);
	fa4->SetParameter(0,0.025);

	fa3->Draw("same");
	fa4->Draw("same");

}
