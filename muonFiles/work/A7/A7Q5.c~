/*

KE = 3/2*(1.38066E-23)*T

KE = 1/2*mv^2

v= sqrt((3/m)*(1.38066E-23)*T)

where T is poisson distributed - do this for Vx, Vy, Vz

*/
const double m = 2.32586705*pow(10,-26);
const double kb = 1.38*pow(10,-23);
const double T = 300.0;

double VfromT(){
	return pow(-1, round(gRandom->Uniform(1,2)))*sqrt((3/m)*1.38066*pow(10,-23)*gRandom->Poisson(300));
}

TF1 *drawGauss(double mu, double sigma, int col=kRed) {
  TF1 *gaus = new TF1("","TMath::Gaus(x,[0],[1],1)",-2000,2000);
  gaus->SetLineColor(col);
  gaus->SetParameters(mu,sigma);
  gaus->Draw("same");
  return gaus;
}



void A7Q5(){
	double norm = 1.0;
	gRandom->SetSeed(123456);
	double v, vx, vy, vz;
	TH1D *histVx = new TH1D("Vx","Normalized Velocity in X", 100, -850, 850);
	TH1D *histV = new TH1D("V","Normalized Speed", 100, 0, 1400);
	double var=0;

	for(int trial=0;trial<1000;trial++){
		vx = gRandom->Gaus(0,sqrt(kb*T/(2*m)));
		cout << vx << endl;
		vy = gRandom->Gaus(0,sqrt(kb*T/(2*m)));
		vz = gRandom->Gaus(0,sqrt(kb*T/(2*m)));
		v = sqrt(pow(vx,2) + pow(vy,2) + pow(vz,2));
		histVx->Fill(vx);
		histV->Fill(v);
		var += pow(v-476,2);
	}
	var = var/999.0;
	histVx->Scale(norm/histVx->Integral("width"));
	histV->Scale(norm/histV->Integral("width"));

	histVx->GetXaxis()->SetTitle("Vx");
	histVx->GetXaxis()->CenterTitle();
	histVx->GetYaxis()->SetTitle("Normalized number of data in Bins");
	histVx->GetYaxis()->CenterTitle();

	histV->GetXaxis()->SetTitle("V");
	histV->GetXaxis()->CenterTitle();
	histV->GetYaxis()->SetTitle("Normalized number of data in Bins");
	histV->GetYaxis()->CenterTitle();



	TCanvas *can1 = new TCanvas("C1");
	histVx->Draw();
	TF1 *fitvx = drawGauss(0, sqrt(kb*T/(2*m)));

	TCanvas *can2 = new TCanvas("C2");
	histV->Draw();
	
	//TF1 *fit1 = drawGauss(476,sqrt(var));
	TF1 *fitv = new TF1("FitV","pow(2.32586705*pow(10,-26)/(2*1.38*pow(10,-23)*300.0*TMath::Pi()),3/2)*pow(x,3)*4*TMath::Pi()*exp(-2.32586705*pow(10,-26)*pow(x,2)/(2*1.38*pow(10,-23)*300.0))",0,2000);
	fitv->Draw("same");

}

