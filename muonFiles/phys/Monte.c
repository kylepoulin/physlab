#include <stdio.h>
#include <math.h>

//Define a global for the histogram
TH1D *hist;

double chi2_MUON(double N1, double N2, double tau1, double tau2, double b){
  double chi2 = 0;
  for(int i=0; i<hist->GetNbinsX(); i++){
    chi2 += pow(hist->GetBinContent(i) - b - N1*exp(-(hist->GetBinCenter(i))/tau1) - N2*exp(-(hist->GetBinCenter(i))/tau2),2)/(b + N1*exp(-(hist->GetBinCenter(i))/tau1) + N2*exp(-(hist->GetBinCenter(i))/tau2));
  }
  return chi2;
}

void minuitChi2(int &nDim, double* gout, double& result, double *par, int flg) {
  double n1 = par[0], n2 = par[1], t1 = par[2], t2 = par[3], back = par[4];
  result = chi2_MUON(n1,n2,t1,t2,back);
}


int count=0;
int point=0;

void Monte(){
	//Setup
	double TMIN=0.0;
	double TMAX=0.00002;
	double eventProb;
	double tau = 0.000002197;
	double tau2 = 0.0000002;
	double DELT = 0.000003;
	int NBINS = round(512);
	gRandom->SetSeed(12345678);

	double numEvents = 1000000;


	TFitter *fitter = new TFitter(10);
  	fitter->SetFCN(minuitChi2);
  	
	TGraphErrors *gr1 = new TGraphErrors();
	TGraphErrors *gr2 = new TGraphErrors();
	
	
	int count=0;
	int point=0;
	double arg=0;
	int countB=0;
	int countSpont=0;
	int countCap=0;

	hist = new TH1D("Time hist", "Histogram of observed decay time",NBINS,0,TMAX-TMIN);
	cout << "NBINS: " << NBINS << endl;
	
	for(int i=0; i<numEvents; i++){
		count++;
		eventProb = gRandom->Uniform(0,1);
		if(eventProb>=0.1){
			//simulate a spontaneous u decay
			double SMIN = tau*(1-exp(-TMIN/tau));
			double SMAX = tau*(1-exp(-TMAX/tau));
			double s = SMIN+(gRandom->Uniform(0,1)*(SMAX-SMIN));
			double t = -tau*log(1-s/tau);
			//cout << "t: " << t-TMIN << endl;
			hist->Fill((t-TMIN));
			countSpont++;
		} else {
			//choose between a background event or captured u event
			eventProb = gRandom->Uniform(0,1);
			if(eventProb<0.1){
				//simulate captured u		     
				double SMIN = tau2*(1-exp(-TMIN/tau2));
				double SMAX = tau2*(1-exp(-TMAX/tau2));
				double s = SMIN+(gRandom->Uniform(0,1)*(SMAX-SMIN));
				double t = -tau2*log(1-s/tau2);
				hist->Fill((t-TMIN));
				countCap++;
			} else {
				//simulate background event
				countB++;
				double t =gRandom->Uniform(0,1)*(TMAX-TMIN);				
				hist->Fill(t);
			}
		}

		if(count==50000){

			fitter->SetParameter(0, "N1", 50000, 0.1, 1, 100000);
  			fitter->SetParameter(1, "N2", 1500, 0.1, 1, 100000);
			fitter->SetParameter(2, "tau1", 0.0000002, 0.0000000001, 0.00000005, 0.000001);
			fitter->SetParameter(3, "tau2", 0.000001, 0.0000000001, 0.0000018, 0.000003);
  			fitter->SetParameter(4, "b", 10, 0.01, 0, 1000);
			fitter->ExecuteCommand("MIGRAD",&arg,0);
			double Bhat = fitter->GetParameter(4);
  			double N1hat = fitter->GetParameter(0);
  			double N2hat = fitter->GetParameter(1);
  			double tau1hat = fitter->GetParameter(2);
  			double tau2hat = fitter->GetParameter(3);

			//Calc Tau1 error and store point
			TF1 *chi2Tau1= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],x,[2],[3])",0.1*pow(10,-9),10*tau1hat);
  			chi2Tau1->SetParameters(N1hat,N2hat,tau2hat,Bhat);
  			double stdTau11 = tau1hat-chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,0,tau1hat);
  			double stdTau12 = chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,tau1hat,10*tau1hat)-tau1hat;
			double stdTau1 = (stdTau11+stdTau12)/2;
			gr1->SetPoint(point,i,tau1hat);
			gr1->SetPointError(point,0,stdTau1);

			//Calc Tau2 error and store point
			TF1 *chi2Tau2= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],[2],x,[3])",0.1*pow(10,-6),10*tau2hat);
  			chi2Tau2->SetParameters(N1hat,N2hat,tau1hat,Bhat);
  			double stdTau21 = tau2hat-chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,0,tau2hat);
  			double stdTau22 = chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,tau2hat,10*tau2hat)-tau2hat;
			double stdTau2 = (stdTau21+stdTau22)/2;
			gr2->SetPoint(point,i,tau2hat);
			gr2->SetPointError(point,0, stdTau2);
			printf("\033[22;32mTau1: %e\033[22;30m\n", tau1hat);
			printf("\033[22;32mstdTau11: %e , stdTau12: %e\033[22;30m\n", stdTau11, stdTau12);
			printf("\033[22;32mTau2: %e\033[22;30m\n", tau2hat);
			printf("\033[22;32mstdTau21: %e , stdTau22: %e\033[22;30m", stdTau21, stdTau22);
			point++;
			count=0;
		}
	}
	//hist->Scale(1/hist->Integral("width"));
	TCanvas *can1 = new TCanvas("C1");
	fitter->SetParameter(0, "N1", 50000, 0.1, 11000, 1000000);
  	fitter->SetParameter(1, "N2", 1500, 0.1, 1, 1000000);
	fitter->SetParameter(2, "tau1", 0.0000002, 0.0000000001, 0.00000005, 0.000001);
	fitter->SetParameter(3, "tau2", 0.000002, 0.0000000001, 0.0000018, 0.000003);
  	fitter->SetParameter(4, "b", 10, 0.01, 0, 1000);
	fitter->ExecuteCommand("MIGRAD",&arg,0);
	double Bhat = fitter->GetParameter(4);
  	double N1hat = fitter->GetParameter(0);
  	double N2hat = fitter->GetParameter(1);
  	double tau1hat = fitter->GetParameter(2);
  	double tau2hat = fitter->GetParameter(3);

	TF1 *theFit = new TF1("Exponential Fit", "[4] + [0]*exp(-(x)/[2]) + [1]*exp(-(x)/[3])",0,0.00006);
  	theFit->SetParameters(N1hat,N2hat,tau1hat,tau2hat,Bhat);
  	
	hist->Draw();
	theFit->Draw("SAME");
	printf("Chi2: %f\n", chi2_MUON(N1hat,N2hat,tau1hat,tau2hat,Bhat));


	//Calculate the std of tau1 and tau2
  TCanvas *can2 = new TCanvas("C2");
  TF1 *chi2Tau1= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],x,[2],[3])",0.1*pow(10,-9),10*tau1hat);
  
  chi2Tau1->SetParameters(N1hat,N2hat,tau2hat,Bhat);
  chi2Tau1->Draw();
  double stdTau11 = tau1hat-chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,0,tau1hat);
  double stdTau12 = chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,tau1hat,10*tau1hat)-tau1hat;

  printf("Tau1: %e\n", tau1hat);
  printf("Tau1 std1: %e\n", stdTau11);
  printf("Tau1 std2: %e\n", stdTau12);
  printf("Tau1 std: %e\n", (stdTau12+stdTau11)/2);
  printf("Tau1 std from MIGRAD: %e\n", fitter->GetParError(2));
 
  TCanvas *can3 = new TCanvas("C3");
  TF1 *chi2Tau2= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],[2],x,[3])",0.1*pow(10,-9),5*tau2hat);
  chi2Tau2->SetParameters(N1hat,N2hat,tau1hat,Bhat);
  printf("Tau2: %e\n", tau2hat);
  chi2Tau2->Draw();
  double stdTau21 = tau2hat-chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,0,tau2hat);
  double stdTau22 = chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,tau2hat,5*tau2hat)-tau2hat;

  printf("Tau2 std1: %e\n", stdTau21);
  printf("Tau2 std2: %e\n", stdTau22);
  printf("Tau2 std: %e\n", (stdTau22+stdTau21)/2);
  printf("Tau2 std from MIGRAD: %e\n", fitter->GetParError(3));

	TCanvas *can4 = new TCanvas("C4");
	gr1->SetMarkerStyle(20);
  	gr1->Draw("AP");

	TCanvas *can5 = new TCanvas("C5");
	gr2->SetMarkerStyle(20);
	gr2->Draw("AP");

  printf("Spontaneous: %d\n, Captured: %d\n, Background: %d\n", countSpont, countCap, countB);
}
