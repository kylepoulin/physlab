#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>


vector<double> time1;
vector<double> counts;

double tau0 = 0.0000002;

double chi2_MUON(double N1, double N2, double tau1, double tau2, double b){
  double chi2 = 0;
  for(int i=0; i<time1.size(); i++){
    chi2 += pow(counts[i] - b - N1*exp(-(time1[i])/tau1) - N2*exp(-(time1[i])/tau2),2)/(b + N1*exp(-(time1[i])/tau1) + N2*exp(-(time1[i])/tau2));
  }
  return chi2;
}

void minuitChi2(int &nDim, double* gout, double& result, double *par, int flg) {
  double n1 = par[0], n2 = par[1], t1 = par[2], t2 = par[3], back = par[4];
  result = chi2_MUON(n1,n2,t1,t2,back);
}



void Fitter()
{
  //Open the files and fill the data vectors
  ifstream infile;
  infile.open("28time2.txt");
  
  double value;
  if(infile.is_open()){
    while(infile>>value){
      time1.push_back(value);
    }
  }
  infile.close();
  
  infile.open("28counts2.txt");
  
  if(infile.is_open()){
    while(infile>>value){
      counts.push_back(value);
    }
  }
  infile.close();

  double arg = 0;
  
  TFitter *fitter = new TFitter(10);
  fitter->SetFCN(minuitChi2);
  //Set the fit parameters
  fitter->SetParameter(0, "N1", 75000, 0.001, 1000, 1000000);
  fitter->SetParameter(1, "N2", 1000, 0.001, 1000, 1000000);
  fitter->SetParameter(2, "tau1", 0.0000001, 0.0000000001, 0, 0.000001);
  fitter->SetParameter(3, "tau2", 0.000001, 0.0000000001, 0, 0.000003);
  fitter->SetParameter(4, "b", 10, 0.1, 0, 1000);
  
  fitter->ExecuteCommand("MIGRAD",&arg,0);
  

  double Bhat = fitter->GetParameter(4);
  double N1hat = fitter->GetParameter(0);
  double N2hat = fitter->GetParameter(1);
  double tau1hat = fitter->GetParameter(2);
  double tau2hat = fitter->GetParameter(3);
  

  printf("Background: %.3f\n", fitter->GetParameter(4));
  printf("N1: %f\n", fitter->GetParameter(0));
  printf("N2: %f\n", fitter->GetParameter(1));
  printf("Tau1: %e\n", fitter->GetParameter(2));
  printf("Tau2: %e\n", fitter->GetParameter(3));
  printf("Chi2: %f\n", chi2_MUON(fitter->GetParameter(0),fitter->GetParameter(1),fitter->GetParameter(2),fitter->GetParameter(3),fitter->GetParameter(4)));
  
  TCanvas *can1 = new TCanvas("C1");
  
  TH1D *hist = new TH1D("Muon Times", "Muon Times", time1.size(), 0, 0.00001);
  for(int i=0; i<time1.size(); i++){
    for(int y=0; y<counts[i]; y++){
      hist->Fill(time1[i]);
    }
  }
  hist->Draw();
  TF1 *theFit = new TF1("Exponential Fit", "[4] + [0]*exp(-(x)/[2]) + [1]*exp(-(x)/[3])",0,0.00001);

  theFit->SetParameters(N1hat,N2hat,tau1hat,tau2hat,Bhat);
  theFit->Draw("SAME");

  //Calculate the std of tau1 and tau2
  TGraph *gr1 = new TGraph();
  TCanvas *can2 = new TCanvas("C2");
  TF1 *chi2Tau1= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],x,[2],[3])",0.1*pow(10,-9),2*tau1hat);
  
  chi2Tau1->SetParameters(N1hat,N2hat,tau2hat,Bhat);
  chi2Tau1->Draw();
  double stdTau11 = tau1hat-chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,0,tau1hat);
  double stdTau12 = chi2Tau1->GetX(chi2Tau1->Eval(tau1hat)+1,tau1hat,3*tau1hat)-tau1hat;

  gr1->SetPoint(0, tau1hat-stdTau11, chi2Tau1->Eval(tau1hat)+1);
  gr1->SetPoint(1, tau1hat, chi2Tau1->Eval(tau1hat));
  gr1->SetPoint(2, tau1hat+stdTau12, chi2Tau1->Eval(tau1hat)+1);
  gr1->SetMarkerStyle(20);
  gr1->Draw("C* SAME");

  printf("Tau1 chi minimum: %f\n", chi2Tau1->GetMinimum(0.1*pow(10,-9),2*tau1hat));
  printf("Tau1 std1: %e\n", stdTau11);
  printf("Tau1 std2: %e\n", stdTau12);
 
  TGraph *gr2 = new TGraph();
  TCanvas *can3 = new TCanvas("C3");
  TF1 *chi2Tau2= new TF1("chi2 dist of tau1", "chi2_MUON([0],[1],[2],x,[3])",0.1*pow(10,-6),2*tau2hat);
  chi2Tau2->SetParameters(N1hat,N2hat,tau1hat,Bhat);
  printf("Tau2 chi minimum: %f\n", chi2Tau1->GetMinimum(0.1*pow(10,-6),2*tau2hat));
  chi2Tau2->Draw();
  double stdTau21 = tau2hat-chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,0,tau2hat);
  double stdTau22 = chi2Tau2->GetX(chi2Tau2->Eval(tau2hat)+1,tau2hat,3*tau2hat)-tau2hat;

  gr2->SetPoint(0, tau2hat-stdTau21, chi2Tau2->Eval(tau2hat)+1);
  gr2->SetPoint(1, tau2hat, chi2Tau2->Eval(tau2hat));
  gr2->SetPoint(2, tau2hat+stdTau22, chi2Tau2->Eval(tau2hat)+1);
  gr2->SetMarkerStyle(20);
  gr2->Draw("C* SAME");

  printf("Tau2 std1: %e\n", stdTau21);
  printf("Tau2 std2: %e\n", stdTau22);

}
