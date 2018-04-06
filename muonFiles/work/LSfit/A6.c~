// make the data globally available
TGraphErrors *data;

// calculate the chiSquare quantity assuming a linear function
// with the parameters alpha and beta
double chi2_linear(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha - beta*data->GetX()[i],2)/pow(data->GetEY()[i],2);
  return chi2;
}

double chi2_zeroOrder(double alpha){
	double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow(data->GetY()[i] - alpha,2)/pow(data->GetEY()[i],2);
  return chi2;
}

double chi2_secondOrder(double alpha, double beta, double gamma) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha - beta*data->GetX()[i] - gamma*pow(data->GetX()[i],2),2)/pow(data->GetEY()[i],2);
  return chi2;
}

// Minuit function used by fitter
// Minuit minimize, so this function should return either the chi^2 or -2 ln L
// Error will be given by +/- 1 of this quantity
void minuitChi20(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0];
  result = chi2_zeroOrder(alpha);
}

void minuitChi21(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_linear(alpha,beta);
}

void minuitChi22(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1], gamma = par[2];
  result = chi2_secondOrder(alpha,beta,gamma);
}

void A6() {
  // true parameters
  double a = 1.0, b=0.5, stdDev = 1.0;
  
  // 1. Generate the data
  gRandom->SetSeed(12345);
  data = new TGraphErrors();
  for (int i=1;i<=10;++i) {
    data->SetPoint(i-1,i,gRandom->Gaus(a+b*i,stdDev));
    data->SetPointError(i-1,0,stdDev);
  }
  data->SetMarkerStyle(20);
  data->Draw("AP");
  gPad->Print("toy_data.pdf");
  
  TFitter *fitter0 = new TFitter(10);
  fitter0->SetFCN(minuitChi20);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter0->SetParameter(0,"alpha0",1,err,low_val,high_val);
 
  
  double arg = 0;
  fitter0->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat0 = fitter0->GetParameter(0), err_alphaHat0 = fitter0->GetParError(0);
  

 //fitter for the pol1
  TFitter *fitter1 = new TFitter(10);
  fitter1->SetFCN(minuitChi21);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  //double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter1->SetParameter(0,"alpha1",1,err,low_val,high_val);
  fitter1->SetParameter(1,"beta1",1,err,low_val,high_val);
  
  
  //double arg = 0;
  fitter1->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat1 = fitter1->GetParameter(0), err_alphaHat1 = fitter1->GetParError(0);
  double betaHat1  = fitter1->GetParameter(1), err_betaHat1  = fitter1->GetParError(1);


//fitter for the pol2
  TFitter *fitter2 = new TFitter(10);
  fitter2->SetFCN(minuitChi22);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  //double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter2->SetParameter(0,"alpha2",1,err,low_val,high_val);
  fitter2->SetParameter(1,"beta2",1,err,low_val,high_val);
  fitter2->SetParameter(2, "gamma2",1,err,low_val,high_val);
  
  //double arg = 0;
  fitter2->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat2 = fitter2->GetParameter(0), err_alphaHat2 = fitter2->GetParError(0);
  double betaHat2  = fitter2->GetParameter(1), err_betaHat2  = fitter2->GetParError(1);
  double gammaHat2  = fitter2->GetParameter(2), err_gammaHat2  = fitter2->GetParError(2);
  
 
  printf("\nDegrees of freedom = %.3d\n",data->GetN()-1);
  printf("Zero-th order Polynomial parameters:\n");
  printf("alpha0: %.3f +/- %.3f\n\n",alphaHat0,err_alphaHat0);

  printf("First order Polynomial parameters:\n");
  printf("alpha1: %.3f +/- %.3f\n",alphaHat1,err_alphaHat1);
  printf("beta1:  %.3f +/- %.3f\n\n",betaHat1,err_betaHat1);

  printf("Second order Polynomial parameters:\n");
  printf("alpha2: %.3f +/- %.3f\n",alphaHat2,err_alphaHat2);
  printf("beta2:  %.3f +/- %.3f\n",betaHat2,err_betaHat2);
  printf("gamma2: %.3f +/- %.3f\n\n",gammaHat2,err_gammaHat2);
  
  printf("Zero order chi2: %.3f\n",chi2_zeroOrder(alphaHat0));
  printf("First order chi2: %.3f\n",chi2_linear(alphaHat1,betaHat1));
  printf("Second order chi2: %.3f\n\n",chi2_secondOrder(alphaHat2,betaHat2,gammaHat2));

  printf("Pval of Zero Order: %.3f\n", TMath::Prob(chi2_zeroOrder(alphaHat0),data->GetN()-1));
  printf("Pval of First Order: %.3f\n", TMath::Prob(chi2_linear(alphaHat1,betaHat1),data->GetN()-1));
  printf("Pval of Second Order: %.3f\n\n", TMath::Prob(chi2_secondOrder(alphaHat2,betaHat2,gammaHat2),data->GetN()-1));

  printf("First order Chi2 with true parameters: %.3f\n",chi2_linear(1,0.5));
  printf("This value is smaller than the minimized chi2 we calculated, since the line we fitted is based on the actual data points compared to the 'true' line which the points are based on.\n");
  printf("It doesn't make sense to do this with the second order, since the data is linear (we don't have a true gamma value.\n");

  printf("First order chi2 modified: %.3f\n", chi2_linear(alphaHat1+err_alphaHat1,betaHat1));

	TF1 *fit0 = new TF1("","pol0",0,11); // pol0, pol1 or pol2 ...
	fit0->SetParameters(alphaHat0,0); // plug in the fitted values
	fit0->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

	TF1 *fit1 = new TF1("","pol1",0,11); // pol0, pol1 or pol2 ...
	fit1->SetParameters(alphaHat1,betaHat1); // plug in the fitted values
	fit1->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

	TF1 *fit2 = new TF1("","pol2",0,11); // pol0, pol1 or pol2 ...
	fit2->SetParameters(alphaHat2,betaHat2,gammaHat2); // plug in the fitted values
	fit2->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

}
