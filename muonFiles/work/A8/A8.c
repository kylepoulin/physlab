// make the data globally available
TGraphErrors *data;

double chi2_secondOrder(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha*data->GetX()[i] - beta*pow(data->GetX()[i],2),2)/pow(data->GetEY()[i],2);
  return chi2;
}

double chi2_linear(double alpha) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha*data->GetX()[i],2)/pow(data->GetEY()[i],2);
  return chi2;
}

void minuitChi21(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0];
  result = chi2_linear(alpha);
}

void minuitChi22(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_secondOrder(alpha,beta);
}

double chi2_exponential(double alpha, double beta){
	double chi2=0;

  for (int i=0; i<data->GetN();i++)
	chi2 += pow( data->GetY()[i] - alpha*pow(data->GetX()[i],beta),2)/pow(data->GetEY()[i],2);
  return chi2;
}

void minuitChi2exp(int &nDim, double* gout, double& result, double *par, int flg){
	double alpha = par[0], beta = par[1];
	result = chi2_exponential(alpha,beta);
}

double chi2_sqrt(double alpha) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha*sqrt(data->GetX()[i]),2)/pow(data->GetEY()[i],2);
  return chi2;
}
void minuitChi2sqrt(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0];
  result = chi2_sqrt(alpha);
}



double h[5] = {1000,828, 800, 600, 300};
double d[5] = {1500, 1340, 1328, 1172, 800};

void A8(){

  double arg=0;

  data = new TGraphErrors();
  for (int i=0;i<5;++i) {
    data->SetPoint(i,h[i],d[i]);
    data->SetPointError(i,0,15);
  }
  data->SetMarkerStyle(20);
  data->Draw("AP");
  gPad->Print("toy_data.pdf");

  //fitter for the pol1
  TFitter *fitter1 = new TFitter(5);
  fitter1->SetFCN(minuitChi21);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter1->SetParameter(0,"alpha1",1,err,low_val,high_val);
  //fitter1->SetParameter(1,"beta1",1,err,low_val,high_val);
  
  
  //double arg = 0;
  fitter1->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat1 = fitter1->GetParameter(0), err_alphaHat1 = fitter1->GetParError(0);
  //double betaHat1  = fitter1->GetParameter(1), err_betaHat1  = fitter1->GetParError(1);


//fitter for the pol2
  TFitter *fitter2 = new TFitter(5);
  fitter2->SetFCN(minuitChi22);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  //double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter2->SetParameter(0,"alpha2",1,err,low_val,high_val);
  fitter2->SetParameter(1,"beta2",1,err,low_val,high_val);

  
  //double arg = 0;
  fitter2->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat2 = fitter2->GetParameter(0), err_alphaHat2 = fitter2->GetParError(0);
  double betaHat2  = fitter2->GetParameter(1), err_betaHat2  = fitter2->GetParError(1);


  //Fitter for the exponential function
  TFitter *fitter3 = new TFitter(5);
  fitter3->SetFCN(minuitChi2exp);

  fitter3->SetParameter(0,"alpha2",1,err,low_val,high_val);
  fitter3->SetParameter(1,"beta2",1,err,low_val,high_val);

  fitter3->ExecuteCommand("MIGRAD", &arg, 0);

  double alphaHat3 = fitter3->GetParameter(0), err_alphaHat3 = fitter3->GetParError(0);
  double betaHat3  = fitter3->GetParameter(1), err_betaHat3  = fitter3->GetParError(1);

  //fitter for the sqrt thing
  TFitter *fitter4 = new TFitter(5);
  fitter4->SetFCN(minuitChi2sqrt);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  //double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter4->SetParameter(0,"alpha1",1,err,low_val,high_val);
  //fitter1->SetParameter(1,"beta1",1,err,low_val,high_val);
  
  
  //double arg = 0;
  fitter4->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat4 = fitter4->GetParameter(0), err_alphaHat4 = fitter4->GetParError(0);
  //double betaHat1  = fitter1->GetParameter(1), err_betaHat1  = fitter1->GetParError(1);

  printf("First order Polynomial parameters:\n");
  printf("alpha1: %.3f +/- %.3f\n",alphaHat1,err_alphaHat1);
  //printf("beta1:  %.3f +/- %.3f\n\n",betaHat1,err_betaHat1);

  printf("Second order Polynomial parameters:\n");
  printf("alpha2: %.3f +/- %.3f\n",alphaHat2,err_alphaHat2);
  printf("beta2:  %.3f +/- %.3f\n",betaHat2,err_betaHat2);
  //printf("gamma2: %.3f +/- %.3f\n\n",gammaHat2,err_gammaHat2);

  printf("Exponential parameters:\n");
  printf("alpha2: %.3f +/- %.3f\n",alphaHat3,err_alphaHat3);
  printf("beta2:  %.3f +/- %.3f\n",betaHat3,err_betaHat3);

  printf("Sqrt parameters:\n");
  printf("alpha2: %.3f +/- %.3f\n",alphaHat4,err_alphaHat4);
  

  printf("First order chi2: %.3f\n",chi2_linear(alphaHat1));
  printf("Second order chi2: %.3f\n\n",chi2_secondOrder(alphaHat2,betaHat2));
  printf("Exponential chi2: %.3f\n\n",chi2_exponential(alphaHat3,betaHat3));
  printf("Sqrt chi2: %.3f\n\n",chi2_sqrt(alphaHat4));

  printf("Pval of First Order: %.3f\n", TMath::Prob(chi2_linear(alphaHat1),data->GetN()-1));
  printf("Pval of Second Order: %.3f\n\n", TMath::Prob(chi2_secondOrder(alphaHat2,betaHat2),data->GetN()-2));
  printf("Pval of Exponential: %.3f\n\n", TMath::Prob(chi2_exponential(alphaHat3,betaHat3),data->GetN()-2));
  printf("Pval of Sqrt: %.3f\n\n", TMath::Prob(chi2_sqrt(alphaHat4),data->GetN()-1));


  TF1 *fit1 = new TF1("","pol1",0,2000); 
  fit1->SetParameters(0,alphaHat1); // plug in the fitted values
  fit1->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

  TF1 *fit2 = new TF1("","pol2",0,2000); // pol0, pol1 or pol2 ...
  fit2->SetParameters(0,alphaHat2,betaHat2); // plug in the fitted values
  fit2->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

  TF1 *fit3 = new TF1("","[0]*pow(x,[1])",0,2000); // pol0, pol1 or pol2 ...
  fit3->SetParameter(0,alphaHat3); // plug in the fitted values
  fit3->SetParameter(1,betaHat3); // plug in the fitted values
  fit3->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)

  TF1 *fit4 = new TF1("","[0]*sqrt(x)",0,2000); // pol0, pol1 or pol2 ...
  fit4->SetParameter(0,alphaHat4); // plug in the fitted values
  fit4->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)


}
