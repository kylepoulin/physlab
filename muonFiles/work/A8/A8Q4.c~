// make the data globally available
TGraphErrors *data;

double t[6] = {0.5, 1.5, 2.5, 3.0, 3.6, 4.2};
double vel[6] = {8.1, 17.6, 27.5, 32.8, 38.8, 44.8};


double chi2_linear(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha - beta*data->GetX()[i],2)/pow(data->GetEY()[i],2);
  return chi2;
}

void minuitChi21(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_linear(alpha,beta);
}


void A8Q4(){

  double arg=0;


  data = new TGraphErrors();
  for (int i=0;i<6;++i) {
    data->SetPoint(i,t[i],vel[i]);
    data->SetPointError(i,0,0.2);
  }
  data->SetMarkerStyle(20);
  data->Draw("AP");
  gPad->Print("toy_data.pdf");


  //fitter for the pol2
  TFitter *fitter2 = new TFitter(10);
  fitter2->SetFCN(minuitChi21);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter2->SetParameter(0,"alpha2",1,err,low_val,high_val);
  fitter2->SetParameter(1,"beta2",1,err,low_val,high_val);
  
  
  //double arg = 0;
  fitter2->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat2 = fitter2->GetParameter(0), err_alphaHat2 = fitter2->GetParError(0);
  double betaHat2  = fitter2->GetParameter(1), err_betaHat2  = fitter2->GetParError(1);
  

  printf("Second order Polynomial parameters:\n");
  printf("alpha2: %.3f +/- %.3f\n",alphaHat2,err_alphaHat2);
  printf("beta2:  %.3f +/- %.3f\n",betaHat2,err_betaHat2);
  


  printf("Second order chi2: %.3f\n\n",chi2_linear(alphaHat2,betaHat2));

 double chi2value = chi2_linear(alphaHat2,betaHat2);

  printf("Pval of Second Order: %.3f\n\n", TMath::Prob(chi2_linear(alphaHat2,betaHat2),data->GetN()-2));

  TF1 *fit2 = new TF1("","pol2",0,11); // pol0, pol1 or pol2 ...
  fit2->SetParameters(alphaHat2,betaHat2); // plug in the fitted values
  fit2->Draw("same"); // draw it, "same" means on same canvas (i.e. on top)


  TCanvas *can1 = new TCanvas("C1");

  fitter2->GetMinuit()->SetErrorDef(chi2value+1);
  TGraph *conf_1 = (TGraph*)fitter2->GetMinuit()->Contour(120);
  conf_1->SetFillColor(kGray);


  fitter2->GetMinuit()->SetErrorDef(chi2value+1.645);
  TGraph *conf_2 = (TGraph*)fitter2->GetMinuit()->Contour(80);
  conf_2->SetFillColor(kAzure);


  conf_2->Draw("af"); // only fill "f" (draw on top, i.e. no axis)
  conf_1->Draw("f"); // A: draw with axis and F: fill
  


}
