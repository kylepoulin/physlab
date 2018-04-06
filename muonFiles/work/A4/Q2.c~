#include <stdio.h>
#include <math.h>

#define _USE_MATH_DEFINES


void Q2(){
	double norm = 1.0;
	gRandom->SetSeed(12345678);
	//5 steps
	TH1D *h_5sX = new TH1D("5 steps X","Five Steps - X", 50, -10,10);
	TH1D *h_5sR = new TH1D("5 steps R","Five Steps - Total Length", 50, -2,8);
	h_5sX->GetXaxis()->SetTitle("X - Position");
	h_5sX->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_5sR->GetXaxis()->SetTitle("Distance from Origin");
	h_5sR->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	

	//20 steps
	TH1D *h_20sX = new TH1D("20 steps X","Twenty Steps - X", 50, -15,15);
	TH1D *h_20sR = new TH1D("20 steps R","Twenty Steps - Total Length", 50, -5,15);
	h_20sX->GetXaxis()->SetTitle("X - Position");
	h_20sX->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_20sR->GetXaxis()->SetTitle("Distance from Origin");
	h_20sR->GetYaxis()->SetTitle("Normalized Number of Points in Bins");


	//100 steps
	TH1D *h_100sX = new TH1D("100 steps X","One Hundred Steps - X", 50, -25,25);
	TH1D *h_100sR = new TH1D("100 steps R","One Hundred Steps - Total Length", 50, -5,30);
	h_100sX->GetXaxis()->SetTitle("X - Position");
	h_100sX->GetYaxis()->SetTitle("Normalized Number of Points in Bins");
	h_100sR->GetXaxis()->SetTitle("Distance from Origin");
	h_100sR->GetYaxis()->SetTitle("Normalized Number of Points in Bins");

	//5 steps
	for(int k=0; k<10000; k++){

		double x=0, y=0;
		
		for(int i=0; i<5; i++){
			double theta = gRandom->Uniform(0,2)*M_PI;
			y += sin(theta);	
			x += cos(theta);
		}
		h_5sX->Fill(x);
		h_5sR->Fill(sqrt(pow(x,2)+pow(y,2)));
		//cout << y << " " << x << endl;
	
	}
	h_5sX->Scale(norm/h_5sX->Integral("width"));
	h_5sR->Scale(norm/h_5sR->Integral("width"));
	TCanvas *cans1 = new TCanvas();
	cans1->Divide(1,2);
	cans1->cd(1);
	h_5sX->Draw();
	cans1->cd(2);
	h_5sR->Draw();

	//20 steps
	for(int k=0; k<10000; k++){

		double x=0, y=0;
		
		for(int i=0; i<20; i++){
			double theta = gRandom->Uniform(0,2)*M_PI;
			y += sin(theta);	
			x += cos(theta);
		}
		h_20sX->Fill(x);
		h_20sR->Fill(sqrt(pow(x,2)+pow(y,2)));
		//cout << y << " " << x << endl;
	
	}
	h_20sX->Scale(norm/h_20sX->Integral("width"));
	h_20sR->Scale(norm/h_20sR->Integral("width"));
	TCanvas *cans2 = new TCanvas();
	cans2->Divide(1,2);
	cans2->cd(1);
	h_20sX->Draw();
	cans2->cd(2);
	h_20sR->Draw();

	//100 steps
	for(int k=0; k<10000; k++){

		double x=0, y=0;
		
		for(int i=0; i<100; i++){
			double theta = gRandom->Uniform(0,2)*M_PI;
			y += sin(theta);	
			x += cos(theta);
		}
		h_100sX->Fill(x);
		h_100sR->Fill(sqrt(pow(x,2)+pow(y,2)));
		//cout << y << " " << x << endl;
	
	}
	h_100sX->Scale(norm/h_100sX->Integral("width"));
	h_100sR->Scale(norm/h_100sR->Integral("width"));
	TCanvas *cans3 = new TCanvas();
	cans3->Divide(1,2);
	cans3->cd(1);
	h_100sX->Draw();
	cans3->cd(2);
	h_100sR->Draw();

	

}
