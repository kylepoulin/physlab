//#include "Utils.h"
#include "stdio.h"

void volumeCalc(){

	int nSpace = 5;
	
	double Ndarts = 1000;
	int Trials = 100;
	TH1D *h_volume = new TH1D("",";Volume;Number of trials",50,0,10);
	TH1D *h_efficiency = new TH1D("",";Efficiency;Number of trials",50,0,1);
	
	//double space[4];
	
	
	for(int i=0; i < Trials; i++){
		double Nhits = 0;
		for (int g=0; g<Ndarts; g++){
			double sum=0;
			for(int j=0; j<nSpace; j++){
				double dart;
				//space[j] = gRandom->Uniform(-1,1);
				dart = gRandom->Uniform(-1,1);
				sum = sum + pow(dart,2);
			}
			bool hit = sum < 1;
			if(hit) Nhits++;
		}
		double V = pow(2,nSpace)*Nhits/Ndarts;
		cout << Nhits << "   " << Ndarts << endl;
		//cout << V << endl;
		h_volume->Fill(V);
		
		h_efficiency->Fill(Nhits/Ndarts);
	}

	h_volume->Draw();
	TCanvas *c2 = new TCanvas();
	h_efficiency->Draw();
}
