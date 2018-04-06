/*



*/



#include <math.h>

void A7(){
	int norm=1;
	int N = 10*40*40;
	int site[40][40];
	
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			site[i][j]=1;
		}
	}

	gRandom->SetSeed(123456);
	for(int trial=0;trial<N;trial++){
		int i = round(gRandom->Uniform(0,39));
		int j = round(gRandom->Uniform(0,39));
		
		if(site[i][j]>0){
			site[i][j]--;
			int i = round(gRandom->Uniform(0,39));
			int j = round(gRandom->Uniform(0,39));
			site[i][j]++;
		}
	}
	TH1D *NormHist1 = new TH1D("Energy Distribution","Normalized Energy Dist.", 11, -1, 10);
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			NormHist1->Fill(site[i][j]);
		}
	}
	//Reset site data for 2N
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			site[i][j]=2;
		}
	}
	for(int trial=0;trial<2*N;trial++){
		int i = round(gRandom->Uniform(0,39));
		int j = round(gRandom->Uniform(0,39));
		
		if(site[i][j]>0){
			site[i][j]--;
			int i = round(gRandom->Uniform(0,39));
			int j = round(gRandom->Uniform(0,39));
			site[i][j]++;
		}
	}
	TH1D *NormHist2 = new TH1D("Energy Distribution 2N","Normalized Energy Dist. for 2N", 21, -1, 20);
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			NormHist2->Fill(site[i][j]);
		}
	}

	//Reset Site Data for 5N
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			site[i][j]=5;
		}
	}
	for(int trial=0;trial<5*N;trial++){
		int i = round(gRandom->Uniform(0,39));
		int j = round(gRandom->Uniform(0,39));
		
		if(site[i][j]>0){
			site[i][j]--;
			int i = round(gRandom->Uniform(0,39));
			int j = round(gRandom->Uniform(0,39));
			site[i][j]++;
		}
	}
	TH1D *NormHist3 = new TH1D("Energy Distribution 5N","Normalized Energy Dist. for 5N", 26, -1, 25);
	for(int i=0;i<40;i++){
		for(int j=0;j<40;j++){
			NormHist3->Fill(site[i][j]);
		}
	}
	/*
	//Reset Site Data for N=400 400quanta in first spot
	int site2[20][20];
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			site2[i][j]=0;
		}
	}
	site2[0][0]=400;
	for(int trial=0;trial<5*N;trial++){
		int i = round(gRandom->Uniform(0,19));
		int j = round(gRandom->Uniform(0,19));
		
		if(site[i][j]>0){
			site[i][j]--;
			int i = round(gRandom->Uniform(0,19));
			int j = round(gRandom->Uniform(0,19));
			site[i][j]++;
		}
	}
	TH1D *NormHist4 = new TH1D("Energy Distribution","Normalized Energy Dist. for 400Quanta in first spot", 20, -1, 405);
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			NormHist4->Fill(site[i][j]);
		}
	}
	*/
	

	NormHist1->Scale(norm/NormHist1->Integral("width"));
	NormHist2->Scale(norm/NormHist2->Integral("width"));
	NormHist3->Scale(norm/NormHist3->Integral("width"));
	//NormHist4->Scale(norm/NormHist4->Integral("width"));

	TCanvas *can1 = new TCanvas("C1");
	NormHist1->Draw();

	TCanvas *can2 = new TCanvas("C2");
	NormHist2->Draw();

	TCanvas *can3 = new TCanvas("C3");
	NormHist3->Draw();
	
	//TCanvas *can4 = new TCanvas("C4");
	//NormHist4->Draw();
}
