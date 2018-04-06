#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

//define variables

double sampleMean(vector<double> data);
double sampleVariance(vector<double> data);
double sampleStandardDeviation(vector<double> data);
double median(vector<double> data);
double covariance(vector<double> data1, vector<double> data2);
double max(vector<double> data);
double min(vector<double> data);

//main function
void exercise1(){
	vector<double> data1;
	vector<double> data2;
	
//open H_pT and place values in the vector data1
	ifstream infile;
	infile.open("H_pT.txt");

	double value;
	if (infile.is_open()){
		while (infile>>value){
		data1.push_back(value);
		}
	}
	infile.close();

//open jet_pT and place values in the vector data2
	ifstream infile2;
	infile2.open("jet1_pT.txt");

	double value2;
	if (infile2.is_open()){
		while (infile2>>value2){
		data2.push_back(value2);
		}
	}
	infile2.close();
	
	
//set variables for the vector sizes
	int Ndata1 = data1.size();
	int Ndata2 = data2.size();	
	cout.precision(3);

//creating the canvas and setting 
	TCanvas *c3=new TCanvas();
//	double maxdata1=max(data1);
//	double maxdata2=max(data2);
//	double mindata1=min(data1);
//	double mindata2=min(data2);


//draw the 2D histogram, since we want the vectors to be unsorted	
	TH2D *h3 = new TH2D("2D Hist", "2D Hist", 200, min(data1)-10, max(data1)+20, 200, min(data2)-10, max(data2)+20);
	
	for(int i=0; i<data1.size(); i++){
		h3->Fill(data1[i],data2[i]);
	}
	h3->Draw();
	cout << "The correlation coefficient from the TH2D hist function is: " << h3->GetCorrelationFactor() << endl;
//calculate the covariance of the 2 data sets
	double covDD12 = covariance(data1,data2);
//output statistics
	cout << "Read " << Ndata1 << " numbers from " << "H_pT.txt" << endl;

	cout << "The sample mean of the dataset is " << sampleMean(data1) << endl;
	cout << "The sample variance of the dataset is " << sampleVariance(data1) << endl;
	cout << "The sample std dev is " << sampleStandardDeviation(data1) << endl;
	cout << "The median is " << median(data1) << endl;

	cout << "Read " << Ndata2 << " numbers from " << "jet_pT.txt" << endl;

	cout << "The sample mean of the dataset is " << sampleMean(data2) << endl;
	cout << "The sample variance of the dataset is " << sampleVariance(data2) << endl;
	cout << "The sample std dev is " << sampleStandardDeviation(data2) << endl;
	cout << "The median is " << median(data2) << endl;
	cout << "The covariance of these 2 data sets is " << covDD12 << " !" << endl;
	cout << "The correlation coefficient for these 2 data sets is: " << covDD12/(sampleStandardDeviation(data1)*sampleStandardDeviation(data2)) << endl;
	
//create new canvas and put in the sorted data1 values into a histogram
	TCanvas *c1 = new TCanvas();
	sort( data1.begin(), data1.end());
	TH1D *h1 = new TH1D("H-pT", "H_pT Hist", 300, -30, data1[data1.size()-1]+20);
	
	for(int i=0; i<data1.size(); i++){
		h1->Fill(data1[i]);
	}
	h1->Draw();
	
//create a new canvas and put in the sorted data2 values into a histogram
	TCanvas *c2=new TCanvas();
	sort( data2.begin(), data2.end());
	TH1D *h2 = new TH1D("jet-pT", "jet_pT Hist", 300, -30, data2[data2.size()-1]+20);
	
	for(int i=0; i<data2.size(); i++){
		h2->Fill(data2[i]);
	}
	h2->Draw();

	
	
}

//various functions for calculating the statistics of the data sets
//these functions are based off of the equations shown in the lecture slides

double covariance(vector<double> data1, vector<double> data2){
	int size = data1.size();
	double mean1 = sampleMean(data1);
	double mean2 = sampleMean(data2);
	double sum = 0;
	int count = 0;
	for (double num : data1){
		sum = sum + (num - mean1)*(data2[count] - mean2);
		count++;
	}
	sum = sum/(size-1);
	return sum;
}

double max(vector<double> data){
	double max = data[0];
	for (double num : data){
		if(num>max){
			max = num;
		}
	}
	return max;
}

double min(vector<double> data){
	double min = data[0];

	for (double num : data){
		if(num<min){
			min = num;
		}
	}
	return min;
}

double sampleMean(vector<double> data) {
	double sum = 0;
	for(double num : data){
		sum = sum + num;
	}
	if(data.size()!=0){
		return sum/data.size();
	} else {
		return 0;
	}
}

double sampleVariance(vector<double> data) {
	double mean = sampleMean(data);
	double squareSum = 0;
	for(double num : data){
		squareSum = squareSum + pow((num - mean),2);
	}
	if(data.size()!=0){
		return squareSum/(data.size()-1);
	} else {
        	return 0;
	}
}

double sampleStandardDeviation(vector<double> data) {
	double variance = sampleVariance(data);
	if(variance>0){
		return sqrt(variance);
	} else {
		return 0;
	}
}

static bool greater_than_sort(double u, double v){
        return u > v;
}


double median(vector<double> data) {
	sort(data.begin(), data.end(), greater_than_sort);
        return (data[199]+data[200])/2;
}

