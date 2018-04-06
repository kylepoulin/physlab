#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>


double sampleMean(vector<double> data);
double sampleVariance(vector<double> data);
double sampleStandardDeviation(vector<double> data);
double median(vector<double> data);

void exercise1a(){
	vector<double> data;
	
	ifstream infile;
	infile.open("jet1_pT.txt");

	double value;
	if (infile.is_open()){
		while (infile>>value){
		data.push_back(value);
		}
	}
	infile.close();
	int Ndata = data.size();
	cout.precision(3);
	cout << "Read " << Ndata << " numbers from the file!" << endl;

	cout << "The sample mean of the dataset is " << sampleMean(data) << endl;
	cout << "The sample variance of the dataset is " << sampleVariance(data) << endl;
	cout << "The sample std dev is " << sampleStandardDeviation(data) << endl;
	cout << "The median is " << median(data) << endl;
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

