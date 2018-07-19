#pragma warning( disable : 4819)
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <ql/quantlib.hpp>
#include <string>
#include <climits>
#include <cstdint>
#include <typeinfo>

#include <fstream> 
#include <string>
#include <iostream>
#include <streambuf> 

#include "VaR.h"
#include "Swaption_MC.h"

using namespace std;
using namespace QuantLib;

void VaR::Historical() {
	cout << "Calculating VaR by Historical method" << endl;
	/*
	This function is used for caculating VaR using historical method,
	it assume future change is same with historical price change
	*/
	vector<double> PL; // vector to store PnL of option value
	double OldPrice = stod(data[0]["OPTION_PRICE"]);

	for (int i = 1; i < data.rowCount() - 1; i++) {
		double NewPrice = stod(data[i]["OPTION_PRICE"]);
		double SimulatedPrice = NewPrice - OldPrice + SpotPrice;
		OldPrice = NewPrice;
		double pl = SimulatedPrice * quantity - SpotPrice * quantity;
		PL.push_back(pl);
	}
	vector<double> v = PL;
	sort(v.begin(), v.end()); // sort P&L from loss to profit
	int ind = floor(percentile / 100 * v.size()); // find the index according to percentile
	vector<double>::const_iterator first = v.begin();
	vector<double>::const_iterator last = v.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	His_VaR = v[ind];
	His_ES = ES;
	cout << endl;
}


void VaR::Analytical() {
	//////////////////// Analytical Method ////////////////////////////////
	/*
	This function is used for caculating VaR using Analytical method,
	it assume the return of option price follows normal distribution
	*/
	// Caculate mean and sample sd of option return
	double Mu = 0, Sigma = 0;
	for (int i = 0; i < data.rowCount(); i++) Mu += stod(data[i]["OPTION_RETURN"]);
	Mu /= data.rowCount(); // mean value of the return of option
	for (int i = 0; i < data.rowCount(); i++) Sigma += pow((stod(data[i]["OPTION_RETURN"]) - Mu), 2);
	Sigma = sqrt(Sigma / (data.rowCount() - 1)); // standard deviation of the return

	cout << "Calculating VaR by Analytical method (mu-alpha*sigma)" << endl;
	// Obtain z statistics at 95%
	boost::math::normal_distribution<> d(0, 1);
	double z_stat = quantile(d, 1 - percentile / 100.0);
	// Calculate the one day 95% VaR
	Anal_VaR = SpotPrice * quantity * (Mu - z_stat * Sigma);
	Anal_ES = - Anal_VaR / z_stat * (-exp(-z_stat * z_stat / 2) / ((percentile / 100)*sqrt(2 * atan(1) * 4)));
	cout << endl;
}

void VaR::MonteCarlo() {
	//////////////////// MonteCarlo Simulation Method ////////////////////////////////
	/*
	This function is used for caculating VaR using Monte Carlo simulation method,
	first use Hull-White One Factor model in Swaption_MC.h to simulate forward swap rate, then apply delta-gamma normal
	approch to calculate.
	Then use historical implied volatility and Vega to calculate the VaR and ES contributed by implied volatility with
	Historical method.
	Finally we combine them into a delta-gamma-vega normal approach in which we only simulate forward swap rate.
	*/
	cout << "Calculating VaR by Monte Carlo method (Hull-White One Factor Model)" << endl;
	vector<double> v = Swaption_MC(r0, strike, Delta, Gamma, Maturity, Tenor, quantity, SpotPrice); // get PnL of swaption value contributed only by simulated swap rate, i.e. delta-gamma normal approach
	sort(v.begin(), v.end()); // sort P&L from loss to profit
	int ind = floor(percentile / 100 * v.size()); // find the index according to percentile
	vector<double>::const_iterator first = v.begin();
	vector<double>::const_iterator last = v.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	MC_VaR = v[ind];
	MC_ES = ES;

	vector<double> PL; // vector to store PnL of swaption value contributed only by historical implied volatility
	double OldIV = stod(data[0]["IMLIED_VOLATILITY"]);
	for (int i = 1; i < data.rowCount() - 1; i++) {
		double NewIV = stod(data[i]["IMLIED_VOLATILITY"]);
		double DeltaIV = NewIV - OldIV;
		OldIV = NewIV;
		double pl = DeltaIV * quantity * Vega;
		PL.push_back(pl);
	}
	vector<double> v2 = PL;
	sort(v2.begin(), v2.end()); // sort P&L from loss to profit
	int ind2 = floor(percentile / 100 * v2.size()); // find the index according to percentile
	vector<double>::const_iterator first2 = v2.begin();
	vector<double>::const_iterator last2 = v2.begin() + ind2 - 1;
	vector<double> ES_vector2(first2, last2); // use a new vector to store losses that are beyond VaR
	double ES2 = std::accumulate(ES_vector2.begin(), ES_vector2.end(), 0.0) / ES_vector2.size(); // caculate expected shortfall (average of ES_vector)
	MC_VaR += v2[ind2]; // sum the VaR from delta-normal approach and the VaR from vega-normal approach
	MC_ES += ES2;// sum the ES from delta-normal approach and the VaR from vega-normal approach

	cout << endl;
}

void::VaR::PrintResult() {
	/*
	This function is used for printing the result to screen.
	*/
	cout << Maturity << "x" << Tenor << " Swaption" << endl
		<< "Confidence level: " << 100 - percentile << "%" << endl
		<< "Current value: " << SpotPrice*quantity << endl;
	if (His_VaR*His_ES) {
		cout << "Historical VaR = " << His_VaR << endl
			 << "Historical ES  = " << His_ES << endl;
	}
	else cout << "Haven't calculate by Historical method yet." << endl;
	if (Anal_VaR*Anal_ES) {
		cout << "Analytical VaR = " << Anal_VaR << endl
			 << "Analytical ES  = " << Anal_ES << endl;
	}
	else cout << "Haven't calculate by Analytical method yet." << endl;
	if (MC_VaR*MC_ES) {
		cout << "MonteCarlo VaR = " << MC_VaR << endl
			 << "MonteCarlo ES  = " << MC_ES << endl;
	}
	else cout << "Haven't calculate by Monte Carlo method yet." << endl;
}

void::VaR::ExportResult() {
	/*
	This function is used for exporting the result to csv file.
	*/
	ofstream oFile;

	oFile.open("result.csv", ios::out | ios::trunc);

	oFile << Maturity << "x" << Tenor << " Swaption" << endl
		<< "Confidence level: " << 100 - percentile << "%" << endl
		<< "Current value: " << SpotPrice * quantity << endl;
	if (His_VaR*His_ES) {
		oFile << "Historical VaR = " << His_VaR << endl
			<< "Historical ES  = " << His_ES << endl;
	}
	else oFile << "Haven't calculate by Historical method yet." << endl;
	if (Anal_VaR*Anal_ES) {
		oFile << "Analytical VaR = " << Anal_VaR << endl
			<< "Analytical ES  = " << Anal_ES << endl;
	}
	else oFile << "Haven't calculate by Analytical method yet." << endl;
	if (MC_VaR*MC_ES) {
		oFile << "MonteCarlo VaR = " << MC_VaR << endl
			<< "MonteCarlo ES  = " << MC_ES << endl;
	}
	else oFile << "Haven't calculate by Monte Carlo method yet." << endl;
}