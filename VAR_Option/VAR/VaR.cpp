#pragma warning( disable : 4819)
#include <iostream>
#include <vector>
#include <cmath>
#include "VaR.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <ql/quantlib.hpp>
#include <string>
#include "BSModel.h"
#include <climits>
#include <cstdint>
#include <typeinfo>


using namespace std;
using namespace QuantLib;
double VaR::risk_VaR(const vector<double> & pl, double percentile)
{
	/*
	this function returns the Value-at-Risk at certain percentile,
	and print out the expected shortfall, given the vector of return distribution
	typically, when there is a loss, VaR is reported as a negative number
	input:
	pl: profit&loss distribution
	percentile: the 0-100 percentile Domer chooses
	(eg: 99% VaR: percentile=1; 95% VaR: percentile=5)
	*/
	vector<double> v = pl;
	sort(v.begin(), v.end()); // sort P&L from loss to profit
	int ind = floor(percentile / 100 * pl.size()); // find the index according to percentile
	vector<double>::const_iterator first = v.begin();
	vector<double>::const_iterator last = v.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	cout << "Expected shortfall in " << (100 - percentile) << "% is " << ES << endl; // print out ES
	return v[ind]; // return the value of VaR
}


void VaR::historical() {
	cout << "Calculate VaR by historical method:" << endl;
	/*
	This function is used for caculating VaR using historical method,
	it assume future change is same with historical price change
	*/
	vector<double> PL;
	double SpotPrice = stod(data[data.rowCount() - 1]["OPTION_PRICE"]); // Price at today
	double OldPrice = stod(data[0]["OPTION_PRICE"]);
	double oldPV = SpotPrice * quantity;

	cout << setiosflags(ios::fixed) << setprecision(2);
	cout << "Option" << setw(15) << "Portfolio" << setw(15) << "P&L" << endl;
	for (int i = 1; i < data.rowCount() - 1; i++) {
		double NewPrice = stod(data[i]["OPTION_PRICE"]);
		double SimulatedPrice = NewPrice - OldPrice + SpotPrice;
		OldPrice = NewPrice;
		double port = SimulatedPrice * quantity;
		double pl = port - oldPV;
		PL.push_back(pl);
		cout << SimulatedPrice << setw(15) << port << setw(15) << pl << endl;
	}
	cout << "The one day 95% VaR based on Historical method is: " << risk_VaR(PL, 5) << endl << endl;
}

void VaR::MC_GBM() {
	cout << "Calculate VaR by Monte Carlo without simulating volatility (GBM only)" << endl;
	vector<double> PL;
	double SpotStockPrice = stod(data[data.rowCount() - 1]["STOCK_PRICE"]);
	double SpotOptionPrice = stod(data[data.rowCount() - 1]["OPTION_PRICE"]);
	double oldPV = SpotOptionPrice * quantity;
	//GBM
	double delta_T = 1.0 / 252; 
	long m = data.rowCount(); // Simulation times

	// Caculate sample sd of stock return
	double Mu = 0, Sigma = 0;
	for (int i = 0; i < data.rowCount(); i++) Mu += stod(data[i]["STOCK_RETURN"]);
	Mu /= data.rowCount();
	for (int i = 0; i < data.rowCount(); i++) Sigma += pow((stod(data[i]["STOCK_RETURN"]) - Mu), 2);
	Sigma = sqrt(Sigma / (data.rowCount() - 1));

	double s0 = SpotStockPrice, sigma_GBM = Sigma * sqrt(data.rowCount());
	BSModel Stock(s0, rf, sigma_GBM);
	vector<double> SimulatedPricePath(m);
	Stock.MC_Simulation(delta_T, m, SimulatedPricePath);
	for (int i = 1; i < m; i++) {
		double dS = SimulatedPricePath[i] - SimulatedPricePath[i-1];
		double dSigma = stod(data[i]["IMPLIED_VOLATILITY"]) - stod(data[i - 1]["IMPLIED_VOLATILITY"]);
		double dF = stod(data[i]["DELTA"])*dS + 0.5*stod(data[i]["GAMMA"])*dS*dS + stod(data[i]["VEGA"])*dSigma;
		double pl = dF * quantity;
		PL.push_back(pl);
	}
	cout << "The one day 95% VaR based on Monte Carlo with only GBM model is: " << risk_VaR(PL, 5) << endl << endl;
}

void VaR::analytical() {
	//////////////////// Analytical Method ////////////////////////////////
	double SpotPrice = stod(data[data.rowCount() - 1]["OPTION_PRICE"]);
	double oldPV = SpotPrice * quantity;

	// Caculate mean and sample sd of option return
	double Mu = 0, Sigma = 0;
	for (int i = 0; i < data.rowCount(); i++) Mu += stod(data[i]["OPTION_RETURN"]);
	Mu /= data.rowCount();
	for (int i = 0; i < data.rowCount(); i++) Sigma += pow((stod(data[i]["OPTION_RETURN"])-Mu),2);
	Sigma = sqrt(Sigma / (data.rowCount() - 1));

	cout << "Calculate VaR by analytical method (mu-alpha*sigma)" <<endl;
	// Obtain z statistics at 95%
	boost::math::normal_distribution<> d(0, 1);
	double z_stat = quantile(d, 0.95);
	// Calculate the one day 95% VaR
	double VaR = oldPV *(Mu - z_stat * Sigma);
	cout << "The one day 95% VaR based on Analytical method is: " << VaR << endl <<endl;
}


void VaR::analytical_Greeks() {
	cout << "Calculate VaR by using greeks letter:" << endl;
	vector<double> PL;
	double SpotPrice = stod(data[data.rowCount() - 1]["OPTION_PRICE"]); // Price at today
	double OldPrice = stod(data[0]["OPTION_PRICE"]);
	double oldPV = SpotPrice * quantity;

	/*cout << setiosflags(ios::fixed) << setprecision(2);
	cout << "dS" << setw(15)<< "Delta"<<setw(15) << "Gamma" << setw(15) << "dSigma" << setw(15) <<  "Vega" << setw(15)<< "P&L" << endl;*/
	for (int i = 1; i < data.rowCount(); i++) {
		double dS = stod(data[i]["STOCK_PRICE"]) - stod(data[i - 1]["STOCK_PRICE"]);
		double dSigma = stod(data[i]["IMPLIED_VOLATILITY"]) - stod(data[i - 1]["IMPLIED_VOLATILITY"]);
		double dF = stod(data[i]["DELTA"])*dS + 0.5*stod(data[i]["GAMMA"])*dS*dS + stod(data[i]["VEGA"])*dSigma;
		double pl = dF * quantity;
		PL.push_back(pl);
		//cout << dS << setw(15) << data[i]["DELTA"] << setw(15) << data[i]["GAMMA"] << setw(15) << dSigma << data[i]["VEGA"] << setw(15)<<setw(15) << pl << endl;
	}
	cout << "The one day 95% VaR based on Analytical method with Greeks letter is: " << risk_VaR(PL, 5) << endl << endl;
}