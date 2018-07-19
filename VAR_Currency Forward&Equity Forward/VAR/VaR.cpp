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
#include "VasicekModel.h"
#include <climits>
#include <cstdint>
#include <typeinfo>

#pragma warning( disable : 4819)
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
	cout << endl <<"Expected shortfall in " << (100 - percentile) << "% is " << ES <<endl; // print out ES
	return v[ind]; // return the value of VaR
}

void VaR::historical() {
	/*
	This function is used for caculating VaR using historical method,
	it assume future change is same with historical price change
	*/

	cout << "Calculate VaR by historical method:" << endl;
		vector<double> pl; // store P&L value
		// find current domestic rate, foreign rate and spot rate
		double DomRate = stod(data[data.rowCount() - 1]["DomRate"]);
		double ForRate = stod(data[data.rowCount() - 1]["ForRate"]);
		double spotRate = stod(data[data.rowCount() - 1]["spotRate"]);
		// caculate today's portfolio value
		double oldPV = spotRate * quantity / (1 + ForRate / 100 * TTM) - k * quantity / (1 + DomRate / 100 * TTM);
		
		// caculate future rate based on historical rate change
		double oldDomRate = stod(data[0]["DomRate"]);
		double oldForRate = stod(data[0]["ForRate"]);
		double oldspotRate = stod(data[0]["spotRate"]);
		cout << setiosflags(ios::fixed) << setprecision(6);
		cout << "Exhange" << setw(15) << "PV (Dom)" << setw(15) << "PV(For)" << setw(15) << "Portfolio" << setw(15) << "P&L" << endl;
		for (int i = 1; i < data.rowCount() - 1; i++) {
			double newDomRate = stod(data[i]["DomRate"]);
			double newForRate = stod(data[i]["ForRate"]);
			double newspotRate = stod(data[i]["spotRate"]);
			// simulated rate = historical rate change + current rate
			double monteForRate = newForRate - oldForRate + ForRate;
			double monteDomRate = newDomRate - oldDomRate + DomRate;
			double monteSpotRate = newspotRate / oldspotRate * spotRate;
			// update histrocial rate at each step
			oldDomRate = newDomRate;
			oldForRate = newForRate;
			oldspotRate = newspotRate;
			double PV_For = 1 / (1 + monteForRate / 100 * TTM);
			double PV_Dom = 1 / (1 + monteDomRate / 100 * TTM);
			double port = (monteSpotRate * quantity * PV_For - k * quantity * PV_Dom); // caculate simulated portfolio value
			double retn = port - oldPV; // caculate P&L based on difference between simulated portfolio value and current portfolio value
			pl.push_back(retn); // store P&L in a vector
			cout << monteSpotRate << setw(15) << PV_Dom << setw(15) << PV_For << setw(15) << port << setw(15) << retn << endl;
		}
		cout << "The one day 95% VaR based on Historical method is: " << risk_VaR(pl, percentile) << endl << endl;
	}

void VaR::MC_GBM() {
	/*
	This function is used for caculating VaR using Monte Carlo Simulation,
	note we only simulate spotrate (exchange rate) based on GBM model,
	and we use historical method to caculate simulated Domestic Rate and Foreign Rate,
	so we didn't use any interest rate model
	*/
	cout << "calculate VaR by the combination of historical and Monte Carlo" << endl;
	vector<double> pl;
	double DomRate = stod(data[data.rowCount() - 1]["DomRate"]);
	double ForRate = stod(data[data.rowCount() - 1]["ForRate"]);
	double spotRate = stod(data[data.rowCount() - 1]["spotRate"]);
	double oldPV = spotRate * quantity / (1 + ForRate / 100 * TTM) - k * quantity / (1 + DomRate / 100 * TTM);
	double oldDomRate = stod(data[0]["DomRate"]);
	double oldForRate = stod(data[0]["ForRate"]);
	double oldspotRate = stod(data[0]["spotRate"]);
	// simulate future spotrate based on GBM
	double delta_T = 1.0 / 360; 
	long m = data.rowCount() - 2; // simulation time, note that because we only simulate spot rate, so simulation time is equal to the length of historical data
	double s0 = stod(data[data.rowCount() - 1]["spotRate"]); // current spotrate
	double rf = stod(data[data.rowCount() - 1]["DomRate"]); // risk-free rate, which happens to be current domestic rate in currency forward
	double sigma_GBM = Sigma * sqrt(360); // historical volatility for spot rate
	BSModel Spot_Rate(s0, rf, sigma_GBM);
	vector<double> Spot_Rate_s(m); // vector to store simulated spot rate
	Spot_Rate.MC_Simulation(delta_T, m, Spot_Rate_s); // update simulated spot rate vector value
	// following part is similar with Historical Method
	cout << setiosflags(ios::fixed) << setprecision(6);
	cout << "Exhange" << setw(15) << "PV (Dom)" << setw(15) << "PV(For)" << setw(15) << "Portfolio" << setw(15) << "P&L" << endl;
	for (int i = 1; i < data.rowCount() - 1; i++) {
		double newDomRate = stod(data[i]["DomRate"]);
		double newForRate = stod(data[i]["ForRate"]);
		double newspotRate = stod(data[i]["spotRate"]);
		double monteForRate = newForRate - oldForRate + ForRate;
		double monteDomRate = newDomRate - oldDomRate + DomRate;
		double monteSpotRate = Spot_Rate_s[i - 1];
		oldDomRate = newDomRate;
		oldForRate = newForRate;
		oldspotRate = newspotRate;
		double PV_For = 1 / (1 + monteForRate / 100 * TTM);
		double PV_Dom = 1 / (1 + monteDomRate / 100 * TTM);
		double port = (monteSpotRate * quantity * PV_For - k * quantity * PV_Dom);
		double retn = port - oldPV;
		pl.push_back(retn);
		cout << monteSpotRate << setw(15) << PV_Dom << setw(15) << PV_For << setw(15) << port << setw(15) << retn << endl;
	}
	cout << "The one day 95% VaR based on Monte Carlo with only GBM model is: " << risk_VaR(pl, percentile) << endl;
}


void VaR::MC_Vasicek(){
	cout << "calculate VaR by Monte Carlo method" << endl;
	/*
	This function is used for caculating VaR using Monte Carlo Simulation,
	note we not only simulate spotrate (exchange rate) based on GBM model,
	but also we use Vasicek Model to caculate simulated Domestic Rate and Foreign Rate.
	However the result is very bad, 
	so this is just to show Vasicek Model is not suitable to model interest rate
	*/
	double delta_T = 1.0 / 360; long m = 1e2;

// US Rate
double r0 = 4.468, a = 15.3537283085531, r_bar = 4.38165807849567, sigma  = 1.63896341410723;
VasicekModel US_Rate(r0, a, r_bar, sigma);
vector<double> US_Rate_r(m);
US_Rate.MC_Simulation(delta_T, m, US_Rate_r);

// UK Rate
r0 = 5.96, a = 19.5079234274301, r_bar = 6.60027905715864, sigma = 3.60951651166847;
VasicekModel UK_Rate(r0, a, r_bar, sigma);
vector<double> UK_Rate_r(m);
UK_Rate.MC_Simulation(delta_T, m, UK_Rate_r);

// Spot Rate
double s0 = 1.6542, rf = 4.658, sigma_GBM = 0.1143513;
sigma_GBM = 0.002928519;
BSModel Spot_Rate(s0, rf, sigma_GBM);
vector<double> Spot_Rate_s(m);
Spot_Rate.MC_Simulation(delta_T, m, Spot_Rate_s);

// Portfolio Value
vector<double> Portfolio(m);
vector<double> Portfolio_Loss(m);
cout << setiosflags(ios::fixed) << setprecision(6);
cout << "Exhange" << setw(15) << "PV (Dom)" << setw(15) << "PV(For)" << setw(15) << "Portfolio" << setw(15) << "P&L" << endl;

for (int i = 0;i < m;i++)
{
	Portfolio[i] = quantity* Spot_Rate_s[i] * (1 / (1 + UK_Rate_r[i] / 100 * TTM)) - k * quantity *(1 / (1 + US_Rate_r[i] / 100 * TTM));
	Portfolio_Loss[i] = Portfolio[i] - Portfolio[0];
	double PV_For = 1 / (1 + UK_Rate_r[i] / 100 * TTM);
	double PV_Dom = 1 / (1 + US_Rate_r[i] / 100 * TTM);
	cout << Spot_Rate_s[i] << setw(15) << PV_Dom << setw(15) << PV_For << setw(15) << Portfolio[i] << setw(15) << Portfolio_Loss[i] << endl;

}
sort(Portfolio_Loss.begin(), Portfolio_Loss.end());
sort(Portfolio.begin(), Portfolio.end());

vector<double>::iterator itr = Portfolio_Loss.begin();
double VaR = *(itr + m * 0.05);
cout << "The one day 95% VaR based on Monte Carlo with Vasciek model is:" << VaR << endl << endl;
cout << "Simply show Vasciek Model is not suitable to use. Ignore these parts." << endl;
}

 

void VaR::analytical() {
	//////////////////// Analytical Method ////////////////////////////////
	cout << "calculate VaR by analytical method" << endl;
	/*
	This function is used for caculating VaR using Analytical Method,
	note we first decompose the currency forward into three individual assets,
	that is, the forward contract is equivalent to
	1. A long position of (SP*) on the spot rate
	2. A long position of (SP*) in the foreign bill
	3. A short position of (KP) in the domestic bill (borrowing)
	where S = current spot price, P/P* = current domestic/foreign zero coupon bond price,
	K = strike price (exchange rate) in the contract.

	Then we assume the returns of these three assets follow normal distribution,
	and we caculate mean and sd of portfolio return based on these three assets and covariance matrix
	VaR = mu - alpha*sigma

	Note we cannot caculate Expected Shortfall in Analytical Method.
	*/

	// Collect the daily return data, note the data is pre-processed in Excel 
	ifstream filein(filename);
	int t;
	double r_S, r_PVf, r_PVd, weight1, weight2;
	// Construct the sequencestatistics from Quantlib libary
	SequenceStatistics ss;
	vector<Real> daily_return_indexes; // vector to store daily return of three assets in the requirement of Quantlib
	for (string line; getline(filein, line);)
	{
		stringstream lineStream(line);
		string cell;
		int cell_idx = 0;

		while (getline(lineStream, cell, ',')) {
			// Classify the cell into different fields
			switch (cell_idx) {
			case 0:
				t = atoi(cell.c_str());
				if (t<1 || t>500) t = -1; // Mark the invalid entry as -1
				break;
			case 9:      r_S = atof(cell.c_str());      break; // get return of spot rate in the 9th column
			case 10:      r_PVf = atof(cell.c_str());      break; // get return of foreign bill in the 10th column
			case 11:      r_PVd = atof(cell.c_str());      break; // get return of domestic bill
			case 13:	  weight1 = atof(cell.c_str());      break; // get weight1, i.e. (SP*)
			case 14:	  weight2 = -k * atof(cell.c_str());      break; // get weight1, i.e. (KP)
			}
			cell_idx++;
		}
		// Skip the invalid data
		if (t == -1) continue;
		daily_return_indexes.clear();
		daily_return_indexes.push_back(r_S);
		daily_return_indexes.push_back(r_PVf);
		daily_return_indexes.push_back(r_PVd);
		ss.add(daily_return_indexes); // Input data for sequencestatistics from Quantlib libary
	}

	// Calculate the covariance matrix
	cout << "Covariance Matrix with equally weighting of the daily returns" <<endl;
	cout.setf(ios_base::right, ios_base::adjustfield);
	cout.setf(ios_base::showpoint);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout.precision(10);
	cout << ss.correlation() << "\n";

	//// Calculate the portfolio variance using normalized weights (sum of weights is 1)
	//// The result is not as good as using unnormalized weights
	//Matrix omega(3, 1); // Investment weighting
	//omega[0][0] = quantity * weight1 / (2*weight1 + weight2); // S
	//omega[1][0] = quantity * weight1 / (2*weight1 + weight2); // PVf
	//omega[2][0] = quantity * weight2 / (2*weight1 + weight2); // PVd

	// Calculate the portfolio variance using unnormalized weights
	Matrix omega(3, 1); // Investment weighting
	omega[0][0] = quantity * weight1; // S
	omega[1][0] = quantity * weight1; // PVf
	omega[2][0] = quantity * weight2; // PVd

	Matrix R(3, 1); // mean of the returns of three assets
	R[0][0] = ss.mean()[0];
	R[1][0] = ss.mean()[1];
	R[2][0] = ss.mean()[2];

	Matrix mean_portfolio(1, 1);
	mean_portfolio = transpose(omega) * R; // mean of portfolio return = weight*(individual asset return)
	cout << "Portfolio mean: " << mean_portfolio[0][0] << "\n";

	Matrix variance_portfolio(1, 1);
	variance_portfolio = transpose(omega) * ss.covariance() * omega; // variance of portfolio return = weight*(covariance matrix)*weight
	cout << "Portfolio variance: " << variance_portfolio[0][0] << "\n";

	// Obtain z statistics at 95%
	boost::math::normal_distribution<> d(0, 1);
	double z_stat = quantile(d, 1-percentile/100);

	// Calculate the one day 95% VaR
	double VaR = mean_portfolio[0][0] - z_stat * sqrt(variance_portfolio[0][0]); // caculate VaR based on portflio return and variance
	cout << "The one day 95% VaR based on Analytical method is: " << VaR << endl <<endl;
}



