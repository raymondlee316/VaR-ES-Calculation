#pragma warning( disable : 4819)
#ifndef VaR_h
#define VaR_h
#include "CSVparser.hpp"
#include <vector>
#include <iostream>
#include <math.h> 
#include "Bonds.h"
#include <ql/quantlib.hpp>


using namespace std;
using namespace QuantLib;

double VaR_BondPorfolio_Analytical_Yield(vector<Bonds*> IndividualCouponBonds) {
	int N = IndividualCouponBonds.size(); // number of zero coupon bonds in portfolio
	double Undiversified_VaR = 0.0;
	SequenceStatistics ss;
	vector<Real> logReturn_yield; // vector to store daily return of three assets in the requirement of Quantlib

	for (int i = 0; i < IndividualCouponBonds[0]->getLength_logReturn(); ++i) {
		logReturn_yield.clear();
		for (auto itr = IndividualCouponBonds.begin(); itr != IndividualCouponBonds.end(); ++itr){
			logReturn_yield.push_back((*itr)->getlogReturn_Yield(i));
		}
		ss.add(logReturn_yield); // Input data of log-return of yield for sequencestatistics from Quantlib libary
	}

	//// Display the correlation matrix
	//cout << "Correlation Matrix of Return" << endl;
	//cout.setf(ios_base::right, ios_base::adjustfield);
	//cout.setf(ios_base::showpoint);
	//cout.setf(ios_base::fixed, ios_base::floatfield);
	//cout.precision(10);
	//cout << ss.correlation() << "\n";

	// Caculate bond portfolio VaR
	Matrix omega(N, 1); // Investment weighting, i.e. VaR of each individual bond
	for (int i = 0; i < N; ++i) {
		omega[i][0] = IndividualCouponBonds[i]->getVaR_Analytical();
		Undiversified_VaR += omega[i][0];
	}
	Matrix Portfolio_VaR(1, 1);
	Portfolio_VaR = transpose(omega) * ss.correlation() * omega;

	double Anal_VaR = -sqrt(Portfolio_VaR[0][0]);
	double Anal_ES = -Anal_VaR / IndividualCouponBonds[0]->getZstat() * (-exp(-IndividualCouponBonds[0]->getZstat() * IndividualCouponBonds[0]->getZstat() / 2) / ((IndividualCouponBonds[0]->getAlpha() / 100)*sqrt(2 * atan(1) * 4)));
	cout << "Diversified porfolio ES  based on yield and Analytical method is: " << Anal_ES << endl;
	cout << "Diversified porfolio VaR based on yield and Analytical method is: " << Anal_VaR <<endl;
	//cout << "Undiversified porfolio VaR based on yield and Analytical method is: " << Undiversified_VaR << endl;

	return Anal_VaR;
}

double VaR_BondPorfolio_Historical_Yield(vector<Bonds*> IndividualCouponBonds) {
	int N = IndividualCouponBonds.size(); // number of zero coupon bonds in portfolio
	double Undiversified_VaR = 0.0, daily_PortfolioPL = 0.0;
	vector<double> PortfolioPL;

	for (int i = 0; i < IndividualCouponBonds[0]->getLength_logReturn(); ++i) {
		daily_PortfolioPL = 0.0;
		for (auto itr = IndividualCouponBonds.begin(); itr != IndividualCouponBonds.end(); ++itr) {
			daily_PortfolioPL += (*itr)->getYield_Change(i)*((*itr)->getMarketValue())*((*itr)->getModifiedDuration());
		}
		PortfolioPL.push_back(daily_PortfolioPL);
	}

	sort(PortfolioPL.begin(), PortfolioPL.end()); // sort P&L from loss to profit
	int ind = floor(IndividualCouponBonds[0]->getAlpha() / 100 * PortfolioPL.size()); // find the index according to percentile
	double His_VaR = PortfolioPL[ind];
	vector<double>::const_iterator first = PortfolioPL.begin();
	vector<double>::const_iterator last = PortfolioPL.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double His_ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	cout << "Porfolio ES  based on yield and Historical method is:" << His_ES << endl;
	cout << "Porfolio VaR based on yield and Historical method is:" << His_VaR << endl;
	return His_VaR;
}

double VaR_BondPorfolio_Analytical_Price(vector<Bonds*> IndividualCouponBonds) {
	int N = IndividualCouponBonds.size(); // number of zero coupon bonds in portfolio
	double Undiversified_VaR = 0.0, PortfolioValue = 0.0;
	SequenceStatistics ss;
	vector<Real> logReturn_price; // vector to store daily return of three assets in the requirement of Quantlib

	for (int i = 0; i < IndividualCouponBonds[0]->getLength_logReturn(); ++i) {
		logReturn_price.clear();
		for (auto itr = IndividualCouponBonds.begin(); itr != IndividualCouponBonds.end(); ++itr) {
			logReturn_price.push_back((*itr)->getlogReturn_Price(i));
		}
		ss.add(logReturn_price); // Input data for sequencestatistics from Quantlib libary
	}

	// Caculate bond portfolio VaR
	Matrix omega(N, 1); // Investment weighting, i.e. VaR of each individual bond
	for (int i = 0; i < N; ++i) {
		omega[i][0] = IndividualCouponBonds[i]->getMarketValue();
		Undiversified_VaR += IndividualCouponBonds[i]->getVaR_Analytical();
	}
	Matrix Portfolio_VaR(1, 1);
	Portfolio_VaR = transpose(omega) * ss.covariance() * omega;
	double Anal_VaR = (-IndividualCouponBonds[0]->getZstat()*sqrt(Portfolio_VaR[0][0]));
	double Anal_ES = -Anal_VaR / IndividualCouponBonds[0]->getZstat() * (-exp(-IndividualCouponBonds[0]->getZstat() * IndividualCouponBonds[0]->getZstat() / 2) / ((IndividualCouponBonds[0]->getAlpha() / 100)*sqrt(2 * atan(1) * 4)));
	cout << "Diversified porfolio ES  based on price and Analytical method is: " << Anal_ES << endl;
	cout << "Diversified porfolio VaR based on price and Analytical method is: " << Anal_VaR << endl;
	//cout << "Undiversified porfolio VaR based on price and Analytical method is: " << Undiversified_VaR << endl;
	return Anal_VaR;
}

double VaR_BondPorfolio_Historical_Price(vector<Bonds*> IndividualCouponBonds) {
	int N = IndividualCouponBonds.size(); // number of zero coupon bonds in portfolio
	double Undiversified_VaR = 0.0, daily_PortfolioPL = 0.0;
	vector<double> PortfolioPL;
	
	for (int i = 0; i < IndividualCouponBonds[0]->getLength_logReturn(); ++i) {
		daily_PortfolioPL = 0.0;
		for (auto itr = IndividualCouponBonds.begin(); itr != IndividualCouponBonds.end(); ++itr) {
			daily_PortfolioPL += (*itr)->getPrice_Change(i)*((*itr)->getQuantity());
		}
		PortfolioPL.push_back(daily_PortfolioPL);
	}
	
	sort(PortfolioPL.begin(), PortfolioPL.end()); // sort P&L from loss to profit
	int ind = floor(IndividualCouponBonds[0]->getAlpha() / 100 * PortfolioPL.size()); // find the index according to percentile
	double His_VaR = PortfolioPL[ind];
	vector<double>::const_iterator first = PortfolioPL.begin();
	vector<double>::const_iterator last = PortfolioPL.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double His_ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	cout << "Porfolio ES  based on price and Historical method is:" << His_ES << endl;
	cout << "Porfolio VaR based on price and Historical method is:" << His_VaR << endl;
	return His_VaR;
}

double VaR_BondPorfolio_MonteCarlo(vector<Bonds*> IndividualCouponBonds) {
	int N = IndividualCouponBonds.size(); // number of zero coupon bonds in portfolio
	double daily_PortfolioPL = 0.0;
	vector<double> PortfolioPL;

	for (int i = 0; i < IndividualCouponBonds[0]->getLength_MCPrice(); ++i) {
		daily_PortfolioPL = 0.0;
		for (auto itr = IndividualCouponBonds.begin(); itr != IndividualCouponBonds.end(); ++itr) {
			daily_PortfolioPL += (*itr)->getMCPrice_Change(i)*((*itr)->getQuantity());
		}
		PortfolioPL.push_back(daily_PortfolioPL);
	}

	sort(PortfolioPL.begin(), PortfolioPL.end()); // sort P&L from loss to profit
	int ind = floor(IndividualCouponBonds[0]->getAlpha() / 100 * PortfolioPL.size()); // find the index according to percentile
	double MC_VaR = PortfolioPL[ind];
	vector<double>::const_iterator first = PortfolioPL.begin();
	vector<double>::const_iterator last = PortfolioPL.begin() + ind - 1;
	vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
	double MC_ES = std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
	cout << "Porfolio ES  based on MonteCarlo Simulation method is:" << MC_ES << endl;
	cout << "Porfolio VaR based on MonteCarlo Simulation method is:" << MC_VaR << endl;
	return MC_VaR;
}
#endif
#pragma once
