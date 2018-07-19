#pragma warning( disable : 4819)
#ifndef Bonds_h
#define Bonds_h
#include "CSVparser.hpp"
#include <vector>
#include <iostream>
#include <math.h> 
#include "Solver.h"
#include "MonteCarlo.h"
#include <ql/quantlib.hpp>



class Bonds {
protected:
	Parser data; // data set
	string filename, maturity;
	long quantity;
	double alpha, vol_yield, vol_price, MD, MV, yield0, price0, z_stat;
	int frequency;
	// alpha: confidence level; vol_yield/vol_price: volatility; MD: modified duration; MV: market value;
	// yield0/price0: current yield/price; z_stat: z statistics at alpha%
	vector<double> yield, price, logReturn_yield, logReturn_price; // vector of daily data for bond yield, bond price, log-return of yield and log-return of price
	vector<double> yield_change, price_change; // vector of data for daily change of yield and price
	double Anal_VaR, His_VaR, Anal_ES, His_ES; // VaR and ES

	vector<double> MC_price, MC_price_change;
	double MC_ES, MC_VaR;
public:
	Bonds(const string &filename_, long quantity_, string maturity_, double alpha_, int frequency_) :data(Parser(filename_)){
		filename.assign(filename_);
		quantity = quantity_;
		maturity = maturity_;
		alpha = alpha_;
		 //Obtain z statistics at alpha%
		boost::math::normal_distribution<> d(0, 1);
		z_stat = quantile(d, 1-alpha/100);
		//z_stat = 1.65;
		frequency = frequency_;
	}
	Bonds(const Bonds &bonds_){
		filename.assign(bonds_.filename);
	}
	Bonds() {};
	~Bonds() {};

	void setYield() {
		string Bond_Type = maturity + " Yr"; // find column name in data set
		for (int i = 0; i < data.rowCount(); ++i) yield.push_back((double)stod(data[i][Bond_Type])); // store data into vector
		yield0 = yield.back(); // current yield
		for (int i = 1; i < yield.size(); ++i) {
			logReturn_yield.push_back(log(yield[i] / yield[i - 1])); // calculate log-return of yield
			yield_change.push_back(yield[i] - yield[i - 1]); // calcualte daily yield change
		}
		vol_yield = StandardDeviation(logReturn_yield); // calculate daily volatility of yield
	}
	void setPrice() {
		string Bond_Type = maturity + " Yr";
		for (int i = 0; i < data.rowCount(); ++i) price.push_back((double)stod(data[i][Bond_Type]));
		price0 = price.back(); // current price
		for (int i = 1; i < price.size(); ++i) {
			logReturn_price.push_back(log(price[i] / price[i - 1])); // calculate daily log-return of price
			price_change.push_back(price[i] - price[i - 1]); // calcualte daily price change
		}
		vol_price = StandardDeviation(logReturn_price); // calculate daily volatility of price
	}

	virtual void Calculation_PriceBased() = 0;
	virtual void Calculation_YieldBased() = 0;

	double VaR_Analytical_Yield() {
		Anal_VaR = - MV * MD * yield0 * vol_yield * z_stat; // calculate VaR based on Analytical method and yield data
		Anal_ES = - Anal_VaR / z_stat * (-exp(-z_stat * z_stat / 2) / ((alpha/100)*sqrt(2 * atan(1) * 4)));// calculate ES based on Analytical method and yield data

		cout << "Analytical ES  of " << maturity << " year bond based on yield is: " << Anal_ES << endl;
		return Anal_VaR;
	}
	double VaR_Analytical_Price() {
		Anal_VaR = - MV * vol_price * z_stat;
		Anal_ES = - Anal_VaR / z_stat * (-exp(-z_stat * z_stat / 2) / ((alpha/100)*sqrt(2 * atan(1) * 4)));
		cout << "Analytical ES  of " << maturity << " year bond based on price is: " << Anal_ES << endl;
		return Anal_VaR;
	}

	double VaR_Historical_Yield() {
		sort(yield_change.rbegin(), yield_change.rend()); // sort yield from positive to negative (positive yield means loss)
		int ind = floor(alpha / 100 * yield_change.size()); // find the index according to percentile
		vector<double>::const_iterator first = yield_change.begin();
		vector<double>::const_iterator last = yield_change.begin() + ind - 1;
		vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
		His_ES = -MV * MD * std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
		cout << "Historical ES  of " << maturity << " year bond based on yield is: " << His_ES << endl;
		His_VaR = - MV * MD * yield_change[ind];
		return His_VaR; // return the value of VaR
	}
	double VaR_Historical_Price() {
		sort(price_change.begin(), price_change.end()); // sort P&L from loss to profit
		int ind = floor(alpha / 100 * price_change.size()); // find the index according to percentile
		vector<double>::const_iterator first = price_change.begin();
		vector<double>::const_iterator last = price_change.begin() + ind - 1;
		vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
		His_ES = quantity * std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
		cout << "Historical ES  of " << maturity << " year bond based on price is: " << His_ES << endl;
		His_VaR = price_change[ind] * quantity;
		return His_VaR; // return the value of VaR
	}

	////////////// Following functions are used to get private data of class Bond ///////////////////// 
	double getPrice(int i) {
		return price[i];
	}
	double getYield(int i) {
		return yield[i];
	}
	double getQuantity() {
		return quantity;
	}

	double getlogReturn_Yield(int i) {
		return logReturn_yield[i];
	}
	double getlogReturn_Price(int i) {
		return logReturn_price[i];
	}

	double getPrice_Change(int i) {
		return price_change[i];
	}
	double getMCPrice_Change(int i) {
		return MC_price_change[i];
	}
	double getYield_Change(int i) {
		return yield_change[i];
	}

	int getLength_logReturn() {
		return logReturn_yield.size();
	}
	int getLength_MCPrice() {
		return MC_price.size();
	}
	double getVol_Price() {
		return vol_price;
	}
	double getMarketValue() {
		return MV;
	}
	double getModifiedDuration() {
		return MD;
	}
	double getZstat(){
		return z_stat;
	}
	double getAlpha() {
		return alpha;
	}
	double getVaR_Analytical() {
		return Anal_VaR;
	}
	double getVaR_Historical() {
		return His_VaR;
	}
	double getES_Analytical() {
		return Anal_ES;
	}
	double getES_Historical() {
		return His_ES;
	}
};

class ZeroCouponBonds : public Bonds {

public:
	ZeroCouponBonds(const string &filename_, long quantity_, string maturity_, double alpha_, double frequency_) : Bonds(filename_, quantity_, maturity_, alpha_, frequency_) {}

	void Calculation_YieldBased() {
		MD = stod(maturity) / (1 + yield0 / frequency); // calculate modified duration of zero-coupon bond (T/(1+y))
		for (auto itr = yield.begin();itr != yield.end(); ++itr) {
			double value = 1 / pow((1 + *itr / frequency), (frequency*(double)stod(maturity)));
			price.push_back(value); // calculate bond price based on yield data
		}
		price0 = price.back(); // current price
		MV = quantity * price0; // current market value of zero-coupon bond
		for (int i = 1; i < price.size(); ++i) {
			logReturn_price.push_back(log(price[i] / price[i - 1]));
			price_change.push_back(price[i] - price[i - 1]);
		}
		vol_price = StandardDeviation(logReturn_price);
	}

	void Calculation_PriceBased() {
		price0 = price.back();
		for (auto itr = price.begin();itr != price.end(); ++itr) {
			double value = (exp(log(1 / (*itr)) / (frequency* (double)stod(maturity))) - 1)*frequency; // calculate yield based on price
			yield.push_back(value);
		}
		yield0 = yield.back();
		MD = (double)stod(maturity) / (1 + yield0 / frequency);
		MV = quantity * price0;
		for (int i = 1; i < yield.size(); ++i) {
			logReturn_yield.push_back(log(yield[i] / yield[i - 1]));
			yield_change.push_back(yield[i] - yield[i - 1]);
		}
		vol_yield = StandardDeviation(logReturn_yield);
	}

	double VaR_MonteCarlo(int numPath, const string filename_df, double r0) {
		Parser DF = Parser(filename_df);
		vector<double> dfs(1,1.0);
		for (int i = 1; i < DF.columnCount(); i++) {
			dfs.push_back(stod(DF[DF.rowCount()-1][i]));
		}
		MC_price = MonteCarlo_ZeroCouponBond(numPath, stod(maturity), r0, dfs);
		for (int i = 0; i < MC_price.size(); ++i) {
			MC_price_change.push_back(MC_price[i] - price0); // calcualte daily price change
		}
		sort(MC_price_change.begin(), MC_price_change.end()); // sort P&L from loss to profit
		int ind = floor(alpha / 100 * MC_price_change.size()); // find the index according to percentile
		vector<double>::const_iterator first = MC_price_change.begin();
		vector<double>::const_iterator last = MC_price_change.begin() + ind - 1;
		vector<double> ES_vector(first, last); // use a new vector to store losses that are beyond VaR
		MC_ES = quantity * std::accumulate(ES_vector.begin(), ES_vector.end(), 0.0) / ES_vector.size(); // caculate expected shortfall (average of ES_vector)
		cout << "Simulated  ES  of " << maturity  << " year bond by Monte Carlo is: " << MC_ES << endl;
		MC_VaR = MC_price_change[ind] * quantity;
		return MC_VaR; // return the value of VaR
	}

};

class CouponBonds : public Bonds {
private:
	double coupon; // coupon at each period based (instead of annually) on $1 face value
	double F = 1 ; int T = stoi(maturity)*frequency;
	
public:
	CouponBonds(const string &filename_, long quantity_, string maturity_, double alpha_, double frequency_, double coupon_) : Bonds(filename_, quantity_, maturity_, alpha_, frequency_) {
		coupon = coupon_ / frequency_;
	}

	double Value(double y) {
		double result = 0;
		for (int n = 1;n <= T;n++)
			result += coupon * 1 / pow((1 + y / frequency), n);
		result += F * 1 / pow((1 + y / frequency), T);
		return result;
	}

	double ModifiedDuration(double y, double p) {
		double result = 0;
		for (int n = 1;n <= T;n++)
			result += coupon * n / pow((1 + y / frequency), (n+1));
		result += F * T / pow((1 + y / frequency), (T+1));
		return result/p/frequency;
	}

	void Calculation_YieldBased() {
		for (auto itr = yield.begin();itr != yield.end(); ++itr) {
			double priceValue = Value(*itr);
			price.push_back(priceValue); // calculate bond price based on yield data
		}
		price0 = price.back(); // current price
		MD = ModifiedDuration(yield0, price0);
		MV = quantity * price0; // current market value of zero-coupon bond
		for (int i = 1; i < price.size(); ++i) {
			logReturn_price.push_back(log(price[i] / price[i - 1]));
			price_change.push_back(price[i] - price[i - 1]);
		}
		vol_price = StandardDeviation(logReturn_price);
	}



	void Calculation_PriceBased() {
		price0 = price.back();
		double Acc = 0.00001;
		double LEnd = 0.0;
		double REnd = 1.0;
		double Guess = 0.04;

		for (auto itr = price.begin();itr != price.end(); ++itr) {
			double yieldValue = SolveByBisect<CouponBonds>(this, *itr, LEnd, REnd, Acc); 
			yield.push_back(yieldValue);
		}
		yield0 = yield.back();
		MD = ModifiedDuration(yield0, price0);
		MV = quantity * price0;
		for (int i = 1; i < yield.size(); ++i) {
			logReturn_yield.push_back(log(yield[i] / yield[i - 1]));
			yield_change.push_back(yield[i] - yield[i - 1]);
		}
		vol_yield = StandardDeviation(logReturn_yield);
	}

};
#endif

