#ifndef VaR_h
#define VaR_h
#include "CSVparser.hpp"

class VaR {
private:
	Parser data;
	string filename;
	double k;
	double TTM;
	long quantity;
	double Mu = 0, Sigma = 0; // mean and sample sd of return of histrocial spot rate
	double percentile;
public:
	VaR(const string &filename_, double k_, double TTM_, long quantity_, double percentile_) :data(Parser(filename_))
	{
		filename.assign(filename_);
		k = k_;
		// time to maturity as portion of 360 days
		TTM = TTM_;
		quantity = quantity_;
		percentile = percentile_;
		for (int i = 1; i < data.rowCount(); i++) Mu += stod(data[i]["r_S"]);
		Mu /= data.rowCount();
		for (int i = 1; i < data.rowCount(); i++) Sigma += pow((stod(data[i]["r_S"]) - Mu), 2);
		Sigma = sqrt(Sigma / (data.rowCount() - 1));
	}

	//given return vector and percentile, get VaR
	double risk_VaR(const vector<double> & pl, double percentile);
	//calculate VaR by historical method
	void historical();
	//calculate VaR by Monte Carlo method
	void MC_Vasicek();
	//calculate VaR by analytical method
	void analytical();
	//calculate VaR by the combination of historical and Monte Carlo
	void MC_GBM();
	
};

#endif
#pragma once
