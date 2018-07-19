#ifndef VaR_h
#define VaR_h
#include "CSVparser.hpp"



class VaR {
private:
	Parser data;
	string filename;
	long quantity;
	double rf; // rf is 3-month T-bill rate of current date
public:
	VaR(const string &filename_, long quantity_, double rf_) :data(Parser(filename_))
	{
		filename.assign(filename_);
		quantity = quantity_;
		rf = rf_;
	}
	//given return vector and percentile, get VaR
	double risk_VaR(const vector<double> & PL, double percentile);
	//calculate VaR by historical method
	void historical();
	//calculate VaR by analytical method
	void analytical();
	void analytical_Greeks();
	//caculate VaR by monte carlo simulation
	void MC_GBM();
};

#endif
#pragma once
