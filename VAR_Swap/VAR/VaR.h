#ifndef VaR_h
#define VaR_h
#include "CSVparser.hpp"

class VaR {
private:
	Parser data;
	string filename;
	long quantity;
	double r0, fixedRate;
	int Maturity, Tenor;
	double percentile; 
	double SpotPrice;
	double His_VaR = 0, Anal_VaR = 0, MC_VaR = 0, His_ES = 0, Anal_ES = 0 , MC_ES = 0;
public:
	VaR(const string &filename_, long quantity_, double r0_, double fixedRate_, double percentile_, int Maturity_, int Tenor_) :data(Parser(filename_))
	{
		filename.assign(filename_);
		quantity = quantity_;
		r0 = r0_;
		fixedRate = fixedRate_;
		percentile = percentile_;
		Maturity = Maturity_;
		Tenor = Tenor_;
		SpotPrice = stod(data[data.rowCount() - 1]["NPV"]);
	}
	//calculate VaR by historical method
	void Historical();
	//calculate VaR by analytical method
	void Analytical();
	//caculate VaR by monte carlo simulation
	void MonteCarlo();
	//print calculation result
	void PrintResult();
	//Export result to csv file
	void ExportResult();
};

#endif
#pragma once
