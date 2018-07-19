#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cmath>
#include <string>
#include "CSVparser.hpp"
#include <typeinfo>
#include "VaR.h""
using namespace std;

int main()
{
	// Try different data sets in both following lines
	//const string filename = "Data_1y.csv"; 
	// const string filename = "Data.csv";
	const string filename = "Data_Equity Forward.csv"; // Equity Forward
	double k = 1.65; // stike price for exchange rate
	double TTM = 3.0 / 12; // time to maturity
	long quantity = 1e7; // quantity of position
	double percentile = 5; // 95% confidence level
	VaR VaR(filename, k, TTM, quantity, percentile);
	VaR.historical();
	VaR.analytical();
	VaR.MC_Vasicek();
	VaR.MC_GBM();
	system("pause");
}

