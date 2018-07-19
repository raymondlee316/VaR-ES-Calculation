#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cmath>
#include <string>
#include <map>
#include <typeinfo>
#include "VaR.h"
#include "CSVparser.hpp"

using namespace std;

int main()
{
	const string filename = "SwaptionIV_8x3.csv";
	long quantity = 1; // quantity of position
	double Delta = 0.4951, Gamma = 5.44/1e4, Vega = 11921.27/1e4; // Greeks data are based on 1 bps
	int Maturity = 8, Tenor = 3; // Maturity and tenor of swaption, e.g. 7x6 swaption has maturity 7 years, tenor 6 years.
	double r0 = 0.12550 / 100; // r0 is current short rate, we use Libor overnight rate in example
	double strike = 5.1016/ 100; // strike of swaption
	double percentile = 5; // 	percentile: the 0-100 percentile Domer chooses (eg: 99 % VaR : percentile = 1; 95 % VaR: percentile = 5)
	
	VaR VaR(filename, quantity, r0, strike, percentile, Delta, Gamma, Vega, Maturity, Tenor);
	VaR.Historical();
	VaR.Analytical();
	VaR.MonteCarlo();
	VaR.PrintResult(); // print result
	VaR.ExportResult(); // export result to csv file

	system("pause");
	return 0;
}