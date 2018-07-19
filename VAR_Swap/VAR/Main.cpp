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
	const string filename = "Swap5Y.csv";
	long quantity = 100; // quantity of position, i.e. notional
	int Maturity = 0, Tenor = 5; // Maturity and tenor of swap, e.g. 7x6 swap means the swap will begin in 7 years and last 6 years
	double r0 = 0.18400 / 100; // r0 is current short rate, we use Libor overnight rate in example
	double fixedRate = 2.9161/100; // fixed rate of swap
	double percentile = 5; // 	percentile: the 0-100 percentile Domer chooses (eg: 99 % VaR : percentile = 1; 95 % VaR: percentile = 5)
	
	VaR VaR(filename, quantity, r0, fixedRate, percentile, Maturity, Tenor);
	VaR.Historical();
	VaR.Analytical();
	VaR.MonteCarlo();
	VaR.PrintResult(); // print result
	//VaR.ExportResult(); // export result to csv file

	system("pause");
	return 0;
}