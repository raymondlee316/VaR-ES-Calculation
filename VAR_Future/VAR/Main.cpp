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
	const string filename = "SP 500 Futures Jun 18.csv"; 
	long quantity = 1;
	double rf = 1.92; // rf is 3-month T-bill rate, i.e. current date risk-free rate today
	VaR VaR(filename, quantity,rf);
	VaR.historical();
	VaR.MC_GBM();
	VaR.analytical();
	system("pause");
}

