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
	const string filename = "StrikeApple145_2015.csv"; 
	long quantity = 1e2; // quantity of position
	double rf = 0.22; // rf is 3-month T-bill rate of current date (here we use 2015.12.31)
	VaR VaR(filename, quantity, rf);
	VaR.historical();
	VaR.MC_GBM();
	VaR.analytical();
	VaR.analytical_Greeks();

	system("pause");
	return 1;
}

