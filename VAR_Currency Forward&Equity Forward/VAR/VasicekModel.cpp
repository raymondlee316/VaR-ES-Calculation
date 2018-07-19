#include "VasicekModel.h"
#include "BSModel.h"
#include <cmath>
#include <algorithm>
#include <iostream>

const double pi = 4.0*atan(1.0);
double Gauss() // used for generating normal random number
{
	double U1 = (rand() + 1.0) / (RAND_MAX + 1.0);
	double U2 = (rand() + 1.0) / (RAND_MAX + 1.0);
	return sqrt(-2.0*log(U1)) * cos(2.0*pi*U2);
}

void VasicekModel::MC_Simulation(double deltaT, long m, vector<double>& r) // simulate future interest rate by Vasicek Model
{
	r[0] = r0;
	for (int k = 1; k<m; k++)
	{
		r[k] = r_bar*(1-exp(-a*deltaT)) + r0 * exp(-a*deltaT) + sigma * sqrt((1-exp(-2*a*deltaT))/(2*a)) *Gauss();
	}
}
