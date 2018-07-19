#ifndef VasicekModel_h
#define VasicekModel_h
using namespace std;
#include <vector>
#include <cstdlib> // for typedef
#include <ctime>		// for srand()


class VasicekModel
{
public:
	double r0;
	double a, r_bar, sigma;
	VasicekModel(double r0_, double a_, double r_bar_, double sigma_)
	{
		r0 = r0_; a = a_; r_bar = r_bar_;sigma = sigma_;
		srand(time(NULL));
	}
	void MC_Simulation(double deltaT, long m, vector<double>& r);
};
#endif

#pragma once
