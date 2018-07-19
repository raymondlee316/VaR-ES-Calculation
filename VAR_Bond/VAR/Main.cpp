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
#include "CSVparser.hpp"
#include <typeinfo>
#include "Bonds.h"
#include "VaR.h"
#include "MonteCarlo.h"

using namespace std;

int main()
{
	const string filename1 = "Daily Treasury Yield Curve Rates(2010-2011).csv";  // data set
	const string filename2 = "Discount Factor(2010-2011).csv";
	const string filename3 = "BAC yield.csv";
	const string filename4 = "BAC price.csv";
	const string filename_df = "Discount Factor(2010-2011).csv";
	double alpha = 5; // confidence level, 5 for 95% confidence level
	int frequency = 2; // coupon frequency, 1 for annually, 2 for semi-annually
	int numPath = 1e4; // number of paths for monte carlo simulation
	double Libor_overnight = 0.12550/100; // overnight Libor rate, which is current short rate
	//Libor_overnight = 2.54125/100;
	vector<Bonds*> ZeroCouponBonds1; // special vector to store zero-coupon bonds to create a zero-coupon bond portfolio
	vector<Bonds*> ZeroCouponBonds2; // special vector to store zero-coupon bonds to create a zero-coupon bond portfolio
	vector<Bonds*> IndividualCouponBonds; // special vector to store individual coupon bonds to create a bond portfolio

	//////////example 1: zero-coupon bond calculation based on bond yield instead of bond price//////////////////////
	cout << "example 1: zero-coupon bond calculation based on bond yield instead of bond price." << endl;
	cout << "(one 3-year-maturity $100 par value bonds with 5% annually coupon rates)" << endl << endl;
	double quantity1 = 5; // par value instead of market value
	string maturity1 = "1"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_1y(filename1, quantity1, maturity1, alpha, frequency); // create one object for each zero coupon bond
	bond_1y.setYield(); // use setyield() function to initialize if use data set of yield; use setprice() for data of price
	bond_1y.Calculation_YieldBased(); // use zerocouponcalculation_yield() function to caculate if use data set of yield; use zerocouponcalculation_price if use data of price
	cout << "analytical VaR of " << maturity1 << " year bond based on yield is: " << bond_1y.VaR_Analytical_Yield() << endl; // caculate VaR of individual zero-coupon bond by analytical method
	cout << "historical VaR of " << maturity1 << " year bond based on yield is: " << bond_1y.VaR_Historical_Yield() << endl; // caculate VaR of individual zero-coupon bond by historical method
	cout << "Simulated  VaR of " << maturity1 << " year bond by Monte Carlo is: " << bond_1y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds1.push_back(&bond_1y); // store in vector zerocouponbonds for later caculation
	cout << endl;

	double quantity2 = 5; // par value instead of market value
	string maturity2 = "2"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_2y(filename1, quantity2, maturity2, alpha, frequency);
	bond_2y.setYield();
	bond_2y.Calculation_YieldBased();
	cout << "analytical VaR of " << maturity2 << " year bond based on yield is: " << bond_2y.VaR_Analytical_Yield() << endl;
	cout << "historical VaR of " << maturity2 << " year bond based on yield is: " << bond_2y.VaR_Historical_Yield() << endl;
	cout << "Simulated  VaR of " << maturity2 << " year bond by Monte Carlo is: " << bond_2y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds1.push_back(&bond_2y);
	cout << endl<< endl;

	double quantity3 = 105; // par value instead of market value
	string maturity3 = "3"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_3y(filename1, quantity3, maturity3, alpha, frequency);
	bond_3y.setYield();
	bond_3y.Calculation_YieldBased();
	cout << "analytical VaR of " << maturity3 << " year bond based on yield is: " << bond_3y.VaR_Analytical_Yield() << endl;
	cout << "historical VaR of " << maturity3 << " year bond based on yield is: " << bond_3y.VaR_Historical_Yield() << endl;
	cout << "Simulated  VaR of " << maturity3 << " year bond by Monte Carlo is: " << bond_3y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds1.push_back(&bond_3y);
	cout << endl;

	////////////example 2: zero-coupon bond calculation based on bond price instead of yield (note different function name)//////////////////////
	cout << "example 2: zero-coupon bond calculation based on bond price instead of yield." << endl << endl;
	double quantity4 = 5; // par value instead of market value
	string maturity4 = "5"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_4y(filename2, quantity4, maturity4, alpha, frequency);
	bond_4y.setPrice();
	bond_4y.Calculation_PriceBased();
	cout << "analytical VaR of " << maturity4 << " year bond based on price is: " << bond_4y.VaR_Analytical_Price() << endl;
	cout << "historical VaR of " << maturity4 << " year bond based on price is: " << bond_4y.VaR_Historical_Price() << endl;
	cout << "Simulated  VaR of " << maturity4 << " year bond by Monte Carlo is: " << bond_4y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds2.push_back(&bond_4y);
	cout << endl;

	double quantity5 = 5; // par value instead of market value
	string maturity5 = "7"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_5y(filename2, quantity5, maturity5, alpha, frequency);
	bond_5y.setPrice();
	bond_5y.Calculation_PriceBased();
	cout << "analytical VaR of " << maturity5 << " year bond based on price is: " << bond_5y.VaR_Analytical_Price() << endl;
	cout << "historical VaR of " << maturity5 << " year bond based on price is: " << bond_5y.VaR_Historical_Price() << endl;
	cout << "Simulated  VaR of " << maturity5 << " year bond by Monte Carlo is: " << bond_5y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds2.push_back(&bond_5y);
	cout << endl;
	
	double quantity6 = 105; // par value instead of market value
	string maturity6 = "10"; // maturity in year, eg:"1/12" means 1 month
	ZeroCouponBonds bond_6y(filename2, quantity6, maturity6, alpha, frequency);
	bond_6y.setPrice();
	bond_6y.Calculation_PriceBased();
	cout << "analytical VaR of " << maturity6 << " year bond based on price is: " << bond_6y.VaR_Analytical_Price() << endl;
	cout << "historical VaR of " << maturity6 << " year bond based on price is: " << bond_6y.VaR_Historical_Price() << endl;
	cout << "Simulated  VaR of " << maturity6 << " year bond by Monte Carlo is: " << bond_6y.VaR_MonteCarlo(numPath, filename_df, Libor_overnight) << endl; // caculate VaR of individual zero-coupon bond by historical method
	ZeroCouponBonds2.push_back(&bond_6y);
	cout << endl;

	//////////example 3: bond portolio consists of zero-coupon bonds//////////////////////
	cout << "example 3: bond portolio consists of zero-coupon bonds (based on yield data)" << endl;
	VaR_BondPorfolio_Analytical_Yield(ZeroCouponBonds1); // caculate bond portfolio which consists of zero-coupon bonds by analytical method
	VaR_BondPorfolio_Historical_Yield(ZeroCouponBonds1); // caculate bond portfolio which consists of zero-coupon bonds by historical method
	VaR_BondPorfolio_MonteCarlo(ZeroCouponBonds1);
	cout << endl;
	
	cout << "bond portolio consists of zero-coupon bonds (based on price data)" << endl;
	VaR_BondPorfolio_Analytical_Price(ZeroCouponBonds2);
	VaR_BondPorfolio_Historical_Price(ZeroCouponBonds2);
	VaR_BondPorfolio_MonteCarlo(ZeroCouponBonds2);
	cout << endl << endl;


	////////////Example 4: Coupon Bond Calculation Based on Yield Data//////////////////////
	cout << "Example 4: Coupon Bond Calculation Based on Yield Data" << endl << endl;
	double quantity7 = 100; // par value instead of market value
	double coupon7 = 4.2 /100; // coupon rate
	string maturity7 = "10"; // maturity in year, eg:"1/12" means 1 month
	CouponBonds bond_7y(filename3, quantity7, maturity7, alpha, frequency, coupon7);
	bond_7y.setYield();
	bond_7y.Calculation_YieldBased();
	cout << "Analytical VaR of " << maturity7 << " year bond based on yield is: " << bond_7y.VaR_Analytical_Yield() << endl;
	cout << "Historical VaR of " << maturity7 << " year bond based on yield is: " << bond_7y.VaR_Historical_Yield() << endl;
	IndividualCouponBonds.push_back(&bond_7y);
	cout << endl << endl;

	////////////Example 5: Coupon Bond Calculation Based on Price Data (Notice different function name)//////////////////////
	cout << "Example 5: Coupon Bond Calculation Based on Price Data (Notice different function name)" << endl << endl;
	double quantity8 = 100; // par value instead of market value
	double coupon8 = 4.2 / 100; // coupon rate
	string maturity8 = "10"; // maturity in year, eg:"1/12" means 1 month
	CouponBonds bond_8y(filename4, quantity8, maturity8, alpha, frequency, coupon8);
	bond_8y.setPrice();
	bond_8y.Calculation_PriceBased();
	cout << "Analytical VaR of " << maturity8 << " year bond based on price is: " << bond_8y.VaR_Analytical_Price() << endl;
	cout << "Historical VaR of " << maturity8 << " year bond based on price is: " << bond_8y.VaR_Historical_Price() << endl;
	IndividualCouponBonds.push_back(&bond_8y);
	cout << endl << endl;

	//////////Example 6: Bond Portfolio//////////////////////
	cout << "Example 6: Bond Portolio Consists of Bonds (based on yield data)" << endl;
	VaR_BondPorfolio_Analytical_Yield(IndividualCouponBonds); // Caculate bond portfolio which consists of zero-coupon bonds by Analytical method
	cout << endl;
	VaR_BondPorfolio_Historical_Yield(IndividualCouponBonds); // Caculate bond portfolio which consists of zero-coupon bonds by Historical method
	cout << endl;

	cout << "Bond Portolio Consists of Bonds (based on price data)" << endl;
	VaR_BondPorfolio_Analytical_Price(IndividualCouponBonds);
	cout << endl;
	VaR_BondPorfolio_Historical_Price(IndividualCouponBonds);


	system("pause");
	return 0;
}

