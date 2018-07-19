#pragma warning( disable : 4819)
#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/quantlib.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/treeswaptionengine.hpp>
#include <ql/pricingengines/swaption/jamshidianswaptionengine.hpp>
#include <ql/pricingengines/swaption/g2swaptionengine.hpp>
#include <ql/pricingengines/swaption/fdhullwhiteswaptionengine.hpp>
#include <ql/pricingengines/swaption/fdg2swaptionengine.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/models/shortrate/onefactormodels/blackkarasinski.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/utilities/dataformatters.hpp>

#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>

using namespace QuantLib;
using namespace std;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

	Integer sessionId() { return 0; }

}
#endif

/*
This file is used for simulating forward swap rate using Hull-White One Factor model, and then calculating the PnL of swaption by delta-gamma normal approach.
We first need to calibrate the parameters by current swaption volatility and discount factors.
Then simulate forward swap rate, calculate the PnL of swaption with delta, gamma.
This file can also be used for swaption pricing by Monte Carlo.

Note: current swaption volatility, date, discount factors need to be input manully.
This file is only for ATM swaption because it can only simulate forward swap rate for ATM swaption
*/

//Number of swaptions to be calibrated
Size numRows = 10; // Maturity
Size numCols = 10; // Tenor
Date todaysDate(01, July, 2011); // current date

Integer swapLengths[] = { 1,2,3,4,5,6,7,8,9,10 }; // tenor

Volatility swaptionVols[] = {
	0.719,0.598,0.5,0.428,0.3878,0.353,0.3297,0.3155,0.3048,0.2915,
	0.5442,0.4475,0.3835,0.3453,0.3215,0.3033,0.2888,0.281,0.273,0.2643,
	0.3808,0.3358,0.3075,0.29,0.2783,0.2688,0.2595,0.2577,0.252,0.2448,
	0.2965,0.2788,0.2655,0.2573,0.25,0.2442,0.2385,0.2375,0.2348,0.2288,
	0.258,0.252,0.2455,0.237,0.233,0.2285,0.2242,0.224,0.2215,0.217,
	0.2385,0.234,0.229,0.225,0.222,0.2188,0.2148,0.2135,0.2125,0.2095,
	0.225,0.2227,0.217,0.2125,0.208,0.208,0.204,0.2055,0.205,0.2008,
	0.212,0.211,0.2073,0.2035,0.2018,0.2005,0.198,0.198,0.1977,0.197,
	0.2045,0.2023,0.2002,0.1968,0.1943,0.193,0.1905,0.1902,0.1905,0.1908,
	0.1968,0.193,0.1907,0.1885,0.1868,0.186,0.1848,0.1865,0.1865,0.184 }; // current Swaption Volatility (2011/7/1), row stands for maturity, column for tenor

void calibrateModel(
	/*
	This function is for calibration by optimization method, we try to minimize the difference between model swaption price and market price
	*/
	const boost::shared_ptr<ShortRateModel>& model,
	const std::vector<boost::shared_ptr<CalibrationHelper> >& helpers) {

	LevenbergMarquardt om;
	model->calibrate(helpers, om,
		EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8));
}


vector<double> Swaption_MC(double r0, double strike, double Delta, double Gamma, int Maturity, int Tenor, double quantity, double SpotPrice) {

	Calendar calendar = TARGET();
	Settings::instance().evaluationDate() = todaysDate;

	// start building current yield curve by discount factors
	std::vector <Date > dates; // Dates of each discount factor point
	std::vector < DiscountFactor > dfs; // Discount Factor
	Calendar cal = TARGET();
	EURLibor3M libor; // only for day count
	DayCounter dc = libor.dayCounter();
	Natural settlementDays = 0;
	Date settlement = cal.advance(todaysDate, settlementDays, Days);

	dates.push_back(settlement); dates.push_back(settlement + 3 * Months);
	dates.push_back(settlement + 5 * Months); dates.push_back(settlement + 8 * Months);
	dates.push_back(settlement + 11 * Months); dates.push_back(settlement + 14 * Months);
	dates.push_back(settlement + 17 * Months); dates.push_back(settlement + 20 * Months);
	dates.push_back(settlement + 2 * Years); dates.push_back(settlement + 3 * Years);
	dates.push_back(settlement + 4 * Years); dates.push_back(settlement + 5 * Years);
	dates.push_back(settlement + 6 * Years); dates.push_back(settlement + 7 * Years);
	dates.push_back(settlement + 8 * Years); dates.push_back(settlement + 9 * Years);
	dates.push_back(settlement + 10 * Years); dates.push_back(settlement + 11 * Years);
	dates.push_back(settlement + 12 * Years); dates.push_back(settlement + 15 * Years);
	dates.push_back(settlement + 20 * Years); dates.push_back(settlement + 25 * Years);
	dates.push_back(settlement + 30 * Years); dates.push_back(settlement + 40 * Years);

	dfs.push_back(1.0); dfs.push_back(0.999372);
	dfs.push_back(0.998612); dfs.push_back(0.997533);
	dfs.push_back(0.99626); dfs.push_back(0.994644);
	dfs.push_back(0.992501); dfs.push_back(0.989723);
	dfs.push_back(0.98572); dfs.push_back(0.96545);
	dfs.push_back(0.936598); dfs.push_back(0.90133);
	dfs.push_back(0.862887); dfs.push_back(0.82346);
	dfs.push_back(0.784508); dfs.push_back(0.746472);
	dfs.push_back(0.709697); dfs.push_back(0.674411);
	dfs.push_back(0.640635); dfs.push_back(0.549721);
	dfs.push_back(0.432787); dfs.push_back(0.343264);
	dfs.push_back(0.273871); dfs.push_back(0.180376);

	Handle<YieldTermStructure> rhTermStructure(
		boost::shared_ptr<InterpolatedDiscountCurve<Linear> >(
			new InterpolatedDiscountCurve <Linear>(dates, dfs, dc, cal))); // building yield curve by QuantLib

	// Define the ATM/OTM/ITM swaps for calibration
	Frequency fixedLegFrequency = Annual;
	BusinessDayConvention fixedLegConvention = Unadjusted;
	BusinessDayConvention floatingLegConvention = ModifiedFollowing;
	DayCounter fixedLegDayCounter = Thirty360(Thirty360::European);
	Frequency floatingLegFrequency = Annual;
	VanillaSwap::Type type = VanillaSwap::Payer;
	Rate dummyFixedRate = 0.03;
	boost::shared_ptr<IborIndex> indexSixMonths(new
		USDLibor(Period(3, Months), rhTermStructure));

	// defining the swaptions to be used in model calibration
	std::vector<Period> swaptionMaturities;
	for (int i = 1; i <= 10; i++) {
		swaptionMaturities.push_back(Period(i, Years)); // swaption maturity
	}

	std::vector<boost::shared_ptr<CalibrationHelper> > swaptions; // used for calibration

	// List of times that have to be included in the timegrid
	std::list<Time> times;

	Size i, j, k = 0;
	for (i = 1; i <= numRows; i++) {
		for (j = 1; j <= numCols; j++) {
			boost::shared_ptr<Quote> vol(new SimpleQuote(swaptionVols[k]));
			swaptions.push_back(boost::shared_ptr<CalibrationHelper>(new
				SwaptionHelper(swaptionMaturities[i - 1],  // maturity
					Period(swapLengths[j - 1], Years),  // tenor
					Handle<Quote>(vol), //const Handle<Quote>& volatility,
					indexSixMonths, // boost::shared_ptr<IborIndex>& index
					indexSixMonths->tenor(), //Period& fixedLegTenor,
					indexSixMonths->dayCounter(), //const DayCounter& fixedLegDayCounter,
					indexSixMonths->dayCounter(), // const DayCounter& floatingLegDayCounter,
					rhTermStructure))); // const Handle<YieldTermStructure>& termStructure,
			swaptions.back()->addTimesTo(times);
			k++;
		}
	} // build swaptions used for calibration based on volatility surface

	  // Building time-grid
	TimeGrid grid(times.begin(), times.end(), 30); // used for numerical/tree

	boost::shared_ptr<HullWhite> modelHW(new HullWhite(rhTermStructure));// defining the models

	std::cout << "Hull-White (analytic formulae) calibration" << std::endl;
	for (i = 0; i < swaptions.size(); i++) {
		swaptions[i]->setPricingEngine(boost::shared_ptr<PricingEngine>(
			new JamshidianSwaptionEngine(modelHW)));
	}

	calibrateModel(modelHW, swaptions); // calibrate

	// print parameters that were calibrated
	std::cout << "calibrated to:\n"
		<< "a = " << modelHW->params()[0] << ", "
		<< "sigma = " << modelHW->params()[1]
		<< std::endl << std::endl;

	// Start simulating forward swap rate
	int Length = Maturity + Tenor;

	Real x0 = r0; // current short rate, e.g. current libor overnight rate
	Real x;
	Real a = modelHW->params()[0]; // parameter from HW model
	Real sigma = modelHW->params()[1];
	BigInteger seed = SeedGenerator::instance().get(); // seed for generating normal random number
	MersenneTwisterUniformRng unifMt(seed);
	BoxMullerGaussianRng < MersenneTwisterUniformRng > bmGauss(unifMt);
	boost::shared_ptr < HullWhiteProcess > shortRateProces(
		new HullWhiteProcess(rhTermStructure, a, sigma));
	Time dt = Maturity, t = 0.0;
	Real dw, temp = 0, swapRateChange = 0, swaptionChange = 0;
	Size numVals = 1e4; // number of simulation
	Rate swapRate = 0, swaptionNPV = 0;

	vector<double> swaptionPnL;

	for (Size j = 1; j <= numVals; ++j) {
		dw = bmGauss.next().value; // brownian motion
		x = shortRateProces->evolve(t, x0, dt, dw); // simulated short rate
		for (Size i = (dt + 1); i <= Length; i++) {
			temp += modelHW->discountBond(dt, i, x);
		}
		swapRate = (1.0 - modelHW->discountBond(dt, Length, x)) / temp; // simulated swap rate
		swaptionNPV += rhTermStructure->discount(settlement + Maturity * Years)*max((swapRate - strike), 0.0)*temp; // simulated swaption price
		temp = 0;

		swapRateChange = swapRate - strike; // change of forward swap rate
		swaptionChange = swapRateChange * Delta - 0.5*swapRateChange*swapRateChange*Gamma; // calculate PnL of swaption by delta-gamma normal approach

		swaptionPnL.push_back(swaptionChange*quantity); // vector to store simulated swaption PnL
	}
	cout << "Swaption NPV by Monte Carlo =" << 100 * quantity* swaptionNPV / numVals << endl; // Print Swaption Price by Monte Carlo simulation
	return swaptionPnL;


	//for (int k = 1; k <= numVals; ++k) {
	//	for (Size j = 1; j <= numVals; ++j) {
	//		dw = bmGauss.next().value;
	//		x = shortRateProces->evolve(t, x0, dt, dw);
	//		for (Size i = (dt + 1); i <= Length; i++) {
	//			temp += modelHW->discountBond(dt, i, x); 
	//		}
	//		swapRate = (1.0 - modelHW->discountBond(dt, Length, x)) / temp;
	//		swaptionNPV += rhTermStructure->discount(settlement + Maturity * Years)*max((swapRate - fixedATMRate), 0.0)*temp;
	//		temp = 0;
	//	}
	//	swaptionNPV = 100 * quantity* swaptionNPV / numVals;
	//	swaptionPnL.push_back(swaptionNPV - SpotPrice *quantity);
	//	std::cout << swaptionNPV << endl;
	//	swaptionNPV = 0;
	//}
	//return swaptionPnL;

}

