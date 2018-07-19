#ifndef MonteCarlo_h
#define MonteCarlo_h
#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
# include <ql/auto_link.hpp>
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

///////////// Swaption Data Used for Calibration of Hull-White Model //////////////
//Number of swaptions to be calibrated to...

Size numRows = 10;
Size numCols = 10;

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
	0.1968,0.193,0.1907,0.1885,0.1868,0.186,0.1848,0.1865,0.1865,0.184 }; // current Swaption Volatility (2011/7/1)

//Volatility swaptionVols[] = {
//	0.367,0.34,0.317,0.3,0.2833,0.272,0.261,0.2465,0.24,0.229,
//	0.2855,0.2725,0.2645,0.2548,0.2445,0.2365,0.2295,0.226,0.2215,0.2115,
//	0.254,0.243,0.2385,0.2318,0.227,0.2185,0.217,0.213,0.203,0.2008,
//	0.233,0.2255,0.2195,0.2145,0.2105,0.2063,0.2018,0.2005,0.197,0.191,
//	0.217,0.211,0.207,0.2002,0.1977,0.195,0.194,0.19,0.1875,0.185,
//	0.203,0.195,0.1945,0.1915,0.188,0.1855,0.182,0.183,0.1785,0.175,
//	0.191,0.185,0.183,0.1808,0.182,0.18,0.177,0.174,0.172,0.17,
//	0.1815,0.177,0.175,0.173,0.1705,0.169,0.1675,0.1655,0.164,0.1615,
//	0.172,0.169,0.1665,0.165,0.163,0.161,0.161,0.1585,0.157,0.1555,
//	0.165,0.162,0.164,0.163,0.161,0.1568,0.158,0.1555,0.154,0.1495 }; // current Swaption Volatility (2008/7/1)


void calibrateModel(
	const boost::shared_ptr<ShortRateModel>& model,
	const std::vector<boost::shared_ptr<CalibrationHelper> >& helpers) {

	LevenbergMarquardt om; // Optimization method
	model->calibrate(helpers, om,
		EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8));
	//EndCriteria(Size maxIterations,
	//	Size maxStationaryStateIterations,
	//	Real rootEpsilon,
	//	Real functionEpsilon,
	//	Real gradientNormEpsilon);

	//for (Size i = 0; i<numRows; i++) {
	//	Size j = numCols - i - 1; // 1x10, 2x9, 3x8, 4x7, 5x6, 6x5, 7x4, 8x3, 9x2, 10x1
	//	Size k = i * numCols + j;
	//	Real npv = helpers[i]->modelValue();
	//	// i = 0,1,2,3,4,5,6,7,8,9
	//	// j = 9,8,7,6,5,4,3,2,1,0
	//	// k = 9,18,27,36,45,54,63,72,81,90
	//	Volatility implied = helpers[i]->impliedVolatility(npv, 1e-4,
	//		1000, 0.05, 0.50);
	//	//Volatility impliedVolatility(Real targetValue,
	//	//	Real accuracy,
	//	//	Size maxEvaluations,
	//	//	Volatility minVol,
	//	//	Volatility maxVol) const;
	//	Volatility diff = implied - swaptionVols[k]; // difference between implied vol and actual vol

	//	std::cout << i + 1 << "x" << swapLengths[j]
	//		<< std::setprecision(5) << std::noshowpos
	//		<< ": model " << std::setw(7) << io::volatility(implied)
	//		<< ", market " << std::setw(7)
	//		<< io::volatility(swaptionVols[k])
	//		<< " (" << std::setw(7) << std::showpos
	//		<< io::volatility(diff) << std::noshowpos << ")\n";
	//}
}

vector<double> MonteCarlo_ZeroCouponBond(int numPath, double maturity, double r0, const vector<double> &DFS) {
	boost::timer timer;
	std::cout << std::endl;

	Date todaysDate(01, July, 2011);
	Calendar calendar = TARGET();
	Settings::instance().evaluationDate() = todaysDate;
	
	std::vector <Date > dates; // Dates of each discount factor point
	std::vector < DiscountFactor > dfs = DFS; // Discount Factor
	Calendar cal = TARGET();
	EURLibor1M libor;
	DayCounter dc = libor.dayCounter();
	Natural settlementDays = 0;
	Date settlement = cal.advance(todaysDate, settlementDays, Days);

	dates.push_back(settlement); dates.push_back(settlement + 1 * Months);
	dates.push_back(settlement + 3 * Months); dates.push_back(settlement + 6 * Months);
	dates.push_back(settlement + 1 * Years); dates.push_back(settlement + 2 * Years);
	dates.push_back(settlement + 3 * Years); dates.push_back(settlement + 5 * Years);
	dates.push_back(settlement + 7 * Years); dates.push_back(settlement + 10 * Years);
	dates.push_back(settlement + 20 * Years); dates.push_back(settlement + 30 * Years);

	//dates.push_back(settlement); dates.push_back(settlement + 3 * Months);
	//dates.push_back(settlement + 2 * Years); dates.push_back(settlement + 3 * Years);
	//dates.push_back(settlement + 4 * Years); dates.push_back(settlement + 5 * Years);
	//dates.push_back(settlement + 6 * Years); dates.push_back(settlement + 7 * Years);
	//dates.push_back(settlement + 8 * Years); dates.push_back(settlement + 9 * Years);
	//dates.push_back(settlement + 10 * Years); dates.push_back(settlement + 11 * Years);
	//dates.push_back(settlement + 12 * Years); dates.push_back(settlement + 15 * Years);
	//dates.push_back(settlement + 20 * Years); dates.push_back(settlement + 25 * Years);
	//dfs.clear();
	//dfs.push_back(1.0); dfs.push_back(0.999372);
	//dfs.push_back(0.98572); dfs.push_back(0.96545);
	//dfs.push_back(0.936598); dfs.push_back(0.90133);
	//dfs.push_back(0.862887); dfs.push_back(0.82346);
	//dfs.push_back(0.784508); dfs.push_back(0.746472);
	//dfs.push_back(0.709697); dfs.push_back(0.674411);
	//dfs.push_back(0.640635); dfs.push_back(0.549721);
	//dfs.push_back(0.432787); dfs.push_back(0.343264);


	Handle<YieldTermStructure> rhTermStructure(
		boost::shared_ptr<InterpolatedDiscountCurve<Linear> >(
			new InterpolatedDiscountCurve < Linear >(dates, dfs, dc, cal))); // Build discount factor curve

	boost::shared_ptr<IborIndex> indexSixMonths(new
		Euribor6M(rhTermStructure));

	// defining the swaptions to be used in model calibration
	std::vector<Period> swaptionMaturities;
	for (int i = 1; i <= 10; i++) {
		swaptionMaturities.push_back(Period(i, Years)); // swaption maturity
	}


	std::vector<boost::shared_ptr<CalibrationHelper> > swaptions;

	// List of times that have to be included in the timegrid
	std::list<Time> times;

	Size i;
	for (i = 0; i<numRows; i++) {
		Size j = numCols - i - 1; 
		Size k = i * numCols + j;

		boost::shared_ptr<Quote> vol(new SimpleQuote(swaptionVols[k]));
		swaptions.push_back(boost::shared_ptr<CalibrationHelper>(new
			SwaptionHelper(swaptionMaturities[i], 
				Period(swapLengths[j], Years), 
				Handle<Quote>(vol), //const Handle<Quote>& volatility,
				indexSixMonths, // boost::shared_ptr<IborIndex>& index
				indexSixMonths->tenor(), //Period& fixedLegTenor,
				indexSixMonths->dayCounter(), //const DayCounter& fixedLegDayCounter,
				indexSixMonths->dayCounter(),
				rhTermStructure))); // const Handle<YieldTermStructure>& termStructure,
		swaptions.back()->addTimesTo(times);
	}

	// Building time-grid
	TimeGrid grid(times.begin(), times.end(), 30); // used for numerical/tree
	
												   // defining the models

	boost::shared_ptr<G2> modelG2(new G2(rhTermStructure));
	boost::shared_ptr<HullWhite> modelHW(new HullWhite(rhTermStructure));
	boost::shared_ptr<HullWhite> modelHW2(new HullWhite(rhTermStructure));
	boost::shared_ptr<BlackKarasinski> modelBK(
		new BlackKarasinski(rhTermStructure));




	//std::cout << "Hull-White (analytic formulae) calibration" << std::endl;
	for (i = 0; i<swaptions.size(); i++)
		swaptions[i]->setPricingEngine(boost::shared_ptr<PricingEngine>(
			new JamshidianSwaptionEngine(modelHW)));

	calibrateModel(modelHW, swaptions);
	//std::cout << "calibrated to:\n"
	//	<< "a = " << modelHW->params()[0] << ", "
	//	<< "sigma = " << modelHW->params()[1]
	//	<< std::endl << std::endl;

	Real x0 = r0; // current short rate
	Real x; // simulated short rate
	Real a = modelHW->params()[0];
	Real sigma = modelHW->params()[1];
	BigInteger seed = SeedGenerator::instance().get(); // seed for normal variable simulation
	MersenneTwisterUniformRng unifMt(seed);
	BoxMullerGaussianRng < MersenneTwisterUniformRng > bmGauss(unifMt);
	boost::shared_ptr < HullWhiteProcess > shortRateProces(
		new HullWhiteProcess(rhTermStructure, a, sigma));
	Time dt = 1 / 360.0, t = 0.0; // simulate short rate for one day
	Real dw;
	Size numVals = numPath; // number of valuations/paths
	vector<double> MC_price; // simulated bond price
	for (Size j = 1; j <= numVals; ++j) {
		dw = bmGauss.next().value;
		x = shortRateProces->evolve(t, x0, dt, dw); // simulate short rate
		MC_price.push_back(modelHW->discountBond(dt, maturity, x)); // price bond based on simulated short rate
		//std::cout << std::fixed;
		//std::cout << std::setprecision(6);
		//cout << "Price = " << modelHW->discountBond(0, maturity, x) <<endl;
	}
	return MC_price;
	
}
#endif


