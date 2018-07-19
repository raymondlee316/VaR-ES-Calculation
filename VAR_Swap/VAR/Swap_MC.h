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
#include "CSVparser.hpp"

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
This file is used for simulating future swap NPV using Hull-White One Factor model, and then calculating the PnL of swap NPV.
We first need to calibrate the parameters by current swaption volatility and discount factors.
Then simulate forward zero-coupon bond price, calculate the PnL of swap NPV.
This file can also be used for swap pricing by Monte Carlo.

Note: current swaption volatility, date, discount factors need to be input manully.
*/


vector <double> DiscountFactorVec(const string &filename) {
	Parser data = Parser(filename);
	vector <double> result;
	result.push_back(1.0);
	for (int i = 0; i < data.rowCount(); i++)
		result.push_back(stod(data[i]["Discount"])); // read to 10Yr
	return result;
}

vector <double> ImpliedVolatilityVec(const string &filename) {
	Parser data = Parser(filename);
	vector <double> result;
	for (int i = 4; i <= 13; i++) {
		for (int j = 1; j <= 10; j++)
			result.push_back(stod(data[i][j]) / 100);
	}
	return result;
}

vector<double> Swap_MC(double r0, double fixedRate, int Maturity, int Tenor, double quantity, double SpotPrice) {
	//Number of swaptions to be calibrated
	Size numRows = 10; // Maturity
	Size numCols = 10; // Tenor
	Integer swapLengths[] = { 1,2,3,4,5,6,7,8,9,10 }; // tenor

	Date todaysDate(01, April, 2011); // current date
	Calendar calendar = TARGET();
	Settings::instance().evaluationDate() = todaysDate;

	// start building current yield curve by discount factors
	vector < DiscountFactor > dfs = DiscountFactorVec("DF_20110401.csv"); // Discount Factor
	vector < Volatility > swaptionVols = ImpliedVolatilityVec("SwaptionIV_20110401.csv");
	std::vector <Date > dates; // Dates of each discount factor point
	Calendar cal = TARGET();
	EURLibor3M libor; // only for day count
	DayCounter dc = libor.dayCounter();
	Natural settlementDays = 2;
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
	dates.push_back(settlement + 50 * Years);

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
	boost::shared_ptr<IborIndex> indexSixMonths(new
		USDLibor(Period(3, Months), rhTermStructure));

	// defining the swaptions to be used in model calibration
	std::vector<Period> swaptionMaturities;
	for (int i = 1; i <= 10; i++) {
		swaptionMaturities.push_back(Period(i, Years)); // swaption maturity
	}
	std::vector<boost::shared_ptr<CalibrationHelper> > swaptions; // used for calibration

	Size i;
	for (i = 0; i<numRows; i++) {
		Size j = numCols - i - 1;
		Size k = i * numCols + j;
		boost::shared_ptr<Quote> vol(new SimpleQuote(swaptionVols[k]));
		swaptions.push_back(boost::shared_ptr<CalibrationHelper>(new
			SwaptionHelper(swaptionMaturities[i], // 1,2,3,4,5
				Period(swapLengths[j], Years), // 5,4,3,2,1
				Handle<Quote>(vol), //const Handle<Quote>& volatility,
				indexSixMonths, // boost::shared_ptr<IborIndex>& index
				indexSixMonths->tenor(), //Period& fixedLegTenor,
				indexSixMonths->dayCounter(), //const DayCounter& fixedLegDayCounter,
				indexSixMonths->dayCounter(),
				rhTermStructure))); // const Handle<YieldTermStructure>& termStructure,
	}  // calibrate parts of volatility
	
	boost::shared_ptr<HullWhite> modelHW(new HullWhite(rhTermStructure));// defining the models
	std::cout << "Hull-White (analytic formulae) calibration" << std::endl;
	for (i = 0; i < swaptions.size(); i++) {
		swaptions[i]->setPricingEngine(boost::shared_ptr<PricingEngine>(
			new JamshidianSwaptionEngine(modelHW)));
	}

	LevenbergMarquardt om;
	modelHW->calibrate(swaptions, om,
		EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8));

	// print parameters that were calibrated
	std::cout << "calibrated to:\n"
		<< "a = " << modelHW->params()[0] << ", "
		<< "sigma = " << modelHW->params()[1]
		<< std::endl << std::endl;

	//////////////// Start swap pricing ////////////////////////
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
	Time dt = 1/365.0, t = 0;
	Real dw;
	Size numVals = 1e4; // number of simulation
	Real swapNPV = 0, swapSum = 0;

	vector<double> swapPnL;

	for (Size j = 0; j < numVals; ++j) {
		dw = bmGauss.next().value; // brownian motion
		x = shortRateProces->evolve(t, x0, dt, dw); // simulated short rate
		Real temp = 0;
		for (Size i = (Maturity + 1); i <= Length; i++) {
			temp += modelHW->discountBond(dt, i, x);
		}
		swapNPV = - modelHW->discountBond(dt, Maturity, x) + modelHW->discountBond(dt, Length, x) + fixedRate * temp;
		swapSum += swapNPV;
		swapPnL.push_back((swapNPV-SpotPrice)*quantity); // vector to store simulated swap PnL
	}
	cout << "Swap NPV by Monte Carlo =" << quantity * swapSum / numVals << endl; // Print Swap NPV by Monte Carlo simulation
	return swapPnL;

}

