////////////////////////////////////////////////////////////////////////////////////////
/// \file weathergen.cpp
/// \brief Global Weather GENerator 
///
/// \author Lars Nieradzik
/// $Date: 2017-11-24 15:04:09 +0200 (Fri, 24 Nov 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
// 
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.GLOBFIRM

/* This code has been translated from the original FORTRAN-90 code provided with the article
 * describing GWGEN. Therefore, the commenting has been taken from the original with some additions.
 * For description and details please refer to:
 *   Sommer, P. S. and Kaplan, J. O.: A globally calibrated scheme for generating daily meteorology
 *     from monthly statistics: Global-WGEN (GWGEN) v1.0, Geosci. Model Dev., 10, 3771-3791,
 *     doi:10.5194/gmd-10-3771-2017, 2017.
 *   Original Code in Fortran available at: https://arve-research.github.io/gwgen/
 */

#include "config.h"
#include "weathergen.h"
#include <limits>

const double DHUGE     = std::numeric_limits<double>::max();
const long LHUGE       = std::numeric_limits<long>::max()  ;
const int IHUGE	       = std::numeric_limits<int>::max()   ;
const double D_EPSILON = std::numeric_limits<double>::min();
const float R_EPSILON  = std::numeric_limits<float>::min() ;

// Parameters used witin GWGEN
// Freezing temperature of freshwater (K)
const double TFREEZE = 273.15;      

// -----------------------------------------------------------------------------
// ------------------- Defaults for the namelist parameters --------------------
// -----------------------------------------------------------------------------

// This parameter can be used to tune between fast and fine: 
// min maxcount=20; max maxcount = 10000000
const int MAXITER = 20;

//Tmin break-off criterion
const double E_TMIN = 2.5;

// Dampening factors for Tmax and Tmin
// (use these to increase/reduce amplitude of daily Tmin/Tmax variability)
const double TMIN_DAMP = 0.25;
const double TMAX_DAMP = 0.8;

//whether to reset residuals at beginning of each year
const bool LRESET = false;

// Correction for circular conditions when using climatology
// (do not use for transient simulations)
const bool CIRCULAR_CLIM = false;

const int QSIZ  = 10 ;
const int CMUL  = 69609;
const int COFFS =   123;

const double  RNG1 = 1. / (2. * (double)IHUGE);  //scales the random integer to -0.5,0.5
const double  RNG2 = 1. / (double)IHUGE;	 //scales the random integer to -1,1

const double  HALF = 0.5;

class MetVariables{

public:
	// MEMBER VARIABLES

	int month;

	// Derived datatype for the monthly weather generator input

	double mprec; // monthly total precipitation amount (mm)
	double mwetd; // number of days in month with precipitation
	double mwetf; // fraction of days in month with precipitation

	double mtmin; // minumum temperture (C)
	double mtmax; // maximum temperture (C)
	double mcldf; // cloud fraction (0=clear sky, 1=overcast) (fraction)
	double mwind; // wind speed (m/s)
	
	
	bool   pday[2];  // precipitation status: true if the day was a rain day
	double resid[4]; //previous day's weather residuals

	// Derived datatype for the daily weather generator output

	double dprec; // 24 hour total precipitation (mm)
	double dtmin; // 24 hour mean minimum temperature (degC)
	double dtmax; // 24 hour mean maximum temperature (degC)
	double dcldf; // 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
	double dwind; // wind speed (m s-1)
	double drhum; // relative humidity (%)
	double unorm[4];

	double tmn;
	double tmx;
	double wnd;
	double cld;

	// Derived datatype (F90: type_daymetvars) for monthly climate variables
	double dmtmax_mn; // maximum temperature monthly mean (degC)
	double dmtmin_mn; // minimum temperature mothly mean (degC)
	double dmcldf_mn; // mean cloud fraction (fraction)
	double dmwind_mn; // wind speed
	// Standard deviations of corresponding variable above
	double dmtmax_sd; 
	double dmtmin_sd; 
	double dmcldf_sd; 
	double dmwind_sd; 

	// the following parameters are computed by the cloud_params subroutine
	double cldf_w1, cldf_w2, cldf_w3, cldf_w4, cldf_d1, cldf_d2, cldf_d3, cldf_d4;
	double cldf_sd_w, cldf_sd_d;

} metvars;

// Threshold for transition from gamma to gp distribution
double thresh = 5.0; 

// Interpret the thresh as percentile
bool thresh_pctl = false; 

// Coefficient to estimate the gamma scale parameter via
// g_scale = g_scale_coeff * mean_monthly_precip / number_of_wet_days
// following Geng et al., 1986.
// coefficient to esimate the gamma scale parameter
double g_scale_coeff = 1.268022;
 
// Shape parameter for the Generalized Pareto distribution
double gp_shape = 1.5; 

// A matrix used for cross correlation following Richardson_1984 equation (4)
double  A[4][4] = { 
	{0.913437,  0.032532,  -0.020658, 0.000573},
	{0.488761,  0.137304,  -0.07251, -0.046058},
	{-0.00199, -0.04573,    0.591761, 0.026439}, 
	{0.010905, -0.044171,  -0.018568, 0.666672}};
// B matrix used for cross correlation following Richardson_1984 equation (4)
double B[4][4] = { 
	{0.361854, 0.0,      0.0,      0.0}, 
	{0.11441,  0.802987, 0.0,      0.0}, 
	{0.144862,-0.060622, 0.782791, 0.0}, 
	{0.080593,-0.015829, 0.066186, 0.736713}};

// Transition probability correlations
double p11_1  = 0.254877   ; // intercept of p11 best line fit
double p11_2  = 0.745123   ; // slope of p11 best line fit
double p101_1 = 0.0       ; // intercept of p101 best line fit
double p101_2 = 0.846326  ; // slope of p101 best line fit
double p001_1 = 0.0       ; // intercept of p001 best line fit
double p001_2 = 0.724019  ; // slope of p001 best line fit

/* Temperature and cloud correlation parameters corresponding to wet or dry day
 * minimum temperature regression results.
 */
double tmin_w1 = 1.164653; // intercept of best line fit of tmin on wet days (see `meansd`)
double tmin_w2 = 0.955787; // slope of best line fit of tmin on wet days (see `meansd`)
double tmin_d1 =-0.528308; // intercept of best line fit of tmin on dry days (see `meansd`)
double tmin_d2 = 1.020964; // slope of best line fit of tmin on dry days (see `meansd`)
double tmin_sd_breaks[3] = { -40., 0.0, 25. };  // breaks of tmin sd correlation

// Polynomial coefficients for correlating tmin sd on wet days
double tmin_sd_w[6][4]  = { 
	// < -40       -40 - 0     0 - 25	> 25
	{9.72715668, 3.05498827, 3.21874237, 0.55707042}, 
	{0.1010504, -0.21158825,-0.04507634, 0.02443123}, 
	{0.0,        0.01374948, 0.02094482, 0.0}, 
	{0.0,        0.00140538,-0.00264577, 0.0}, 
	{0.0,        3.686e-05,   9.818e-05, 0.0},
	{0.0,        3.2e-07,     -1.13e-06, 0.0}} ;

// Polynomial coefficients for correlating tmin sd on dry days
double tmin_sd_d[6][4] = {
	//  < -40       -40 - 0     0 - 25        > 25
	{10.89900605, 3.56755661,  3.79411755, -4.61943457}, 
	{0.12709893, -0.11544588,  0.03298697,  0.22605603}, 
	{0.0,         0.02824401, -0.01504554,  0.0},
	{0.0,         0.00195612,  0.00190346,  0.0},
	{0.0,         4.314e-05,  -0.00011362,  0.0},
	{0.0,         3.2e-07,     2.13e-06,    0.0} };

// Maximum temperature regression results
double tmax_w1 = -0.586296 ; // intercept of best line fit of tmax on wet days (see `meansd`)
double tmax_w2 = 0.948669  ; // slope of best line fit of tmax on wet days (see `meansd`)
double tmax_d1 = 0.386508  ; // intercept of best line fit of tmax on dry days (see `meansd`)
double tmax_d2 = 1.0061    ; // slope of best line fit of tmax on dry days (see `meansd`)
double tmax_sd_breaks[3] = { -30., 0.0, 35. };  // polynomial coefficients for the breaks of tmax sd correlation

// Polynomial coefficients for correlating tmax sd on wet days
double tmax_sd_w[6][4] = {
	//   < -30       -30 - 0     0 - 35        > 35
	{6.67200351,  3.86010858,  3.79193207,  5.55292835}, 
	{0.03643908, -0.21861197, -0.03126021, -0.09734715}, 
	{0.0,         0.00388465,  0.01611473,  0.0},
	{0.0,         0.00146174, -0.00120298,  0.0},
	{0.0,         6.059e-05,   2.912e-05,   0.0},
	{0.0,         7.4e-07,    -2.4e-07,     0.0}}; 

// Polynomial coefficients for correlating tmax sd on dry days
double tmax_sd_d[6][4] = {
	//   < -30       -30 - 0     0 - 35	> 35
	{7.37455165,  4.61701866,  4.74550991,  3.25541815}, 
	{0.01535526, -0.33872824, -0.07609816, -0.02178605}, 
	{0.0,        -0.0187566,   0.01893058,  0.0}, 
	{0.0,        -0.0003185,  -0.00134943,  0.0}, 
	{0.0,         3.5e-06,     3.209e-05,   0.0}, 
	{0.0,         1.1e-07,    -2.5e-07,     0.0}};  

// Cloud regression results
double cldf_w    = -0.738271; // parameter for cloud fit on wet days (see `meansd`)
double cldf_d    =  0.420534; // parameter for cloud fit on dry days (see `meansd`)
double cldf_sd_w = 0.981917;  // parameter for std. dev. of cloud fit on wet days (see `meansd`)
double cldf_sd_d = 1.041732;  // parameter for std. dev. of cloud fit on dry days (see `meansd`)

// Wind regression results
double wind_w1 = 0.0      ;   // intercept of best line fit of wind on wet days (see `meansd`)
double wind_w2 = 1.092938 ;   // slope of best line fit of wind on wet days (see `meansd`)
double wind_d1 = 0.0      ;   // intercept of best line fit of wind on dry days (see `meansd`)
double wind_d2 = 0.945229 ;   // slope of best line fit of wind on wet days (see `meansd`)

// Polygon coefficients for wind standard deviation on wet days
double wind_sd_w[6] = { 0.0, 0.81840997, -0.12633931, 0.00933591, 0.0, 0.0};

// Polygon coefficients for wind standard deviation on dry days
double wind_sd_d[6] = { 0.0, 1.08596114, -0.24073323, 0.02216454, 0.0, 0.0};

// Wind bias correction (Note: Default is no correction)
// min and max range for bias correction (1st and 99th percentile)
double wind_bias_min =-2.3263478740;
double wind_bias_max = 2.3263478740; 

// Parameters for the exponential intercept correction
double wind_intercept_bias_a = 1.1582245720322826;  // slope in the exponent
double wind_intercept_bias_b =-1.3358916953022832;  // intercept in the exponent

// Parameters of the slope - unorm best fit line
// Coefficients for the bias correction of wind speed
double wind_bias_coeffs[6] = {0.995353879899162,   0.8507947091050573, 0.027799823700343333, 
			      -0.06710144300871658, 0.0, 0.0};
double wind_intercept_bias_coeffs[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

// Alternative slope bias correction using a logistic function
double wind_slope_bias_L  = -9999.; // maximum value of logistic function of wind bias correction
double wind_slope_bias_k  = -9999.; // steepness of logistic function of wind bias correction
double wind_slope_bias_x0 = -9999.; // x-value of sigmoid's midpoint of logistic function of wind bias co

// Coefficients for the bias correction of minimum temperature
// (Note: Default is no correction)
double tmin_bias_coeffs[6] = {0., 0., 0., 0., 0., 0.}; // coefficients for the bias correction of minimum temperature

// min. and max range for bias correction (1st and 99th percentile)
double tmin_bias_min =-2.3263478740;
double tmin_bias_max = 2.3263478740; 

// Matrix multiplication
void matrixmult(double AA[4][4], double B[4], double CC[4]) {
	// 
	int xy = 4;
	for(int j=0;j<xy;j++){ 
		double res = 0.;
		for(int i=0;i<xy;i++){ 
			res += AA[i][j] * B[i]; 
		} 
		CC[j] = res; 
	}
}

// Generate seed depending on geolocation
unsigned int geohash(double lat, double lon) {
	const double SCALE  = 120.0; // scale factor for geohash, the larger the number the more unique values
	const double OFFSET = 0.5;   // offset to calculate pixel number assuming gridcell center coordinates
	const unsigned int ROWLEN = (int)roundoff(SCALE * 360.,0);
	
	unsigned int i,j;
	
	i = (int)roundoff(OFFSET + SCALE * (lon + 180.),0);
	j = (int)roundoff(OFFSET + SCALE * (lat +  90.),0);
	unsigned int geohash = i + ROWLEN * ( j-1 ) - IHUGE;
	return geohash;
}

// Set the seed for the random distribution
void get_seed_by_location(double lat, double lon, WeatherGenState& state) {

	state.xs = geohash(lat,lon);

	for (int i=0; i<QSIZ; i++) { 
		
	 	state.xcng = state.xcng * CMUL + COFFS;
	
	 	state.xs = state.xs^(state.xs<<13); 
	 	state.xs = state.xs^(state.xs>>17); 
	 	state.xs = state.xs^(state.xs>>5); 
	 	state.q[i] = state.xcng + state.xs;
	 }
}

int refill(WeatherGenState& state) {
	
	// Reset random state

	int s;
	int z;
	int h;
	for (int i=0; i<QSIZ; i++) {
		h = state.carry & 1;
		z = ((((unsigned int)state.q[i])<<9)>>1) + ((((unsigned int)state.q[i])<<7)>>1) + ((unsigned int)state.carry>>1);
		state.carry = ((unsigned int)state.q[i]>>23) + (((unsigned int)state.q[i])>>25) + ((unsigned int)z>>31);

		state.q[i] = ~(((unsigned int)z<<1)+h); 					
	}

	state.indx = 1;
	s = state.q[0];
	return s;
}

/* Generates a uniformly distributed random 4 byte integer with the range (-huge(i4),+huge(i4))
 * based on the 32-bit super KISS random number generator by George Marsaglia, published online
 * and translated to Fortran 90 by user "mecej4" and Marsaglia, 
 * http://forums.silverfrost.com/viewtopic.php?t=1480
 * Further modifications to pass the complete state of the generator as an argument by J.O. Kaplan, 2011
 * 
 * Do not use this function outside of this weathergen.cpp module. 
 * Instead use the standard random number generator of LPJ-GUESS.
 */
int ranu(WeatherGenState& state) {

	int supr, ranu;
	if (state.indx < QSIZ) {
		supr = state.q[state.indx];
		state.indx += 1;
	}
	else {
		supr = refill(state); //reset the generator
	}
	state.xcng = state.xcng * CMUL + COFFS;

	state.xs = state.xs^(state.xs<<13); 
	state.xs = state.xs^(state.xs>>17); 
	state.xs = state.xs^(state.xs>> 5); 

	ranu = state.xcng + state.xs + supr;
	
	return ranu;
}

// Generate a random number in the range (0,1)
/* Do not use this function outside of this weathergen.cpp module. 
 * Instead use the standard random number generator of LPJ-GUESS.
 */
double ranur(WeatherGenState& state) {

	double ranur;
	
	ranur = (double)ranu(state) * RNG1 + HALF;
	return ranur ;
}

/* Calculate the parameters used for the first approximation in "meansd"
 * This subroutine calculates the necessary parameters for the adjustment of
 * the monthly cloud fraction mean depending on the wet/dry state
 */
void calc_cloud_params(MetVariables& metvars) {
	

	metvars.cldf_w1   = -cldf_w - 1.0;
	metvars.cldf_w2   = cldf_w * cldf_w;
	metvars.cldf_w3   = -(cldf_w * cldf_w) - cldf_w;
	metvars.cldf_w4   = - 1.0/cldf_w;
	metvars.cldf_sd_w = cldf_sd_w * cldf_sd_w;

	metvars.cldf_d1   = -cldf_d - 1.0;
	metvars.cldf_d2   = cldf_d * cldf_d;
	metvars.cldf_d3   = -(cldf_d * cldf_d) - cldf_d;
	metvars.cldf_d4   = - 1.0/cldf_d;
	metvars.cldf_sd_d = cldf_sd_d * cldf_sd_d;

} 

void temp_sd(MetVariables& metvars) {

	double dmtmin_sd = 0.;
	double dmtmax_sd = 0.;

	bool rainday = metvars.pday[0];

	for (int i=0; i<4; i++) {
		if (i == 3 || metvars.dmtmin_mn <= tmin_sd_breaks[i]) {
			for (int x=0; x<6; x++) {
				if (rainday) {
					dmtmin_sd += tmin_sd_w[x][i] * pow(metvars.dmtmin_mn, (double)x);
				}
				else{
					dmtmin_sd += tmin_sd_d[x][i] * pow(metvars.dmtmin_mn, (double)x);
				}
			}
			break;
		}
	}

	for (int i=0; i<4; i++) {
		if (i == 3 || metvars.dmtmax_mn <= tmax_sd_breaks[i]) {
			for (int x=0; x<6; x++) {
				if (rainday) {
					dmtmax_sd += tmax_sd_w[x][i] * pow(metvars.dmtmax_mn,(double)x);
				}
				else{
					dmtmax_sd += tmax_sd_d[x][i] * pow(metvars.dmtmax_mn,(double)x);
				}
			}
			break;
		}
	}

	metvars.dmtmin_sd = dmtmin_sd;
	metvars.dmtmax_sd = dmtmax_sd;	
} 

/* Adjust the monthly means of temperature, cloud and wind corresponding to the wet/dry state
 *
 * This routine makes the first approximation inside the weather generator to adjust the monthly
 * mean according to the wet/dry state using the best fit lines from the parameterization.
 *
 * Min. and max. temperature, as well as the wind speed, are calculated via
 *
 * .. math::
 *
 *     x_{w/d} = x_{w/d1} + x_{w/d2} \cdot \bar{x}
 *
 * Where :math:`x` stands either for the :math:`T_{min}, T_{max}, T_{min, sd}, T_{max, sd}, wind`
 * or :math:`wind_{sd}`. :math:`w/d` stands for the wet dry state deterimined by `pday`.
 *
 * The cloud fraction is calculated via
 *
 * .. math::
 *
 *     c_{w/d} = \frac{-a_{w/d} - 1}{a_{w/d}^2 * \bar{c} - a_{w/d}^2 - a_{w/d}}  - \frac{1}{a_{w/d}}
 *
 * and it's standard deviation via
 *
 * .. math::
 *
 *     c_{sd, w/d} = a_{sd, w/d}^2 \cdot c_{w/d} \cdot (1 - c_{w/d})
 *
 * Input
 *  pday      : precipitation status (mm/day)
 *  tmn       : smooth interpolation of monthly minimum temperature (degC)
 *  tmx       : smooth interpolation of monthly maximum temperature (degC)
 *  cld       : fraction (0-1)
 *  wind      : wind speed (m/s)
 * Output
 *  metvars.dm* : first guess for the first daily approximation
 */
void meansd(MetVariables& metvars) {
	
	// calculate mean and SD for a wet day
	if (metvars.pday[0]) {  

		metvars.dmtmin_mn = tmin_w1 + tmin_w2 * metvars.tmn;
		metvars.dmtmax_mn = tmax_w1 + tmax_w2 * metvars.tmx;
		metvars.dmwind_mn = wind_w1 + wind_w2 * metvars.wnd;
		metvars.dmcldf_mn = metvars.cldf_w1 / (metvars.cldf_w2 * metvars.cld + metvars.cldf_w3)
			+ metvars.cldf_w4;
		metvars.dmwind_sd = 0.;
		for (int i=0; i<6; i++) {
			metvars.dmwind_sd += wind_sd_w[i] * pow(metvars.dmwind_mn,(double)i);
		}
		metvars.dmcldf_sd = metvars.cldf_sd_w * metvars.dmcldf_mn * (1. - metvars.dmcldf_mn);
	}
	// calculate mean and SD for a dry day
	else {

		metvars.dmtmin_mn = tmin_d1 + tmin_d2 * metvars.tmn;
		metvars.dmtmax_mn = tmax_d1 + tmax_d2 * metvars.tmx;
		metvars.dmwind_mn = wind_d1 + wind_d2 * metvars.wnd;
		metvars.dmcldf_mn = metvars.cldf_d1 / (metvars.cldf_d2 * metvars.cld + metvars.cldf_d3)
			+ metvars.cldf_d4;
		metvars.dmwind_sd = 0.;
		for (int i=0; i<6; i++) {
			metvars.dmwind_sd += wind_sd_d[i] * pow(metvars.dmwind_mn, (double)i);
		}
		metvars.dmcldf_sd = metvars.cldf_sd_d * metvars.dmcldf_mn * (1. - metvars.dmcldf_mn);
	}
	temp_sd(metvars);

} 

/* Sampler for the normal distribution centered at 0 with std. dev. of unity,
 * based on Marsaglia polar method
 * state : state of the uniform random number generator
 * nval  : output: The random number from the normal distribution
 */
double ran_normal(WeatherGenState& state) {
	
	double vals[2], v[2];

	int u[2];

	double s,a,nval;

	vals[0] = state.gamma_vals[0];
	vals[1] = state.gamma_vals[1];
	
	if (state.have) {

		state.have = false;
		nval = vals[1];
	}
	else {
		int b = 0;
		do {

			u[0] = ranu(state);
			u[1] = ranu(state);

			v[0] = (double)u[0] * RNG2;   //convert integer (-huge,+huge) to (-1,+1)
			v[1] = (double)u[1] * RNG2;   //convert integer (-huge,+huge) to (-1,+1)
			s = pow(v[0],2) + pow(v[1],2) ;
				
			if (s < 1.)
				b=1;
		} while (b==0);
 
		a = sqrt(-2. * log(s) / s);

		vals[0] = v[0] * a;
		vals[1] = v[1] * a;
		state.gamma_vals[0] = vals[0];
		state.gamma_vals[1] = vals[1];

		nval = vals[0];

		state.have = true;
	}
 	return nval;
}

/* Select a random number from a Gamma distribution
 *
 * adapted from the cpp adaptation of the Marsaglia & Tsang random gamma algorithm in:
 * http://www.johndcook.com/SimpleRNG.cpp
 *
 * Uses the algorithm in
 * Marsaglia, G. and Tsang, W.W. (2000), *A simple method for generating
 * gamma variables*, Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.
 *
 * state  : state of the uniform random number generator
 * first  : flag if this is the first call to the distribution with this shape
 * shape  : shape parameter of the Gamma distribution (k or alpha, unitless) > 0
 * scale  : scale parameter of the Gamma distribution (theta = 1/beta)
 */
double ran_gamma(WeatherGenState& state,bool first, double shape, double scale) {
	
	static double c;
	static double d;
	double u;
	double v;
	double x;
	double ret; 
	if (shape <= 0.) {
		printf("shape parameter value must be positive \n");
		exit(-1);
	}
	else if (shape >= 1.) {
		if (first) {
			d = shape - 1. / 3.;
			c = 1. / sqrt(9. * d);
		}
		while ( 1==1 ) {
			while  ( 1==1 ) { //generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.
				x=ran_normal(state);
				v = pow((1. + c * x),3);
				if (v > 0.)
					break;
			}
			u = ranur(state);  // generate uniform variable u in the range (0,1)
			if (u < 1. - 0.0331 * pow(x,4) || log(u) < HALF * pow(x,2) + d*(1. - v + log(v))) {
				ret = scale * d * v  ;
				break;
			}
		}
	}
	else {
		// Order is important w.r.t.  state -> reversed as in FORTRAN!!!
		ret = scale * pow(ranur(state),(1. / shape)) * ran_gamma(state,first,shape + 1.,1.) ;
	}
	return ret;
}

/* Select a random number from a generalized pareto (GP) distribution
 *
 * state  : state of the uniform random number generator
 * shape  : shape parameter of the GP distribution (k or alpha, unitless) > 0
 * scale  : scale parameter of the GP distribution (theta = 1/beta)
 * loc    : the location of the GP distribution
 */
double ran_gp(WeatherGenState& state,double shape,double scale, double loc) {
	
	double u, rangp;

	u = ranur(state); // generate uniform variable u in the range (0,1)

	if (shape == 0.0) { 
		rangp = loc - scale * log(u);
	}
	else {
		rangp = loc + scale * (pow(u,(-shape)) - 1) / shape;
	}
	return rangp;
}

/* Select a random number from a hybrid Gamma-GP distribution
 * Variables
 * state     : state of the uniform random number generator
 * first     : flag if this is the first call to the distribution with this shape
 * shape     : shape parameter of the Gamma distribution (k or alpha, unitless) > 0
 * scale     : scale parameter of the Gamma distribution (theta = 1/beta)
 * thresh    : the threshold above which to choose the GP distribution
 * shape_gp  : shape parameter of the GP distribution
 * scale_gp  : scale parameter of the GP distribution
 */
double ran_gamma_gp(WeatherGenState& state,bool first,double shape,double scale,double thresh,double shape_gp,double scale_gp) {
	
	double ret;
	
	ret = ran_gamma(state,first,shape,scale);
	if (ret > thresh) {
		ret = ran_gp(state,shape_gp,scale_gp,thresh);
	}
	return ret;
}

/* Calculate the natural logarithm of GAMMA ( X ).
 *
 * Computation is based on an algorithm outlined in references 1 and 2.
 * The program uses rational functions that theoretically approximate
 * :math:`\log(\Gamma(X))` to at least 18 significant decimal digits.  The
 * approximation for 12 < X is from Hart et al, while approximations
 * for X < 12.0E+00 are similar to those in Cody and Hillstrom,
 * but are unpublished.
 *
 * The accuracy achieved depends on the arithmetic system, the compiler,
 * intrinsic functions, and proper selection of the machine dependent
 * constants.
 *
 *  Licensing:
 *    This code is distributed under the GNU LGPL license.
 *  Modified:
 *    - 16 June 1999
 *    - Extracted June, 2016
 *  Author:
 *    - Original FORTRAN77 version by William Cody, Laura Stoltz.
 *    - FORTRAN90 version by John Burkardt.
 *    - Extracted by Philipp Sommer
 *  Reference:
 *      - William Cody, Kenneth Hillstrom,
 *        Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
 *        Mathematics of Computation,
 *        Volume 21, 1967, pages 198-203.
 *      - Kenneth Hillstrom,
 *        ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
 *        May 1969.
 *      - John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
 *        Charles Mesztenyi, John Rice, Henry Thacher, Christoph Witzgall,
 *        Computer Approximations, Wiley, 1968.
 *
 *  Local Parameters:
 *
 *    Local, real ( kind = 8 ) BETA, the radix for the floating-point
 *    representation.
 *
 *    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
 *
 *    Local, real ( kind = 8 ) XBIG, the largest argument for which
 *    LN(GAMMA(X)) is representable in the machine, the solution to the equation
 *      LN(GAMMA(XBIG)) = BETA**MAXEXP.
 *
 *    Local, real ( kind = 8 ) FRTBIG, a rough estimate of the fourth root
 *    of XBIG.
 *
 *  Approximate values for some important machines are:
 *
 *                            BETA      MAXEXP         XBIG     FRTBIG
 *
 *  CRAY-1        (S.P.)        2        8191       9.62E+2461  3.13E+615
 *  Cyber 180/855 (S.P.)        2        1070       1.72E+319   6.44E+79
 *  IEEE (IBM/XT) (S.P.)        2         128       4.08E+36    1.42E+9
 *  IEEE (IBM/XT) (D.P.)        2        1024       2.55E+305   2.25E+76
 *  IBM 3033      (D.P.)       16          63       4.29E+73    2.56E+18
 *  VAX D-Format  (D.P.)        2         127       2.05E+36    1.20E+9
 *  VAX G-Format  (D.P.)        2        1023       1.28E+305   1.89E+76
 *
 * Input
 *  x : the argument of the Gamma function (> 0.0)
 */
double gamma_log( double x ) { 

	double c[7] = { 
		-1.910444077728E-03, 
		8.4171387781295E-04, 
		-5.952379913043012E-04, 
		7.93650793500350248E-04, 
		-2.777777777777681622553E-03, 
		8.333333333333333331554247E-02, 
		5.7083835261E-03 };
	double corr;
	const double D1 = -5.772156649015328605195174E-01;
	const double D2 =  4.227843350984671393993777E-01;
	const double D4 =  1.791759469228055000094023E+00;

	const double FRTBIG = 1.42E+09;
	const double P1[8] = { 
		4.945235359296727046734888E+00, 
		2.018112620856775083915565E+02, 
		2.290838373831346393026739E+03, 
		1.131967205903380828685045E+04, 
		2.855724635671635335736389E+04, 
		3.848496228443793359990269E+04, 
		2.637748787624195437963534E+04, 
		7.225813979700288197698961E+03 };
	const double P2[8] = {
		4.974607845568932035012064E+00,
		5.424138599891070494101986E+02,
		1.550693864978364947665077E+04,
		1.847932904445632425417223E+05,
		1.088204769468828767498470E+06,
		3.338152967987029735917223E+06,
		5.106661678927352456275255E+06,
		3.074109054850539556250927E+06 };
	const double P4[8] = {
		1.474502166059939948905062E+04,
		2.426813369486704502836312E+06,
		1.214755574045093227939592E+08,
		2.663432449630976949898078E+09,
		2.940378956634553899906876E+10,
		1.702665737765398868392998E+11,
		4.926125793377430887588120E+11,
		5.606251856223951465078242E+11 };
	const double PNT68 = 0.6796875E+00;
	const double Q1[8] = {
		6.748212550303777196073036E+01,
		1.113332393857199323513008E+03,
		7.738757056935398733233834E+03,
		2.763987074403340708898585E+04,
		5.499310206226157329794414E+04,
		6.161122180066002127833352E+04,
		3.635127591501940507276287E+04,
		8.785536302431013170870835E+03 };
	const double Q2[8] = {
		1.830328399370592604055942E+02,
		7.765049321445005871323047E+03,
		1.331903827966074194402448E+05,
		1.136705821321969608938755E+06,
		5.267964117437946917577538E+06,
		1.346701454311101692290052E+07,
		1.782736530353274213975932E+07,
		9.533095591844353613395747E+06 };
	const double Q4[8] = {
		2.690530175870899333379843E+03,
		6.393885654300092398984238E+05,
		4.135599930241388052042842E+07,
		1.120872109616147941376570E+09,
		1.488613728678813811542398E+10,
		1.016803586272438228077304E+11,
		3.417476345507377132798597E+11,
		4.463158187419713286462081E+11 };
	double res;
	const double SQRTPI = 0.9189385332046727417803297E+00;
	const double XBIG   = 4.08E+36;
	double xden;
	double xm1;
	double xm2;
	double xm4;
	double xnum;
	double xsq;

	//  Return immediately if the argument is out of range.
	if ( x <= 0.0E+00 || XBIG < x ) {
		return DHUGE;
	}
	if ( x <= D_EPSILON ) {
		res = -log ( x );
	}
	else if ( x <= 1.5E+00 ) {

		if ( x < PNT68 ) {
			corr = - log ( x );
			xm1 = x;
		}
		else {
			corr = 0.0E+00;
			xm1 = ( x - 0.5E+00 ) - 0.5E+00;
		}
 
		if ( x <= 0.5E+00 || PNT68 <= x ) {

			xden = 1.0E+00;
			xnum = 0.0E+00;

			for(int i=0; i<8; i++) {
				xnum = xnum * xm1 + P1[i];
				xden = xden * xm1 + Q1[i];
			}
			res = corr + ( xm1 * ( D1 + xm1 * ( xnum / xden ) ) );
		}
		else {

			xm2 = ( x - 0.5E+00 ) - 0.5E+00;
			xden = 1.0E+00;
			xnum = 0.0E+00;
 			for(int i=0; i<8; i++) {
				xnum = xnum * xm2 + P2[i];
				xden = xden * xm2 + Q2[i];
			}
			res = corr + xm2 * ( D2 + xm2 * ( xnum / xden ) );
		}
	}
	else if ( x <= 4.0E+00 ) {

		xm2 = x - 2.0E+00;
		xden = 1.0E+00;
		xnum = 0.0E+00;
		for(int i=0; i<8; i++) {
			xnum = xnum * xm2 + P2[i];
			xden = xden * xm2 + Q2[i];
		}
		res = xm2 * ( D2 + xm2 * ( xnum / xden ) );
	}       
	else if ( x <= 12.0E+00 ) {

		xm4 = x - 4.0E+00;
		xden = - 1.0E+00;
		xnum = 0.0E+00;
		for(int i=0; i<8; i++) {
			xnum = xnum * xm4 + P4[i];
			xden = xden * xm4 + Q4[i];
		}
		res = D4 + xm4 * ( xnum / xden );
	}
	else {

		res = 0.0E+00;

		if ( x <= FRTBIG ) {

			res = c[6];
			xsq = x * x;
			for(int i=0; i<6; i++) {
				res = res / xsq + c[i];
			}
		}
		res = res / x;
		corr = log ( x );
		res = res + SQRTPI - 0.5E+00 * corr;
		res = res + x * ( corr - 1.0E+00 );
	}
	return res;
}

/* Evaluate a polynomial using Horner's method.
 *
 * The polynomial
 *
 * .. math::
 *
 *     p(x) = c_0 + c_1 * x + c_2 * x^2 + ... + c_m * x^m
 *
 * is to be evaluated at the value X.
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *     - 02 January 2014
 *     - Extracted: November, 2016
 * Author:
 *     - John Burkardt
 *     - Extracted by Philipp Sommer
 * Input
 *  m      : the degree
 *  c(0:m) : the polynomial coefficients. C(I) is the coefficient of  :math:`X^I`
 *  x      : the polynomial value
 */
double r8poly_value_horner ( int m, double *c, double x ) {
	
	double value;

	value = c[m]; 
	for(int i=m-1; i>=0; i--) {
		value = value * x + c[i];
	}
	return value;
}

/* Invert the standard normal CDF.
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *    - 05 June 2007
 *    - Extracted: November, 2016
 * Author:
 *     - Original FORTRAN77 version by Michael Wichura.
 *     - FORTRAN90 version by John Burkardt.
 *     - Extracted by Philipp Sommer
 * Reference:
 *     Michael Wichura,
 *     Algorithm AS241:
 *     The Percentage Points of the Normal Distribution,
 *     Applied Statistics,
 *     Volume 37, Number 3, pages 477-484, 1988.
 *
 * .. note::
 *
 *     The result is accurate to about 1 part in 10^16.
 *
 * Input
 *  p : the value of the cumulative probability densitity function.  0 < P < 1.
 *      If P is outside this range, an "infinite" value will be returned.
 *  x : the normal deviate value with the property that the probability of a
 *      standard normal deviate being less than or equal to the value is P.
 */
void normal_01_cdf_inv (double p,double x) {
	
	double a[8] = {
		3.3871328727963666080E+00,
		1.3314166789178437745E+02,
		1.9715909503065514427E+03,
		1.3731693765509461125E+04,
		4.5921953931549871457E+04,
		6.7265770927008700853E+04,
		3.3430575583588128105E+04,
		2.5090809287301226727E+03 };
	double b[8] = {
		1.0E+00, 
		4.2313330701600911252E+01,
		6.8718700749205790830E+02,
		5.3941960214247511077E+03,
		2.1213794301586595867E+04,
		3.9307895800092710610E+04,
		2.8729085735721942674E+04,
		5.2264952788528545610E+03 };
	double c[8] = {
		1.42343711074968357734E+00,
		4.63033784615654529590E+00,
		5.76949722146069140550E+00,
		3.64784832476320460504E+00,
		1.27045825245236838258E+00,
		2.41780725177450611770E-01,
		2.27238449892691845833E-02,
		7.74545014278341407640E-04 };
	const double CONST1 = 0.180625E+00;
	const double CONST2 = 1.6E+00;
	double d[8] = {
		1.0E+00, 
		2.05319162663775882187E+00,
		1.67638483018380384940E+00,
		6.89767334985100004550E-01,
		1.48103976427480074590E-01,
		1.51986665636164571966E-02,
		5.47593808499534494600E-04,
		1.05075007164441684324E-09 };
	double e[8] = {
		6.65790464350110377720E+00,
		5.46378491116411436990E+00,
		1.78482653991729133580E+00,
		2.96560571828504891230E-01,
		2.65321895265761230930E-02,
		1.24266094738807843860E-03,
		2.71155556874348757815E-05,
		2.01033439929228813265E-07 };
	double f[8] = {
		1.0E+00, 
		5.99832206555887937690E-01,
		1.36929880922735805310E-01,
		1.48753612908506148525E-02,
		7.86869131145613259100E-04,
		1.84631831751005468180E-05,
		1.42151175831644588870E-07,
		2.04426310338993978564E-15 };
	double q;
	double r;
	const double SPLIT1 = 0.425E+00;
	const double SPLIT2 = 5.0E+00;

	if ( p <= 0.0E+00 ) {
		x = - DHUGE ;
		return;
	}

	if ( 1.0E+00 <= p ) {
		x = DHUGE;
		return;
	}
	q = p - 0.5E+00;

	if ( fabs ( q ) <= SPLIT1 ) {
		r = CONST1 - q * q;
		x = q * r8poly_value_horner( 7, a, r ) / r8poly_value_horner( 7, b, r );
	}
	else {
		if ( q < 0.0E+00 ) {
			r = p;
		}
		else {
			r = 1.0E+00 - p;
		}

		if ( r <= 0.0E+00 ) {
			x = DHUGE;
		}
		else {
			r = sqrt ( - log ( r ) );

			if ( r <= SPLIT2 ) {

				r = r - CONST2;
				x = r8poly_value_horner ( 7, c, r ) / r8poly_value_horner ( 7, d, r );
			}
			else {

				r = r - SPLIT2;
				x = r8poly_value_horner ( 7, e, r ) / r8poly_value_horner ( 7, f, r );
			}
		}
 
		if ( q < 0.0E+00 ) {
			x = -x;
		}
	}
}

/* Evaluate the Normal 01 CDF.
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *     - 10 February 1999
 *     - Extracted: June, 2016
 * Author:
 *     - John Burkardt
 *     - Extracted by Philipp Sommer
 * Reference:
 *     AG Adams,
 *     Algorithm 39,
 *     Areas Under the Normal Curve,
 *     Computer Journal,
 *     Volume 12, pages 197-198, 1969.
 * Input
 *  x   : the argument of the CDF.
 *  cdf : the value of the CDF.
 */
void normal_01_cdf ( double x, double cdf ) {
	
	const double A1 = 0.398942280444E+00 ;
	const double A2 = 0.399903438504E+00 ;
	const double A3 = 5.75885480458E+00  ;
	const double A4 = 29.8213557808E+00  ;
	const double A5 = 2.62433121679E+00  ;
	const double A6 = 48.6959930692E+00  ;
	const double A7 = 5.92885724438E+00  ;
	const double B0 = 0.398942280385E+00 ;
	const double B1 = 3.8052E-08	     ;
	const double B2 = 1.00000615302E+00  ;
	const double B3 = 3.98064794E-04     ;
	const double B4 = 1.98615381364E+00  ;
	const double B5 = 0.151679116635E+00 ;
	const double B6 = 5.29330324926E+00  ;
	const double B7 = 4.8385912808E+00   ;
	const double B8 = 15.1508972451E+00  ;
	const double B9 = 0.742380924027E+00 ;
	const double B10 = 30.789933034E+00  ;
	const double B11 = 3.99019417011E+00 ;
	double q;
	double y;

	//  |X| <= 1.28.
	if ( fabs ( x ) <= 1.28E+00 ) {

		y = 0.5E+00 * x * x;

		q = 0.5E+00 - fabs ( x ) * ( A1 - A2 * y / ( y + A3 - A4 /
				( y + A5 + A6 / ( y + A7 ) ) ) );
	}

	//  1.28 < |X| <= 12.7
	else if ( fabs ( x ) <= 12.7E+00 ) {

		y = 0.5E+00 * x * x;

		q = exp ( - y ) * B0 / ( fabs ( x ) - B1 
					 + B2 / ( fabs ( x ) + B3 
					 + B4 / ( fabs ( x ) - B5 
					 + B6 / ( fabs ( x ) + B7 
					 - B8 / ( fabs ( x ) + B9 
					 + B10 /( fabs ( x ) + B11 ) ) ) ) ) );
	}

	//  12.7 < |X|
	else {
		q = 0.0E+00;
	}


	//  Take account of negative X.
	if ( x < 0.0E+00 ) {
		cdf = q;
	}
	else {
		cdf = 1.0E+00 - q;
	}
}

/* Invert the Normal CDF.
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *     - 23 February 1999
 *     - Extracted: November, 2016
 * Author:
 *     - John Burkardt
 *     - Extracted by Philipp Sommer
 *
 * Input
 *  cdf : the value of the CDF. 0.0 <= CDF <= 1.0.
 *  a   : the mean of the pdf
 *  b   : the standard deviation of the pdf
 *  x   : the corresponding argument
 */
void normal_cdf_inv ( double cdf, double a, double b, double x ) {
	
	double x2 = 0.;

	if ( cdf < 0.0E+00 || 1.0E+00 < cdf ) {
		printf("\n ");
		printf("NORMAL_CDF_INV - Fatal error!\n");
		printf("  CDF < 0 or 1 < CDF.\n");
		exit(-1);
	}
	normal_01_cdf_inv ( cdf, x2 );
	x = a + b * x2;
}

/* chi-square approximation for the :f:func:`gamma_cdf_inv` function
 *
 * p   : the quantile
 * nu  : twice the gamma shape
 * g   : the logarithm of the gamma function at the gamma shape
 * tol : the tolerance for the approximation
 */
double qchisq_appr(double p, double nu, double g, double tol) {
	
	double alpha, a, c, ch, p1;
	double p2, t, lgam1pa;

	alpha = 0.5 * nu;
	c = alpha - 1.0;

	p1 = log(p);

	if (nu < (-1.24) * p1) {
		// For small chi-squared 
		//    log(alpha) + g = log(alpha) + log(gamma(alpha)) =
		//       = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
		//    catastrophic cancellation when alpha << 1
		if (alpha < 0.5) 
			lgam1pa = gamma_log(alpha + 1.0);
		else
			lgam1pa = log(alpha) + g;
	    
		ch = exp((lgam1pa + p1)/alpha + log(2.0));
	}
	else if (nu > 0.32) { //  using Wilson and Hilferty estimate
		double x = 0.;
		normal_cdf_inv(p, (double) 0., (double) 1., x);
		p1 = 2. / (9.0 * nu);
		ch = nu * pow((x * sqrt(p1) + 1.0 - p1),3);
		// Approximation for p tending to 1:
		if (ch > 2.2 * nu + 6)
			ch = -2.0 * (log(1 - p) - c * log(0.5 * ch) + g);
	}
	else {
		double q = 0.;
		ch = 0.4;
		a = log(1 - p) + g + c * log(2.0);
		while (fabs(q - ch) > tol * fabs(ch)) {
			q = ch;
			p1 = 1. / (1 + ch * (4.67 + ch));
			p2 = ch * (6.73 + ch * (6.66 + ch));
			t = -0.5 + (4.67 + 2 * ch) * p1 - (6.73 + ch*(13.32 + 3 * ch)) / p2;
			ch = ch - (1 - exp(a + 0.5 * ch) * p2 * p1) / t;
		}
	}
	return ch;
}

/* Compute the incomplete Gamma function.
 *
 * Formulas:
 *     .. math::
 *         \Gamma_{inc}(P, 0) = 0
 *
 *     .. math::
 *         \Gamma_{inc}(P, \infty) = 1.
 *
 *     .. math::
 *         \Gamma_{inc}(P,X) = \int_0^x{T^{P-1} \exp{(-T)} \mathrm{d}t} / \Gamma(P)
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *     - 01 May 2001
 *     - Extracted: June, 2016
 * Author:
 *     - Original FORTRAN77 version by B L Shea.
 *     - FORTRAN90 version by John Burkardt
 *     - Extracted by Philipp Sommer
 * Reference:
 *    BL Shea,
 *    Chi-squared and Incomplete Gamma Integral,
 *    Algorithm AS239,
 *    Applied Statistics,
 *    Volume 37, Number 3, 1988, pages 466-473.
 *
 * Input
 *  p  : the exponent parameter (0.0 < P)
 *  x  : the integral limit parameter. If X is less than or equal to 0, GAMMA_INC is returned as 0.
 */
double gamma_inc ( double p, double x ) {
	
	double  a;
	double  arg;
	double  b;
	double  c;
	double  cdf = 0.;
	const double EXP_ARG_MIN = -88.0E+00;
	const double L_OVERFLOW = 1.0E+37;
	const double PLIMIT = 1000.0E+00;
	double  pn1;
	double  pn2;
	double  pn3;
	double  pn4;
	double  pn5;
	double  pn6;
	double  rn;
	const double TOL = 1.0E-07;
	const double XBIG = 1.0E+08;

	double gamma_inc = 0.0;

	if ( p <= 0.0E+00 ) {
		printf("\n GAMMA_INC - Fatal error! \n");
		printf("  Parameter P <= 0.");
		exit(-1);
	}

	if ( x <= 0.0E+00 ) {
		return 0.0E+00;
	}

	//  Use a normal approximation if PLIMIT < P.
	if ( PLIMIT < p ) {
		pn1 = 3.0E+00 * sqrt ( p ) * ( pow(( x / p ), ( 1.0E+00 / 3.0E+00 )) 
					       + 1.0E+00 / ( 9.0E+00 * p ) - 1.0E+00 );
		normal_01_cdf ( pn1, cdf );
		return cdf;
	} 

	//  Is X extremely large compared to P?
	if ( XBIG < x ) {
		return 1.0E+00;
	}
	//  Use Pearson's series expansion.
	//  (P is not large enough to force overflow in the log of Gamma.
	if ( x <= 1.0E+00 || x < p ) {

		arg = p * log ( x ) - x - gamma_log( p + 1.0E+00 );
		c = 1.0E+00;
		gamma_inc = 1.0E+00;
		a = p;

		for (int j=0; j<1; j+=0) {
			a = a + 1.0E+00;
			c = c * x / a;
			gamma_inc = gamma_inc + c;
			if ( c <= TOL ) {
				break;
			}
		}
		
		arg = arg + log ( gamma_inc );

		if ( EXP_ARG_MIN <= arg ) {
			gamma_inc = exp ( arg );
		} 
		else {
			gamma_inc = 0.0E+00;
		}
	}
	else {
		//  Use a continued fraction expansion.
		arg = p * log ( x ) - x - gamma_log ( p );
		a = 1.0E+00 - p;
		b = a + x + 1.0E+00;
		c = 0.0E+00;
		pn1 = 1.0E+00;
		pn2 = x;
		pn3 = x + 1.0E+00;
		pn4 = x * b;
		gamma_inc = pn3 / pn4;

		for (int j=0; j<1; j+=0) {

			a = a + 1.0E+00;
			b = b + 2.0E+00;
			c = c + 1.0E+00;
			pn5 = b * pn3 - a * c * pn1;
			pn6 = b * pn4 - a * c * pn2;

			if ( 0.0E+00 < fabs ( pn6 ) ) {

				rn = pn5 / pn6;
				
				if ( fabs ( gamma_inc - rn ) <= min ( TOL, TOL * rn ) ) {

					arg = arg + log ( gamma_inc );

					if ( EXP_ARG_MIN <= arg ) {
						gamma_inc = 1.0E+00 - exp ( arg );
					}
					else {
						gamma_inc = 1.0E+00;
					}
					return gamma_inc;
				}
				gamma_inc = rn;
			}

			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;

			//  Rescale terms in continued fraction if terms are large.
			if ( L_OVERFLOW <= fabs ( pn5 ) ) {
				pn1 = pn1 / L_OVERFLOW;
				pn2 = pn2 / L_OVERFLOW;
				pn3 = pn3 / L_OVERFLOW;
				pn4 = pn4 / L_OVERFLOW;
			}
		} 
	}
	return gamma_inc;
}

/* Evaluate the Gamma CDF.
 * 
 *  Licensing:
 *    This code is distributed under the GNU LGPL license.
 *  Modified:
 *    02 January 2000
 *    Extracted: June, 2016
 *  Author:
 *    John Burkardt
 *    Extracted by Philipp Sommer
 * Input
 *  x    : the input value for which to compute the CDF
 *  a    : the location (< `x`) of the gamma distribution (usually 0)
 *  b    : the shape (> 0.0) of the distribution
 *  c    : the scale (>0.0) of the distribution
 * Output
 *  cdf  : the returned value of the CDF
 */
double gamma_cdf( double x, double a, double b, double c ) {

	double p2, x2;

	x2 = ( x - a ) / b;
	p2 = c;
	double cdf = gamma_inc( p2, x2 );
	return cdf;
}

/* Compute the quantile function of the gamma distribution.
 *
 * This function is based on the Applied Statistics Algorithm AS 91
 * ("ppchi2") and via pgamma(.) AS 239.
 *
 * References
 *	    Best, D. J. and D. E. Roberts (1975).
 *	    Percentage Points of the Chi-Squared Distribution.
 *	    Applied Statistics 24, page 385.
 *
 * .. note::
 *
 *     Compared to the original R function, we do not use the final
 *     newton step which might lead to values going to infinity for
 *     quantiles close to 1
 * p:       the quantile between 0 and 1
 * alpha:   the shape of the gamma distribution
 * scale:   the scale of the gamma distribution
 */
double gamma_cdf_inv(double p, double alpha, double scale) {
	
	double a, b, c, g, ch, ch0, p1;
	double p2, q, s1, s2, s3, s4, s5, s6, t;

	const double EPS1    = 1.0e-2;
	const double EPS2    = 5.0e-7;
	const double PMIN    = 1.0e-25;
	const double PMAX    = (1-1e-14) ;
	const double PNEGINF = -10e34;
	const int MAXIT = 1000;
	
	int i;
	
	double gamma_cdf_inv = 0.;

	if (alpha == 0) {
		return gamma_cdf_inv = 0;
	}
	
	g = gamma_log(alpha);
	// Phase I : Starting Approximation
	ch = qchisq_appr(p, 2.0 * alpha, g, EPS1);

	if ((ch < EPS2) || (p > PMAX) || (p < PMIN) )
		return gamma_cdf_inv = 0;
	
	// Phase II: Iteration
	// Call pgamma() [AS 239] and calculate seven term taylor series
	c  = alpha - 1.0;
	s6 = (120.0 + c * (346.0 + 127.0 * c)) / 5040.0;

	ch0 = ch;  // save initial approx.
	for ( i=0; i<MAXIT; i++) {

		q = ch;
		p1 = 0.5 * ch;

		p2 = gamma_cdf(p1 * scale, (double)0., scale, alpha);
		p2 = p - p2;

		if ((p2 < PNEGINF) || ch <= 0.0) {
			ch = ch0;
			break;
		}

		t = p2*exp(alpha * log(2.0) + g + p1 - c * log(ch));
		b = t / ch;
		a = 0.5 * t - b * c;
		s1 = (210.0 + a * (140.0 + a * (105.0 + a * (84.0 + a * (70.0 + 60.0 * a))))) / 420.0;
		s2 = (420.0 +  a * (735.0 + a * (966.0 + a * (1141.0 + 1278.0 * a)))) / 2520.;
		s3 = (210.0 + a * (462.0 + a * (707.0 + 932.0 * a))) / 2520.0;
		s4 = (252.0 + a * (672.0 + 1182.0 * a) + c * (294.0 +a * ( 889.0 + 1740.0 * a))) / 5040.0;
		s5 = (84.0 + 2264.0 * a + c*(1175.0 + 606.0 * a)) / 2520.0;

		ch = ch +  t * (1.0 + 0.5 * t * s1 - b * c * ( 
			s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));

		if (fabs(q - ch) < EPS2 * ch)
			break;

		if (fabs(q - ch) > 0.1 * ch) {
			if (ch < q) 
				ch = 0.9 * q;
			else
				ch = 1.1 * q;
		}
	}

	gamma_cdf_inv = 0.5  * scale * ch;
	return gamma_cdf_inv;
}
 
/* Evaluate Gamma(X) for a real argument.
 *
 * This routine calculates the gamma function for a real argument X.
 *
 * Computation is based on an algorithm outlined in reference 1.
 * The program uses rational functions that approximate the gamma
 * function to at least 20 significant decimal digits.  Coefficients
 * for the approximation over the interval (1,2) are unpublished.
 * Those for the approximation for 12 <= X are from reference 2.
 *
 * Modified:
 *     - 11 February 2008
 *     - Extracted: June, 2016
 *
 * Author:
 *     - Original FORTRAN77 version by William Cody, Laura Stoltz.
 *     - FORTRAN90 version by John Burkardt.
 *     - Extracted by Philipp Sommer
 * Reference:
 *     - William Cody,
 *       An Overview of Software Development for Special Functions,
 *       in Numerical Analysis Dundee, 1975,
 *       edited by GA Watson, Lecture Notes in Mathematics 506,
 *       Springer, 1976.
 *     - John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
 *       Charles Mesztenyi, John Rice, Henry Thatcher,
 *       Christoph Witzgall, Computer Approximations, Wiley, 1968,
 *       LC: QA297.C64.
 *
 * Input
 *  x : the argument of the function.
 *
 *  Coefficients for minimax approximation over (12, INF).
 */
double r8_gamma ( double x ) {
	
	double c[7] = {
		-1.910444077728E-03,
		8.4171387781295E-04, 
		-5.952379913043012E-04, 
		7.93650793500350248E-04, 
		-2.777777777777681622553E-03, 
		8.333333333333333331554247E-02, 
		5.7083835261E-03 };
	const double L_EPS = 2.22E-16;
	double fact;
	const double HALF = 0.5E+00;
	int n;
	const double ONE = 1.0E+00;
	double p[8] = {
		-1.71618513886549492533811E+00,
		2.47656508055759199108314E+01,
		-3.79804256470945635097577E+02,
		6.29331155312818442661052E+02,
		8.66966202790413211295064E+02,
		-3.14512729688483675254357E+04,
		-3.61444134186911729807069E+04,
		6.64561438202405440627855E+04 };
	bool parity = false;
	double q[8] = {
		-3.08402300119738975254353E+01,
		3.15350626979604161529144E+02,
		-1.01515636749021914166146E+03,
		-3.10777167157231109440444E+03,
		2.25381184209801510330112E+04,
		4.75584627752788110767815E+03,
		-1.34659959864969306392456E+05,
		-1.15132259675553483497211E+05 };
	const double PI = 3.1415926535897932384626434E+00;
	const double SQRTPI = 0.9189385332046727417803297E+00;
	const double TWELVE = 12.0E+00;
	const double TWO    = 2.0E+00;
	const double XBIG   = 171.624E+00;
	const double XINF   = 1.0E+30;
	const double XMININ = 2.23E-308;
	double res;
	double sum;
	double xden;
	double xnum;
	double y;
	double y1;
	double ysq;
	double z;
	const double ZERO = 0.0E+00;

	double r8_gamma = 0.;
	
	fact = ONE;
	n = 0;
	y = x;

	//  Argument is negative.
	if ( y <= ZERO ) {

		y = - x;
		y1 = (int) y;
		res = y - y1;
		if ( res != ZERO ){

			if ( y1 != (int) ( y1 * HALF ) * TWO ) {
				parity = true;
			}
			fact = - PI / sin ( PI * res );
			y = y + ONE;
		}
		else {
			res = XINF;
			r8_gamma = res;
			return r8_gamma;
		}
	}

	//  Argument is positive.
	if ( y < L_EPS ) {

		//  Argument < EPS.
		if ( XMININ <= y ) {
			res = ONE / y;
		}
		else {
			res = XINF;
			r8_gamma = res;
			return r8_gamma;
		}
	}
	else if ( y < TWELVE ) {

		y1 = y;

		//  0.0 < argument < 1.0.
		if ( y < ONE ) {
			z = y;
			y = y + ONE;
		}
		//  1.0 < argument < 12.0.
		//  Reduce argument if necessary.
		else {

			n = int( y ) - 1;
			y = y - (double) n;
			z = y - ONE;
		}

		//  Evaluate approximation for 1.0 < argument < 2.0.
		xnum = ZERO;
		xden = ONE;
		for(int i=0; i<8; i++) {
			xnum = ( xnum + p[i] ) * z;
			xden = xden * z + q[i];
		}

		res = xnum / xden + ONE;

		//  Adjust result for case  0.0 < argument < 1.0.
		if ( y1 < y ) {
			res = res / y1;
		}
		//  Adjust result for case 2.0 < argument < 12.0.
		else if ( y < y1 ) {
			for(int i=0; i<n; i++) {
				res = res * y;
				y = y + ONE;
			}
		}
	}
	else {
		//  Evaluate for 12.0 <= argument.
		if ( y <= XBIG ) {

			ysq = y * y;
			sum = c[6];
 			for(int i=0; i<6; i++) {
				sum = sum / ysq + c[i];
			}
			sum = sum / y - y + SQRTPI;
			sum = sum + ( y - HALF ) * log ( y );
			res = exp ( sum );
		}
		else {
			res = XINF;
			r8_gamma = res;
			return r8_gamma;
		}
	}

	//  Final adjustments and return.
	if ( parity ) {
		res = - res;
	}

	if ( fact != ONE ) {
		res = fact / res;
	}
	r8_gamma = res;

	return r8_gamma ;
}

/* Evaluate the Gamma PDF.
 *
 * .. math::
 *
 *     PDF(a,b,c;x) = \exp({-(x-a)/b}) \cdot ((x-a)/b)^{c-1} / (b \cdot \Gamma(c))
 *
 * - GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
 * - GAMMA_PDF(A,B,1;X) is the Exponential PDF.
 * - GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.
 *
 * Licensing:
 *     This code is distributed under the GNU LGPL license.
 * Modified:
 *     - 02 January 2000
 *     - Extracted: June, 2016
 * Author:
 *     - John Burkardt
 *     - Extracted by Philipp Sommer
 *
 * Input
 *  x    : the argument of the PDF. A <= X
 *  a    : the location of the peak;  A is often chosen to be 0.0.
 *  b    : the "scale" parameter; 0.0 < B, and is often 1.0.
 *  c    : the "shape" parameter; 0.0 < C, and is often 1.0.
 * Output
 *  pdf  : the returned value of the PDF.
 */
double gamma_pdf ( double x, double a, double b, double c ) {
	
	double y;
	double pdf = 0.;

	if ( x <= a ) {
		pdf = 0.0E+00;
	}
	else {
		y = ( x - a ) / b;
		pdf = pow(y,( c - 1.0E+00 )) / ( b * r8_gamma ( c ) * exp ( y ) );
	}
	return pdf;
}
									     
/* Iterative, mean preserving method to smoothly interpolate mean data to pseudo-sub-timestep values
 * From Rymes, M.D. and D.R. Myers, 2001. Solar Energy (71) 4, 225-231
 * Input
 *  lm     : left margin month
 *  rm     : right margin month
 *  m      : vector of mean values at super-time step (e.g., monthly), minimum three values
 *  dmonth : vector of number of intervals for the time step (e.g., days per month)
 *  bcond  : boundary conditions for the result vector (1=left side, 2=right side)
 *  r      : result vector of values at chosen time step
 */
void rmsmooth(int lm,int rm, double *m,int *dmonth,double bcond[2], double *m_curr) {
	
	// Parameters
	double const OT = 1. / 3.;

	// Local variables
	int g[100];	
	double r[100];

	double ck;

	int ni = 0;
	for (int i=lm; i<=rm; i++)
		ni +=dmonth[i];	
	
	double bc[2];
	bc[0] = bcond[0];
	bc[1] = bcond[1];

	// Initialize the result vector
	int i = 0;
	int j = 0;
	int a = 0;
	int b = 0;
	for ( a=lm; a<=rm; a++) {  
		j=i;
		for ( b=0; b<dmonth[a]; b++) {
			r[i] = m[a];
			g[i] = j;
			i++;
		}
	}
	
	// Iteratively smooth and correct the result to preserve the mean
	for (i=0;i<ni;i++) {
		for (j=1;j<ni-1;j++) {
			r[j] = OT * (r[j-1] + r[j] + r[j+1]);   //Eqn. 1
		}
		r[0]    = OT * (bc[0]  + r[0]  +  r[1]);        //Eqns.2
		r[ni-1] = OT * (r[ni-2] + r[ni-1] + bc[1]);

		j=0;
		for (int k=lm;k<=rm;k++) {      // calculate one correction factor per super-timestep
			a = g[j];               // index of the first timestep value of the super-timestep
			b = g[j] + dmonth[k] ;  // index of the last timestep value of the super-timestep
			
			ck = 0.;
			for (int x=a; x<b; x++) {
				ck += m[k] - r[x];
			}
			ck /= (double)ni;       // !Eqn. 4

			// Apply the correction to all timestep values in the super-timestep		
			for(int l=0; l<dmonth[k]; l++) { 
				r[j] += ck;
				j++;
			}
			// Correction for circular conditions when using climatology 
			// (do not use for transient simulations)
			if (CIRCULAR_CLIM) {
				bc[0] = r[ni-1];
				bc[1] = r[0];
			}
		}
	}
	for (int k=0; k<dmonth[1]; k++) {
		if (lm==0) {
			m_curr[k] = r[k+dmonth[0]];
		}
		else {
			m_curr[k] = r[k];
		}
	}
}

// Initialize the weather generator	
void init_weathergen(MetVariables& metvars, WeatherGenState& rndst) {

	metvars.pday[0] = false;
	metvars.pday[1] = false;
	
	for (int i=0;i<4;i++) {
		metvars.resid[i] = 0.;
	}
	
	for (int i=0;i<QSIZ;i++) {
		rndst.q[i] = 0;
	}
	
	rndst.carry =       362;
	rndst.xcng  =   1236789;
	rndst.xs    = 521288629; //default seed
	rndst.indx  = QSIZ+1;
	rndst.have  = false;
	
	for (int i=0;i<2;i++) {
		rndst.gamma_vals[i] = 0.; 
	}
}

/* Calculation of daylength, insolation, and equilibrium evapotranspiration
 * for each day, given mean daily temperature, insolation (as percentage
 * of full sunshine or mean daily instantaneous downward shortwave
 * radiation flux, W/m2), latitude and day of year
 *
 * INPUT AND OUTPUT PARAMETER
 * climate = gridcell climate
 */
double cldf2rad(double input, double lat, int doy, bool cldf2rad) {

	const double QOO = 1360.0;
	const double BETA = 0.17;
	const double C = 0.25;
	const double D = 0.5;
	const double K = 13750.98708;
	const double PI = 3.14159269;
	const double DEGTORAD    = PI / 180.;
	const double DAYLENGTHS  = 86400.;
	const double YEARLENGTHD = date.year_length();

	double w,u,v,qo,delta,hh,sinehh;
	double rad, cldfr;

	//	CALCULATION OF NET DOWNWARD SHORT-WAVE RADIATION FLUX
	//	Refs: Prentice et al 1993, Monteith & Unsworth 1990,
	//	      Henderson-Sellers & Robinson 1986

	//	 (1) rs = (c + d*ni) * (1 - beta) * Qo * cos Z * k
	//	       (Eqn 7, Prentice et al 1993)
	//	 (2) Qo = Qoo * ( 1 + 2*0.01675 * cos ( 2*pi*(i+0.5)/365) )
	//	       (Eqn 8, P8rentice et al 1993; angle in radians)
	//	 (3) cos Z = sin(lat) * sin(delta) + cos(lat) * cos(delta) * cos h
	//	       (Eqn 9, Prentice et al 1993)
	//	 (4) delta = -23.4 * pi / 180 * cos ( 2*pi*(i+10.5)/365 )
	//	       (Eqn 10, Prentice et al 1993, angle in radians)
	//	 (5) h = 2 * pi * t / 24 = pi * t / 12

	//	     where rs    = instantaneous net downward shortwave radiation
	//	                   flux, including correction for terrestrial shortwave albedo
	//	                   (W/m2 = J/m2/s)
	//	           c, d  = empirical constants (c+d = clear sky
	//	                   transmissivity)
	//	           ni    = proportion of bright sunshine
	//	           beta  = average 'global' value for shortwave albedo
	//	                   (not associated with any particular vegetation)
	//	           i     = julian day, (0-364, 0=1 Jan)
	//	           Qoo   = solar constant, 1360 W/m2
	//	           Z     = solar zenith angle (angular distance between the
	//	                   sun's rays and the local vertical)
	//	           k     = conversion factor from solar angular units to
	//	                   seconds, 12 / pi * 3600
	//	           lat   = latitude (+=N, -=S, in radians)
	//	           delta = solar declination (angle between the orbital
	//	                   plane and the Earth's equatorial plane) varying
	//	                   between +23.4 degrees in northern hemisphere
	//	                   midsummer and -23.4 degrees in N hemisphere
	//	                   midwinter
	//	           h     = hour angle, the fraction of 2*pi (radians) which
	//	                   the earth has turned since the local solar noon
	//	           t     = local time in hours from solar noon

	//	From (1) and (3), shortwave radiation flux at any hour during the
	//	day, any day of the year and any latitude given by
	//	 (6) rs = (c + d*ni) * (1 - beta) * Qo * ( sin(lat) * sin(delta) +
	//	          cos(lat) * cos(delta) * cos h ) * k
	//	Solar zenith angle equal to -pi/2 (radians) at sunrise and pi/2 at
	//	sunset.  For Z=pi/2 or Z=-pi/2,
	//	 (7) cos Z = 0
	//	From (3) and (7),
	//	 (8)  cos hh = - sin(lat) * sin(delta) / ( cos(lat) * cos(delta) )
	//	      where hh = half-day length in angular units
	//	Define
	//	 (9) u = sin(lat) * sin(delta)
	//	(10) v = cos(lat) * cos(delta)
	//	Thus
	//	(11) hh = acos (-u/v)
	//	To obtain the daily net downward short-wave radiation sum, integrate
	//	equation (6) from -hh to hh with respect to h,
	//	(12) rad = 2 * (c + d*ni) * (1 - beta) * Qo *
	//	              ( u*hh + v*sin(hh) )
	//	Define
	//	(13) w = (c + d*ni) * (1 - beta) * Qo
	//	From (12) & (13), and converting from angular units to seconds
	//	(14) rad = 2 * w * ( u*hh + v*sin(hh) ) * k

	// Calculate values of saved parameters for this day
	qo = QOO * (1.0 + 2.0 * 0.01675 * cos(2.0 * PI * ((double)doy + 0.5) / YEARLENGTHD)); // Eqn 2
	delta = -23.4 * DEGTORAD * cos(2.0 * PI * ((double)doy + 10.5) / YEARLENGTHD);
				// Eqn 4, solar declination angle (radians)
	u = sin(lat * DEGTORAD) * sin(delta); // Eqn 9
	v = cos(lat * DEGTORAD) * cos(delta); // Eqn 10

	if (u >= v) {
		hh = PI; // polar day
	}
	else if (u <= -v) {
		hh = 0.0; // polar night
	}
	else { 
		hh = acos(-u / v); // Eqn 11
	}

	sinehh = sin(hh);
	
	if ( cldf2rad ) {
		cldfr = input;
		w    = (C+D * (1.-cldfr)) * (1.0 - BETA) * qo; // Eqn 13
		rad  = 2.0 * w * (u * hh + v * sinehh) * K; // Eqn 14
		rad /= DAYLENGTHS;
		return rad;
	}
	else {
		rad   = input*DAYLENGTHS;
		if ( hh > 0. ) {
			w     = rad / (2. *(u * hh + v * sinehh) * K);
			cldfr = 1.-((w/((1.0 - BETA) * qo) -C)/D);
			cldfr = max(0.,min(1.,cldfr));
		} 
		else {
			cldfr = 0.;
		}
		return cldfr;
	}
}

// Compute daily weather data from monthly input
void weathergen_get_daily_met(MetVariables& metvars, WeatherGenState& rndst) {

	// Local variables
	int i = 0;

	// Monthly total precipitation amount (mm)
	double pre = metvars.mprec; 
	// Number of days in month with precipitation (fraction)
	double wetd= metvars.mwetd;  
	// Fraction of days in month with precipitation (fraction) 
	double wetf= metvars.mwetf;
	// Minumum temperture (C)
	double tmn = metvars.dtmin; 
	// Maximum temperture (C)
	double tmx = metvars.dtmax;  
	// Cloud fraction (0=clear sky, 1=overcast) (fraction)
	double cld = metvars.dcldf;  
	// Wind (m/s)
	double wnd = metvars.dwind;  

	double prec;
	double tmin;
	double tmax;
	double cldf;
	double wind;

	double pbar = 0.;     // mean amount of precipitation per wet day (mm)
	double pwet = 0.;     // probability that today will be wet
	double u    = 0.;     // uniformly distributed random number (0-1)

	double g_shape     = 0.;
	double g_scale	   = 0.;
	double gp_scale	   = 0.;
	double thresh2use  = 0.;

	// Bias correction
	double slopecorr       = 0.;  // slope correction for wind
	double intercept_corr  = 0.;  // intercept correction for wind
	double tmin_bias       = 0.;  // intercept correction for tmin

	double cdf_thresh   = 0.;  // gamma cdf at the threshold
	double pdf_thresh   = 0.;  // gamma pdf at the threshold

	double unorm[4];  // vector of uniformly distributed random numbers (0-1)
	
	// Precipitation occurrence
	// If there is precipitation this month, calculate the precipitation state for today
	if (wetf > 0. && pre > 0.) {
		
		// Calculate transitional probabilities for dry to wet and wet to wet days
		// relationships from Geng & Auburn, 1986, Weather simulation models 
		// based on summaries of long-term data
		
		// If yesterday was raining, use p11
		if (metvars.pday[0]) { 
			pwet = p11_1 + p11_2 * wetf;
		}
		// If yesterday was not raining but the day before yesterday was raining, use p101
		else if (metvars.pday[1]) { 
			pwet = p101_1 + p101_2 * wetf;
		}
		// Both yesterday and the day before were dry, use p001
		else {  
			pwet = p001_1 + p001_2 * wetf;
		}

		// Determine the precipitation state of the current day 
		// using the Markov chain approach
		u = ranur(rndst);
		
		metvars.pday[1] = metvars.pday[0];
		if (u <= pwet) { // today is a rain day
			metvars.pday[0] = true;
		}
		else  { //today is dry
			metvars.pday[0] = false;
		}

		// Precipitation amount
		
		if (metvars.pday[0]) { //today is a wet day, calculate the rain amount
			// Calculate parameters for the distribution function of precipitation amount
			pbar = pre / wetd;

			g_scale = g_scale_coeff * pbar;
			g_shape = pbar / g_scale;
			
			if (thresh_pctl) {
				thresh2use = gamma_cdf_inv(thresh, g_shape, g_scale);
			}
			else {
				thresh2use = thresh;
			}
			
			cdf_thresh = gamma_cdf(thresh2use, (double)0.0, g_scale, g_shape);
			pdf_thresh = gamma_pdf(thresh2use, (double)0.0, g_scale, g_shape);

			gp_scale = (1.0 - cdf_thresh)/ pdf_thresh;
			
			for  (i=0; i<1000; i++) { // enforce positive precipitation
				
				// Today's precipitation
				prec = ran_gamma_gp(rndst,true,g_shape,g_scale,thresh2use,gp_shape,gp_scale);
				
				if (prec > 0. && prec <= 1.05 * pre) {
					// Simulated precipitation should have no more precision than the input (0.1mm)
					prec = roundoff(prec,1);
					break;
				}
				if (i == 1000) { 
					printf("Could not find good precipitation with %f mm and %f wet days...",pre,wetd);
					exit(-1);
				}
			}
		}
		else {
			prec = 0.;
		}
	}
	else {
		metvars.pday[0] = false;
		metvars.pday[1] = false;
		prec = 0.;
	}

	// Temperature min and max, cloud fraction, wind.
	// Calculate a baseline mean and SD for today's weather depending on precip status.
	metvars.tmn = tmn;
	metvars.tmx = tmx;
	metvars.cld = cld;
	metvars.wnd = wnd;

	meansd(metvars);

	// Use random number generator for the normal distribution
	for (i=0;i<4;i++) {
		unorm[i] = ran_normal(rndst);
	}

	// Calculate today's residuals for weather variables
	double CC[4],DD[4];
	matrixmult(A,metvars.resid,CC);
	matrixmult(B,unorm,DD);

	for (int j=0; j<4; j++) {
		metvars.resid[j] = CC[j]+DD[j];
	}

	tmin = roundoff(TMIN_DAMP * metvars.resid[0] * metvars.dmtmin_sd + metvars.dmtmin_mn,1);
	tmax = roundoff(TMAX_DAMP * metvars.resid[1] * metvars.dmtmax_sd + metvars.dmtmax_mn,1);
	cldf = metvars.resid[2] * metvars.dmcldf_sd + metvars.dmcldf_mn;
	wind = max(0.0, metvars.resid[3] * sqrt(max(0.0, metvars.dmwind_sd)) + sqrt(max(0.0, metvars.dmwind_mn)));
	wind = roundoff(wind * wind, 1);
	
	// wind bias correction
	slopecorr = 0.;
	if (wind_slope_bias_L > 0.0) {
		slopecorr = wind_slope_bias_L / ( 1 + exp( - wind_slope_bias_k * 
							   ( metvars.resid[3]) - wind_slope_bias_x0));
	}
	else {
		for (i=0; i<6; i++) {
			slopecorr += wind_bias_coeffs[i] * 
				(pow(max(wind_bias_min, min(wind_bias_max, metvars.resid[3])), i));
		}
	}

	intercept_corr = 0.;
	if (fabs(wind_intercept_bias_a + 9999.) > 1e-7) {
		intercept_corr = exp(wind_intercept_bias_b + 
				 wind_intercept_bias_a * max(wind_bias_min, min(wind_bias_max, metvars.resid[3])));
	} 
	else {
		for (i=0; i<6; i++) {
			intercept_corr += wind_intercept_bias_coeffs[i] * 
				(pow(max(wind_bias_min, min(wind_bias_max, metvars.resid[3])),i));
		}
	}

	wind = (wind - intercept_corr) / max(slopecorr, 9e-4);

	// tmin bias correction
	for (int i=0; i<6; i++) {
		tmin_bias += tmin_bias_coeffs[i] * 
				 (pow(max(tmin_bias_min, min(tmin_bias_max, metvars.resid[0])),i));
	}
	tmin = tmin - roundoff(tmin_bias, 1);

	// Add checks for invalid values here
	if (cldf>1.) {
		cldf = 1.0;
	}
	else if (cldf < 0.0) {
		// Bugfix for negative cldf allows for redistr on initial vals
		cldf = metvars.dcldf * 0.001;
	}

	if (wind<=0.) {
		wind = metvars.dwind * 0.1;
	}

	if (tmin+TFREEZE < 0.) {
		printf("Unphysical min. temperature with %f K from a monthly mean %f degC with bias correction %f K for residual %f", tmin, metvars.dmtmin_mn,tmin_bias,metvars.resid[0]);
		exit(-1);
	}
	else if (tmax+TFREEZE < 0.) {
		printf("Unphysical max. temperature with %f K from a monthly mean %f degC",tmax, metvars.dmtmax_mn);
		exit(-1);
	}

	// Repopulate daily arrays
 	metvars.dprec  = prec;
	metvars.dtmin  = tmin;
	metvars.dtmax  = tmax;
	metvars.dcldf  = cldf;
	metvars.dwind  = wind;

}

/* Redistribute daily values when there is a max,min limit or both
 * while ensuring conservation and relative distribution (in terms of <=,>=)
 * like cloud-fraction or relative humidity ([0,1]).
 */
void redist_restricted_vals(double *inval, int ll, double scalval, double *limit, double *wght) {   
	
	int cnt       = 0;
	double rest   = 0.;
	bool flag[31] ;		
	for (int i=0; i<ll; i++) {
		flag[i] = true ;
	}
	double corfac;

	if ( limit[0] != 0. ) {
		fail("Error in weathergen.cpp redist_restricted_vals()\n");
	}

	bool go = true;
	while ( go ) {
		double remsum = 0.; 
		double gonsum = 0.; 
		for (int i=0; i<ll; i++) {
			if ( flag[i] )
				remsum += inval[i] * wght[i] / (double)ll;
			else 
				gonsum += inval[i] * wght[i] / (double)ll;
		}
		corfac  = remsum / ( scalval - gonsum );

		// Now check if active values have exceeded upper limit
		rest = 0.;
		for (int i=0; i<ll; i++) {
			if ( flag[i] ) {
				inval[i] *= wght[i] / corfac;
				if (inval[i] > limit[1]) {
					rest    += inval[i] - limit[1];
					inval[i] = limit[1];
					flag [i] = false;
				}
			}
		}
		cnt++;
		if ( rest<0.000001 * limit[1] || cnt>30 )
			go = false;
	}
	if ( rest > 0.000001 * limit[1] ) {
		fail ("Redistribution in weathergen.cpp failed!");
	}
}

/* Compute relative humidity dericed from Buck 1981
 * input
 * T_avg: daily mean temperature[°C]
 * T_dew: dew-point temperature[°C]
 * output
 * relhum: estimated daily mean relative humidity [fract.]
 */
double get_arden_rh(double T_avg, double T_dew) {
	
	const double B = 18.678;
	const double C = 257.14; // °C
	const double D = 234.5 ; // °C

	double rh1 = B * T_dew / (C + T_dew);
	double rh2 = ( B - T_avg/D )* T_avg / (C + T_avg);
	double relhum = min(1.,max(0., exp(rh1 - rh2)));
	return relhum;
}

double correlation(int len, double *xarr, double *yarr) {

	double correlation = 0.;
	double osum, xsum2, xmean, ysum2, ymean;
	osum = xsum2 = xmean = ysum2 = ymean = 0.;

	for (int i=0; i<len; i++) {
		xmean += xarr[i]; 
		ymean += yarr[i];
	}
	xmean /= (double) len;
	ymean /= (double) len;

	for (int i=0; i<len; i++) {
		osum  += (xarr[i]-xmean)*(yarr[i]-ymean); 
		xsum2 += (xarr[i]-xmean)*(xarr[i]-xmean);
		ysum2 += (yarr[i]-ymean)*(yarr[i]-ymean);
	}
	correlation = osum / sqrt(xsum2*ysum2);
	return correlation;
}

/* The driver routine of GWGEN.
 * Computes one year's daily met data from monthly averages.
 */
void weathergen_get_met(Gridcell& gridcell, double* in_mtemp, double* in_mprec, double* in_mwetd, 
		   double* in_msol, double* in_mdtr, double* in_mwind, double* in_mrhum, 
		   double* out_dtemp, double* out_dprec, double* out_dsol,double* out_ddtr,
		   double* out_dwind, double* out_drhum) {

	bool is_first_day = ( date.day == 0 && ( date.year == 0 || 
		( restart && date.year == state_year ) ) );

	WeatherGenState& rndst = gridcell.climate.weathergenstate;

	double lat = gridcell.get_lat();
	double lon = gridcell.get_lon();

	// Meteorological vars derived from input vars mtemp,mdtr,msol
	double in_mtmin[12];
	double in_mtmax[12];
	double in_mcldf[12];

	int doy = 1;  
	for (int m=0; m<12; m++) {
		// Use mid-of-month length-of-day 
		int ndaymon = date.ndaymonth[m];
		in_mtmin[m] = in_mtemp[m] - 0.5 * in_mdtr[m];
		in_mtmax[m] = in_mtemp[m] + 0.5 * in_mdtr[m];

		in_mcldf[m] = 0.;
		for ( int day=0; day<ndaymon; day++) {
			in_mcldf[m] += cldf2rad(in_msol[m], lat, doy, false);	
			doy++;
		}
		in_mcldf[m] /= (double)ndaymon;
		// Have a min cldf of 1% to introduce a monthly variability 
		// to fit lower sol vals with rainfall
		in_mcldf[m] = max(0.01,in_mcldf[m]);

		if (in_mwetd[m]==0 && in_mprec[m]>1.) {
			in_mwetd[m] = 1.;
		}
	}

	int accumday = 0;

	for (int mon=0;mon<12;mon++) {
 
		const int NDAYMONTH = 31;
		int ndaymon = date.ndaymonth[mon]; 
		//
		double mtmin_curr [NDAYMONTH];
		double mtmax_curr [NDAYMONTH];
		double mcloud_curr[NDAYMONTH];
		double mwind_curr [NDAYMONTH];

		// Intermediate daily values
		double dprec[NDAYMONTH];
		double dtmin[NDAYMONTH];
		double dtmax[NDAYMONTH];
		double dcldf[NDAYMONTH];
		double dwind[NDAYMONTH];
		double drhum[NDAYMONTH];
		double dsol [NDAYMONTH];
		double cldwght[NDAYMONTH];
		double dprec_sav[NDAYMONTH];
		double dtmin_sav[NDAYMONTH];
		double dtmax_sav[NDAYMONTH];
		double dcldf_sav[NDAYMONTH];
		double dwind_sav[NDAYMONTH];
		
		// Break-off parameters
		int pdaydiff    = 0;
		double precdiff = 0.; 
		double tmindiff = 0.;

		int i_count = 1;
		// Initially populate cloud params
		if ( is_first_day && mon==0 ) {

			calc_cloud_params(metvars);

			// Set initial vals if spinning up
			if ( ! restart ) {
				init_weathergen(metvars, rndst);
				get_seed_by_location(lat,lon ,rndst);
				i_count = 0;
			}
			else {

				// Get restart values from WeatherGen-class
				metvars.pday[0] = rndst.pday[0];
				metvars.pday[1] = rndst.pday[1];
				for (int i=0; i<4;i++) {
					metvars.resid[i] = rndst.resid[i];
				}
			}	
		
		}

		// At beginning of month:
		if ( mon > 0 )  
			accumday += date.ndaymonth[mon-1];

		// Dummy weighting array
		double dum[NDAYMONTH];
		for (int day=0; day<ndaymon; day++)
			dum[day] = 1.;
		
		// Index for annual arrays
		int lm = max(mon-1,0); 
		int rm = min(11,mon+1);

		// Index for rmsmooth
		int ilm = 0;
		int irm = 2;
		if ( mon == 0 ) {
			ilm = 1;
		}
		else if ( mon == 11 ) { 
			irm = 1;
		}

		double bcond[2];
		double tmvals[3];
		int cmdays[3] = 
			{date.ndaymonth[lm],date.ndaymonth[mon],date.ndaymonth[rm]};

		// Smooth tmin
 		bcond[0]  = in_mtmin[lm] ;
		bcond[1]  = in_mtmin[rm] ;
		tmvals[0] = in_mtmin[lm] ;
		tmvals[1] = in_mtmin[mon];
		tmvals[2] = in_mtmin[rm] ;
		rmsmooth( ilm,irm,tmvals,cmdays,bcond, mtmin_curr );

		// Smooth tmax
		bcond[0]  = in_mtmax[lm];
		bcond[1]  = in_mtmax[rm];
		tmvals[0] = in_mtmax[lm];
		tmvals[1] = in_mtmax[mon];
		tmvals[2] = in_mtmax[rm];
		rmsmooth( ilm,irm,tmvals,cmdays,bcond, mtmax_curr );
		
		// Smooth cloud cover
		bcond[0]  = in_mcldf[lm];
		bcond[1]  = in_mcldf[rm];
		tmvals[0] = in_mcldf[lm];
		tmvals[1] = in_mcldf[mon];
		tmvals[2] = in_mcldf[rm];
		rmsmooth( ilm,irm,tmvals,cmdays,bcond, mcloud_curr );
		
		// Ensure positivity for cloud-cover
		for (int day=0;day<ndaymon;day++) {
			mcloud_curr[day] = max(0.01,mcloud_curr[day]); 
		}

		// Smooth wind
		bcond[0]  = in_mwind[lm];
		bcond[1]  = in_mwind[rm];
		tmvals[0] = in_mwind[lm];
		tmvals[1] = in_mwind[mon];
		tmvals[2] = in_mwind[rm];
		rmsmooth( ilm,irm,tmvals,cmdays,bcond, mwind_curr );

		// Ensure positivity for wind
		for (int day=0;day<ndaymon;day++) {
			mwind_curr[day] = max(0.1,mwind_curr[day]); 
		}
		
		// Reset residuals at beginning of year if desired
		if (LRESET && mon==0) {
			for (int i=0;i<4;i++)
				metvars.resid[i] = 0.;
			i_count = 0;
		}
		
		// Set quality threshold for preciptation amount
		double prec_t = max(0.5,0.05 * in_mprec[mon]);
		
		metvars.mprec = in_mprec[mon];

		MetVariables metvar_sav = metvars;
		
		// Here, a bugfix for CRU data is applied when there is non-zero rain
		// but no wet days.
		if ( metvars.mprec > 0. ) {
			metvars.mwetd = max(1.,in_mwetd[mon]);
		} else {
			metvars.mwetd = 0.;
		}
		metvars.mwetf = metvars.mwetd/(double)ndaymon;;
		
		double metric_sav = 99999.;
		
		double chk_dtemp = 0.; 
		double chk_ddtr  = 0.; 
		double chk_dprec = 0.; 
		double chk_dsol  = 0.; 
		double chk_dwind = 0.; 
		double chk_drhum = 0.;

		// Set breakoff-threshold for raindays according to
		// total amount of raindays in month.
		int pday_thresh;

		if(metvars.mwetd <= 5.1) {
			pday_thresh = 0;
		}
		else if (metvars.mwetd <= 10.1) {
			pday_thresh = 1;
		}
		else if (metvars.mwetd <= 20.1) {
			pday_thresh = 2;
		}
		else {
			pday_thresh = 3;
		}
		
		do {
			int mwetd_sim    = 0;
			double mprec_sim = 0.0;
			double tmin_acc  = 0.;
	
			//Reload residuals and pday for start of month
			metvars.pday[0] = metvar_sav.pday[0];
			metvars.pday[1] = metvar_sav.pday[1];
			for (int i=0;i<4;i++) {
				metvars.resid[i] = metvar_sav.resid[i];
			}

			// Dayloop
			for (int day=0; day<ndaymon; day++) {
				
				metvars.dtmin   = mtmin_curr[day] ;
				metvars.dtmax   = mtmax_curr[day] ;
				metvars.dcldf   = mcloud_curr[day];
				metvars.dwind   = mwind_curr[day] ;

				// Now get day's weather
				weathergen_get_daily_met(metvars, rndst);
				
				dprec[day]= metvars.dprec;
				dtmin[day]= metvars.dtmin;
				dtmax[day]= metvars.dtmax;
				dcldf[day]= metvars.dcldf;
				dwind[day]= metvars.dwind;
				
				if ( metvars.dprec > 0. ) {
					mwetd_sim++; 
					mprec_sim += metvars.dprec;
				}
				tmin_acc += metvars.dtmin;
			} 
			// Break off criterium
			tmindiff = fabs(in_mtmin[mon] - tmin_acc / (double)ndaymon);
			
			// Reset met_out_save after initialization
			if (i_count == 0) {
				metvar_sav = metvars;
			}

			if (metvars.mprec <= 0.01 && tmindiff < E_TMIN) {
				pdaydiff = 0;
				precdiff = 0.;
				break;
			}
			// Enforce at least two times over the month to get initial values ok
			else if (i_count >= 1) {
				
				pdaydiff = (int)metvars.mwetd - mwetd_sim;
				precdiff = metvars.mprec - mprec_sim;

				// Breakoff-criteria for sufficient skill
				if ( (abs(pdaydiff) <= pday_thresh && fabs(precdiff) <= prec_t && tmindiff < E_TMIN) ||
				     (pdaydiff == 0 && fabs(precdiff) <= 1.5*prec_t  && tmindiff < E_TMIN))  {
					break;
				}

				double metric = abs(pdaydiff)*20/((double)pday_thresh + 1.0)  + fabs(precdiff) + 2.*tmindiff;
				
				// Save state if better w.r.t. metric 
				if ( metric < metric_sav ) {
					for ( int day=0; day<ndaymon; day++) {
						dprec_sav[day]= dprec[day];
						dtmin_sav[day]= dtmin[day];
						dtmax_sav[day]= dtmax[day];
						dcldf_sav[day]= dcldf[day];
						dwind_sav[day]= dwind[day];
					}
					metric_sav = metric;
					rndst.pday[0]  = metvars.pday[0];
					rndst.pday[1]  = metvars.pday[1];
					for (int i=0; i<4; i++) {
						rndst.resid[i] = metvars.resid[i];
					}
				}

				// After max amount of iterations is reached, use 
				// best set of data so far.
				if (i_count==MAXITER) {
					for (int day=0; day<ndaymon; day++) {
						dprec[day]= dprec_sav[day];
						dtmin[day]= dtmin_sav[day];
						dtmax[day]= dtmax_sav[day];
						dcldf[day]= dcldf_sav[day];
						dwind[day]= dwind_sav[day];
					}
					metvars.pday[0]= rndst.pday[0];
					metvars.pday[1]= rndst.pday[1];
					for (int i=0; i<4; i++) {
						metvars.resid[i] = rndst.resid[i];
					}
					break;
				}
			}
			i_count++;
		} while ( i_count <= MAXITER ); // 10000000 );

		// Write current settings for restart
		rndst.pday[0]  = metvars.pday[0];
		rndst.pday[1]  = metvars.pday[1];
		for (int i=0; i<4; i++) {
			rndst.resid[i] = metvars.resid[i];
		}

		// Enforce conservation by scaling with monthly averages

		// Correct Temperature biases by shifting
		double tmeancor= 0.;
		double dtrcor  = 0.;
		double preccor = 0.;
		double windcor = 0.;
		double cldfcor = 0.;
		double solcor  = 0.;
		double tot_cldwght = 0.;
		doy = accumday;
		for (int day=0; day<ndaymon;day++) {
			// Sometimes negative dtr can occur at values around 0. -> swap min,max
			if ( dtmin[day] > dtmax[day]) {
				double dummy = dtmin[day];
				dtmin[day]   = dtmax[day];
				dtmax[day]   = dummy;
			}

			tmeancor+= 0.5*(dtmax[day]+dtmin[day])/(double)ndaymon;
			dtrcor  += (dtmax[day]-dtmin[day])/(double)ndaymon;
			preccor += dprec[day];
			windcor += dwind[day];
			cldfcor += dcldf[day];
			doy++;

			// Compute days max rad (i.e. cldfr=0.) for weighting
			cldwght[day] = max(0.01,cldf2rad(0.0,lat,doy,true));
			tot_cldwght += cldwght[day];
		}

		tmeancor-= in_mtemp[mon];
		if (in_mdtr[mon] > 0.) {
			dtrcor  /= in_mdtr[mon];
		}
		else {
			dtrcor = 1.;
		}
		if (in_mprec[mon] > 0.){
			preccor /= in_mprec[mon];
		}
		else{
			preccor = 1.;
		}
		windcor /= (in_mwind[mon]*(double)ndaymon);

		tot_cldwght /= (double)ndaymon;

		for (int day=0; day<ndaymon;day++) {

			// Correct temp by shifting
			double ttmean = 0.5 * (dtmin[day]+dtmax[day]) - tmeancor;
			double tdtr   = (dtmax[day]-dtmin[day]) / dtrcor;
			dtmin[day] = ttmean - 0.5 * tdtr;
			dtmax[day] = ttmean + 0.5 * tdtr;
			
			// Correct wind by factor
			dwind[day] /= windcor;

			// Correct precip by factor 
			if ( preccor > 0. ) {
				dprec[day] /= preccor;
			} else {
				dprec[day] = 0.;
			}

			// Correct cldfr by factor
			cldwght[day] /= tot_cldwght;

			// Compute relative humidity. 
			// Use daylight avg temp following Running et al. 1987.
			double tdavg = 0.606*dtmax[day] + 0.394*dtmin[day];
			tdavg        = -1.14 + 1.12*tdavg;
			drhum[day]   = get_arden_rh(tdavg,dtmin[day]);
		}

		// Redistribute limited parameters like relhum and cloud-fraction		
		double limit[2] = {0.,1.};
		if ( in_msol[mon] > 0. ) {
	
			if (in_mcldf[mon] > 0.)
				redist_restricted_vals(dcldf, ndaymon, in_mcldf[mon], limit, cldwght);
			
			solcor = 0.;
			for (int day=0;day<ndaymon;day++) {
				doy = accumday+day+1;
				dsol[day] = max(0.001,cldf2rad(dcldf[day],lat,doy,true));
				solcor += dsol[day];
			}
			solcor /= (in_msol[mon]*(double)ndaymon);
			for (int day=0;day<ndaymon;day++) 
				dsol[day] /= solcor;
		}
		
		if ( in_mrhum[mon] > 0. ) {
			redist_restricted_vals(drhum, ndaymon, in_mrhum[mon], limit, dum);
		}

		// Compute solar radiation from cloud-fraction
		chk_dtemp = 0.; 
		chk_ddtr  = 0.; 
		chk_dprec = 0.; 
		chk_dsol  = 0.; 
		chk_dwind = 0.; 
		chk_drhum = 0.;
		for (int day=0; day<ndaymon;day++) {
			out_dtemp[day+accumday] =(dtmax[day] + dtmin[day]) / 2.;
			out_ddtr [day+accumday] = dtmax[day] - dtmin[day];
			out_dprec[day+accumday] = dprec[day];
			out_dsol [day+accumday] = dsol [day];
			out_dwind[day+accumday] = dwind[day];
			out_drhum[day+accumday] = drhum[day];

			chk_dtemp += out_dtemp[day+accumday]/(double)ndaymon; 
			chk_ddtr  += out_ddtr [day+accumday]/(double)ndaymon; 
			chk_dprec += out_dprec[day+accumday]; 
			chk_dsol  += out_dsol [day+accumday]/(double)ndaymon; 
			chk_dwind += out_dwind[day+accumday]/(double)ndaymon; 
			chk_drhum += out_drhum[day+accumday]/(double)ndaymon;
		}

		// If GWGen doesn't find a day for precipitation add it at the first third
		// of the month.
		if (metvars.mwetd > 0 && fabs(chk_dprec - in_mprec[mon]) > in_mprec[mon]-0.01) {
			out_dprec[accumday+9] =  in_mprec[mon];
		}

	} // month loop
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Buck, A. L., New equations for computing vapor pressure and enhancement factor, J. Appl.
//   Meteorol., 20, 1527-1532, 1981
// Prentice IC, Sykes MT, Cramer W, 1993. A simulation model of the transient effects of
//   climate change on forest landscapes. Ecological Modelling, 65, 51-70.
// Sommer, P. S. and Kaplan, J. O.: A globally calibrated scheme for generating daily meteorology
//   from monthly statistics: Global-WGEN (GWGEN) v1.0, Geosci. Model Dev., 10, 3771-3791,
//   doi:10.5194/gmd-10-3771-2017, 2017.
// Original Code available in Fortran at: https://arve-research.github.io/gwgen/
