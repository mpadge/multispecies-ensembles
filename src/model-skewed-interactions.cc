/*
 * model_skewed_interactions.cc
 * 
 * Adapted from model.cc to compare population trajectories with normally-distributed
 * interactions to those with squared values that increase the proportion of lower
 * values. Correlations are calculated between population trajectories only; the
 * global noise signal is irrelevant in this case. Differences between the two
 * increase with ensemble size (of course), so a fixed size is used throughout, and
 * the effects of different values of alpha are examined.
 * 
 * The squared values have to maintain the same distributional variance. Squaring an
 * N(0,1) distribution and dividing values by 1.75 produces once again a distribution 
 * with a variance of 1, so this is what is implemented below. The whole things
 * makes bugger all difference.
 * 
 * A plotting routine is included at the bottom. It loads up the two files,
 * "results_skewed_interactions50sp.txt" and "...100sp.txt". These were both
 * generated with only 10 trials in each case. The value for alpha=0 and 100 species
 * when averaged over 100 trials is 0.78 +/1 0.04.
 * 
 * The point of this is that with no shared environment, then different distributions
 * of interaction strengths make quite a bit of difference, but any appreciable
 * degree of environmental overlap reduces these differences very profoundly. That
 * in itself is an important result.
 * 
 * One further issues arises, however, in that simply squaring the interaction strengths
 * changes the entire distribution. If it is subsequently standardised to the same
 * integral as before, by multiplying all transformed values by the sum of absolute
 * strengths prior to transformation, then values are increased by quite a bit. This
 * then makes a considerable difference. Even more importantly, once more than about 70
 * species are included, correlations start to become strongly negative, and for 100 
 * species end up perfectly inverted from what they "should" be. This is very strange,
 * but it is genuine --- and maybe worthy of further investigation. For the moment,
 * however, it's probably fairly safe to conclude that this happens because with such
 * strong competition in large ensembles, the effects of predation become so fierce
 * that populations effectively become perfectly desynchronised, just as the correlation
 * lines in the manuscript start decreasing in large ensembles if predation becomes too 
 * strong. So this may be neglected, meaning that rescaling only works for smaller
 * ensembles up to 50-60 species.
 * 
 * The plotting routine at the bottom nevertheless includes all results just to
 * demonstrate this effect.
 * 
 * NOTE however ... after all that ... that the "proper" way to adjust the distribution
 * is actually to ensure that the variance remains constant, because that is the
 * assumption of the analyses. This is
 */

#include <stdlib.h> // has abs function
#include <math.h>
#include <iostream>
#include <fstream> // file in & out
//#include <string.h>
#include <stdio.h>
#include <dirent.h> // For reading directories
#include <time.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "boost/multi_array.hpp"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
//using namespace std;

#define PI 3.1415926535897932384626433832795

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;
typedef boost::numeric::ublas::vector<int> ivec;
typedef boost::numeric::ublas::vector<double> dvec;
typedef boost::numeric::ublas::matrix<double> dmat;
typedef boost::numeric::ublas::matrix<double> darr;

const int inix = -9999, len_t = 1000;
const double dnix = -9999.0, bmn = -0.001, bsd = 0.05, sigma_a = 0.5, sigma_alpha = 0.5;
/*
 * NOTE that setting len_t = 100 produces notably noisier results, so it
 * needs to be greater, hence the value of 1,000, even through computation
 * then takes some time.
 * 
 * Also NOTE that the sigma values are actually SDs, not variances!
 */

struct SpeciesPars{
	double ar [4], alpha;
};
/* SpeciesPars is a vector of (nSpecies) with:
 * 	a		= AR coefficients around means of [1,0,0], with
 * 				SDs of sigma_a * [1, 1/2. 1/4]
 * 	alpha	= degree of environment sharing (mean = alpha_0)
 */
typedef boost::numeric::ublas::vector<SpeciesPars> species_vec;
struct SpeciesData{
	species_vec spvec;
	dmat compMat;
};

struct RegrResults {
	double r2, slope, intercept, SS, tval;	};

void makeCommunity (int nSpecies, double sync, SpeciesData & speciesData, base_generator_type & generator);
void makeNoise (int nSpecies, darr & noise, base_generator_type & generator);
void runPop (int nSpecies, darr noise, SpeciesData speciesData, dvec & p_vec);
RegrResults regression(dvec x, dvec y);

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[])
{
	int count;
	double tempd, sync, sumr, sumr2;
	SpeciesData speciesData;
	RegrResults regr;
	std::ifstream in_file;
	std::string linetxt, fname;
	std::stringstream ss_code;
	std::ofstream out_file;
	clock_t timer[2];
	time_t seed;
	base_generator_type generator(42u);

	//std::cout.setf (std::ios::fixed, std::ios::floatfield);   // floatfield set to fixed
	//std::cout.precision(2);
	int nTrials = atoi (argv [2]);
	if (nTrials < 1) { nTrials = 10;	}
	int nSpecies = atoi (argv [3]);
	if (nSpecies < 2) { nSpecies = 50;	}
	std::cout << "nTrials = " << nTrials << "; nSpecies = " << nSpecies << std::endl;

	time (&seed);
	generator.seed(static_cast<unsigned int>(seed));
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif(generator, uni_dist);
	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);
	// Burn both generators in
	for (int i=0; i<20; i++) {
		tempd = runif();
		tempd = rnorm();
	}

	// Set up data to fixed size of nSpecies_max
	boost::numeric::ublas::vector<SpeciesPars> tempvec (nSpecies);
	boost::numeric::ublas::matrix<double> tempdmat (nSpecies, nSpecies);
	speciesData.spvec = tempvec;
	speciesData.compMat = tempdmat;
	boost::numeric::ublas::matrix<double> noise (len_t, nSpecies + 1);
	boost::numeric::ublas::vector<double> pop_t1 (len_t); // normal distribution
	boost::numeric::ublas::vector<double> pop_t2 (len_t); // squared values

	out_file.open ("aaajunk.txt");
	out_file << "alpha,\tr2mn,\tr2sd" << std::endl;
	for (int s=0; s<50; s++) {
		sync = (double) s / 100.0;
		sumr = 0.0;
		sumr2 = 0.0;
		count = 0;
		for (int i=0; i<nTrials; i++) {
			makeCommunity (nSpecies, sync, speciesData, generator);
			makeNoise (nSpecies, noise, generator);
			runPop (nSpecies, noise, speciesData, pop_t1);
			// Then square the species interactions to skew the distribution. Note
			// that adjusting the distribution to maintain constant variance
			// requires re-scaling the rnorm() values to (0,1), dividing by
			// the value of 1.75 --- i just got this from R, not by analysing what
			// it should really be! --- and then re-scaling again by bsd.
			for (int j=0; j<(nSpecies - 1); j++) {
				for (int k=(j + 1); k<nSpecies; k++) {
					tempd = (speciesData.compMat (j, k) - bmn) / bsd;
					if (tempd > 0.0) {
						speciesData.compMat (j, k) = bmn + bsd * tempd * tempd / 1.75;
						speciesData.compMat (k, j) = bmn - bsd * tempd * tempd / 1.75;
					}
					else {
						speciesData.compMat (j, k) = bmn - bsd * tempd * tempd / 1.75;
						speciesData.compMat (k, j) = bmn + bsd * tempd * tempd / 1.75;
					}
			}	} // end for k & j
			runPop (nSpecies, noise, speciesData, pop_t2);
			regr = regression (pop_t1, pop_t2);
			// For some reason, correlations are sometimes very close to zero (dnix's?),
			// so this if just removes those ones.
			if (regr.r2 > 0.5) {
				count++;
				sumr += regr.r2;
				sumr2 += regr.r2 * regr.r2;
			}
		}
		sumr = sumr / (double) count;
		sumr2 = sumr2 / (double) count - sumr * sumr;
		sumr2 = sqrt (sumr2);
		out_file << sync << ",\t" << sumr << ",\t" << sumr2 << std::endl;
		std::cout << s << "."; std::cout.flush();
	} // end for s
	out_file.close();
	std::cout << std::endl;

	return 0;
}; // end main

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SUBFUNCTIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


void makeCommunity (int nSpecies, double sync, SpeciesData & speciesData, base_generator_type & generator)
{
	double tempd;

	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);

	for (int i=0; i<nSpecies; i++) {
		speciesData.spvec (i).ar [0] = 1.0 + sigma_a * rnorm ();
		speciesData.spvec (i).ar [1] = sigma_a * rnorm () / 2.0;
		speciesData.spvec (i).ar [2] = sigma_a * rnorm () / 4.0;
		speciesData.spvec (i).ar [3] = sigma_a * rnorm () / 8.0;
		tempd = 0.0;
		for (int j=0; j<4; j++) { tempd += fabs (speciesData.spvec (i).ar [j]);	}
		for (int j=0; j<4; j++) {
			speciesData.spvec (i).ar [j] = speciesData.spvec (i).ar [j] / tempd;
		}

		/* NOTE the importance of the following condition. The upper limit has to be
		 * the opposite of the lower limit, so the distribution stays at all times
		 * symmetrical. Without this, the distribution becomes biased with increasing
		 * sigma_alpha, producing greater *mean* values, and skewing everything.*/
		if (sync == 0.0) { speciesData.spvec (i).alpha = 0.0;	}
		else {
			speciesData.spvec (i).alpha = sync + sigma_alpha * rnorm ();
			while (speciesData.spvec (i).alpha <= 0.0 || speciesData.spvec (i).alpha >= (2.0 * sync)) {
				speciesData.spvec (i).alpha = sync + sigma_alpha * rnorm ();	}
		}
	}

	for (int i=0; i<nSpecies; i++) {
		speciesData.compMat (i, i) = 0.0;	}
	for (int i=0; i<(nSpecies - 1); i++) {
		for (int j=(i + 1); j<nSpecies; j++) {
			tempd = rnorm ();
			speciesData.compMat (i, j) = bmn + bsd * tempd;
			speciesData.compMat (j, i) = bmn - bsd * tempd;
	}	}
} // end SpeciesData makeCommunity


void makeNoise (int nSpecies, darr & noise, base_generator_type & generator)
{
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif(generator, uni_dist);
	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);

	for (int t=0; t<len_t; t++) {
		for (int i=0; i<=nSpecies; i++) {
			noise (t, i) = rnorm();
		} // end for j over nSpecies
	} // end for t over len_t
} // end darr makeNoise

void runPop (int nSpecies, darr noise, SpeciesData speciesData, dvec & pop_t)
{
	double tempd;

	boost::numeric::ublas::vector<double> spvec0 (nSpecies), spvec1 (nSpecies), spvec2 (nSpecies), spvec3 (nSpecies), tempvec (nSpecies);

	// Add shared signal to species noise
	for (int t=0; t<len_t; t++) {
		pop_t (t) = 0.0;
		for (int i=0; i<nSpecies; i++) {
			noise (t, i) = speciesData.spvec (i).alpha * noise (t, nSpecies) +
				(1.0 - speciesData.spvec (i).alpha) * noise (t, i);
		}
	}

	// Run in first three AR steps
	for (int i=0; i<nSpecies; i++) {
		spvec0 (i) = noise (0, i);
		for (int j=0; j<nSpecies; j++) {
			spvec0 (i) += noise (0, j) * speciesData.compMat (j, i);
		}
		pop_t (0) += spvec0 (i);
		spvec1 (i) = speciesData.spvec (i).ar [0] * noise (1, i) + speciesData.spvec (i).ar [1] * spvec0 (i);
	}
	tempvec = spvec1;
	for (int i=0; i<nSpecies; i++) {
		for (int j=0; j<nSpecies; j++) {
			tempvec (i) += spvec1 (j) * speciesData.compMat (j, i);
	}	}
	spvec1 = tempvec;
	for (int i=0; i<nSpecies; i++) {
		pop_t (1) += spvec1 (i);
		spvec2 (i) = speciesData.spvec (i).ar [0] * noise (2, i) + speciesData.spvec (i).ar [1] * spvec1 (i) +
			speciesData.spvec (i).ar [2] * spvec0 (i);
	}
	tempvec = spvec2;
	for (int i=0; i<nSpecies; i++) {
		for (int j=0; j<nSpecies; j++) {
			tempvec (i) += spvec2 (j) * speciesData.compMat (j, i);
	}	}
	spvec2 = tempvec;
	for (int i=0; i<nSpecies; i++) { pop_t (2) += spvec2 (i);	}

	// Then the full run
	for (int t=3; t<len_t; t++) {
		for (int i=0; i<nSpecies; i++) {
			spvec3 (i) = speciesData.spvec (i).ar [0] * noise (t, i) + speciesData.spvec (i).ar [1] * spvec2 (i) +
				speciesData.spvec (i).ar [2] * spvec1 (i) + speciesData.spvec (i).ar [3] * spvec0 (i);
		}
		tempvec = spvec3;
		for (int i=0; i<nSpecies; i++) {
			for (int j=0; j<nSpecies; j++) {
				tempvec (i) += spvec3 (j) * speciesData.compMat (j, i);
		}	}
		for (int i=0; i<nSpecies; i++) { pop_t (t) += tempvec (i);	}
		spvec0 = spvec1;
		spvec1 = spvec2;
		spvec2 = tempvec;
	} // end for i over len_t
} // end dev runPop



RegrResults regression (dvec x, dvec y)
{
	double sx, sx2, sy, sy2, sxy, t1, t2, xmn, ymn;
	RegrResults regr_results;

	sx = 0.0; sx2 = 0.0;
	sy = 0.0; sy2 = 0.0; sxy = 0.0;
	int count = 0, n = x.size();
	for (int i=0; i<n; i++) {
		if (x(i) > dnix && y(i) > dnix) {
			count++;
			sx += x(i);
			sx2 += x(i) * x(i);
			sy += y(i);
			sy2 += y(i) * y(i);
			sxy += x(i) * y(i);	}	}
	xmn = sx / (double) count;
	ymn = sy / (double) count;
	if (count > 0) {
		t1 = (sxy - sx * sy / (double) count);
		t2 = (sx2 - sx * sx / (double) count) * (sy2 - sy * sy / (double) count);
		regr_results.r2 = t1 / sqrt(t2); // the R-value
		regr_results.slope = t1 / (sx2 - sx * sx / (double) count); // Slope
		regr_results.intercept = sy / (double) count - regr_results.slope * sx / (double) count; // Intercept

		regr_results.r2 = regr_results.r2 * regr_results.r2;
		if (regr_results.slope < 0) { regr_results.r2 = -regr_results.r2;	}

		// Then calculate SS and tval
		sy2 = 0.0; sx2 = 0.0; count = 0;
		regr_results.SS = 0.0;
		for (int i=0; i<n; i++) {
			if (x(i) > dnix && y(i) > dnix) {
				count++;
				t1 = regr_results.slope * x(i) + regr_results.intercept;
				regr_results.SS += (y(i) - t1) * (y(i) - t1);
				sx2 += (x(i) - xmn) * (x(i) - xmn);
				sy2 += (y(i) - ymn) * (y(i) - ymn);	}
			} // end for i
		if (count > 0) { // tval calculation
			regr_results.SS = regr_results.SS / (double) count;
			regr_results.tval = sqrt(sy2 / ((double) count - 2)) / sqrt(sx2);
			regr_results.tval = regr_results.slope / regr_results.tval;
		}
		else {
			regr_results.SS = dnix; regr_results.tval = dnix;	}
	} // end if count > 0
	else {
		regr_results.r2 = dnix;
		regr_results.slope = dnix;
		regr_results.intercept = dnix;
		regr_results.SS = dnix;
		regr_results.tval = dnix;
	}

	return regr_results;
} // end function regression


/* R script to plot results
junk <- function()
{
	ylims=c(0.95, 1)

	setwd("/data/Documents/analyses/birds & climate/multispecies model/")
	dat <- read.csv ("results_skewed_interactions50sp.txt", header=TRUE)
	indx <- which (dat$alpha != 0.29)
	dat <- dat[indx,]
	plot (dat$alpha, dat$r2mn, "l", col="lawngreen", lwd=2, ylim=ylims, xlab="alpha", ylab="R2")
	lines (dat$alpha, dat$r2mn - dat$r2sd, col="lawngreen", lwd=2, lty=2)
	lines (dat$alpha, dat$r2mn + dat$r2sd, col="lawngreen", lwd=2, lty=2)

	#dat <- read.csv ("results_skewed_interactions100sp.txt", header=TRUE)
	dat <- read.csv ("aaajunk.txt", header=TRUE)
	#indx <- which (dat$alpha != 0.48)
	#dat <- dat [indx, ]
	lines (dat$alpha, dat$r2mn, col="red", lwd=2)
	lines (dat$alpha, dat$r2mn - dat$r2sd, col="red", lwd=2, lty=2)
	lines (dat$alpha, dat$r2mn + dat$r2sd, col="red", lwd=2, lty=2)

	ypos <- ylims[1] + 0.75 * diff (ylims)
	legend (0.5 * max (dat$alpha), ypos, lwd=2, bty="n",
		col=c("lawngreen", "red"), legend=c("50 species", "100 species"))
}

# And demonstration that scaling by 1.75 converts squared normal distribution
# to the same variance as the untransformed version
x1 <- rnorm (1e5, 0, 1)
h1 <- hist(x1,plot=FALSE)
h1$counts <- h1$counts / sum (h1$counts)

x2 <- sign (x1) * x1 ^ 2 / 1.75
h2 <- hist(x2,plot=FALSE)
h2$counts <- h2$counts / sum (h2$counts)

xlims <- range (c (h1$mids, h2$mids))
ylims <- range (c (h1$counts, h2$counts))
plot (h1$mids, h1$counts, "l", col="red", xlim=xlims, ylim=ylims)
lines (h2$mids, h2$counts, col="blue")
legend (min(xlims), max(ylims), lwd=1, col=c("red","blue"), bty="n",
	legend=c("normal","norm ^ 2"))

v1 <- formatC( var (x1), format="f", digits=4)
v2 <- formatC( var (x2), format="f", digits=4)
title (main=paste("Variances = (", v1, ", ", v2, ")", sep=""))
*/