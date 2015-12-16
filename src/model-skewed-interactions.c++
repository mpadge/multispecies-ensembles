/*
 * model_skewed_interactions.cc
 * 
 * Adapted from model.cc to compare population trajectories with
 * normally-distributed interactions to those with squared values that increase
 * the proportion of lower values. Correlations are calculated between
 * population trajectories only; the global noise signal is irrelevant in this
 * case. Differences between the two increase with ensemble size (of course), so
 * a fixed size is used throughout, and the effects of different values of alpha
 * are examined.
 * 
 * The squared values have to maintain the same distributional variance.
 * Squaring an N(0,1) distribution and dividing values by 1.75 produces once
 * again a distribution with a variance of 1, so this is what is implemented
 * below. The whole thing makes bugger all difference.
 * 
 * A plotting routine is included at the bottom. It loads up the two files,
 * "results_skewed_interactions50sp.txt" and "...100sp.txt". These were both
 * generated with only 10 trials in each case. The value for alpha=0 and 100
 * species when averaged over 100 trials is 0.78 +/1 0.04.
 * 
 * The point of this is that with no shared environment, then different
 * distributions of interaction strengths make quite a bit of difference, but
 * any appreciable degree of environmental overlap reduces these differences
 * very profoundly. That in itself is an important result.
 * 
 * One further issues arises, however, in that simply squaring the interaction
 * strengths changes the entire distribution. If it is subsequently standardised
 * to the same integral as before, by multiplying all transformed values by the
 * sum of absolute strengths prior to transformation, then values are increased
 * by quite a bit. This then makes a considerable difference. Even more
 * importantly, once more than about 70 species are included, correlations start
 * to become strongly negative, and for 100 species end up perfectly inverted
 * from what they "should" be. This is very strange, but it is genuine --- and
 * maybe worthy of further investigation. For the moment, however, it's probably
 * fairly safe to conclude that this happens because with such strong
 * competition in large ensembles, the effects of predation become so fierce
 * that populations effectively become perfectly desynchronised, just as the
 * correlation lines in the manuscript start decreasing in large ensembles if
 * predation becomes too strong. So this may be neglected, meaning that
 * rescaling only works for smaller ensembles up to 50-60 species.
 * 
 * The plotting routine at the bottom nevertheless includes all results just to
 * demonstrate this effect.
 * 
 * NOTE however ... after all that ... that the "proper" way to adjust the
 * distribution is actually to ensure that the variance remains constant,
 * because that is the assumption of the analyses. 
 */

#include "model-skewed-interactions.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[])
{
	int nSpecies, nTrials, count;
	double sync, sumr, sumr2, tempd;
	SpeciesData speciesData;
	RegrResults regr;
	std::ifstream in_file;
	std::string linetxt, fname;
	std::stringstream ss_code;
	std::ofstream out_file;
	clock_t timer[2];
	time_t seed;
	base_generator_type generator(42u);

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "version 1.0.1")
            ("help", "this is not helpful")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("nTrials,n", boost::program_options::value <int> 
                (&nTrials)->default_value (10), "nTrials")
            ("nSpecies,s", boost::program_options::value <double>
                (&sync)->default_value (20), "nSpecies")
            ;

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        boost::program_options::options_description visible("Allowed options");
        visible.add(generic).add(config);

        boost::program_options::variables_map vm;
        store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).run(), vm);

        notify(vm);

        if (vm.count("help")) {
            std::cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "model_dependent_xi, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    

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
	dvec noise_global (len_t);
	dvec pop_t (len_t);
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
			makeCommunity (nSpecies, sync, speciesData, generator, 0);
			makeNoise (nSpecies, noise, generator);
			runPop (nSpecies, noise, speciesData, pop_t1);
            // Then square the species interactions to skew the distribution.
            // Note that adjusting the distribution to maintain constant
            // variance requires re-scaling the rnorm() values to (0,1),
            // dividing by the value of 1.75 --- i just got this from R, not by
            // analysing what it should really be! --- and then re-scaling again
            // by bsd.
			for (int j=0; j<(nSpecies - 1); j++) 
				for (int k=(j + 1); k<nSpecies; k++) 
                {
					tempd = (speciesData.compMat (j, k) - bmn) / bsd;
					if (tempd > 0.0) 
                    {
						speciesData.compMat (j, k) = bmn + 
                                bsd * tempd * tempd / 1.75;
						speciesData.compMat (k, j) = bmn - 
                                bsd * tempd * tempd / 1.75;
                    } else 
                    {
						speciesData.compMat (j, k) = bmn - 
                                bsd * tempd * tempd / 1.75;
						speciesData.compMat (k, j) = bmn + 
                                bsd * tempd * tempd / 1.75;
					}
                }
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
