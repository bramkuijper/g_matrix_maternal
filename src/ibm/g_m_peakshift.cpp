//
// evolving G matrix in environments characterized by peak-shifts


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>

// random number generation
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// various functions, such as unique filename creation
#include "bramauxiliary.h"

//#define NDEBUG
//
// the compilation sign should only be turned on when one wants
// to assess the complete distribution of phenotypes

using namespace std;


///////////////////     PARAMS      ///////////////////

// number of generations
const int NumGen = 50000;

// population size
const int Npop = 1366; 

// number of generations to skip when outputting data
const int skip = 20;

// number of loci
const int n_loci_g = 50;

// mutational variance trait 1, 2
double a1 = 0;
double a2 = 0;

// mutational correlation
double rmu = 0;

// mutation rate conventional gene loci
double mu = 0;

// mutation rate maternal effects loci
double mu_m = 0;
double sdmu_m = 0;

// initial number of generations without change 
int burnin = 5000;

// strengths of selection
double omega[2][2] = {{0,0},{0,0}};

// strength of correlated selection
double r_omega = 0;

// optima 
double theta1 = 0;
double theta2 = 0;

// deterministic change in optima
double delta_t1 = 0;
double delta_t2 = 0;

// interval in which no change occurs
int interval = 0;

// variance to mimick brownian motion
double sigma_theta1 = 0;
double sigma_theta2 = 0;

// fecundity per individual
double B = 0;

///////////////////     STATS       ///////////////////

// track number of individuals 
size_t Nm = 0;
size_t Nf = 0;
size_t NKids = 0;

// printing the covariance between z1 and z2 for offspring
// to have a gist of what is wrong with this model
double meancov = 0;

// mean fitness
double meanw = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
int generation = 0;

// random seed
unsigned seed = 0;

// gnu scientific library random number generator initialization
// http://www.gnu.org/software/gsl/ 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    // the trait alleles are given by z_{ijk}
    // 0<i<=n_loci_g: the number of unlinked gene loci 
    // j: the number of independent traits (2 in this case)
    // k: ploidy (in this case diploid)
    double g[2][n_loci_g][2]; 

    double m[2][2][2];

    double phen[2];
    double gen[2];

};

// allocate a population and a population of survivors
typedef Individual Population[Npop];
typedef Individual NewPopulation[Npop*2*10];
Population Males, Females;
NewPopulation NewPop;

// generate a unique filename for the output file
string filename("sim_Gmatrix");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

#ifdef DISTRIBUTION
// generate a filename for the phenotype distribution file
string filename_new2(create_filename("sim_evolving_m_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION


// bivariate Gaussian fitness function as in Jones et al 2003
double v(double const z1, double const z2)
{
    return(
         exp(-.5 * ( 
             (z1 - theta1) * (z1 - theta1) * omega[1][1] 
             - 2 * omega[0][1] * (z1 - theta1) * (z2 - theta2)
             + (z2 - theta2) * (z2 - theta2) * omega[0][0] 
            )/ (omega[0][0] * omega[1][1] - omega[1][0] * omega[1][0])
         ));
}


// initialize simulations from command line arguments
void initArguments(int argc, char *argv[])
{
    a1 = atof(argv[1]);
    a2 = atof(argv[2]);
    rmu = atof(argv[3]);
    mu = atof(argv[4]);
    mu_m = atof(argv[5]);
    sdmu_m = atof(argv[6]);
    omega[0][0] = atof(argv[7]);
    omega[1][1] = atof(argv[8]);
    r_omega = atof(argv[9]);
    B = atof(argv[10]);
    delta_t1 = atof(argv[11]);
    delta_t2 = atof(argv[12]);
    interval = atoi(argv[13]);
    sigma_theta1 = atof(argv[14]);
    sigma_theta2 = atof(argv[15]);

    omega[1][0] = omega[0][1] = r_omega * sqrt(omega[0][0] * omega[1][1]);

    assert(B < 10);
}


// mutation according to a continuum of alleles model
// see description in Jones et al 2012 J Evol Biol 
//
// "During the progeny-production phase, we choose
// gametes at random to be affected by mutation. We
// assume a uniform per-locus mutation rate of l. We
// draw a pseudorandom number between 0 and 1 from a
// uniform distribution and assume that a gamete carries a
// new mutation if the random number is less than nloci * mu 
void MutateG(Individual &ind)
{
    // random deviates from bivariate gaussian dist
    double a = 0;
    double b = 0;

    for (size_t locus_i = 0; locus_i < n_loci_g; ++locus_i)
    {
        // mutation in each locus occurs with probability mu
        if (gsl_rng_uniform(r) < mu)
        {
            // generate new allelic increments a,b by drawing them
            // from a bivariate gaussian distribution with std deviations a1, a2
            // and 
            gsl_ran_bivariate_gaussian(r, a1, a2, rmu, &a, &b);

            //cout << a << ";" << b << endl;

            // add the new allelic effect to trait one
            ind.g[0][locus_i][0] += a;
            // add the new allelic effect to trait two
            ind.g[1][locus_i][0] += b;

            ind.gen[0] += a;
            ind.gen[1] += b;
        }

        // the other genome copy (diploidy)
        if (gsl_rng_uniform(r) < mu)
        {
            // generate new allelic increments a,b by drawing them
            // from a bivariate gaussian distribution with std deviations a1, a2
            // and 
            gsl_ran_bivariate_gaussian(r, a1, a2, rmu, &a, &b);

            // add the new allelic effect to trait one
            ind.g[0][locus_i][1] += a;
            // add the new allelic effect to trait two
            ind.g[1][locus_i][1] += b;

            ind.gen[0] += a;
            ind.gen[1] += b;
        }
    }
}

// mutate the maternal effects matrix
double MutateM(double m)
{
    if (gsl_rng_uniform(r) < mu_m)
    {
        m += gsl_ran_gaussian(r, sdmu_m);
    }

    return(m);
}

// write the parameters (typically at the end of the output file)
void WriteParameters()
{
	DataFile << endl
		<< endl 
        << "npop;" << Npop << endl
        << "nloci_g;" << n_loci_g << endl
        << "a1;" << a1 << endl 
        << "a2;" << a2 << endl
        << "rmu;" << rmu << endl
        << "seed;" << seed << endl
        << "r_omega;" << r_omega << endl
        << "mu;" << mu << endl
        << "mu_m;" << mu_m << endl
        << "sdmu_m;" << sdmu_m << endl
        << "omega_11;" << omega[0][0] << endl
        << "omega_12;" << omega[0][1] << endl
        << "omega_21;" << omega[1][0] << endl
        << "omega_22;" << omega[1][1] << endl
        << "delta_t1;" << delta_t1 << endl
        << "delta_t2;" << delta_t2 << endl
        << "interval;" << interval << endl
        << "sigma_theta1;" << sigma_theta1 << endl
        << "sigma_theta2;" << sigma_theta2 << endl
        << "B;" << B << endl;
}

// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void Init()
{
    // get the timestamp (with nanosecs)
    // to initialize the seed
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);


	// initialize the whole population
	for (size_t i = 0; i < Npop/2; ++i)
	{
        // loop through the different characters
        for (size_t trait_i = 0; trait_i < 2; ++trait_i)
        {
            Males[i].phen[trait_i] = 0;
            Males[i].gen[trait_i] = 0;
            Females[i].phen[trait_i] = 0;
            Females[i].gen[trait_i] = 0;
                        
            // loop through each of the loci
            for (size_t j = 0; j < n_loci_g; ++j)
            {
                // loop through each of the alleles
                for (size_t k = 0; k < 2; ++k)
                {
                    Males[i].g[trait_i][j][k] = 0;
                    Females[i].g[trait_i][j][k] = 0;
                }
            }
        }
	}

    Nf = Npop/2;
    Nm = Npop/2;
}

// create an offspring
void Create_Kid(size_t const mother, size_t const father, Individual &kid)
{
    // copy mother and father, preventing many array lookups
    Individual Mother = Females[mother];
    Individual Father = Males[father];

    // reset the total genetic values corresponding to each trait to 0
    kid.gen[0] = 0;
    kid.gen[1] = 0;

    // inherit 'normal' gene loci
    // loop through all the loci
    for (size_t i = 0; i < n_loci_g; ++i)
    {
        if (gsl_rng_uniform(r) < 0.5)
        {
            kid.g[0][i][0] = Mother.g[0][i][0];
            kid.gen[0] += kid.g[0][i][0];
            
            kid.g[1][i][0] = Mother.g[1][i][0];
            kid.gen[1] += kid.g[1][i][0];
        }
        else
        {
            kid.g[0][i][0] = Mother.g[0][i][1];
            kid.gen[0] += kid.g[0][i][0];
            
            kid.g[1][i][0] = Mother.g[1][i][1];
            kid.gen[1] += kid.g[1][i][0];
        }
        
        
        if (gsl_rng_uniform(r) < 0.5)
        {
            kid.g[0][i][1] = Father.g[0][i][0];
            kid.gen[0] += kid.g[0][i][1];
            
            kid.g[1][i][1] = Father.g[1][i][0];
            kid.gen[1] += kid.g[1][i][1];
        }
        else
        {
            kid.g[0][i][1] = Father.g[0][i][1];
            kid.gen[0] += kid.g[0][i][1];
            
            kid.g[1][i][1] = Father.g[1][i][1];
            kid.gen[1] += kid.g[1][i][1];
        }
    }

    // inherit maternal effect loci
    for (size_t trait_i = 0; trait_i < 2; ++trait_i)
    {
        for (size_t trait_j = 0; trait_j < 2; ++trait_j)
        {
            kid.m[trait_i][trait_j][0] = MutateM(Mother.m[trait_i][trait_j][gsl_rng_uniform_int(r,2)]);
            kid.m[trait_i][trait_j][1] = MutateM(Father.m[trait_i][trait_j][gsl_rng_uniform_int(r,2)]);
        }
    }

    MutateG(kid);
   
    // add environmental variance to each trait by adding a random number
    // drawn from a normal distribution to each phenotype
    kid.phen[0] = kid.gen[0] + gsl_ran_gaussian(r,1.0) + (kid.m[0][0][0] + kid.m[0][0][1]) * Mother.phen[0] + (kid.m[0][1][0] + kid.m[0][1][1]) * Mother.phen[1];
    kid.phen[1] = kid.gen[1] + gsl_ran_gaussian(r,1.0) + (kid.m[1][0][0] + kid.m[1][0][1]) * Mother.phen[0] + (kid.m[1][1][0] + kid.m[1][1][1]) * Mother.phen[1];
}


// Survival of juveniles to reproductive adults
void Reproduce_Survive()
{
    // set kids counter to 0 prior to reproduction
    NKids = 0;

    // stats for average fitness
    meanw = 0;

    double w;

    // stats for genetic covariance within offspring
    meancov = 0;

    for (size_t i = 0; i < Nf; ++i)
    {
        // random mating
        size_t father = gsl_rng_uniform_int(r, Nm);

        // produce kids and let them survive
        for (size_t j = 0; j < 2 * B; ++j)
        {
            Individual Kid;

            // create a kid from maternal and paternal genes
            Create_Kid(i, father, Kid);

            // calculate survival
            w = v(Kid.phen[0], Kid.phen[1]);

            //cout << i << " " << father << " " << Kid.phen[0] << " " << Kid.phen[1] << " " << w << endl;


            assert(w >= 0 && w <= 1.0);

            meanw += w;

            // individual survives; add to stack
            if (gsl_rng_uniform(r) < w)
            {
                NewPop[NKids++] = Kid;
                assert(NKids < Npop * 2 * 10);
            }
        }

        // remove dad
        Males[father] = Males[--Nm];

        if (Nm == 0)
        {
            break;
        }
    }

    //cout << meancov / (Npop * 2 * B) << endl;

    meanw /= Nf * 2 * B;

    //cout << NKids << endl;
    
    if (NKids < 100)
    {
        cout << "extinct " << NKids << endl;
        WriteParameters();
        exit(1);
    }

    size_t random_kid;

    Nm = 0;
    Nf = 0;

    // sample new generation from kids
    for (size_t i = 0; i < Npop; ++i)
    {
        assert(NKids >= 1);

        // get a randomly sampled kid
        random_kid = gsl_rng_uniform_int(r, NKids);

        if (gsl_rng_uniform(r) < 0.5)
        {
            Males[Nm++] = NewPop[random_kid];
        }
        else
        {
            Females[Nf++] = NewPop[random_kid];
        }
    }

    // change the environment
    if (generation > burnin && generation % interval == 0)
    {
        theta1 += delta_t1 + gsl_ran_gaussian(r, sigma_theta1);
        theta2 += delta_t2 + gsl_ran_gaussian(r, sigma_theta2);
    }
}


// write down summary statistics
void WriteData()
{
    // genotypic and phenotypic means
    double meangen[2] = {0,0};
    double meanphen[2] = {0,0};

    double meanm[2][2] = {{0,0},{0,0}};
    double ssm[2][2] = {{0,0},{0,0}};

    // genotypic and phenotypic sums of squares 
    double ssgen[2][2] = {{0,0},{0,0}};
    double ssphen[2][2] = {{0,0},{0,0}};

    // get stats from the population
    for (size_t i = 0; i < Nf; ++i)
    {
        // loop through the different traits
        for (size_t j1 = 0; j1 < 2; ++j1)
        {
            meanphen[j1] += Females[i].phen[j1];
            meangen[j1] += Females[i].gen[j1];

            for (size_t j2 = 0; j2 < 2; ++j2)
            {
                ssgen[j1][j2] += Females[i].gen[j1] * Females[i].gen[j2];
                ssphen[j1][j2] += Females[i].phen[j1] * Females[i].phen[j2];

                meanm[j1][j2] += Females[i].m[j1][j2][0] + Females[i].m[j1][j2][1];
                ssm[j1][j2] += (Females[i].m[j1][j2][0] + Females[i].m[j1][j2][1]) * 
                    (Females[i].m[j1][j2][0] + Females[i].m[j1][j2][1]);
            }
        }
    }


    for (size_t i = 0; i < Nm; ++i)
    {
        // loop through the different traits
        for (size_t j1 = 0; j1 < 2; ++j1)
        {
            meanphen[j1] += Males[i].phen[j1];
            meangen[j1] += Males[i].gen[j1];

            for (size_t j2 = 0; j2 < 2; ++j2)
            {
                ssgen[j1][j2] += Males[i].gen[j1] * Males[i].gen[j2];
                ssphen[j1][j2] += Males[i].phen[j1] * Males[i].phen[j2];

                meanm[j1][j2] += Males[i].m[j1][j2][0] + Males[i].m[j1][j2][1];
                ssm[j1][j2] += (Males[i].m[j1][j2][0] + Males[i].m[j1][j2][1]) * 
                    (Males[i].m[j1][j2][0] + Males[i].m[j1][j2][1]);
            }
        }
    }

    DataFile << generation << ";";

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        meanphen[j1] /= Nf + Nm;
        meangen[j1] /= Nf + Nm; 
        DataFile << meanphen[j1] << ";";

        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            meanm[j1][j2] /= Nf + Nm;
            DataFile << meanm[j1][j2] << ";";
        }

    }

    double G[2][2];
    double P[2][2];
    double Gm[2][2];

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            G[j1][j2] = ssgen[j1][j2] / (Nf + Nm) 
                    - meangen[j1] * meangen[j2];

            P[j1][j2] = ssphen[j1][j2] / (Nf + Nm)
                    - meanphen[j1] * meanphen[j2];

            Gm[j1][j2] = ssm[j1][j2] / (Nf + Nm)
                    - meanm[j1][j2] * meanm[j1][j2];


            DataFile << G[j1][j2] << ";" << P[j1][j2] << ";" << Gm[j1][j2] << ";";
        }
    }

    double trace = G[0][0] + G[1][1];
    double det = G[0][0] * G[1][1] - G[1][0] * G[0][1];

    // calculate eigenvalues
    double ev1 = .5 * (trace + sqrt(trace*trace - 4 * det));
    double ev2 = .5 * (trace - sqrt(trace*trace - 4 * det));

    DataFile << trace << ";" << det << ";" << ev1 << ";" << ev2 << ";" << meanw << ";" << theta1 << ";" << theta2 << ";" << endl;
}

// write the headers of a datafile
void WriteDataHeaders()
{
    DataFile << "generation;";

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        DataFile << "meanz" << j1 + 1 << ";";

        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            DataFile << "meanm" << j1 + 1 << j2 + 1 << ";";
        }
    }

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            DataFile << "G" << (j1 + 1) << (j2 + 1) << ";"
                        << "P" << (j1 + 1) << (j2 + 1) << ";"
                        << "Gm" << (j1 + 1) << (j2 + 1) << ";";
        }
    }

    DataFile << "trace;det;ev1;ev2;meanw;theta1;theta2;" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
		Reproduce_Survive();

        // output stats every xth generation except for the last 2000 gens
        do_stats = generation % skip == 0;
        if (do_stats || generation > NumGen - 2000)
		{
			WriteData();
		}

	}

	WriteParameters();
}
