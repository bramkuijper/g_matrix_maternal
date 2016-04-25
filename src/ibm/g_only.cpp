//
// evolving G matrix in fluctuating environments


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
const int NumGen = 30000;

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

// mutation rate
double mu = 0;

// strengths of selection
double omega[2][2] = {{0,0},{0,0}};

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
    double z1[n_loci_g][2]; 
    double z2[n_loci_g][2]; 

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


// bivariate Gaussian fitness function as in Revell 2007
double v(double const z1, double const z2)
{
    // Exp[-1/2 * (z^T . Inverse[omega] . z )] can be simplified as
    return(
            exp(
                    -.5 * (z1*z1 * omega[1][1] + z2 * (-omega[0][1] * z1 - omega[1][0] * z1 + z2 * omega[0][0])) / 
                    (-omega[1][0] * omega[1][0] + omega[0][0] * omega[1][1])
               )
          );
}


// initialize simulations from command line arguments
void initArguments(int argc, char *argv[])
{
    a1 = atof(argv[1]);
    a2 = atof(argv[2]);
    rmu = atof(argv[3]);
    mu = atof(argv[4]);
    omega[0][0] = atof(argv[5]);
    omega[0][1] = atof(argv[6]);
    omega[1][0] = atof(argv[7]);
    omega[1][1] = atof(argv[8]);
    B = atof(argv[9]);

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
            ind.z1[locus_i][0] += a;
            // add the new allelic effect to trait two
            ind.z2[locus_i][0] += b;

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
            ind.z1[locus_i][1] += a;
            // add the new allelic effect to trait two
            ind.z2[locus_i][1] += b;

            ind.gen[0] += a;
            ind.gen[1] += b;
        }
    }
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
        << "mu;" << mu << endl
        << "omega_11;" << omega[0][0] << endl
        << "omega_12;" << omega[0][1] << endl
        << "omega_21;" << omega[1][0] << endl
        << "omega_22;" << omega[1][1] << endl
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
    
    // set the seed to the random number generator
    // stupidly enough, for gsl this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);


	// initialize the whole population
	for (size_t i = 0; i < Npop/2; ++i)
	{
        // loop through the different characters
        for (size_t l = 0; l < 2; ++l)
        {
            Males[i].phen[l] = 0;
            Males[i].gen[l] = 0;
            Females[i].phen[l] = 0;
            Females[i].gen[l] = 0;

        }
                    
        // loop through each of the loci
        for (size_t j = 0; j < n_loci_g; ++j)
        {
            // loop through each of the alleles
            for (size_t k = 0; k < 2; ++k)
            {
                Males[i].z1[j][k] = 0;
                Males[i].z2[j][k] = 0;
                Females[i].z1[j][k] = 0;
                Females[i].z2[j][k] = 0;
            }
        }
	}

    Nf = Npop/2;
    Nm = Npop/2;
}

// create an offspring
void Create_Kid(size_t const mother, size_t const father, Individual &kid)
{
    Individual Mother = Females[mother];
    Individual Father = Males[father];

    kid.gen[0] = 0;
    kid.gen[1] = 0;

    // loop through all the loci
    for (size_t i = 0; i < n_loci_g; ++i)
    {
        // // randomly choose one of both maternal alleles to inherit
        // kid.z1[i][0] = Mother.z1[i][gsl_rng_uniform_int(r,2)];
        // kid.gen[0] += kid.z1[i][0];

        // // randomly choose one of both paternal alleles to inherit
        // kid.z1[i][1] = Father.z1[i][gsl_rng_uniform_int(r,2)];
        // kid.gen[0] += kid.z1[i][1];
        // 
        // // randomly choose one of both maternal alleles to inherit
        // kid.z2[i][0] = Mother.z2[i][gsl_rng_uniform_int(r,2)];
        // kid.gen[1] += kid.z2[i][0];

        // // randomly choose one of both paternal alleles to inherit
        // kid.z2[i][1] = Father.z2[i][gsl_rng_uniform_int(r,2)];
        // kid.gen[1] += kid.z2[i][1];
        //
        if (gsl_rng_uniform(r) < 0.5)
        {
            kid.z1[i][0] = Mother.z1[i][0];
            kid.gen[0] += kid.z1[i][0];
            
            kid.z2[i][0] = Mother.z2[i][0];
            kid.gen[1] += kid.z2[i][0];
        }
        else
        {
            kid.z1[i][0] = Mother.z1[i][1];
            kid.gen[0] += kid.z1[i][0];
            
            kid.z2[i][0] = Mother.z2[i][1];
            kid.gen[1] += kid.z2[i][0];
        }
        
        
        if (gsl_rng_uniform(r) < 0.5)
        {
            kid.z1[i][1] = Father.z1[i][0];
            kid.gen[0] += kid.z1[i][1];
            
            kid.z2[i][1] = Father.z2[i][0];
            kid.gen[1] += kid.z2[i][1];
        }
        else
        {
            kid.z1[i][1] = Father.z1[i][1];
            kid.gen[0] += kid.z1[i][1];
            
            kid.z2[i][1] = Father.z2[i][1];
            kid.gen[1] += kid.z2[i][1];
        }
    }

    
    MutateG(kid);
   
    // add environmental variance to each trait by adding a random number
    // drawn from a normal distribution to each phenotype
    kid.phen[0] = kid.gen[0] + gsl_ran_gaussian(r,1.0);
    kid.phen[1] = kid.gen[1] + gsl_ran_gaussian(r,1.0);

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
                if (generation == 11500)
                {
                    double sumz1 = 0;
                    double sumz2 = 0;
                    for (size_t i = 0; i < n_loci_g; ++i)
                    {
                        sumz1 += Kid.z1[i][0] + Kid.z1[i][1];
                        sumz2 += Kid.z2[i][0] + Kid.z2[i][1];
                        //cout << NKids << ";" << Kid.z1[i][0] + Kid.z1[i][1] << ";" << Kid.z1[i][0] + Kid.z1[i][1] << endl;
                    }

                    //cout << NKids << ";" << sumz1 << ";" << sumz2 << endl;
                }

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
    
    if (NKids < Npop)
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

        // delete kid (no resampling possible) by copying kid
        // from the end of the stack and reducing Nkids by one
        NewPop[random_kid] = NewPop[--NKids];
    }
}


// write down summary statistics
void WriteData()
{
    // genotypic and phenotypic means
    double meangen[2] = {0,0};
    double meanphen[2] = {0,0};

    // genotypic and phenotypic sums of squares 
    double ssgen[2][2] = {{0,0},{0,0}};
    double ssphen[2][2] = {{0,0},{0,0}};

    double covar = 0;

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
            }
        }
    }

//    double mean0 = meangen[0] / Nf;
//    double mean1 = meangen[1] / Nf;
//
//    for (size_t i = 0; i < Nf; ++i)
//    {
//        covar += (Females[i].gen[0] - mean0) * (Females[i].gen[1] - mean1);
//
//    }

    //cout << mean0 << ";" << mean1 << ";" << (covar / Nf) << endl;


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
            }
        }
    }

//    // get stats from the population
//    for (size_t i = 0; i < Npop; ++i)
//    {
//        // loop through the different traits
//        for (size_t j1 = 0; j1 < 2; ++j1)
//        {
//            // loop again through the different traits 
//            // as we calculate covs
//            for (size_t j2 = 0; j2 < 2; ++j2)
//            {
//                // phenotypic variance. Should be very similar
//                // to genetic variance calculated above, except that
//                // genetic variance measure also takes disequilibria
//                // between loci into account
//                ssphen[j1][j2] += Pop[i].phen[j1] * Pop[i].phen[j2];
//
//                // loop through the different loci
//                for (size_t k1 = 0; k1 < n_loci_g; ++k1)
//                {
//                    // loop again through the different loci
//                    for (size_t k2 = 0; k2 < n_loci_g; ++k2)
//                    {
//                        // loop through the different alleles per locus
//                        for (size_t l1 = 0; l1 < 2; ++l1)
//                        {
//                            // loop through the different alleles per locus
//                            for (size_t l2 = 0; l2 < 2; ++l2)
//                            {
//                                // calculate G
//                                ssg[j1][j2] += Pop[i].z[k1][j1][l1] * Pop[i].z[k2][j2][l2];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//
    
    DataFile << generation << ";";

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        meanphen[j1] /= Nf + Nm;
        meangen[j1] /= Nf + Nm; 
        DataFile << meanphen[j1] << ";";
    }

    double G[2][2];
    double P[2][2];

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            G[j1][j2] = ssgen[j1][j2] / (Nf + Nm) 
                    - meangen[j1] * meangen[j2];

            P[j1][j2] = ssphen[j1][j2] / (Nf + Nm)
                    - meanphen[j1] * meanphen[j2];

//            cout << "ssg_" << j1 + 1 << "_" << j2 + 1 << " " 
//                << ssg[j1][j2] << " " << meangen[j1] << " " 
//                << meangen[j2] << " " << endl;
//
            DataFile << G[j1][j2] << ";" << P[j1][j2] << ";";
        }
    }

    double trace = G[0][0] + G[1][1];
    double det = G[0][0] * G[1][1] - G[1][0] * G[0][1];

    // calculate eigenvalues
    double ev1 = .5 * (trace + sqrt(trace*trace - 4 * det));
    double ev2 = .5 * (trace - sqrt(trace*trace - 4 * det));

    DataFile << trace << ";" << det << ";" << ev1 << ";" << ev2 << ";" << meanw << ";" << endl;
}

// write the headers of a datafile
void WriteDataHeaders()
{
    DataFile << "generation;meanz1;meanz2;";

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            DataFile << "G" << (j1 + 1) << (j2 + 1) << ";"
                        << "P" << (j1 + 1) << (j2 + 1) << ";";
        }
    }

    DataFile << "trace;det;ev1;ev2;meanw;" << endl;
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

        // some output for testing
//        if (generation == NumGen)
//        {
//            for (int i = 0; i < Npop; ++i)
//            {
//                for (int k = 0; k < n_loci_g; ++k)
//                {
//                    cout << generation << ";" << i << ";" << Pop[i].z[k][0][0] << ";" << Pop[i].z[k][1][0] << ";" << endl
//                        << generation << ";" << i << ";" << Pop[i].z[k][0][1] << ";" << Pop[i].z[k][1][1] << ";" << endl;
//                }
//            }
//        }
	}

	WriteParameters();
}
