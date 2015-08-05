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
const int NumGen = 12000;

// population size
const int Npop = 1000; 

// number of generations to skip when outputting data
const int skip = 10;

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

// optima
double theta[2] = {0,0};

// fecundity per individual
double B = 0;

///////////////////     STATS       ///////////////////

// track number of individuals 
int NPop = 0;
int NKids = 0;


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
    // j: the number of independent traits
    // k: ploidy (in this case diploid), hence k in {0,1}
    double z[n_loci_g][2][2]; 

    double phen[2];

};

// allocate a population and a population of survivors
typedef Individual Population[Npop];
typedef Individual NewPopulation[Npop*2*10];
Population Pop;
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
void MutateG(Individual &ind)
{
    double a, b;
    int random_locus;
    
    /// first mutate the first haploid part of the genome
    if (gsl_rng_uniform(r) < n_loci_g * mu)
    {
        // generate new allelic increments
        gsl_ran_bivariate_gaussian(r, a1, a2, rmu, &a, &b);

        random_locus = gsl_rng_uniform_int(r, n_loci_g);

        // add the new allelic effect to trait one
        ind.z[random_locus][0][0] += a;
        // add the new allelic effect to trait two
        ind.z[random_locus][1][0] += b;

        ind.phen[0] += a;
        ind.phen[1] += b;
    }

    // repeat the same for the second haploid part of the genome
    if (gsl_rng_uniform(r) < n_loci_g * mu)
    {
        // generate new allelic increments
        gsl_ran_bivariate_gaussian(r, a1, a2, rmu, &a, &b);
        
        random_locus = gsl_rng_uniform_int(r, n_loci_g);

        // add the new allelic effect to trait one
        ind.z[random_locus][0][1] += a;
        // add the new allelic effect to trait two
        ind.z[random_locus][1][1] += b;
        
        ind.phen[0] += a;
        ind.phen[1] += b;
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
        << "theta_1;" << theta[0] << endl
        << "theta_2;" << theta[1] << endl
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


	// initialize the whole populatin
	for (size_t i = 0; i < Npop; ++i)
	{
        // loop through the different characters
        for (size_t l = 0; l < 2; ++l)
        {
            Pop[i].phen[l] = 0;

            // loop through each of the loci
            for (size_t j = 0; j < n_loci_g; ++j)
            {
                // loop through each of the alleles
                for (size_t k = 0; k < 2; ++k)
                {
                    Pop[i].z[j][k][l] = 0;
                }
            }
        }
	}
}

// create an offspring
void Create_Kid(int mother, int father, Individual &kid)
{
    Individual Mother = Pop[mother];
    Individual Father = Pop[father];

    assert(mother >= 0 && mother < Npop);
    assert(father >= 0 && father < Npop);

    kid.phen[0] = 0;
    kid.phen[1] = 0;

    for (size_t i = 0; i < n_loci_g; ++i)
    {
        for (size_t j = 0; j < 2; ++j)
        {
            kid.z[i][j][0] = Mother.z[i][j][gsl_rng_uniform_int(r,2)];
            kid.phen[0] += kid.z[i][j][0];

            kid.z[i][j][1] = Father.z[i][j][gsl_rng_uniform_int(r,2)];
            kid.phen[1] += kid.z[i][j][1];
        }
    }

    kid.phen[0] += gsl_ran_gaussian(r,1.0);
    kid.phen[1] += gsl_ran_gaussian(r,1.0);

    MutateG(kid);
}

double w(double const z0, double const z1)
{
    return(
            exp(
                -.5 * (z1 * (
                            z1 * omega[0][0] / (-omega[0][0]*omega[0][0] + omega[0][0] * omega[1][1])
                            - z0 * omega[1][0] / (-omega[1][0]*omega[1][0] + omega[0][0] * omega[1][1])
                            )

                        + z0 *(
                            - z1 * omega[1][0] / (-omega[0][0]*omega[0][0] + omega[0][0] * omega[1][1])
                            - z0 * omega[1][1] / (-omega[1][0]*omega[1][0] + omega[0][0] * omega[1][1])
                            )
                    )
               )
          );
}


// Survival of juveniles to reproductive adults
void Reproduce_Survive()
{
    NKids = 0;

    for (size_t i = 0; i < Npop; ++i)
    {
        // random mating
        size_t father;

        do {
            father = gsl_rng_uniform_int(r, Npop);
        }
        while (father == i);

        // produce kids and let them survive
        for (size_t j = 0; j < 2 * B; ++j)
        {
            Individual Kid;

            Create_Kid(i, father, Kid);

            if (gsl_rng_uniform(r) < w(Kid.phen[0], Kid.phen[1]))
            {
                NewPop[NKids++] = Kid;
            }
        }
    }
    
    if (NKids < Npop)
    {
        WriteParameters();
        exit(1);
    }


    for (size_t i = 0; i < Npop; ++i)
    {
        // let individual survive or not
        Pop[i] = NewPop[gsl_rng_uniform_int(r, NKids)];
    }

}


// write down summary statistics
void WriteData()
{
    double meanphen[2] = {0,0};
    double meangen[2] = {0,0};
    double ssg[2][2] = {{0,0},{0,0}};
    double ssphen[2][2] = {{0,0},{0,0}};


    // get stats from the population
    for (size_t i =  0; i < Npop; ++i)
    {
        // loop through the different traits
        for (size_t j1 = 0; j1 < 2; ++j1)
        {
            meanphen[j1] += Pop[i].phen[j1];

            // loop through the different loci
            for (size_t k1 = 0; k1 < n_loci_g; ++k1)
            {
                meangen[j1] += Pop[i].z[k1][j1][0] + Pop[i].z[k1][j1][1];
            }
        }
    }

    // get stats from the population
    for (size_t i = 0; i < Npop; ++i)
    {
        // loop through the different traits
        for (size_t j1 = 0; j1 < 2; ++j1)
        {
            // loop again through the different traits 
            // as we calculate covs
            for (size_t j2 = 0; j2 < 2; ++j2)
            {
                // phenotypic variance. Should be very similar
                // to genetic variance calculated above, except that
                // genetic variance measure also takes disequilibria
                // between loci into account
                ssphen[j1][j2] += Pop[i].phen[j1] * Pop[i].phen[j2];

                // loop again through the different loci
                for (size_t k1 = 0; k1 < n_loci_g; ++k1)
                {
                    // loop again through the different loci
                    for (size_t k2 = 0; k2 < n_loci_g; ++k2)
                    {
                        // loop through the different alleles per locus
                        for (size_t l1 = 0; l1 < 2; ++l1)
                        {
                            for (size_t l2 = 0; l2 < 2; ++l2)
                            {
                                // stats for m
                                ssg[j1][j2] += Pop[i].z[k1][j1][l1] * Pop[i].z[k2][j2][l2];
                            }
                        }
                    }
                }
            }
        }
    }

    
    DataFile << generation << ";";

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        meanphen[j1] /= Npop;
        meangen[j1] /= Npop; 
        DataFile << meanphen[j1] << ";";
    }

    double G[2][2];
    double P[2][2];

    for (size_t j1 = 0; j1 < 2; ++j1)
    {
        for (size_t j2 = 0; j2 < 2; ++j2)
        {
            G[j1][j2] = (ssg[j1][j2] / Npop 
                    - meangen[j1] * meangen[j2]);

            P[j1][j2] = (ssphen[j1][j2] / Npop
                    - meanphen[j1] * meanphen[j2]);

//            cout << "ssg_" << j1 + 1 << "_" << j2 + 1 << " " 
//                << ssg[j1][j2] << " " << meangen[j1] << " " 
//                << meangen[j2] << " " << endl;
//
            DataFile << G[j1][j2] << ";" << P[j1][j2] << ";";
        }
    }

    double trace = G[0][0] + G[1][1];
    double det = G[0][0] * G[1][1] - G[1][0] * G[0][1];

    double ev1 = .5 * (trace + sqrt(trace*trace - 4 * det));
    double ev2 = .5 * (trace - sqrt(trace*trace - 4 * det));

    DataFile << trace << ";" << det << ";" << ev1 << ";" << ev2 << ";" << endl;
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

    DataFile << "trace;det;ev1;ev2;" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		Reproduce_Survive();

        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
