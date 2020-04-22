// simulation code for virus capsid assembly
// from Johnston et al. "Modelling the self-assembly of virus capsids", J Phys Condens Matt 22 104101 (2010)
// free preprint here https://arxiv.org/pdf/0910.1916.pdf

// a Metropolis Monte Carlo scheme is used to simulate the thermodynamic interactions of different agents in a box with periodic boundary conditions
// agents can be pentamers, hexamers, or crowding agents
// command-line parameters can change the default parameterisation
// by default, the size of bound clusters in the simulation and a current snapshot is output over time
// time-labelled snapshots an also be output to allow animations to be produced

// this code is an odd C/C++ hybrid. agents are objects with member functions, but no other object-orientated approaches are used. compile with g++ (with maths library)

// this is a refactoring and amalgamation of several bits of old code and as such is a work in progress!

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// the physical details of the simulation are brought in from here
#include "virus-all.c"

#define VERBOSE 0

// maximum number of agents for memory allocation
#define MAXPOLY 200

#define RND drand48()

// impose limits on integer val to be between lo and hi
void limiti(int *val, int lo, int hi)
{
  if(*val < lo) *val = lo;
  if(*val > hi) *val = hi;
} 

// impose limits on double val to be between lo and hi
void limitf(double *val, int lo, int hi)
{
  if(*val < lo) *val = lo;
  if(*val > hi) *val = hi;
}

int main(int argc, char *argv[])
{
  Protein *P;
  double *E;
  double *N;
  Vector temp;
  int i, j;
  EnergyParams EP;
  int t;
  double total;
  double Diff;
  double T = 0.05;
  int seed;
  Vector X;
  double alpha;
  int move;
  int acc;
  double movesize[MAXPOLY];
  double add, subtract;
  FILE *fp;
  int sizes[MAXPOLY];
  int types[MAXPOLY];
  int runs;
  int a;
  int crowd, prot;
  char labelstr[1000];
  char fstr[1000];
  int accrate, rejrate;
  int outputstep;

  // default params
  int NUMPOLY = 120;
  int NUMHEX = 0;
  int NUMCROWD = 0;
  double BOXSIZE = 65;
  double HEIGHT = 2.5;
  double TEMP = 0.2;
  int NCYCLES = 1e6;
  int NRUNS = 10;
  int RSEED = 121;
  int CROWD_TYPE = REP_CROWDER_S;
  int VIDEO = 0;
  int BURNIN = 1e5;
    
  // deal with command-line arguments
  for(i = 1; i < argc; i+=2)
    {
      if(strcmp(argv[i], "--npoly\0") == 0) NUMPOLY = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--nhex\0") == 0) NUMHEX = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--ncrowd\0") == 0) NUMCROWD = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--boxsize\0") == 0) BOXSIZE = atof(argv[i+1]);
      else if(strcmp(argv[i], "--height\0") == 0) HEIGHT = atof(argv[i+1]);
      else if(strcmp(argv[i], "--temp\0") == 0) TEMP = atof(argv[i+1]);
      else if(strcmp(argv[i], "--ncycles\0") == 0) NCYCLES = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--nruns\0") == 0) NRUNS = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--rseed\0") == 0) RSEED = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--video\0") == 0) VIDEO = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--help\0") == 0)
	{
	  printf("Options [defaults]:\n\n--npoly N\tnumber of agents [120]\n--nhex N\tnumber of hexamers [0]\n--ncrowd N\tnumber of crowders [0]\n--boxsize X\tbox dimension [65]\n--height X\tcapsomer height [2.5]\n--temp X\ttemperature [0.2]\n--ncycles N\tnumber MC cycles [1e6]\n--nruns N\tnumber simulations [3]\n--rseed N\trandom seed [121]\n--video N\trecord coords for animation [0]\n\n");
	  return 0;
	}
      else printf("Didn't understand argument %s\n", argv[i]);
    }
  limiti(&NUMPOLY, 0, 1000);
  limiti(&NUMHEX, 0, NUMPOLY);
  limiti(&NUMCROWD, 0, NUMPOLY-NUMHEX);
  limitf(&BOXSIZE, 0, 1000);
  limitf(&HEIGHT, 0, 10);
  limitf(&TEMP, 0, 10);
  limiti(&NCYCLES, 0, 1e7);
  limiti(&NRUNS, 0, 1e3);

  outputstep = (NUMPOLY*NCYCLES)/100;
  
  printf("Running with options:\nnumpoly %i\nnumhex %i\nnumcrowd %i\nboxsize %.2f\nheight %.3f\ntemp %.3f\nncycles %i\nnruns %i\nrseed %i\noutput step %i\n\n", NUMPOLY, NUMHEX, NUMCROWD, BOXSIZE, HEIGHT, TEMP, NCYCLES, NRUNS, RSEED, outputstep);
  printf("(run with --help to see command-line options)\n\n");
  
  sprintf(labelstr, "virusout-%i-%i-%i-%.2f-%.3f-%.3f-%i-%i-%i", NUMPOLY, NUMHEX, NUMCROWD, BOXSIZE, HEIGHT, TEMP, NCYCLES, NRUNS, RSEED);
	        
  // allocate memory for all players
  P = (Protein*)malloc(sizeof(Protein)*NUMPOLY);
  E = (double*)malloc(sizeof(double)*NUMPOLY*NUMPOLY);
  N = (double*)malloc(sizeof(double)*NUMPOLY);

  // initialise energy function params
  EP.rho = 3;
  for(i = 0; i < 6; i++)
    {
      for(j = 0; j < 6; j++)
	EP.sigma[i][j] = (i == 3 || j == 3 ? 20 : 10);
    }
  EP.e_r = 0.5; EP.e_att = 0.5;
  EP.basecutoff = 3.5; EP.apexcutoff = 20;

  // initialise random number generator
  srand48(RSEED);

  // initialise protein forms and randomly distribute crowders
  temp.x = temp.y = temp.z = 0;
  crowd = prot = 0;
  for(i = 0; i < NUMPOLY; i++)
    types[i] = 5;
  for(i = 0; i < NUMCROWD; i++)
    {
      a = RND*NUMPOLY;
      if(types[a] == 5)
	types[a] = 0;
      else i--;
    }
  for(i = 0; i < NUMHEX; i++)
    {
      a = RND*NUMPOLY;
      if(types[a] == 5)
	types[a] = 6;
      else i--;
    }

  // initialise agents in simulation
  for(i = 0; i < NUMPOLY; i++)
    {
      switch(types[i])
	{
	case 0:
	  crowd++;
	  P[i].Initialise(0, 0, 0, temp, CROWD_TYPE);
	  break;
	case 5:
	  P[i].Initialise(HEIGHT, 5, 5, temp, CAPSOMER);
	  prot++;
	  break;
	case 6:
	  P[i].Initialise(HEIGHT, 5, 6, temp, CAPSOMER);
	  prot++;
	  break;
	}
    }

  // loop through iterations
  for(runs = 0; runs < NRUNS; runs++)
    {
      sprintf(fstr, "%s-%i.sizes.txt", labelstr, runs);
      fp = fopen(fstr, "w"); fclose(fp);
      
      srand48(runs*5);
      for(i = 0; i < NUMPOLY; i++)
	{
	  movesize[i] = 2;
	}
     
      Grid(P, NUMPOLY, BOXSIZE);

      // calculate initial energy matrix. this is the only time we should need to do the full NxN calculation
      total = 0;
      for(i = 0; i < NUMPOLY; i++)
	{
	  for(j = i+1; j < NUMPOLY; j++)
	    {
	      E[i*NUMPOLY+j] = P[i].CalcPot(P[j], BOXSIZE, EP);
	      E[j*NUMPOLY+i] = E[i*NUMPOLY+j];
	      total += E[j*NUMPOLY+i];
	    }
	}

      accrate = rejrate = 0;
      
      // begin running through time
      for(t = 0; t < NCYCLES*NUMPOLY; t++)
	{
	  // set temperature
	  if(t < BURNIN)
	    T = 10;
	  else
	    T = TEMP;
	  
	  // pick agent and random move to make
	  seed = RND*NUMPOLY;
  	  if(P[seed].numpoints > 0)
	    move = RND*2; 
	  else move = 0;
	  X.x = (RND-0.5)*movesize[seed]; X.y = (RND-0.5)*movesize[seed]; X.z = (RND-0.5)*movesize[seed];
	  alpha = (RND-0.5)*movesize[seed];
  
	  switch(move)
	    {
	    case 0: P[seed].TranslateProtein(X); break;
	    case 1: P[seed].RotateProtein(X, P[seed].BaseCentre, alpha); break;
	    }

	  // compute new interaction energies with N-1 partners, and work out change in total energy
	  add = 0; subtract = 0;
	  for(i = 0; i < NUMPOLY; i++)
	    {
	      if(i != seed)
		{
		  N[i] = P[seed].CalcPot(P[i], BOXSIZE, EP);
		  add += N[i]; subtract += E[i*NUMPOLY+seed];
		}
	    }
	  Diff = add-subtract;

	  // MC acceptance step
	  acc = 1;
	  if(Diff > 0)
	    {
	      if(RND > exp(-Diff/T))
		acc = 0;
	    }

	  if(VERBOSE)
	    printf("Diff %f\n", Diff);
	  
	  if(acc == 1)
	    {
	      // accept this move, so retain new state and update energy matrix
	      accrate++;
	      for(i = 0; i < NUMPOLY; i++)
		{
		  if(i != seed)
		    {
		      E[i*NUMPOLY+seed] = N[i];
		      E[seed*NUMPOLY+i] = N[i];
		    }
		}
	      total += add-subtract;
	    }
	  else
	    {
	      rejrate++;
	      // reject this move, so reverse move
	      switch(move)
		{
		case 0: P[seed].TranslateProtein(X*-1); break;
		case 1: P[seed].RotateProtein(X, P[seed].BaseCentre, alpha*-1); break;
		}
	    }

	  // output state of simulation
	  if(t % outputstep == 0)
	    {
	      printf("run %i step %i total energy %f acc rate %f\n", runs, t/NUMPOLY, total, (double)accrate/(accrate+rejrate));
	      accrate = rejrate = 0;
	      if(runs == 0)
		{
		  // if this is the first simulation, produce a snapshot and possibly a time-labelled snapshot
		  if(VIDEO != 0)
		    {
		      sprintf(fstr, "%s-%i.txt", labelstr, t/outputstep);
		      OutputForMathematica(P, NUMPOLY, fstr, BOXSIZE);
		    }
		  sprintf(fstr, "%s.txt", labelstr);
		  OutputForMathematica(P, NUMPOLY, fstr, BOXSIZE);
		}
	      // compute distribution of cluster sizes and output
	      MaxClusSizeNew(E, NUMPOLY, sizes);
	      sprintf(fstr, "%s-%i.sizes.txt", labelstr, runs);
	      fp = fopen(fstr, "a");
	      fprintf(fp, "%i %f ", t/NUMPOLY, total);
	      for(i = 0; i < NUMPOLY; i++)
		fprintf(fp, "%i ", sizes[i]);
	      fprintf(fp, "\n");
	      fclose(fp);
	    }
	}
    }


  return 0;
}
