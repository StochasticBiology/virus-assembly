// definition of physical structures and functions for virus capsid assembly

#ifndef __VIRUS_C_

#define __VIRUS_C_

// different types of protein agent
#define CAPSOMER 0
#define REP_CROWDER_S 1
#define ATT_CROWDER_S 2
#define REP_CROWDER_L 3
#define ATT_CROWDER_L 4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// get vector functions
#include "vectors-all.c"

// structure for energy function parameters
// V(c_i, c_j) = V_apex(|a_i - a_j|) + sum_{u,v} V_M (|p_iu - p_jv|)
// V_apex(r) = e_r (sigma / r)^12
// V_M(r) = e_att [ exp( rho*(1 - r/r_e) ) - 2 ] exp( rho*(1 - r/r_e) )
//
// apexcutoff is a distance beyond which we don't interact at all
// basecutoff is a distance threshold for base vertex interactions
struct EnergyParams
{
  double rho;
  double sigma[5][5];
  double e_r;
  double e_att;
  double apexcutoff, basecutoff;
  double kappa, theta;
};

// structure for a protein
// this can be a pentamer, hexamer, or crowder
// quasi-object-orientated: member functions do various calculations and perturbations
struct Protein
{
  double Height;
  double Radius;
  int numpoints;
  int type;

  Vector *Vertices;

  Vector Apex;
  Vector FirstPoint;
  Vector BaseCentre;

  void Copy(Protein Src);
  void DoVertices(void);
  double CalcPot(Protein Partner, double BOX, EnergyParams P);
  void RotateProtein(Vector NDir, Vector Origin, double Rotation);
  void TranslateProtein(Vector X);
  void Initialise(double h, double r, int n, Vector V, int t);
  void Free(void);
};

double BaseEnergy(double r, EnergyParams P);
double ApexEnergy(double r, EnergyParams P, int t1, int t2, double sigma);
void OutputForMathematica(Protein *Proteins, int POLYNUM, char *label, double BOX);
void OutputForGnuplot(Protein *Proteins, int POLYNUM, char *label, double BOX);
int ClusterCount(int *BondNeighbours, int *AlreadyCounted, int cluscount, int NumInCluster, int POLYNUM);
void ReadCapsid(Protein *Proteins, int POLYNUM, double HEIGHT, char *filename);
int MaxClusSize(double *Energies, int POLYNUM, int *sizes, int *biggest, int *numcomp);
void Grid(Protein *Proteins, int POLYNUM, double BOX);

#define RND drand48()

// copy Src to self
void Protein::Copy(Protein Src)
{
  int i;

  if(Src.numpoints != numpoints)
    {
      // need to reallocate a different number of points
      free(Vertices);
      Vertices = (Vector*)malloc(sizeof(Vector)*Src.numpoints);
    }
  
  // copy everything
  numpoints = Src.numpoints;
  Height = Src.Height;
  Radius = Src.Radius;
  Apex = Src.Apex;
  FirstPoint = Src.FirstPoint;
  BaseCentre = Src.BaseCentre;
  for(i = 0; i < Src.numpoints; i++)
    Vertices[i] = Src.Vertices[i];
}

// compute individual vertex positions from descriptive vectors Apex, BaseCentre, FirstPoint
void Protein::DoVertices(void)
{
  Vector Corner;
  Vector Direction;
  int counter;
  double alpha = 3.1415926*2/numpoints;

  if(numpoints > 0)
    {
      Direction = Apex - BaseCentre;
      Direction = NormaliseDir(Direction);

      Corner = FirstPoint - BaseCentre;

      Vertices[0] = Corner + BaseCentre;

      for(counter = 1; counter < numpoints; counter++)
	{
	  Corner = Corner * cos(alpha) + Direction * (Direction,Corner) * (1 - cos(alpha)) + (Corner^Direction)*sin(alpha);
	  Vertices[counter] = Corner + BaseCentre;
	}
    }
}

// calculate interaction energy with an agent Partner
double Protein::CalcPot(Protein Partner, double BOX, EnergyParams P)
{
  int counter1, counter2;
  double energy;
  double buf;
  double theta;
  int onebonded = 0;
  double bufdist;
  Vector V1, V2;
  double vv;
  double apexdist;

  // compute apex distance cutoff and return zero if we're out of range
  apexdist = VectorDist(Apex, Partner.Apex, BOX);
  if(P.apexcutoff && apexdist > P.apexcutoff) return 0;

  energy = 0;

  // loop through vertices, adding interaction energy for each pair
  for(counter1 = 0; counter1 < numpoints; counter1++)
    {
      for(counter2 = 0; counter2 < Partner.numpoints; counter2++)
	{
	  bufdist = VectorDist(Vertices[counter1], Partner.Vertices[counter2], BOX);
	  energy += BaseEnergy(bufdist, P); // cutofff and Morse?
	}
    }

  // add apex repulsion
  energy += ApexEnergy(apexdist, P, type, Partner.type, P.sigma[type][Partner.type]);

  return energy;
}

// rotate self Rotation radians in direction NDir given Origin
void Protein::RotateProtein(Vector NDir, Vector Origin, double Rotation)
{
  NDir = NormaliseDir(NDir);

  Apex = Apex - Origin;
  BaseCentre = BaseCentre - Origin;
  FirstPoint = FirstPoint - Origin;

  Apex = Apex * cos(Rotation) + NDir * (NDir, Apex) * (1 - cos(Rotation)) + (Apex^NDir) * sin(Rotation);
  BaseCentre = BaseCentre * cos(Rotation) + NDir * (NDir, BaseCentre) * (1 - cos(Rotation)) + (BaseCentre^NDir) * sin(Rotation);
  FirstPoint = FirstPoint * cos(Rotation) + NDir * (NDir, FirstPoint) * (1 - cos(Rotation)) + (FirstPoint^NDir) * sin(Rotation);

  Apex = Apex + Origin;
  BaseCentre = BaseCentre + Origin;
  FirstPoint = FirstPoint + Origin;

  DoVertices();
}

// translate self by vector X
void Protein::TranslateProtein(Vector X)
{
  Apex = Apex + X;
  BaseCentre = BaseCentre + X;
  FirstPoint = FirstPoint + X;
  DoVertices();
}

// initialise self
void Protein::Initialise(double h, double r, int n, Vector V, int t)
{
  Height = h;
  Radius = r;
  numpoints = n;
  type = t;
  Vertices = (Vector*)malloc(sizeof(Vector)*n);
  BaseCentre = V;
  Apex = V; Apex.z += h;
  FirstPoint = V; FirstPoint.x += r;
  DoVertices();
}

// free allocated memory for self vertices
void Protein::Free(void)
{
  free(Vertices);
}

// interaction energy for distance r
double BaseEnergy(double r, EnergyParams P)
{
  double t;

  if(P.basecutoff && r > P.basecutoff) return 0;
  t = exp(P.rho*(1-r));
  t = (t-2)*t;

  return t;
}

// apex interactions with distance r and agent types t1, t2
double ApexEnergy(double r, EnergyParams P, int t1, int t2, double sigma)
{
  double t;

  if(t1 == CAPSOMER && t2 == CAPSOMER)
    {
      // capsomer-capsomer apex repulsion
      t = sigma/r;
      t = t*t; t = t*t;
      t = t*t*t;
      t *= P.e_r;
    }
  else if ((t1 == ATT_CROWDER_S && t2 == ATT_CROWDER_S) || (t1 == ATT_CROWDER_L && t2 == ATT_CROWDER_L))
    {
      // attractive crowding agents interacting
      t = sigma/r;
      t = t*t*t; 
      t = t*t;
      t = t*t-t;
      t *= P.e_att;
    }
  else
    {
      // repulsive crowding agent
      t = sigma/r;
      t = t*t; t = t*t;
      t = t*t*t;
      t *= P.e_r;
    }



  return t;
}

// arrange a set of agents in a grid in a box 
void Grid(Protein *Proteins, int POLYNUM, double BOX)
{
  int a = ceil(pow((double)POLYNUM, 1.0/3.0));
  int x, y, z;
  int count = 0;

  for(x = 0; x < a; x++)
    {
      for(y = 0; y < a; y++)
	{
	  for(z = 0; z < a; z++)
	    {
	      if(count < POLYNUM)
		{
		  Proteins[count].BaseCentre.x = x*BOX/a;
		  Proteins[count].BaseCentre.y = y*BOX/a;
		  Proteins[count].BaseCentre.z = z*BOX/a;
		  Proteins[count].Apex = Proteins[count].BaseCentre;
		  Proteins[count].Apex.z += Proteins[count].Height;
		  Proteins[count].FirstPoint = Proteins[count].BaseCentre;
		  Proteins[count].FirstPoint.x += Proteins[count].Radius;
		  Proteins[count].DoVertices();
		  count++;
		}
	    }
	}
    }
}


// output current state of the system for visualisation in Mathematica
void OutputForMathematica(Protein *Proteins, int POLYNUM, char *label, double BOX)
{
  int i, j, k;
  FILE *fp1, *fp2;
  char s1[400], s2[400];
  Vector off;
  int prot, crowd;

  prot = crowd = 0;

  // wipe existing files with this label
  for(i = 0; i <= 6; i++)
    {
      sprintf(s1, "%s-%i-sides-faces.dat", label, i);
      sprintf(s2, "%s-%i-sides-bases.dat", label, i);
      fp1 = fopen(s1, "w"); fp2 = fopen(s2, "w");
      fclose(fp1); fclose(fp2);
      if(i == 0) i = 4;
    }

  // go through agents in system
  for(i = 0; i < POLYNUM; i++)
    {
      // open file corresponding to this agent type
      sprintf(s1, "%s-%i-sides-faces.dat", label, Proteins[i].numpoints);
      sprintf(s2, "%s-%i-sides-bases.dat", label, Proteins[i].numpoints);
      fp1 = fopen(s1, "a"); fp2 = fopen(s2, "a");

      off.x = BOX*nint(Proteins[i].BaseCentre.x/BOX);
      off.y = BOX*nint(Proteins[i].BaseCentre.y/BOX);
      off.z = BOX*nint(Proteins[i].BaseCentre.z/BOX);

      if(Proteins[i].numpoints == 0)
	{
	  // output crowding agent
	  fprintf(fp1, "%.4f %.4f %.4f\n", Proteins[i].Apex.x-off.x, Proteins[i].Apex.y-off.y, Proteins[i].Apex.z-off.z);
	  crowd++;
	}
      else prot++;

      // output vertices for protein agent
      for(j = 0; j < Proteins[i].numpoints; j++)
	{
	  if(j == 0) k = Proteins[i].numpoints-1;
	  else k = j-1;

	  fprintf(fp1, "%.4f %.4f %.4f\n", Proteins[i].Apex.x-off.x, Proteins[i].Apex.y-off.y, Proteins[i].Apex.z-off.z);
	  fprintf(fp1, "%.4f %.4f %.4f\n", Proteins[i].Vertices[j].x-off.x, Proteins[i].Vertices[j].y-off.y, Proteins[i].Vertices[j].z-off.z);
	  fprintf(fp1, "%.4f %.4f %.4f\n", Proteins[i].Vertices[k].x-off.x, Proteins[i].Vertices[k].y-off.y, Proteins[i].Vertices[k].z-off.z);

	  fprintf(fp2, "%.4f %.4f %.4f\n", Proteins[i].Vertices[j].x-off.x, Proteins[i].Vertices[j].y-off.y, Proteins[i].Vertices[j].z-off.z);
	}
      fclose(fp1); fclose(fp2);
    }
  printf("Output: %i proteins, %i crowders\n", prot, crowd);
}

void OutputForGnuplot(Protein *Proteins, int POLYNUM, char *label, double BOX)
{
  int i, j, k;
  FILE *fp1;
  char s1[400];
  Vector t;
  
  sprintf(s1, "%s-gnuplot-plot.dat", label);
  fp1 = fopen(s1, "w"); 

  for(i = 0; i < POLYNUM; i++)
    {
      for(j = 0; j < Proteins[i].numpoints; j++)
	{
	  t = Proteins[i].Vertices[j];
	  ImposeBC(&t, BOX);
	  fprintf(fp1, "%.4f %.4f %.4f\n", t.x, t.y, t.z);
	}
    }
  fclose(fp1);
}

// recursive function for cluster counting
int ClusterCount(int *BondNeighbours, int *AlreadyCounted, int cluscount, int NumInCluster, int POLYNUM)
{
  int pyramidcounter;

  // only continue if we haven't already met this agent
  if(AlreadyCounted[cluscount] == 0)
    {
      // mark this agent as met and increment cluster size
      AlreadyCounted[cluscount] = 1;
      NumInCluster++;

      // go through all agents
      for(pyramidcounter = 0; pyramidcounter < POLYNUM; pyramidcounter++)
	{
	  if(pyramidcounter != cluscount && BondNeighbours[cluscount*POLYNUM+pyramidcounter] == 1)
	    {
	      // this is not self (though we shouldn't have self pairs marked as partners anyway) and is a bond partner, so recurse through
	      NumInCluster = ClusterCount(BondNeighbours, AlreadyCounted, pyramidcounter, NumInCluster, POLYNUM);
	    }
	}
    }

  return NumInCluster;
}


// get cluster sizes in the system
// this has been updated in refactoring, must be tested
int MaxClusSizeNew(double *Energies, int POLYNUM, int *sizes)
{
  int maxclussize = 0;
  int cluscount;
  int NumInCluster;
  int *AlreadyCounted;
  int *Clusters;
  int c1, c2;
  int *BondNeighbours;
  int i;
  int count, tcount;

  // allocate memory for structures we need
  BondNeighbours = (int *)malloc(sizeof(int)*POLYNUM*POLYNUM);
  AlreadyCounted = (int *)malloc(sizeof(int)*POLYNUM);
  Clusters = (int*)malloc(sizeof(int)*POLYNUM);

  // initialise all recorded cluster sizes to zero
  for(c1 = 0; c1 < POLYNUM; c1++)
    sizes[c1] = 0;

  // loop through pairs of agents, marking those that are interacting (energy < -1) in the BondNeighbours matrix
  for(c1 = 0; c1 < POLYNUM; c1++)
    {
      for(c2 = c1; c2 < POLYNUM; c2++)
	{
	  if(c2 == c1) BondNeighbours[c1*POLYNUM+c2] = 0;
	  else
	    {
	      BondNeighbours[c1*POLYNUM+c2] = (Energies[c1*POLYNUM+c2] < -1); 
	      BondNeighbours[c2*POLYNUM+c1] = (Energies[c1*POLYNUM+c2] < -1);
	    }
	}
    }

  // initialise AlreadyCounted flag to zero for everyone
  for(c1 = 0; c1 < POLYNUM; c1++)
    AlreadyCounted[c1] = 0;

  // loop through agents
  for(c1 = 0; c1 < POLYNUM; c1++)
    {
      NumInCluster = 0;
      // recursively count this size of the cluster that this agent is in, unless this cluster has already been counted
      NumInCluster = ClusterCount(BondNeighbours, AlreadyCounted, c1, NumInCluster, POLYNUM);
      // record this size of this cluster in Clusters
      Clusters[c1] = NumInCluster;
    }

  // build up a histogram of cluster sizes
  for(cluscount = 0; cluscount <= POLYNUM; cluscount++)
    {
      sizes[cluscount] = 0;
      for(c1 = 0; c1 < POLYNUM; c1++)
	{
	  if(Clusters[c1] == cluscount)
	    sizes[cluscount]++;
	}

      if(sizes[cluscount] != 0)
	maxclussize = cluscount;
    }

  // free up memory
  free(AlreadyCounted);
  free(Clusters);
  free(BondNeighbours);

  return maxclussize;
}

// suspect this is deprecated -- reading a structure in from a file containing coordinates?
void ReadCapsid(Protein *Proteins, int POLYNUM, double HEIGHT, char *filename)
{
  int counter = 0;
  FILE *fp;
  int capcount, numcount, pointcount;
  char inp;
  double num;
  Vector Temp;
  char string[40];

  fp = fopen(filename, "r");

  capcount = 0; numcount = 0; pointcount = 0;
  inp = ' ';
  while(capcount < POLYNUM)
    {
      if(inp == EOF) break;
      inp = fgetc(fp);
      if((inp < '0' || inp > '9') && inp != '.')
	{
	  num = atof(string);
	  sprintf(string, "                    ");
	  counter = 0;
	  numcount++;
	  if(numcount == 1) Temp.x = num-34;
	  if(numcount == 2) Temp.y = num-34;
	  if(numcount == 3)
	    {
	      inp = fgetc(fp);
              Temp.z = num;
	      numcount = 0;
	      Proteins[capcount].Vertices[pointcount] = Temp;
	      pointcount++;
	      if(pointcount == 6)// + (capcount < NUMPENTS ? 0 : 1))
		{
		  Proteins[capcount].Apex = Proteins[capcount].Vertices[0];
		  Proteins[capcount].FirstPoint = Proteins[capcount].Vertices[1];
		  Proteins[capcount].BaseCentre = Proteins[capcount].Vertices[1] + Proteins[capcount].Vertices[2] + Proteins[capcount].Vertices[3] + Proteins[capcount].Vertices[4] + Proteins[capcount].Vertices[5];
		  Proteins[capcount].BaseCentre = Proteins[capcount].BaseCentre * 0.2;
		  Proteins[capcount].Apex = Proteins[capcount].BaseCentre + (NormaliseDir(Proteins[capcount].Apex - Proteins[capcount].BaseCentre)) * HEIGHT;
		  pointcount = 0;
		  capcount++;
		}
	    }
	}
      else if(inp != EOF)
	{
	  string[counter] = inp;
	  counter++;
	}
    }
  //   fclose(fp);
}


#endif


