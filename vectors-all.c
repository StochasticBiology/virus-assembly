// custom processes for 3D vectors in periodic boundary conditions

#ifndef __VECTORS_C_

#define __VECTORS_C_

#include <math.h>

// vector structure
typedef struct
{
  double x;
  double y;
  double z;
} Vector;

// operators for overloading
Vector operator +(Vector V1, Vector V2);
Vector operator -(Vector V1, Vector V2);
Vector operator ^(Vector V1, Vector V2);
double operator ,(Vector V1, Vector V2);
Vector operator *(Vector V1, double alpha);

// useful functions
int nint(double arg);
Vector NormaliseDir(Vector Dir);
void ImposeBC(Vector *V, double BOX);
double VectorDist(Vector A, Vector B, double BOX);

// vector addition
Vector operator +(Vector V1, Vector V2)
{
  register Vector V3;

  V3.x = V1.x + V2.x;
  V3.y = V1.y + V2.y;
  V3.z = V1.z + V2.z;

  return V3;
}

// vector subtraction
Vector operator -(Vector V1, Vector V2)
{
  register Vector V3;

  V3.x = V1.x - V2.x;
  V3.y = V1.y - V2.y;
  V3.z = V1.z - V2.z;

  return V3;
}

// vector cross product
Vector operator ^(Vector V1, Vector V2)
{
  register Vector V3;

  V3.x = V1.y * V2.z - V1.z * V2.y;
  V3.y = V1.z * V2.x - V1.x * V2.z;
  V3.z = V1.x * V2.y - V1.y * V2.x;

  return V3;
}

// vector dot product
double operator ,(Vector V1, Vector V2)
{
  return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
}

// vector multiplication by a scalar alpha
Vector operator *(Vector V1, double alpha)
{
  register Vector V2;

  V2.x = V1.x * alpha;
  V2.y = V1.y * alpha;
  V2.z = V1.z * alpha;

  return V2;
}

// rounding function
int nint(double arg)
{
  if(arg == 0) return 0;
  if(arg > 0)
    {
      if(arg - (int) arg > 0.5) return (int) arg + 1;
      else return (int) arg;
    }
  else
    {
      return -nint(-arg);
    }
}

// normalise a direction vector
Vector NormaliseDir(Vector Dir)
{
  register double Total;
  register double radical;
  Vector X;

  // deal with awkward vectors that shouldn't arise
  if(fabs(Dir.x) > 1000 ||fabs(Dir.y) > 1000 || fabs(Dir.z) > 1000)
    {
      Dir.x = 1; Dir.y = 0; Dir.z = 0;
    }

  // something is wrong here
  if((radical = (pow(Dir.x, 2) + pow(Dir.y, 2) + pow(Dir.z, 2))) <= 0)
    {
      printf("Negative norm?!\n");
    return Dir;
    }

  // do the normalisation
  X.x = Dir.x; X.y = Dir.y; X.z = Dir.z;

  Total = pow(radical, 0.5);

  if(Total == 0) Total = 1;
  X.x /= Total;
  X.y /= Total;
  X.z /= Total;

  return X;
}

// impose periodic boundary conditions on a vector
void ImposeBC(Vector *V, double BOX)
{
  if(V->x < 0) V->x = (floor(-V->x / BOX)+1) * BOX + V->x;
  if(V->y < 0) V->y = (floor(-V->y / BOX)+1) * BOX + V->y;
  if(V->z < 0) V->z = (floor(-V->z / BOX)+1) * BOX + V->z;

  if(V->x > BOX) V->x -= (floor(V->x / BOX)) *BOX;
  if(V->y > BOX) V->y -= (floor(V->y / BOX)) *BOX;
  if(V->z > BOX) V->z -= (floor(V->z / BOX)) *BOX;
}


// compute distance between two vectors given periodic boundary conditions
double VectorDist(Vector A, Vector B, double BOX)
{
  Vector V1 = A, V2 = B;
  double dx, dy, dz;

  dx = V1.x - V2.x;
  dy = V1.y - V2.y;
  dz = V1.z - V2.z;

  if(BOX)
    {
      dx -= BOX*nint(dx/BOX);
      dy -= BOX*nint(dy/BOX);
      dz -= BOX*nint(dz/BOX);
    }

  return sqrt(dx*dx + dy*dy + dz*dz);
}

#endif
