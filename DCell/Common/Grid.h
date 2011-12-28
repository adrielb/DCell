#ifndef GRID_H_
#define GRID_H_

#include "Common.h"

typedef struct _Grid *Grid;

struct _Grid {
  Vec v;
  char name[PETSC_MAX_PATH_LEN]; // TODO: name when saving to disk
  iCoor p; // Position of this local grid within a global grid
  iCoor n; // Number of grid points
  Coor d;  // grid spacing in units of length
  PetscTruth is2D;
  PetscReal *v1;   // Resizable array
  PetscReal *grid; // Rectangular array of vec
  int dof;
  int SIZE;    // Acutal size: n.x * n.y * n.z * dof
  int MAXSIZE; // Size allocated in memory >= SIZE
  PetscViewer filePos, fileSize;
  PetscErrorCode (*Interpolate)( Grid g, Coor X, PetscReal *val);
};

PetscErrorCode GridCreate( Coor dh, iCoor pos, iCoor size, int dof, Grid *grid );
PetscErrorCode GridDestroy( Grid g );
PetscErrorCode GridSetDx( Grid g, Coor d );
PetscErrorCode GridSetName( Grid g, const char *n ); //TODO: Implement GridSetName
PetscErrorCode GridWrite( Grid g, int t );
PetscErrorCode GridResize( Grid g, iCoor pos, iCoor size );
PetscErrorCode GridGet( Grid g, void *grid );
PetscErrorCode GridGetBounds( Grid g, iCoor *p, iCoor *q );
PetscErrorCode GridDrawBorder( Grid g, int borderWidth, PetscReal fill );
PetscErrorCode GridFillRectangle( Grid g, Coor lo, Coor hi, PetscReal fill);
PetscErrorCode GridInterpolate( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate2D( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate3D( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridCopy( Grid g, Grid copy );
PetscErrorCode GridDuplicate( Grid g, Grid *newg );

typedef double (*GridFunction2D)(PetscReal **p, int x, int y, Coor d);
double GridFunction2D_Identity( PetscReal **p, int i, int j, Coor d);
double GridFunction2D_DerivX( PetscReal **p, int i, int j, Coor d);
double GridFunction2D_DerivY( PetscReal **p, int i, int j, Coor d);
double GridFunction2D_Curv( double **p, int i, int j, Coor d );
double Bilinear2D( GridFunction2D gf, PetscReal **v2, Coor dh, double x, double y);

typedef double (*GridFunction3D)(PetscReal ***p, int x, int y, int z, Coor d);
double GridFunction3D_Identity( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivX( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivY( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivZ( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivXX( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivYY( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivZZ( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivXZ( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivYZ( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_DerivXY( PetscReal ***p, int i, int j, int k, Coor d);
double GridFunction3D_Curv( double ***p, int i, int j, int k, Coor d );
double Bilinear3D( GridFunction3D gf, PetscReal ***v3, Coor dh, Coor p );

#endif /*GRID_H_*/
