#ifndef GRID_H_
#define GRID_H_

#include "Common.h"

typedef struct _Grid *Grid;

struct _Grid {
  Vec v;
  char name[PETSC_MAX_PATH_LEN]; // TODO: name when saving to disk
  AABB aabb; // bounding box in world coor
  iCoor p; // Position of this local grid within a global grid
  iCoor q; // q = p + n
  iCoor n; // Number of grid points
  Coor d;  // grid spacing in units of length
  PetscBool is2D;
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
PetscErrorCode GridFillBox( Grid g, Coor lo, Coor hi, PetscReal fill);
PetscErrorCode GridInterpolate( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate_Cubic( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate2D( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate3D( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridInterpolate2( Grid g, Coor X, PetscReal *val);
PetscErrorCode GridCopy( Grid g, Grid copy );
PetscErrorCode GridDuplicate( Grid g, Grid *newg );
inline PetscBool GridIndexInBox( Grid g, iCoor a );
inline void GridIndexToCoor( const Grid g, const Coor i, Coor *p );

typedef double (*GridFunction)(PetscReal *v, iCoor p, Coor d );
double GridFunction2D_DerivX(  PetscReal *v, iCoor p, Coor d );
double GridFunction2D_DerivY(  PetscReal *v, iCoor p, Coor d );
double GridFunction2D_Curv(    PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivX(  PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivY(  PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivZ(  PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivXX( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivYY( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivZZ( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivXZ( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivYZ( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_DerivXY( PetscReal *v, iCoor p, Coor d );
double GridFunction3D_Curv(    PetscReal *v, iCoor p, Coor d );
double GridBilinear( Grid g, GridFunction gf, Coor p );

#endif /*GRID_H_*/
