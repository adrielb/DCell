#include "LevelSetMethod.h"

typedef struct _Image {
  iCoor p,q; // low and high coordinates of img
  iCoor s; // scale
  Coor d;
  PetscReal negThickness, posThickness;
  int* img;
} *Image;

PetscErrorCode ImageCreate( Image *img ) {
  Image i;
  PetscErrorCode ierr;
  ierr = PetscNew(struct _Image,&i); CHKERRQ(ierr);
  *img = i;
  return 0;
}

PetscErrorCode ImageAllocate( Image img ) {
  PetscErrorCode ierr;
  ierr = PetscMalloc(sizeof(int)*n.x*n.y,&img->img); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode LevelSetDrawContour( LevelSet ls, Image img ) {
  int i,j,b;
  iCoor p; // image coor
  iCoor *band;// grid coor
  Coor a;  // grid coor
  iCoor s = img->s;  //scale
  PetscReal phi;
  PetscErrorCode ierr;

  for (b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,i,&band); CHKERRQ(ierr);
    for (j = 0; j < s.y; ++j) {
      for (i = 0; i < s.x; ++i) {
        a.x = i / (double)s.x + band->x;
        a.y = j / (double)s.y + band->y;
        p.x = band->x * s.x + i;
        p.y = band->y * s.y + j;
        phi = Bilinear( phi, a );
        if( img->negThickness < phi && phi < img->posThickness ) {
          img->img[p.x + n.x*p.y] = 255;
        }
      }
    }
  }
}

