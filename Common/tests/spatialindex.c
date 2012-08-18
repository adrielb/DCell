#include "Common.h"

int SphereQueryTest(int argc, char **args);
int BoxQueryTest( );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = BoxQueryTest(); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BoxQueryTest"
int BoxQueryTest( )
{
  const PetscReal eps = .1;
  Coor lo = {0,0,0};
  Coor hi = {4,4,4};
  PetscReal dx = 1;
  int n = 64*100;
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  Coor dh = {dx,dx,dx};
  SpatialIndex sidx;
  ierr = SpatialIndexCreate( "test", &sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( sidx, lo, hi, dh ); CHKERRQ(ierr);

  PetscRandom rnd;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd,lo.x-2,hi.x+2); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND48); CHKERRQ(ierr);

  Coor *pt;
  int i, j;
  Coor center;
  Array items;
  AABB box;

  PetscViewer txtfile;
  char filename[] = "spatialindex.dat";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&txtfile); CHKERRQ(ierr);

  ierr = PetscMalloc(sizeof(Coor)*n,&pt); CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
    ierr = PetscRandomGetValueReal(rnd,&pt[i].x); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&pt[i].y); CHKERRQ(ierr);
//    ierr = PetscRandomGetValueReal(rnd,&pt[i].z); CHKERRQ(ierr);
    pt[i].z = 0;
    ierr = SpatialIndexInsertPoint(sidx,pt[i],&pt[i]); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(txtfile, "%f %f %f ", pt[i].x, pt[i].y, pt[i].z); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(txtfile, "\n", pt[i].x, pt[i].y, pt[i].z); CHKERRQ(ierr);

  // loop to test memory leaks
  for (j = 0; j < 10; ++j) {
    printf("j: %d\n",j);

    ierr = PetscRandomGetValueReal(rnd,&center.x); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&center.y); CHKERRQ(ierr);
//    ierr = PetscRandomGetValueReal(rnd,&center.z); CHKERRQ(ierr);
    center.z = 0;
    box.lo.x = center.x - eps;
    box.lo.y = center.y - 2*eps;
    box.lo.z = center.z - 2*eps;
    box.hi.x = center.x + eps;
    box.hi.y = center.y + 2*eps;
    box.hi.z = center.z + 2*eps;

    ierr = SpatialIndexQueryPointsBox( sidx, box, &items ); CHKERRQ(ierr);

    const int len = ArrayLength(items);
    const Coor **data = ArrayGetData(items);
    ierr = PetscViewerASCIIPrintf(txtfile,
        "%f %f %f  %f %f %f\n", box.lo.x, box.lo.y, box.lo.z,
                                box.hi.x, box.hi.y, box.hi.z ); CHKERRQ(ierr);
    if( len == 0 ) {
      ierr = PetscViewerASCIIPrintf( txtfile, "0 0 0 ");
    }
    for (i = 0; i < len; ++i) {
      ierr = PetscViewerASCIIPrintf( txtfile,
          "%f %f %f ", data[i]->x, data[i]->y, data[i]->z);
    }
    ierr = PetscViewerASCIIPrintf( txtfile, "\n");
  }

  ierr = PetscViewerDestroy( &txtfile ); CHKERRQ(ierr);
  ierr = SpatialIndexClear( sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexDestroy( sidx ); CHKERRQ(ierr);
  ierr = PetscFree(pt); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SphereQueryTest"
int SphereQueryTest(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);
  Coor lo = {0,0,0};
  Coor hi = {4,4,4};
  PetscReal dx = 0.02;
  Coor dh = {dx,dx,dx};
  SpatialIndex sidx;
  ierr = SpatialIndexCreate( "test", &sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( sidx, lo, hi, dh ); CHKERRQ(ierr);

  PetscRandom rnd;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd,lo.x-2,hi.x+2); CHKERRQ(ierr);
//  ierr = PetscRandomSetInterval(rnd,1.2,1.5); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND48); CHKERRQ(ierr);

  int n = 64*1000;
  Coor *pt;
  int i, j, k;
  Coor center = {1.1, 1.1, 1.1};
  PetscReal r = 0.08;
  const int Np = 10;
  int len;
  Coor *list[Np];

  ierr = PetscMalloc(sizeof(Coor)*n,&pt); CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
    ierr = PetscRandomGetValueReal(rnd,&pt[i].x); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&pt[i].y); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&pt[i].z); CHKERRQ(ierr);
    ierr = SpatialIndexInsertPoint(sidx,pt[i],&pt[i]); CHKERRQ(ierr);
  }

  // loop to test memory leaks
  for (j = 0; j < 100; ++j) {
    printf("j: %d\n",j);
    for (k = 0; k < 100; ++k) {
      ierr = PetscRandomGetValueReal(rnd,&center.x); CHKERRQ(ierr);
      ierr = PetscRandomGetValueReal(rnd,&center.y); CHKERRQ(ierr);
      ierr = PetscRandomGetValueReal(rnd,&center.z); CHKERRQ(ierr);
      ierr = SpatialIndexQueryPoints(sidx, center, r, Np, &len, &list); CHKERRQ(ierr);
    }

    for (i = 0; i < len; ++i) {
      //printf("{%f,%f,%f},", list[i]->x, list[i]->y, list[i]->z);
    }
  }

  ierr = SpatialIndexClear( sidx ); CHKERRQ(ierr);

  ierr = SpatialIndexDestroy( sidx ); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
