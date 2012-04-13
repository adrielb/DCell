#include "LevelSetMethod.h"

PetscErrorCode CurvatureFlowVn( LevelSet ls );
PetscErrorCode PrintICoorArray( Array a, int i);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  int tend = 300;
  iCoor CELL_CENTER={0,0,0};
  int d1 = 64;
  Coor d = {1./(d1-1), 1./(d1-1), 1./(d1-1)};
  iCoor size = {d1,d1,0};
  
  Grid c;
  ierr = GridCreate(size,&c); CHKERRQ(ierr);

  LevelSet ls;
  ierr = LevelSetCreate( size, &ls); CHKERRQ(ierr);
  ierr = GridSetDx(ls->g,d); CHKERRQ(ierr);
  ierr = LevelSetInitializeToBall(ls); CHKERRQ(ierr);
//ierr = PrintICoorArray(ls->band,0); CHKERRQ(ierr);
  for (int i = 0; i < tend; ++i)
  {
    printf("i: %d\n", i);
    ierr = VecWrite(ls->g->v, "ls",i); CHKERRQ(ierr);
    ierr = CurvatureWrite(ls,c,i); CHKERRQ(ierr);
    ierr = IrregularNodeListUpdate(CELL_CENTER,ls); CHKERRQ(ierr);
    ierr = CurvatureFlowVn(ls); CHKERRQ(ierr);
    ierr = IrregularNodeListWrite(ls->irregularNodes,i); CHKERRQ(ierr);
    ierr = LevelSetAdvectVn(ls); CHKERRQ(ierr);
    ierr = VecWrite(ls->Vn->v,"Vn",i); CHKERRQ(ierr);
//ierr = PrintICoorArray(ls->band,i); CHKERRQ(ierr);
  }


  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CurvatureWrite"
PetscErrorCode CurvatureWrite( LevelSet ls, Grid c, int t )
{
  int x,y;
  iCoor *b;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(c->v,0.); CHKERRQ(ierr);
  for (int j = 0; j < ArrayLength(ls->band); ++j) {
    ierr = ArrayGet(ls->band,j,(void*)&b); CHKERRQ(ierr);
    x = b->x;
    y = b->y;
    if( PetscAbs(ls->g->v2[y][x]) < 3 )
    {
//      c->v2[b->y][b->x] = GridFunction2D_Curv(ls->g->v2,b->y,b->x,d);
      c->v2[y][x] = (ls->g->v2[y][x-1] - 2.*ls->g->v2[y][x] + ls->g->v2[y][x+1]);
    }
  }
  ierr = VecWrite(c->v,"curv",t); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CurvatureFlowVn"
PetscErrorCode CurvatureFlowVn( LevelSet ls )
{
  IrregularNode *n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  for (int i = 0; i < ArrayLength(ls->irregularNodes); ++i)
  {
    ierr = ArrayGet(ls->irregularNodes, i,(void*)&n); CHKERRQ(ierr);
    double X = n->x + n->op.x,
           Y = n->y + n->op.y,
           Z = n->z + n->op.z;
    
    if( ls->g->is2D )
    {
      n->Vn = -Bilinear2D( GridFunction2D_Curv, ls->g, X, Y );
      n->Vn = n->Vn + 4;

      PetscReal nx = Bilinear2D(GridFunction2D_DerivX, ls->g, n->x+n->ox, n->y+n->oy);
      PetscReal ny = Bilinear2D(GridFunction2D_DerivY, ls->g, n->x+n->ox, n->y+n->oy);
      PetscReal h = sqrt( nx*nx + ny*ny );
      n->nx = nx / h;
      n->ny = ny / h;
      n->Vn = -n->ny;
    } else {
      n->Vn = -Bilinear3D( GridFunction3D_Curv, ls->g, X, Y, Z );
//      n->Vn = -1;
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PrintICoorArray"
PetscErrorCode PrintICoorArray( Array a, int i)
{
  char filename[256];
  iCoor *coor;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  sprintf(filename, "/home/abergman/Research/DCell/temp/a.%d",i);
  FILE *f =  f = fopen(filename, "w");
  
  ierr = ArrayGet(a,0,(void*)&coor); CHKERRQ(ierr);
  //Mathematica ordering starts at 1,..,N rather than C's 0,...,N-1
  fprintf(f,"{{%d,%d,%d}", coor->x, coor->y, coor->z);
  for( int i = 1; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,(void*)&coor); CHKERRQ(ierr);
    fprintf(f,",{%d,%d,%d}", coor->x, coor->y, coor->z);
  }
  fprintf(f,"}");
  
  fclose(f);
  
  PetscFunctionReturn(0);
}
