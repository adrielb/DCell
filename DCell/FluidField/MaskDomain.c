#include "FluidField.h"
#include "FluidField_private.h"

void CreateStencilTemplate(DMDALocalInfo i, MatStencil row, MatStencil *col);

#undef __FUNCT__
#define __FUNCT__ "FluidFieldMaskDomain"
PetscErrorCode FluidFieldMaskDomain(FluidField f)
{
  DMDALocalInfo i;
  MatStencil row;
  //                 dp,dp, D,hx,hx,hy,hy,hz,hz
  PetscReal val[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0}; // [ px, L(u) ]
  MatStencil col[9];
  int len;
  PetscReal *valL  = &val[2];
  MatStencil *colL = &col[2];
  int lenL;
  int c;
  int zstr, zend;
  PetscReal **mask2D, ***mask3D, mp, mu, mv, mw;
  const PetscReal M = 1; // default positive value (outside domain)
  PetscErrorCode ierr=0;

  PetscFunctionBegin;
  ierr = PetscInfo(0,"Applying mask to fluid field\n"); CHKERRQ(ierr);
  if( !f->mask ) PetscFunctionReturn(0);
  ierr = DMDAGetLocalInfo(f->daB, &i); CHKERRQ(ierr);
  ierr = GridGet(f->mask, &mask2D); CHKERRQ(ierr);
  ierr = GridGet(f->mask, &mask3D); CHKERRQ(ierr);

  if( f->is3D )
  {
    len  = 9;
    lenL = 7;
    zstr = i.zs;
    zend = i.zs+i.zm;
  } else {
    len  = 7;
    lenL = 5;
    zstr = 0;
    zend = 1;
  }

  for(     row.k = zstr; row.k <      zend; ++row.k) {
    for(   row.j = i.ys; row.j < i.ys+i.ym; ++row.j) {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i) {
        if( f->is3D ) {
          mp = mask3D[row.k][row.j][row.i];
          mu = row.i == 0 ? M : mask3D[row.k][row.j][row.i-1];
          mv = row.j == 0 ? M : mask3D[row.k][row.j-1][row.i];
          mw = row.k == 0 ? M : mask3D[row.k-1][row.j][row.i];
        } else {
          mp = mask2D[row.j][row.i];
          mu = row.i == 0 ? M : mask2D[row.j][row.i-1];
          mv = row.j == 0 ? M : mask2D[row.j-1][row.i];
          mw = 1;
        }

        // m? < 0 : inside rigid part of domain
        // Set u[i,j,k], v[i,j,k], w[i,j,k], u[i+1,j,k], v[i,j+1,k], w[i,j,k+1]  = 0

        // u[i,j] = 0
        if( mp > 0 || mu > 0 ) {
          row.c = U_FACE;
          CreateStencilTemplate(i,row,(MatStencil*)&col);
          col[1].i--; // pressure gradient cols: p[i,j] and p[i-1,j]
          ierr = MatSetValuesStencil(f->mat,1,&row,len,col,val,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->mat,lenL,colL,1,&row,valL,INSERT_VALUES); CHKERRQ(ierr);
          ierr = FluidField_AppendDBC(f->dirichletBC,row); CHKERRQ(ierr);
        }

        // v[i,j] = 0
        if( mp > 0 || mv > 0 ) {
          row.c = V_FACE;
          CreateStencilTemplate(i,row,(MatStencil*)&col);
          col[1].j--; // pressure gradient cols: p[i,j] and p[i,j-1]
          ierr = MatSetValuesStencil(f->mat,1,&row,len,col,val,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->mat,lenL,colL,1,&row,valL,INSERT_VALUES); CHKERRQ(ierr);
          ierr = FluidField_AppendDBC(f->dirichletBC,row); CHKERRQ(ierr);
        }

        // w[i,j] = 0
        if( f->is3D && (mp > 0 || mw > 0) ) {
          row.c = W_FACE;
          CreateStencilTemplate(i,row,(MatStencil*)&col);
          col[1].k--; // pressure gradient cols: p[i,j,k] and p[i,j,k-1]
          ierr = MatSetValuesStencil(f->mat,1,&row,len,col,val,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->mat,lenL,colL,1,&row,valL,INSERT_VALUES); CHKERRQ(ierr);
          ierr = FluidField_AppendDBC(f->dirichletBC,row); CHKERRQ(ierr);
        }

        // p[i,j] = 0
        if( mp > 0 ) {
          row.c = CELL_CENTER;
          for( c = 0; c < len; ++c)
            col[c] = row;

          col[0].c = CELL_CENTER;                // p[i,j,k]   = 1
          col[1].c = U_FACE; col[1].i = row.i+1; // u[i+1,j,k] = 0
          col[2].c = U_FACE;                     // u[i,j,k]   = 0
          col[3].c = V_FACE; col[3].j = row.j+1; // v[i,j+1,k] = 0
          col[4].c = V_FACE;                     // v[i,j,k]   = 0
          col[5].c = W_FACE; col[5].k = row.k+1; // v[i,j,k+1] = 0
          col[6].c = W_FACE;                     // v[i,j,k]   = 0
          if( col[1].i >= i.mx ) col[1].i = -1;
          if( col[3].j >= i.my ) col[3].j = -1;
          if( col[5].k >= i.mz ) col[5].k = -1;
          ierr = MatSetValuesStencil(f->mat,1,&row,lenL,col,valL,INSERT_VALUES); CHKERRQ(ierr);
//          TODO: Zero pressure RHS????
//          ierr = FluidField_AppendDBC(f->dirichletBC,row); CHKERRQ(ierr);
        }
      } // row.i
    } // row.j
  } // row.k
  ierr = PetscInfo(0,"Finished mask application to fluid field\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished Masking\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidField_PressureBC"
PetscErrorCode FluidField_PressureBC(FluidField f)
{
  DMDALocalInfo i;
  MatStencil row;
  //                  P,ux,ux,vy,vy,wx,wy
  PetscReal val[] = { 1, 0, 0, 0, 0, 0, 0};
  MatStencil col[9];
  int len=9;
  int c;
  int zstr, zend, yend, xend;
  PetscReal **mask2D, ***mask3D, m;
  PetscErrorCode ierr=0;

  PetscFunctionBegin;
  if( !f->mask ) PetscFunctionReturn(0);
  ierr = DMDAGetLocalInfo(f->daB, &i); CHKERRQ(ierr);
  ierr = GridGet(f->mask, &mask2D); CHKERRQ(ierr);
  ierr = GridGet(f->mask, &mask3D); CHKERRQ(ierr);

  if( f->is3D )
  {
    len  = 9;
    zstr = i.zs;
    zend = i.zs+i.zm;
    if( zend == i.mz ) zend--;
  } else {
    len  = 7;
    zstr = 0;
    zend = 1;
  }
  yend = i.ys+i.ym;
  if( yend == i.my ) yend--;
  xend = i.xs+i.xm;
  if( xend == i.mx ) xend--;

  for(     row.k = zstr; row.k < zend; ++row.k) {
    for(   row.j = i.ys; row.j < yend; ++row.j) {
      for( row.i = i.xs; row.i < xend; ++row.i) {
        m = f->is3D ? mask3D[row.k][row.j][row.i] : mask2D[row.j][row.i];
        if( m > -1.5 ) continue;

        // p[i,j] = p(x,y)
        row.c = CELL_CENTER;
        for( c = 0; c < len; ++c)
          col[c] = row;

        col[0].c = CELL_CENTER;                // p[i,j,k]   = 1
        col[1].c = U_FACE; col[1].i = row.i+1; // u[i+1,j,k] = 0
        col[2].c = U_FACE;                     // u[i,j,k]   = 0
        col[3].c = V_FACE; col[3].j = row.j+1; // v[i,j+1,k] = 0
        col[4].c = V_FACE;                     // v[i,j,k]   = 0
        col[5].c = W_FACE; col[5].k = row.k+1; // v[i,j,k+1] = 0
        col[6].c = W_FACE;                     // v[i,j,k]   = 0
        ierr = MatSetValuesStencil(f->mat,1,&row,len,col,val,INSERT_VALUES); CHKERRQ(ierr);
      } // row.i
    } // row.j
  } // row.k
  PetscFunctionReturn(0);
}

void CreateStencilTemplate(DMDALocalInfo i, MatStencil row, MatStencil *col)
{
  int c;
  const int n = 9;

  for( c = 0; c < n; ++c)
    col[c] = row;

  // laplace stencil
  col[3].i--; col[4].i++;  // [-hx, hx]
  col[5].j--; col[6].j++;  // [-hy, hy]
  col[7].k--; col[8].k++;  // [-hz, hz]
  if( col[4].i >= i.mx ) col[4].i = -1;
  if( col[6].j >= i.my ) col[6].j = -1;
  if( col[8].k >= i.mz ) col[8].k = -1;

  // 2-point gradient stencil
  col[0].c = CELL_CENTER;
  col[1].c = CELL_CENTER;
}

#undef __FUNCT__
#define __FUNCT__ "FluidField_EnforceNoSlipBC"
PetscErrorCode FluidField_EnforceNoSlipBC( FluidField f )
{
  int i;
  int len = ArrayLength(f->dirichletBC);
  MatStencil *dbc;
  PetscReal *rhs;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMDAVecGetArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  dbc = ArrayGetData(f->dirichletBC);
  if (f->is3D) {
    PetscReal ****rhs3D = (PetscReal****)rhs;
    for ( i = 0; i < len; ++i) {
      rhs3D[dbc[i].k][dbc[i].j][dbc[i].i][dbc[i].c] = 0;
    }
  } else {
    PetscReal ***rhs2D = (PetscReal***)rhs;
    for ( i = 0; i < len; ++i) {
      rhs2D[dbc[i].j][dbc[i].i][dbc[i].c] = 0;
    }
  }
  ierr = DMDAVecRestoreArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidField_AppendDBC"
PetscErrorCode FluidField_AppendDBC( Array dirichletBC, MatStencil row )
{
  MatStencil *dbc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend(dirichletBC,(void*)&dbc); CHKERRQ(ierr);
  dbc->i = row.i;
  dbc->j = row.j;
  dbc->k = row.k;
  dbc->c = row.c;
  PetscFunctionReturn(0);
}

/* checking file existence
 * parallel loading into vec
#include "sys/stat.h"
  char file[64] = "/tmp/mask"; //TODO: remove hardcoded mask file location
  PetscViewer fd;
  struct stat statbuf;
  ierr = stat(file, &statbuf);
  if( ierr ) {
    ierr = 0;
    PetscFunctionReturn(0);
  }

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd); CHKERRQ(ierr);

  ierr = VecLoad(fd,VECMPI,&f->buf); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(fd); CHKERRQ(ierr);

  ierr = DAVecGetArray(f->daB,f->buf,&mask); CHKERRQ(ierr);
  */
