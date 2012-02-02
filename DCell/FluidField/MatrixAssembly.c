#include "FluidField.h"
#include "FluidField_private.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldMatAssemble"
PetscErrorCode FluidFieldMatAssemble( FluidField f )
{
  int dof;
  iCoor dims = f->dims;
  PetscLogDouble t1,t2;
  int size;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  MPI_Comm_size(f->comm, &size);
  if( f->is3D ) {
    dof = 4; // [u v w p]
    ierr = DMDACreate3d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,
              dims.x,dims.y,dims.z,  PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
              dof,1,  0,0,0, &f->daV); CHKERRQ(ierr);

    dof = 6; // [xx xy xz yy yz zz]
    ierr = DMDACreate3d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,
              dims.x,dims.y,dims.z, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
              dof,1,  0,0,0, &f->daE); CHKERRQ(ierr);

    dof = 1; // buf
    ierr = DMDACreate3d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,
              dims.x,dims.y,dims.z, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
              dof,1,  0,0,0, &f->daB); CHKERRQ(ierr);
//    ierr = FluidFieldMatAssemble_3D( f->daV, f->dh, &f->mat); CHKERRQ(ierr);
  } else {
    dof = 3; // [u v p]
    ierr = DMDACreate2d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,
              dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &f->daV); CHKERRQ(ierr);

    dof = 3; // [xx xy yy]
    ierr = DMDACreate2d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,
              dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &f->daE); CHKERRQ(ierr);

    dof = 1; // buf
    ierr = DMDACreate2d(f->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,
              dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &f->daB); CHKERRQ(ierr);
  }
  if( size == 1 ) {
    ierr = DMGetMatrix(f->daV, MATSEQAIJ, &f->mat); CHKERRQ(ierr);
  } else {
    ierr = DMGetMatrix(f->daV, MATMPIAIJ, &f->mat); CHKERRQ(ierr);
  }
  ierr = MatSetFromOptions(f->mat); CHKERRQ(ierr);

  ierr = FluidField_MatAssemble( f->mu, f->dirichletBC, f->daV, f->dh, f->mat); CHKERRQ(ierr);

  // Apply masking to rectangular domain
  ierr = FluidFieldMaskDomain(f); CHKERRQ(ierr);

//  TODO: need better way to incorporate a pressure BC, perhaps a separate mask?
//  ierr = FluidField_PressureBC(f); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(f->mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished Assembly: %f sec\n", t2-t1); CHKERRQ(ierr);

//  ierr = MatWrite(f->mat,"mat",0); CHKERRQ(ierr);

  ierr = DMSetApplicationContext(f->daV,f); CHKERRQ(ierr);
  ierr = DMSetJacobian(f->daV,FluidField_ComputeMatrix);CHKERRQ(ierr);
  ierr = DMSetInitialGuess(da,ComputeInitialGuess);CHKERRQ(ierr);
  ierr = DMSetFunction(da,ComputeRHS);CHKERRQ(ierr);
  ierr = KSPSetDM(f->ksp,f->daV);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidField_ComputeMatrix"
PetscErrorCode FluidField_ComputeMatrix(DM dm,Vec x,Mat jac,Mat B,MatStructure *stflg)
{
  FluidField fluid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm, &fluid); CHKERRQ(ierr);
  ierr = FluidField_MatAssemble( fluid->mu, fluid->dirichletBC, fluid->daV, fluid->dh, B ); CHKERRQ(ierr);
  *stflg = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}

/* -mu*Laplace(u) + Grad(p) == mu*C{Laplace(u)} - C{Grad(p)}
 *         Div(u)           ==       -C{Div(u)}
 */

#undef __FUNCT__
#define __FUNCT__ "FluidField_MatAssemble"
PetscErrorCode FluidField_MatAssemble( PetscReal mu, Array dbc, DM da, Coor dH, Mat mat )
{
  int c;
  DMDALocalInfo i;
  int n;
  int numdiv = 6;
  int zstr, zend, zenddiv, uzstr;
  VelFace endface;
  Coor H = {1./dH.x, 1./dH.y, 1./dH.z};
  PetscReal *dh = &H.x - U_FACE;
  PetscReal hx2 = mu*H.x*H.x,
            hy2 = mu*H.y*H.y,
            hz2 = mu*H.z*H.z;
  MatStencil row, col[9];// *ucol=&col[2];
  PetscReal val[9] = {
      0, 0,
      2*(hx2 + hy2 + hz2),
      -hx2, -hx2,
      -hy2, -hy2,
      -hz2, -hz2};
  PetscReal val_D[9] = { 0, 0, 1,
                         0, 0, 0,
                         0, 0, 0};
  PetscReal val_div[6] = {H.x,-H.x,
                          H.y,-H.y,
                          H.z,-H.z};
  int end[][3] = {{1,2,2}, // u
                  {2,1,2}, // v
                  {2,2,1}};// w
  int op[][3] = {{1,0,0},  // px
                 {0,1,0},  // py
                 {0,0,1}}; // pz

  PetscErrorCode ierr;
  PetscFunctionBegin;

  DMDAGetLocalInfo(da, &i);

  if( i.dim == 2 ) {
    endface = V_FACE;
    n    = 7;
    zstr = 0;
    zend = 1;
    zenddiv = 1;
    numdiv = 4;
    uzstr = -1;
    end[0][2] = 0;
    end[1][2] = 0;
    val[2] = 2*(hx2 + hy2);
  } else {
    endface = W_FACE;
    n = 9;
    zstr = i.zs;
    zend = i.zs+i.zm;
    zenddiv = i.mz-1;
    numdiv = 6;
    uzstr = 0;
  }

  for(     row.k = zstr; row.k < zend;      ++row.k ) {
    for(   row.j = i.ys; row.j < i.ys+i.ym; ++row.j ) {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i ) {
        for( c = 0; c < n; ++c)
          col[c] = row;
        col[3].i--; col[4].i++;  // [-hx, hx]
        col[5].j--; col[6].j++;  // [-hy, hy]
        col[7].k--; col[8].k++;  // [-hz, hz]
        if( col[4].i >= i.mx ) col[4].i = -1;
        if( col[6].j >= i.my ) col[6].j = -1;
        if( col[8].k >= i.mz ) col[8].k = -1;

        // -mu*Laplace([u,v,w]) + grad([px,py,pz])
        for( row.c = U_FACE; row.c <= endface; row.c++ )  // [u,v,w]
        {
          for( c = 0; c < n; ++c)
            col[c].c = row.c;

          if( 0 < row.i && row.i < i.mx - end[row.c-U_FACE][0] &&
              0 < row.j && row.j < i.my - end[row.c-U_FACE][1] &&
          uzstr < row.k && row.k < i.mz - end[row.c-U_FACE][2] ) {

            // px[i,j,k] = ( p[i,j,k] - p[i-1,j,k] ) / hx
            val[0] =  dh[row.c];
            val[1] = -dh[row.c];
            col[0].c = CELL_CENTER;
            col[1].c = CELL_CENTER;
            col[1].i = row.i - op[row.c-U_FACE][0];
            col[1].j = row.j - op[row.c-U_FACE][1];
            col[1].k = row.k - op[row.c-U_FACE][2];
            ierr = MatSetValuesStencil(mat,1,&row,n,col,val,INSERT_VALUES); CHKERRQ(ierr);
          } else {
            ierr = MatSetValuesStencil(mat,1,&row,n,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
            ierr = FluidField_AppendDBC(dbc,row); CHKERRQ(ierr);
          }
        } // for -L(u,v,w) + px

        // Pressure rows
        row.c = CELL_CENTER;
        for( c = 0; c < n; ++c)
          col[c] = row;

        if( row.j == 0 && (row.i == 0 || row.i == i.mx-2) ) {
          //Neumann BC on lower corners
          val[0] = -1;
          val[1] =  1;

          col[0].i = row.i;
          col[0].j = row.j+1;
          col[0].k = row.k;
          col[1].i = row.i;
          col[1].j = row.j;
          col[1].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,2,col,val,INSERT_VALUES); CHKERRQ(ierr);
        } else if( row.j == i.my-2 && (row.i == 0 || row.i == i.mx-2) ) {
          //Neumann BC on upper corners
          val[0] =  1;
          val[1] = -1;

          col[0].i = row.i;
          col[0].j = row.j;
          col[0].k = row.k;
          col[1].i = row.i;
          col[1].j = row.j-1;
          col[1].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,2,col,val,INSERT_VALUES); CHKERRQ(ierr);
        } else if( 0 <= row.i && row.i < i.mx-1 &&
                   0 <= row.j && row.j < i.my-1 &&
                   0 <= row.k && row.k < zenddiv ) {
          //Divergence equation for interior nodes
          col[0].c = U_FACE;
          col[0].i = row.i+1;
          col[1].c = U_FACE;
          col[1].i = row.i;
          col[2].c = V_FACE;
          col[2].j = row.j+1;
          col[3].c = V_FACE;
          col[3].j = row.j;
          col[4].c = W_FACE;
          col[4].k = row.k+1;
          col[5].c = W_FACE;
          col[5].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,numdiv,col,val_div,INSERT_VALUES); CHKERRQ(ierr);
        } else if( row.i == i.mx-1 ||
                   row.j == i.my-1 ||
                   row.k == zenddiv ) {
          //Overhanging pressure nodes on upper boundaries
          PetscReal one = 1.;
          ierr = MatSetValuesStencil(mat,1,&row,1,&row,&one,INSERT_VALUES); CHKERRQ(ierr);
          ierr = FluidField_AppendDBC(dbc,row); CHKERRQ(ierr);
        } // P
      } // row.i
    } // row.j
    PetscPrintf(PETSC_COMM_WORLD,"*");
  } // row.k
  PetscPrintf(PETSC_COMM_WORLD,"\n");

  ierr = MatAssemblyBegin(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  /*
   *  Make velocity block symmetric by eliminating BCs
   */
  for(     row.k = zstr; row.k < zend;      ++row.k) {
    for(   row.j = i.ys; row.j < i.ys+i.ym; ++row.j) {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i) {
        for( c = 0; c < n; ++c)
          col[c] = row;
        col[3].i--; col[4].i++;  // [-hx, hx]
        col[5].j--; col[6].j++;  // [-hy, hy]
        col[7].k--; col[8].k++;  // [-hz, hz]
        if( col[4].i >= i.mx ) col[4].i = -1;
        if( col[6].j >= i.my ) col[6].j = -1;
        if( col[8].k >= i.mz ) col[8].k = -1;

        for( row.c = U_FACE; row.c <= endface; row.c++ )  // [u,v,w]
        {
          for( c = 0; c < n; ++c)
            col[c].c = row.c;

          if( 0 < row.i && row.i < i.mx - end[row.c-U_FACE][0] &&
              0 < row.j && row.j < i.my - end[row.c-U_FACE][1] &&
          uzstr < row.k && row.k < i.mz - end[row.c-U_FACE][2] ) {
          } else {
            ierr = MatSetValuesStencil(mat,n,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
          } // if BC
        } // row.c
      } // row.i
    } // row.j
    ierr = PetscPrintf(PETSC_COMM_WORLD,"."); CHKERRQ(ierr);
  } // row.k
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
