#include "FluidField.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldAssembleStaggeredGrid_2D"
PetscErrorCode FluidFieldAssembleStaggeredGrid_2D( PetscReal dx, FluidField f )
{
  DALocalInfo i;  
  DAGetLocalInfo(f->da, &i);
  f->d.x = f->d.y = dx;
  f->l.x = f->d.x * (i.mx-2);
  f->l.y = f->d.y * (i.my-2);
  PetscReal hx = 1/(f->d.x*f->d.x),
            hy = 1/(f->d.y*f->d.y);
  MatStencil row, col[5];
  PetscReal one = 1, val_P[5] = {hx,hx,hy,hy,-2*(hx+hy)},
                     val_D[5] = { 0, 0, 0, 0, 1};
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DAGetMatrix(f->da, MATMPIAIJ,&f->matU); CHKERRQ(ierr);
  ierr = DAGetMatrix(f->da, MATMPIAIJ,&f->matV); CHKERRQ(ierr);
  ierr = DAGetMatrix(f->da, MATMPIAIJ,&f->matP); CHKERRQ(ierr);

  int c;
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      for( c = 0; c < 5; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;
      col[2].j--; col[3].j++;
      
      if( row.i==0 || row.j==0 || row.i>=i.mx-2 || row.j>=i.my-2) {
        
      } else {
          ierr = MatSetValuesStencil(f->matP,1,&row,5,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matU,1,&row,5,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matV,1,&row,5,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  
  int r[4], bd[4] = {i.mx-2, 0, i.my-2, 0};
  PetscReal val[5] = {0,0,0,0,-hx};
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      for( c = 0; c < 5; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;
      col[2].j--; col[3].j++;
      if( col[1].i >= i.mx ) col[1].i = -1;
      if( col[3].j >= i.my ) col[3].j = -1;

      if( row.i==0 || row.j==0 || row.i==i.mx-2 || row.j==i.my-2) {
        ierr = MatSetValuesStencil(f->matU,5,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValuesStencil(f->matV,5,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        r[0] = r[1] = row.i;
        r[2] = r[3] = row.j;
        val[0]=val[1]=val[2]=val[3]=0;
        for ( c = 0; c < 4; ++c)
        {
          if( r[c] == bd[c] )
          {
            val[c] = hx;
            ierr = MatSetValuesStencil(f->matP,1,&row,5,col,val,INSERT_VALUES); CHKERRQ(ierr);
            break;
          }
        }
        if( ( row.i==0      && row.j==0      ) || 
            ( row.i==0      && row.j==i.my-2 ) || 
            ( row.i==i.mx-2 && row.j==0      ) ||
            ( row.i==i.mx-2 && row.j==i.my-2 ) || 
            ( row.i==i.mx-1 && row.j==i.my-2 ) ||
            ( row.i==i.my-2 && row.j==i.my-1 ) )
        {
          ierr = MatSetValuesStencil(f->matP,1,&row,5,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
      }
      
      if( row.i == 1 )
      {
        MatSetValuesStencil(f->matU,1,&row,5,col,val_D,INSERT_VALUES);
        MatSetValuesStencil(f->matU,5,col,1,&row,val_D,INSERT_VALUES);
      }
      if( row.j == 1 )
      {
        MatSetValuesStencil(f->matV,1,&row,5,col,val_D,INSERT_VALUES);
        MatSetValuesStencil(f->matV,5,col,1,&row,val_D,INSERT_VALUES);
      }
      if( row.i == i.mx-1 ) 
      {
        MatSetValuesStencil(f->matU,1,&row,5,col,val_D,INSERT_VALUES); 
        MatSetValuesStencil(f->matV,1,&row,5,col,val_D,INSERT_VALUES);
        MatSetValuesStencil(f->matP,5,col,1,&row,val_D,INSERT_VALUES);
      }
      //TODO: set velocity nodes outside domain to Neuman BC
      if( row.j == i.my-1 )
      {
        MatSetValuesStencil(f->matU,1,&row,5,col,val_D,INSERT_VALUES);
        MatSetValuesStencil(f->matV,1,&row,5,col,val_D,INSERT_VALUES);
        MatSetValuesStencil(f->matP,5,col,1,&row,val_D,INSERT_VALUES);
      }
    }
  }
    
  ierr = MatAssemblyBegin(f->matU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(f->matV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(f->matP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldAssembleStaggeredGrid_3D"
PetscErrorCode FluidFieldAssembleStaggeredGrid_3D( PetscReal dx, FluidField f )
{
  PetscErrorCode ierr;
  DALocalInfo i;
  ierr = DAGetLocalInfo(f->da, &i); CHKERRQ(ierr);
  f->d.x = f->d.y = f->d.z = dx;
  f->l.x = f->d.x * (i.mx-2);
  f->l.y = f->d.y * (i.my-2);
  f->l.z = f->d.z * (i.mz-2);
  int g = 2; // gap from end of index
  int xe = i.mx-g-1,
      ye = i.my-g-1,
      ze = i.mz-g-1;
  MatStencil row, col[7];
  PetscReal hx = 1/(f->d.x*f->d.x),
            hy = 1/(f->d.y*f->d.y),
            hz = 1/(f->d.z*f->d.z);
  PetscReal val_P[7] = {hx,hx,hy,hy,hz,hz,-2*(hx+hy+hz)},
            val_D[7] = { 0, 0, 0, 0, 0, 0, 1},
            val_N[7] = { 0, 0, 0, 0, 0, 0, 0}; //nueman values dynamically set
  
  PetscFunctionBegin;
  ierr = DAGetMatrix(f->da, MATSEQAIJ,&f->matU); CHKERRQ(ierr);
  ierr = DAGetMatrix(f->da, MATSEQAIJ,&f->matV); CHKERRQ(ierr);
  ierr = DAGetMatrix(f->da, MATSEQAIJ,&f->matW); CHKERRQ(ierr);
  ierr = DAGetMatrix(f->da, MATSEQAIJ,&f->matP); CHKERRQ(ierr);

  /* eg, Laplace(u) needs:
   *   i: [3---xe]
   *   j: [2---ye]
   *   k: [2---ze]
   */ 
  
  /* First do standard 7-point stencil within the interior region
   * Then take care of the boundaries in a second loop
   * Done in two loops so boundary conditions replace interior equations
   * This is done in order to construct a symmetric matrix
   */ 
  int c;
  for( row.k = i.zs; row.k < i.zs+i.zm; ++row.k)
  {
    for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
    {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
      {
        for( c = 0; c < 7; ++c)
          PetscMemcpy(&col[c],&row,sizeof(MatStencil));
        col[0].i--; col[1].i++;
        col[2].j--; col[3].j++;
        col[4].k--; col[5].k++;

        if( row.i == i.mx-1 ) col[1].i = -1;
        if( row.j == i.my-1 ) col[3].j = -1;
        if( row.k == i.mz-1 ) col[5].k = -1;
        
        if( row.i<2  || row.j<2  || row.k<2 ||
            row.i>xe || row.j>ye || row.k>ze ) {
          // To avoid corner
          ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
        } else { // [2---xe]
          ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matU,1,&row,7,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matV,1,&row,7,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matW,1,&row,7,col,val_P,INSERT_VALUES); CHKERRQ(ierr);
        }
      }
    }
  }

  /* Now take care of the bounday conditions
   * No-slip for velocity
   * Nueman for pressure
   */ 
  
  for( row.k = i.zs; row.k < i.zs+i.zm; ++row.k)
  {
    for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
    {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
      {
        for( c = 0; c < 7; ++c)
          PetscMemcpy(&col[c],&row,sizeof(MatStencil));
        col[0].i--; col[1].i++;
        col[2].j--; col[3].j++;
        col[4].k--; col[5].k++;
         
        if( row.i == i.mx-1 ) col[1].i = -1;
        if( row.j == i.my-1 ) col[3].j = -1;
        if( row.k == i.mz-1 ) col[5].k = -1;
        
        if( row.i<2  || row.j<2  || row.k<2 ||
            row.i>xe || row.j>ye || row.k>ze ) 
        {
          ierr = MatSetValuesStencil(f->matU,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matU,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matV,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matV,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matW,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matW,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
        if( row.i == 2 || row.i == xe )
        {
          ierr = MatSetValuesStencil(f->matU,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matU,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
        if( row.j == 2 || row.j == ye )
        {
          ierr = MatSetValuesStencil(f->matV,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matV,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
        if( row.k == 2 || row.k == ze )
        {
          ierr = MatSetValuesStencil(f->matW,1,&row,7,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValuesStencil(f->matW,7,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
        /* TODO: there has to be a simpler way to set nuemann bc for pressure
         * Boundary condition for dp / dn == 0;
         */
        val_N[0] = val_N[1] = val_N[2] = val_N[3] = val_N[4] = val_N[5] = val_N[6] = 0;    
        if( 2 <= row.i && row.i <= xe && 
            2 <= row.k && row.k <= ze ) // i,k parallel planes 
        {
          val_N[6] = -hy;
          if( row.j == 1 ) // lower plane
          {
            val_N[3] = hy;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[3] = 0;
          }
          if( row.j == ye+1 ) // upper plane
          {
            val_N[2] = hy;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[2] = 0;
          }
        }
        if( 2 <= row.j && row.j <= ye && // j,k parallel planes 
            2 <= row.k && row.k <= ze )
        {
          val_N[6] = -hx;
          if( row.i == 1 )
          {
            val_N[1] = hx;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[1] = 0;
          }
          if( row.i == xe+1 )
          {
            val_N[0] = hx;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[0] = 0;
          }
        }
        if( 2 <= row.i && row.i <= xe && // i,j parallel plane
            2 <= row.j && row.j <= ye )
        {
          val_N[6] = -hz;
          if( row.k == 1 )
          {
            val_N[5] = hz;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[5] = 0;
          }
          if( row.k == ze+1 )
          {
            val_N[4] = hz;
            ierr = MatSetValuesStencil(f->matP,1,&row,7,col,val_N,INSERT_VALUES); CHKERRQ(ierr);
            val_N[4] = 0;
          }
        }
      }
    }
  }
  
  ierr = MatAssemblyBegin(f->matU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(f->matV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(f->matW,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(f->matP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matW,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(f->matP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}  
  
void AssembleLaplacian_BoundaryCheck2D( Mat mat, MatStencil row, int* lens, 
    PetscReal dx2, PetscReal dy2)
{
  const int S = 5;
  PetscReal val[5] = {-1/dx2, -1/dy2, 
                      -1/dx2, -1/dy2, 
                      2*(1/dx2 + 1/dy2)};
  MatStencil col[7];
  PetscInt *sten, m, n;
  
  for(m = 0; m < S; ++m) {
    col[m].i = row.i;
    col[m].j = row.j;
    col[m].k = row.k;
    col[m].c = row.c;
  }

  col[0].i--; col[1].j--;  
  col[2].i++; col[3].j++;
  
//  printf("\n(%d, %d, %d; %d) %f\n", row.i, row.j, row.k, row.c, D );
  for(m = 0; m < S-1; ++m)     // For each non-central node in the stencil
  {
    sten = (PetscInt*)&col[m]; // Convert MatStencil to an array of ints[]
    for( n = 1; n < 3; ++n)    // For each dimension (k, j, i)
    {
      // Does the node extend beyond the domain?
      if( sten[n] < 0 || lens[n] <= sten[n] ) {
        sten[n] = -1;       // Negative indexes are ignored by MatSetValues
        val[S-1] += val[m]; // Change center node's value to satisfy no-flux BC
        break;              // Outside node accounted for
      }
    }
  }
  
  MatSetValuesStencil(mat, 1, &row, S, col, val, INSERT_VALUES);
}