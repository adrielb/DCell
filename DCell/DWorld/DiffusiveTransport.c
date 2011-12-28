#include "petscda.h"

/* 
 * k = D dt / dh^2
 * cn = -k c + (1+6k) cc
 */
void AssembleLaplacianWithBoundaryCheck( MatStencil row, int* lens, 
    MatStencil *col, PetscReal val[7] );

PetscErrorCode AssembleDiffusionCartesian(DALocalInfo i, Mat J, PetscReal *DD,
    PetscReal dx2, PetscReal dy2, PetscReal dz2 )
{
  //TODO: again, domain properties should be define globally, remove local definitions eg. 'i.mx - 2' to 'M'
  PetscInt lens[3] = {i.mz-1,i.my-1,i.mx-1}, ze, ye, xe;
  PetscReal D;
  MatStencil row, col[7];
  
  PetscErrorCode ierr;
  
  ze = i.zs+i.zm;
  ze = ze == i.mz ? lens[0] : ze;
  if( ze == 0 ) { ze = 1; dz2 = 1; lens[0] = 1; }
  ye = i.ys+i.ym;
  ye = ye == i.my ? lens[1] : ye;
  xe = i.xs+i.xm;
  xe = xe == i.mx ? lens[2] : xe;
  
  for (row.k = i.zs; row.k < ze; ++row.k)
  {
    for (row.j = i.ys; row.j < ye; ++row.j)
    {
      for (row.i = i.xs; row.i < xe; ++row.i)
      {
        for( row.c = 0; row.c < i.dof; ++row.c)
        {
          D = DD[row.c];
          PetscReal val[7] = {D/dx2, D/dy2, 0*D/dz2,
                              D/dx2, D/dy2, 0*D/dz2,
                        -2*D*(1/dx2 +1/dy2 +0*1/dz2)};
          
          AssembleLaplacianWithBoundaryCheck(row,lens,col, val);
          
          ierr = MatSetValuesStencil(J, 1, &row, 7, col, val, ADD_VALUES); CHKERRQ(ierr);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

void AssembleLaplacianWithBoundaryCheck( MatStencil row, int* lens, 
    MatStencil *col, PetscReal val[7] )
{
  const int  S = 7;
  PetscInt *sten, m, n;
  
  for(m = 0; m < S; ++m) {
    col[m].i = row.i;
    col[m].j = row.j;
    col[m].k = row.k;
    col[m].c = row.c;
  }

  col[0].i--; col[1].j--; col[2].k--;   
  col[3].i++; col[4].j++; col[5].k++;
  
//  printf("\n(%d, %d, %d; %d) %f\n", row.i, row.j, row.k, row.c, D );
  for(m = 0; m < S-1; ++m)     // For each non-central node in the stencil
  {
    sten = (PetscInt*)&col[m]; // Convert MatStencil to an array of ints[]
    for( n = 0; n < 3; ++n)    // For each dimension (k, j, i)
    {
      // Does the node extend beyond the domain?
      if( sten[n] < 0 || lens[n] <= sten[n] ) {
        sten[n] = -1;       // Negative indexes are ignored by MatSetValues
        val[S-1] += val[m]; // Change center node's value to satisfy no-flux BC
        break;              // Outside node accounted for
      }
    }
  }
}

void AssembleLaplacianCartesianBlockInterior( Mat mat, PetscReal *D, DALocalInfo *info, \
    PetscReal lenX, PetscReal lenY, PetscReal lenZ )
{
  DALocalInfo i;
  PetscMemcpy(&i,info,sizeof(DALocalInfo));
  MatStencil row, col[7];
  int l;
  PetscReal dx, idx2, dy, idy2, dz, idz2;
  dx = lenX / (i.mx-1);
  idx2= 1/(dx*dx);
  dy = lenY / (i.my-1);
  idy2= 1/(dy*dy);
  dz = lenZ / (i.mz-1);
  idz2= 1/(dz*dz);
  PetscReal val[7];
  const PetscReal val_dh[7] = {
      idz2, idy2, idx2, 
      -2*(idx2+idy2+idz2),
      idx2, idy2, idz2};
  
  
  for (row.k = i.zs+1; row.k < i.zs+i.zm-1; ++row.k)
  {
    for (row.j = i.ys+1; row.j < i.ys+i.ym-1; ++row.j)
    {
      for (row.i = i.xs+1; row.i < i.xs+i.xm-1; ++row.i)
      {
        for( row.c = 0; row.c < i.dof; ++row.c)
        {
          for (l = 0; l < 7; ++l)
          {
            col[l].i = row.i;
            col[l].j = row.j;
            col[l].k = row.k;
            col[l].c = row.c;
            val[l] = D[row.c] * val_dh[l]; 
          }
          col[0].k--; col[1].j--; col[2].i--;
          col[4].i++; col[5].j++; col[6].k++;
          MatSetValuesStencil(mat, 1, &row, 7, col, val, INSERT_VALUES);
        }
      }
    }
  }

  /* Interface 
   * 
  PetscReal numNei=0;
  IrregularNode *n;
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i );
    
    row.i = n->x;
    row.j = n->y;
    row.k = n->z;
    
    for (l = 0; l < 7; ++l)
    {
      col[l].i = row.i;
      col[l].j = row.j;
      col[l].k = row.k;
    }
    
    col[0].i--; col[1].j--; col[2].k--;   
    col[3].i++; col[4].j++; col[5].k++;   // DIFFERENT COL ORDERING!!! (keep in mind for val[i])

    numNei = 0;
    for( l = 0; l < 6; ++l)
    {
      if( g->v3[col[l].i][col[l].j][col[l].k] * n->sign < 0. ) {
        val[l] = 0;
      } else {
        val[l] = -K;
        numNei++;
      }
    }
    val[6] = 1 + numNei * K ;
    MatSetValuesStencil(mat, 1, &row, 7, col, val, INSERT_VALUES );
  }
  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat,   MAT_FINAL_ASSEMBLY);
  */
}

void AssembleDiffusionBoundary()
{
/*
 * 
 * PetscInt lens[3] = {0,my,mx};
  PetscReal D;
PetscReal dx2 = dx*dx, dy2=dy*dy;
  for( row.i = xs; row.i < xe; ++row.i)
  {
    for( row.c = 0; row.c < dof; ++row.c )
    {
      D = rxn->D[row.c];
      row.j = 0;
      AssembleDiffusion2DBC( J, row, lens, D, dx2, dy2);
      row.j = ye-1;
      AssembleDiffusion2DBC( J, row, lens, D, dx2, dy2);
    }
  }
  for( row.j = ys; row.j < ye; ++row.j)
  {
    for( row.c = 0; row.c < dof; row.c++ )
    {
      D = rxn->D[row.c];
      row.i = 0;
      AssembleDiffusion2DBC( J, row, lens, D, dx2, dy2);
      row.i = xe-1;
      AssembleDiffusion2DBC( J, row, lens, D, dx2, dy2);
    }
  }
  
  
    
  int m, i, j, k;
  PetscInt dims[3] = {d3,d2,d1}; //Matching MatStencil[k,j,i]

  for( j = 0; j < d2; ++j)
    for( k = 0; k < d3; ++k)
    {
      DCellAssembleDiffusion_BC(mat, 0,    j, k, dims, K );
      DCellAssembleDiffusion_BC(mat, d1-1, j, k, dims, K );
    }

  for( i = 0; i < d2; ++i)
    for( k = 0; k < d3; ++k)
    {
      DCellAssembleDiffusion_BC(mat, i, 0,    k, dims, K );
      DCellAssembleDiffusion_BC(mat, i, d2-1, k, dims, K );
    }
      
  for( i = 0; i < d2; ++i)
    for( j = 0; j < d2; ++j)
    {
      DCellAssembleDiffusion_BC(mat, i, j, 0,    dims, K );
      DCellAssembleDiffusion_BC(mat, i, j, d3-1, dims, K );
    }
    */
}


void DCellAssembleDiffusion_BC( Mat mat, int i, int j, int k, PetscInt *dim, PetscReal K )
{
  PetscReal val[7];
  MatStencil row, col[7];
  PetscReal numNei = 6;
  
  row.i = i;  row.j = j;  row.k = k;
  for (int m = 0; m < 7; ++m) {
    col[m].i = i;
    col[m].j = j;
    col[m].k = k;
  }
  
  col[0].i--; col[1].j--; col[2].k--;   
  col[3].i++; col[4].j++; col[5].k++;
  
  PetscInt *sten, n;
  for (int m = 0; m < 6; ++m)
  {
    val[m]   = -K;
    sten = (PetscInt*)&col[m];
    for( n = 0; n < 3; ++n)
    {
      if( 0 > sten[n] || sten[n] >= dim[n] ) {
        sten[n] = -1;
        val[m]  = 0;
        numNei--;
        break;
      }
    }
  }
  val[6] = numNei * K + 1;
  MatSetValuesStencil(mat, 1, &row, 7, col, val, INSERT_VALUES);
}

/*

void f(LevelSet3D ls, Mat mat, Mat submat)
{
  PetscInt *localToGlobal, *globalToLocal;
  PetscNew(ls->g3d->len*sizeof(PetscInt),&localToGlobal);
  PetscNew(ls->g3d->len*sizeof(PetscInt),&globalToLocal);
  
  PetscInt m = 0;
  for( int i = 0; i < ls->g3d->len; i++ )
  {
    if( ls->g3d->v1[i] < 0 )
    {
      localToGlobal[m] = i;
      globalToLocal[i] = m;
      m++;
    }
  }
  
  IS is;
  ISCreateGeneralNC(PETSC_COMM_SELF,m,idx,&is);
  MatGetSubMatrix(mat, is, is, -1, MAT_INITIAL_MATRIX, &submat );
  PetscReal val[7]={0,0,0,1,0,0,0};
  PetscInt col[7], row;
  
  IrregularNode *n;
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i );
    if( n->sign > 0. ) continue;
    row = From3Dto1D( n->x, n->y, n->z);
    col[0] = From3Dto1D( n->x-1, n->y, n->z);
    col[1] = From3Dto1D( n->x, n->y-1, n->z);
    col[2] = From3Dto1D( n->x, n->y, n->z-1);
    col[3] = row;
    col[4] = From3Dto1D( n->x+1, n->y, n->z);
    col[5] = From3Dto1D( n->x, n->y+1, n->z);
    col[6] = From3Dto1D( n->x, n->y, n->z+1);
    MatSetValues(submat,1,&row,7,col,val,INSERT_VALUES);
  }
}

*/