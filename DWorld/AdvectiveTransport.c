#include "DCell.h"

/*           .-----vij+1---.
 *           |      |      |
 * pij-1----uij----pij----uij+1----pij+1
 *           |      |      |
 *           .-----vij-----.
 */

#undef __FUNCT__
#define __FUNCT__ "AssembleAdvectiveTransport"
PetscInt EVENT_AssembleAdvectiveTransport;
PetscErrorCode AssembleAdvectiveTransport( DALocalInfo info, Mat mat, Coor d, PetscReal **u, PetscReal **v )
{
  int i,j,s,S=5;
  int xe, ye, M, N;
  PetscReal xx = 2 * d.x, yy = 2 * d.y;
  PetscReal U, V;
  PetscReal val[5];
  MatStencil row, col[5];
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleAdvectiveTransport,0,0,0,0);
//  PetscLogEventRegister(&EVENT_AssembleAdvectiveTransport,"AssembleAdvectiveTransport", 0);
  // Cell-centered boundary extends by one grid point beyond domain
  xe = info.xs+info.xm;
  ye = info.ys+info.ym;
  xe = xe == info.mx ? xe - 1 : xe;
  ye = ye == info.my ? ye - 1 : ye;
  M  = info.mx - 1;
  N  = info.my - 1;
  
  for (row.j = info.ys; row.j < ye; ++row.j)
  {
    for (row.i = info.xs; row.i < xe; ++row.i)
    {
      for( s = 0; s < S; s++ )
        PetscMemcpy(&col[s],&row,sizeof(MatStencil));
//TODO: rearrange column indexing to be sorted       
      col[0].i--; col[1].i++;
      col[2].j--; col[3].j++;
      
      /*    3
       *  0-4-1
       *    2
       */
      
      i = row.i; j = row.j;
      
      //No-flux BC
      /* second order flux form
      val[0] = i == 0 ? 0 :    u[i][j] / xx;
      val[1] = i == M ? 0 : -u[i+1][j] / xx;
      val[2] = j == 0 ? 0 :    v[i][j] / yy;
      val[3] = j == N ? 0 : -v[i][j+1] / yy;
      val[4] = ( u[i][j] - u[i+1][j] ) / xx + 
               ( v[i][j] - v[i][j+1] ) / yy;
      */
      
      for( s = 0; s < S; s++ )
        val[s] = 0;
      
      //No-flux BC

      U = u[j][i] / d.x;
      if( U > 0. )  
        val[0] = i == 0 ? 0 : U;
      else
        val[4]+= i == 0 ? 0 : U;
      
      V = v[j][i] / d.y;
      if( V > 0. )
        val[2] = j == 0 ? 0 : V;
      else
        val[4]+= j == 0 ? 0 : V;
      
      U = u[j][i+1] / d.x;
      if( U > 0. )
        val[4]+= i == M ? 0 : -U;
      else
        val[1] = i == M ? 0 : -U;
      
      V = v[j+1][i] / d.y;
      if( V > 0. )
        val[4]+= j == N ? 0 : -V; 
      else
        val[3] = j == N ? 0 : -V;

      for( row.c = 0; row.c < info.dof; ++row.c)
      {
        for( s = 0; s < S; s++ )
          col[s].c = row.c;
        ierr = MatSetValuesStencil(mat, 1, &row, S, col, val, ADD_VALUES); CHKERRQ(ierr);
      }
    }
  }  
  PetscLogEventEnd(EVENT_AssembleAdvectiveTransport,0,0,0,0);
  PetscFunctionReturn(0);
}

/*
AssembleConvectiveTransport( Mat mat, PetscReal **u, PetscReal **v )
{
  int i,j,m;
  PetscReal xx = 2*d.x, yy = 2 * d.y, 2 * d.z;
  MatStencil row, col[7];
  
  PetscErrorCode ierr;
 
  for (row.k = i.zs; row.k < i.zs+i.zm; ++row.k)
  {
    for (row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
    {
      for (row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
      {
        for( m = 0; m < 7; m++ )
          PetscMemcpy(&col[m],row,sizeof(MatStencil));
        
        PetscReal val[7] = {
                 u[i][j][k] / xx,
              -u[i+1][j][k] / xx,
                 v[i][j][k] / yy,
              -v[i][j+1][k] / yy,
                 w[i][j][k] / zz,
              -w[i][j][k+1] / zz,
              ( u[i][j][k] - u[i+1][j][k] ) / xx + 
              ( v[i][j][k] - v[i][j+1][k] ) / yy +
              ( v[i][j][k] - v[i][j][k+1] ) / zz };
        
        for( row.c = 0; row.c < i.dof; ++row.c)
        {
          
          MatSetValuesStencil(mat, 1, &row, S, col, val, ADD_VALUES);
        }
      }
    }
  }
  
  
  col[0].i
  };
}
*/