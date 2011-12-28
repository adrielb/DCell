#include "LevelSetMethod.h"

/*
 * NOTE: An irregular grid point cannot be on the bounday of the grid!!!
 */
#undef __FUNCT__
#define __FUNCT__ "IrregularNodeListUpdate"
PetscInt EVENT_IrregularNodeListUpdate;
PetscErrorCode IrregularNodeListUpdate( int x, int y, LevelSet ls )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IrregularNodeListUpdate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IrregularNodeListUpdate,"IrregularNodeListUpdate", 0);
  
  if( ls->g->is2D )
  {
    IrregularNodeListUpdate_2D( x, y, ls );
  } else {
//    IrregularNodeListUpdate_3D( iCoor s, ls );
  }
  
  PetscLogEventEnd(EVENT_IrregularNodeListUpdate,0,0,0,0);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "IrregularNodeListUpdate_2D"
PetscErrorCode IrregularNodeListUpdate_2D( int x, int y, LevelSet ls )
{
  int i, j, k, I, J, ni, nj, bc;
  const int nei[4][2] = {1,0,0,1,-1,0,0,-1};
  IrregularNode n;
  PetscReal **phi = ls->g->v2;
  PetscReal sten[3][3];
  GArray *g = ls->irregularNodes;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  n.z = 0;
  n.oz = 0;

//  GLib errors are generated when index = 0 < len = 0
  if( g->len != 0 )
  {
    g_array_remove_range(g, 0, g->len);
  }
  
//  Index irregular grid points
  for( j = 2; j < ls->g->n.y-2; ++j)
  {
    for( i = 2; i < ls->g->n.x-2; ++i)
    {
//      printf("%d\t%d\n", i,j);
      for( J = -1; J < 2; ++J)
      {
        for( I = -1; I < 2; ++I)
        {
          sten[J+1][I+1] = (phi[j+J+y][i+I+x] + phi[j+J][i+I]) / 2.;
        }
      }
      for( k = 0; k < 4; ++k)
      {
        ni = 1 + nei[k][0];
        nj = 1 + nei[k][1];
        if( sten[1][1] * sten[ni][nj] <= 0. )
        {
          n.x = i;
          n.y = j;
          n.sign = PetscSign( sten[1][1] );
          OrthogonalProjection2D( sten, &n.ox, &n.oy);
          n.ox += x / 2.; //TODO: need to test OP shift direction
          n.oy += y / 2.;
          g_array_append_val(g, n );
          break;
        }
      }
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IrregularNodeListUpdate_3D"
PetscInt EVENT_IrregularNodeListUpdate_3D;
PetscErrorCode IrregularNodeListUpdate_3D( LevelSet ls )
{
  PetscErrorCode ierr;
  PetscReal ***phi = ls->g->v3;
  iCoor s = ls->g->n;
  IrregularNode *n;
  int i,j,k, l, count = 0, len = ls->irregularNodes->len;
  int ni,nj,nk;
  const int nei[3][6] = {{1,0,0,-1, 0, 0},
                         {0,1,0, 0,-1, 0},
                         {0,0,1, 0, 0,-1}};
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IrregularNodeListUpdate_3D,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IrregularNodeListUpdate_3D,"IrregularNodeListUpdate_3D", 0);
  
  for (k = 1; k < s.z-1; ++k)
  {
    for (j = 1; j < s.y-1; ++j)
    {
      for (i = 1; i < s.x-1; ++i)
      {
        for (l = 0; l < 6; ++l)
        {
          ni = i + nei[0][l];
          nj = j + nei[1][l];
          nk = k + nei[2][l];
          
          if( phi[i][j][k] * phi[ni][nj][nk] < 0 )
          {
            if( count < len ) // Reuse IrregularNode allocations 
            {
              n = &g_array_index(ls->irregularNodes,IrregularNode,count);
              ierr = PetscMemzero(n,sizeof(IrregularNode)); 
              n->x = i;
              n->y = j;
              n->z = k;
              n->sign = PetscSign(phi[i][j][k]);
              OrthogonalProjection3D(phi, i, j, k, &n->ox, &n->oy, &n->oz);
              count++;
            } else { // If not enough, append more nodes to end of list
              IrregularNode n;
              n.x = i;
              n.y = j;
              n.z = k;
              n.sign = PetscSign(phi[i][j][k]);
              OrthogonalProjection3D(phi, i, j, k, &n.ox, &n.oy, &n.oz);
              g_array_append_val(ls->irregularNodes, n);
            }
            break;
          }
        }
      }
    }
  }

  //TODO: test if the node numbering is correct when removing excess length
  if( ls->irregularNodes->len != 0 && count < len )
  {
    g_array_remove_range(ls->irregularNodes, count, len);
  }
  PetscLogEventEnd(EVENT_IrregularNodeListUpdate_3D,0,0,0,0);
  PetscFunctionReturn(0);
}

//TODO: Write out irregular node list directly, instead of copying to array

#undef __FUNCT__
#define __FUNCT__ "IrregularNodeListWrite"
PetscInt EVENT_IrregularNodeListWrite;
PetscErrorCode IrregularNodeListWrite( GArray *irreg, int idx )
{
  PetscErrorCode ierr;
  IrregularNode n;
  PetscReal *list;
  int i,j;
  char outname[25];
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IrregularNodeListWrite,0,0,0,0);
 
  ierr = PetscMalloc(irreg->len * 18 * sizeof(PetscReal), &list); CHKERRQ(ierr);
  
  j = 0;
  for (i = 0; i < irreg->len; ++i)
  {
    n = g_array_index(irreg, IrregularNode, i);
    list[j++] = n.x + n.ox;  // 0
    list[j++] = n.y + n.oy;  // 1
    list[j++] = n.z + n.oz;  // 2
    list[j++] = n.x;   // 3
    list[j++] = n.y;   // 4
    list[j++] = n.z;   // 5
    list[j++] = n.nx;  // 6
    list[j++] = n.ny;  // 7
    list[j++] = n.nz;  // 8
    list[j++] = n.sign;// 9
    list[j++] = n.k;   // 10
    list[j++] = n.f1;  // 11
    list[j++] = n.df1; // 12
    list[j++] = n.ddf1;// 13
    list[j++] = n.f2;  // 14
    list[j++] = n.df2; // 15
    list[j++] = n.ddf2;// 16
    list[j++] = n.c;   // 17
  }
  
  sprintf(outname, "irregNode.%d", idx);
  ierr = WriteVectorArray(outname, j, list); CHKERRQ(ierr);
  PetscFree(list);
  
  PetscLogEventEnd(EVENT_IrregularNodeListWrite,0,0,0,0);
  PetscFunctionReturn(0);
}

void RegisterEvents_IrregularNodes()
{
  
  PetscLogEventRegister(&EVENT_IrregularNodeListWrite,"IrregularNodeListWrite", 0);
}