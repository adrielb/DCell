#include "ImmersedInterfaceMethod.h"
/* 
 * Title: A simple implementation of the immersed interface methods for Stokes flows with singular forces
 * Authors: Lai,Ming-Chih; Tseng,Hsiao-Chieh
 * Source: Computers & Fluids, 2008, 37, 2, 99-106
 */

#undef __FUNCT__
#define __FUNCT__ "IIMSimplePressureGradientCorrection2D"
PetscInt EVENT_IIMSimplePressureGradientCorrection2D;
PetscErrorCode IIMSimplePressureGradientCorrection2D(IIM iim, LevelSet ls, Coor dh, PetscReal **px, PetscReal **py)
{
  IrregularNode *n, **nodes;
  PetscReal pc;
  PetscReal **phi=ls->g->v2;
  PetscReal mid;
  PetscInt len;
  Jump jump;
  int i, j;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IIMSimplePressureGradientCorrection2D,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IIMSimplePressureGradientCorrection2D,"IIMSimplePressureGradientCorrection2D", 0);
  
  for( j = 0; j < ls->irregularNodes->len; j++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, j);

    KDTreeSetEps(iim->kdtree, 4.1);
    ierr = KDTreeRange(iim->kdtree, n, &nodes, &len ); CHKERRQ(ierr);
    
    pc = 0;
    for( i = 0; i < len; i++ )
    {
      if( n->x - nodes[i]->x == 1 && nodes[i]->y == n->y && 
          phi[n->y][nodes[i]->x] * phi[n->y][n->x] < 0) //Correction for px
      {
        mid = PetscSign( ( phi[n->y][n->x] + phi[n->y][n->x-1] ) / 2. );
        if( n->sign > 0 ) { //TODO: this is different from what was published because this cancels the gradient in a normal only jump
          JumpConditionPressure2D( nodes[i], &jump );
          pc = IIMSimpleCorrection( &jump ); // p[i] - ( p[i-1] + pc[i-1] ) 
          px[n->y][n->x] -= pc / dh.x;       // pi-1--|--ui-----pi
        } else {
          JumpConditionPressure2D( n, &jump );
          pc = IIMSimpleCorrection( &jump ); // p[i] + pc[i] - p[i-1]  
          px[n->y][n->x] += pc / dh.x;       // pi-1-----ui--|--pi
        }
      }
      if( n->y - nodes[i]->y == 1 && nodes[i]->x == n->x &&
          phi[nodes[i]->y][n->x] * phi[n->y][n->x] < 0 ) //Correction for py        
      {
        mid = PetscSign( ( phi[n->y][n->x] + phi[n->y-1][n->x] ) / 2. );
//        if( n->sign == mid ) { //how it was in Lai et al
        if( n->sign > 0 ) { //TODO: is the jump based on the current irregular node 'n' or its neighbor, 'nodes[i]' ?
          JumpConditionPressure2D( nodes[i], &jump );
          pc = IIMSimpleCorrection( &jump );
          py[n->y][n->x] -= pc / dh.y;
        } else {
          JumpConditionPressure2D( n, &jump );
          pc = IIMSimpleCorrection( &jump );
          py[n->y][n->x] += pc / dh.y;
        }
      } // if node part of stencil and crosses level set 
    } // for each irreg 'nodes' in KDTreeRange() 
  } // for each irreg 'n' in LevelSet
  
  PetscLogEventEnd(EVENT_IIMSimplePressureGradientCorrection2D,0,0,0,0);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IIMSimpleCorrection2D"
PetscInt EVENT_IIMSimpleCorrection2D;
PetscErrorCode IIMSimpleCorrection2D( IIM iim, LevelSet ls, JumpCondition jc, Coor d, PetscReal **rhs )
{
  PetscReal h;
  PetscInt len;
  IrregularNode *n, **nodes;
  Jump jump;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IIMSimpleCorrection2D,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IIMSimpleCorrection2D,"IIMSimpleCorrection2D", 0);

  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    
    //TODO: Get rid of KDTree and just do brute force search of 3x3 neighborhood of integer index 2D field of irregular nodes 
    KDTreeSetEps(iim->kdtree, 4.1);
    ierr = KDTreeRange(iim->kdtree, n, &nodes, &len ); CHKERRQ(ierr);

    for( int i = 0; i < len; i++ )
    {
      if( n->sign * nodes[i]->sign > 0 ) continue;
      if( ( PetscAbs( nodes[i]->x - n->x ) == 1 && nodes[i]->y == n->y ) ||
          ( PetscAbs( nodes[i]->y - n->y ) == 1 && nodes[i]->x == n->x )  )
      {
        if( PetscAbs( nodes[i]->x - n->x ) == 1 && nodes[i]->y == n->y )
          h = d.x * d.x;
        else 
          h = d.y * d.y;
          
        jc( nodes[i], &jump );
        if( n->sign < 0 ) //TODO: use VecSetValues()
          rhs[n->y][n->x] += IIMSimpleCorrection( &jump ) / h;
        else 
          rhs[n->y][n->x] -= IIMSimpleCorrection( &jump ) / h;
      }
    }
  }
  
  PetscLogEventEnd(EVENT_IIMSimpleCorrection2D,0,0,0,0);
  PetscFunctionReturn(0);
}

double IIMSimpleCorrection( Jump* j )
{
  const PetscReal jpnn = j->jf - j->k * j->jpn - j->jpss;
  return j->jp + j->d * j->jpn + 0.5 * j->d * j->d * jpnn;
}
