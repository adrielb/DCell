#undef __FUNCT__
#define __FUNCT__ "LevelSet2DCheck"
PetscInt EVENT_LevelSet2DCheck;
PetscErrorCode LevelSet2DCheck( LevelSet2D ls, PetscTruth *nan )
{
  int i, j;
  PetscReal **phi = ls->g2d->v2;
  PetscReal d1 = ls->g2d->d1, d2 = ls->g2d->d2;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSet2DCheck,0,0,0,0);
  
  for( j = 0; j < d2; ++j)
  {
    for( i = 0; i < d1; ++i)
    {
      if( phi[i][j] != phi[i][j] )
      {
        printf("{%d, %d},", i, j);
        *nan = PETSC_TRUE;
      }
    }
  }
  
  if( *nan == PETSC_TRUE ) exit(0);
  
  PetscLogEventEnd(EVENT_LevelSet2DCheck,0,0,0,0);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvect2D"
PetscInt EVENT_LevelSetAdvect2D;
PetscErrorCode LevelSetAdvect2D( PetscReal dt, Grid2D gvx, Grid2D gvy, 
  LevelSet2D prev, LevelSet2D new )
{
  PetscErrorCode ierr;
  int i, j, d1 = gvx->d1, d2 = gvx->d2;
  PetscReal **vx = gvx->v2, **fi  = new->g2d->v2,
            **vy = gvy->v2, **phi = prev->g2d->v2, px, py;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSetAdvect2D,0,0,0,0);
  
  for( j = 1; j < d2-1; j++ )
  {
    for( i = 1; i < d1-1; i++ )
    {
      px = vx[i][j] > 0. ? phi[i+1][j]-phi[i][j] : phi[i][j]-phi[i-1][j];
      py = vy[i][j] > 0. ? phi[i][j+1]-phi[i][j] : phi[i][j]-phi[i][j-1];
      fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
    }
  }
  
  for( i = 1; i < d1-1; ++i)
  {
    j = 0;
    px = vx[i][j] > 0. ? phi[i+1][j]-phi[i][j] : phi[i][j]-phi[i-1][j];
    py = phi[i][j+1]-phi[i][j];
    fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
    j = d2-1;
    px = vx[i][j] > 0. ? phi[i+1][j]-phi[i][j] : phi[i][j]-phi[i-1][j];
    py = phi[i][j]-phi[i][j-1];
    fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  }
  
  for( j = 1; j < d2-1; ++j)
  {
    i = 0;
    px = phi[i+1][j]-phi[i][j];
    py = vy[i][j] > 0. ? phi[i][j+1]-phi[i][j] : phi[i][j]-phi[i][j-1];
    fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
    i = d1-1;
    px = phi[i][j]-phi[i-1][j];
    py = vy[i][j] > 0. ? phi[i][j+1]-phi[i][j] : phi[i][j]-phi[i][j-1];
    fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  }
  
  i = j = 0;
  px = phi[i+1][j] - phi[i][j];
  py = phi[i][j+1] - phi[i][j];
  fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  
  i = 0; j = d2-1;
  px = phi[i+1][j] - phi[i][j];
  py = phi[i][j] - phi[i][j-1];
  fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  
  i = d1-1; j = 0;
  px = phi[i][j] - phi[i-1][j];
  py = phi[i][j+1] - phi[i][j];
  fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  
  i = d1-1; j = d2-1;
  px = phi[i][j] - phi[i-1][j];
  py = phi[i][j] - phi[i][j-1];
  fi[i][j] = dt * (vx[i][j] * px + vy[i][j] * py) + phi[i][j];
  
  PetscLogEventEnd(EVENT_LevelSetAdvect2D,0,0,0,0);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "LevelSetNormal2D"
PetscInt EVENT_LevelSetNormal2D;
PetscErrorCode LevelSetNormal2D( LevelSet2D ls, Grid2D gnx, Grid2D gny)
{
  int i,j;
  PetscReal px, py, pn;
  PetscReal **p = ls->g2d->v2;
  PetscReal **nx= gnx->v2, **ny = gny->v2;
  
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSetNormal2D,0,0,0,0);

  for( j = 1; j < ls->g2d->d2 - 1; ++j)
  {
    for( i = 1; i < ls->g2d->d1 - 1; ++i)
    {
      px = (p[i+1][j]-p[i-1][j]) / two;
      py = (p[i][j+1]-p[i][j-1]) / two;
      pn = (px*px + py*py + PETSC_SMALL);
      nx[i][j] = px / pn;
      ny[i][j] = py / pn;
    }
  }
  
  PetscLogEventEnd(EVENT_LevelSetNormal2D,0,0,0,0);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "UpdateIrregularNodeList"
PetscInt EVENT_UpdateIrregularNodeList;
PetscErrorCode UUpdateIrregularNodeList( int x, int y, LevelSet2D ls )
{
  PetscErrorCode ierr;
  int i, j, k, I, J, ni, nj, bc;
  const int nei[4][2] = {1,0,0,1,-1,0,0,-1};
  IrregularNode n;
  PetscReal **phi = ls->g2d->v2;
  PetscReal sten[3][3];
  GArray *g = ls->irregularNodes;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_UpdateIrregularNodeList,0,0,0,0);

  n.z = 0;
  n.oz = 0;

//  GLib errors are generated when index = 0 < len = 0
  if( g->len != 0 )
  {
    g_array_remove_range(g, 0, g->len);
  }
  
//  Index irregular grid points
  for( j = 2; j < ls->g2d->d2-2; ++j)
  {
    for( i = 2; i < ls->g2d->d1-2; ++i)
    {
//      printf("%d\t%d\n", i,j);
      for( I = -1; I < 2; ++I)
      {
        for( J = -1; J < 2; ++J)
        {
          sten[I+1][J+1] = (phi[i+I+x][j+J+y] + phi[i+I][j+J]) / 2.;
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
  
  PetscLogEventEnd(EVENT_UpdateIrregularNodeList,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode LevelSetFindMin(LevelSet ls, iCoor *g)
{
  PetscReal min = 0;
  iCoor n = ls->phi->n;
  int i, j, k;
  //TODO: use VecMin and convert 1D index to coordinate (i,j,k);
  PetscErrorCode ierr;
  if( ls->phi->is2D ) // 2D min search
  {
    PetscReal **phi;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    for( j = 0; j < n.y; ++j) {
      for( i = 0; i < n.x; ++i) {
        if( phi[j][i] < min ) {
          g->x = i;
          g->y = j;
          min = phi[j][i];
        }
      }
    }
  } else {  // 3D min search
    PetscReal ***phi;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    for( k = 0; k < n.z; ++k)
    {
      for( j = 0; j < n.y; ++j)
      {
        for( i = 0; i < n.x; ++i)
        {
          if( phi[k][j][i] < min )
          {
            g->x = i;
            g->y = j;
            g->z = k;
            min = phi[k][j][i];
          }
        }
      }
    }
  }
  return 0;
}
