double ENO1up( double **phi, int i, int j);
double ENO1down(double **phi,int i, int j);
PetscErrorCode BCFromDist( int d1, int d2, double **phi, double **phi1, HEAP *qPos, HEAP *qNeg );
PetscErrorCode PDESignedDistance( int d1, int d2, double **phi, double **phi_buffer );

START_TEST( check_BCFromDist2 )
{
  Vec v, v1;
  int i, j, d1 = 40, d2 = 40, dd = d1*d2;
  double **d;
  PetscErrorCode ierr;
  
  VecCreate(PETSC_COMM_SELF, &v);
  VecSetSizes(v, dd, dd);
  VecSetType(v, VECSEQ);
  VecDuplicate(v, &v1);
  VecSet(v1, INF);
  double *p, **phi, **phi1, **swap;
  VecGetArray(v, &p);

  VecRestoreArray(v, &p);
  
  VecGetArray2d(v, d1,d2,0,0, &phi);
  VecGetArray2d(v1, d1,d2,0,0, &phi1);
  
  for( i = 0; i < d1; i++)
  {
    for( j = 0; j < d2; ++j)
    {
      if( sqrt( SQR( i - d1/2 ) + SQR( j - d2/2 ) ) < d1 / 3 )
      {
        phi[i][j] = -1;
      } else {
        phi[i][j] =  1;
      }
    }
  }

  ierr = PDESignedDistance( d1, d2, phi, phi1 ); CHKERRQ(ierr);  
       
  PrintPhi(d1,d2,phi); 
  VecRestoreArray2d(v,  d1,d2, 0,0, &phi);
  VecRestoreArray2d(v1, d1,d2, 0,0, &phi1);
}
END_TEST

START_TEST( check_BCFromDist )
{
  Vec v, v1;
  double d[64] = {
    -1,1,1, 1, 1,1,1,1,
    -1,-1,1, 1, 1,1,1,1,
    -1,-1,-1, 1, 1,1,1,1,
    -1,-1,-1, -1, 1,1,1,1,
    -1,-1,-1, -1, -1,1,1,1,
    -1,-1,-1, -1, -1,-1,1,1,
    -1,-1,-1, -1, -1,-1,-1,1,
    -1,-1,-1, -1, -1,-1,-1,-1
  };
  PetscErrorCode ierr;
  
  VecCreate(PETSC_COMM_SELF, &v);
  VecSetSizes(v, 64, 64);
  VecSetType(v, VECSEQ);
  VecDuplicate(v, &v1);
  VecSet(v1, INF);
  double *p, **phi, **phi1, **swap;
  VecGetArray(v, &p);
  for( int i = 0; i < 64; i++ )
    p[i] = d[i];
  VecRestoreArray(v, &p);
  
  VecGetArray2d(v, 8,8,0,0, &phi);
  VecGetArray2d(v1, 8,8,0,0, &phi1);
  
  ierr = BCFromDist(8,8, phi, phi1, NULL, NULL); CHKERRQ(ierr);
       
  PrintPhi(8,8,phi); 
  VecRestoreArray2d(v, 8,8, 0,0, &phi);
  VecRestoreArray2d(v1, 8,8,0,0, &phi1);
}
END_TEST

//  tcase_add_test( tc_core, check_BCFromDist );
//  tcase_add_test( tc_core, check_BCFromDist2 );

#undef __FUNCT__
#define __FUNCT__ "PDESignedDistance"
PetscInt EVENT_PDESignedDistance;
PetscErrorCode PDESignedDistance( int d1, int d2, double **phi, double **phi_buffer )
{
	PetscErrorCode ierr;
	int i, j;
	double **swap;
  PetscReal mag, dx = 1, s, d, dt, maxS = 1., dphi, max_dphi;	
	
	PetscFunctionBegin;
	PetscLogEventBegin(EVENT_PDESignedDistance,0,0,0,0);
//	PetscLogEventRegister(&EVENT_PDESignedDistance,"PDESignedDistance", 0);
	
  
  dt = 0.5 * dx / maxS;
  
  for( int t = 0; t < d1; t++ )
  {
    max_dphi = 0;
  	for( i = 1; i < d1-1; i++)
  	{
  		for( j = 1; j < d2-1; j++)
  		{
  			d = phi[i][j];
  			
  			if( d < 0. )
  				mag = ENO1up( phi, i, j);
  			else
  				mag = ENO1down(phi, i, j);
		    
		    s = d / sqrt( SQR(d) + SQR( mag * dx ) );
        dphi = dt / dx * s * ( mag - 1 );
		    phi_buffer[i][j] = phi[i][j] - dphi;
        
        if( ABS( dphi ) > max_dphi )
          max_dphi = dphi;
  		}
  	}
    
    printf("%f \n", max_dphi);
    
    for( i = 0; i < d1; i++)
    {
      j = 0;
      phi_buffer[i][j] = phi_buffer[i][j+1];
      j = d2 - 1;
      phi_buffer[i][j] = phi_buffer[i][j-1];
    }
    
    for( j = 0; j < d2; j++)
    {
      i = 0;
      phi_buffer[i][j] = phi_buffer[i+1][j];
      i = d1 - 1;
      phi_buffer[i][j] = phi_buffer[i-1][j];
    }
  	
    swap = phi;
	  phi = phi_buffer;
	  phi_buffer = swap;
  }
	 
	
	PetscLogEventEnd(EVENT_PDESignedDistance,0,0,0,0);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCFromDist"
PetscInt EVENT_BCFromDist;
PetscErrorCode BCFromDist( int d1, int d2, double **phi, double **phi1, HEAP *qPos, HEAP *qNeg )
{
  PetscErrorCode ierr;
  guint i, j, k;
  int nei[4][2] = {1,0,0,1,-1,0,0,-1};
  int ni, nj, bc;
  GArray *garrayPos = g_array_new( FALSE, FALSE, sizeof(int) ),
         *garrayNeg = g_array_new( FALSE, FALSE, sizeof(int) );
  PetscTruth isRegular = PETSC_FALSE;
  double **swap;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_BCFromDist,0,0,0,0);
  
  

//Index irregular grid points
  for( i = 1; i < d1-1; ++i)
  {
    for( j = 1; j < d2-1; ++j)
    {
      isRegular = PETSC_TRUE;
      for( k = 0; k < 4; ++k)
      {
        ni = i + nei[k][0];
        nj = j + nei[k][1];
        if( phi[i][j] * phi[ni][nj] < 0. )
        {
          isRegular = PETSC_FALSE;
          bc = i + j * d1;
          if( phi[i][j] > 0. )
          {
            g_array_append_val(garrayPos, bc );
            break;
          } else {
            g_array_append_val(garrayNeg, bc );
            break;
          }
        }
      }
      if( isRegular )
      {
//      	phi[i][j] = INF;
      	phi1[i][j]= INF;
      }
    }
  }
  
  for( i = 0; i < d1; i++)
  {
    j = 0;
    phi[i][j] = INF;
    phi1[i][j]= INF;
    j = d2 - 1;
    phi[i][j] = INF;
    phi1[i][j]= INF;
  }
  
  for( j = 0; j < d2; j++)
  {
    i = 0;
    phi[i][j] = INF;
    phi1[i][j]= INF;
    i = d1 - 1;
    phi[i][j] = INF;
    phi1[i][j]= INF;
  }
  int *aPos = (int*)garrayPos->data,
      *aNeg = (int*)garrayNeg->data,
      c1, c2;
  
  PetscReal mag, dx = 1, s, d, dt, maxS = 1.;
  
  dt = 0.5 * dx / maxS;
  
  for( j = 0; j < 50; j++ )
  {
		printf("---\n");
	  for( i = 0; i < garrayNeg->len; ++i)
	  {
	    c2 = aNeg[i] / d1;
	    c1 = aNeg[i] - d1 * c2;
	    mag = ENO1up( phi, c1, c2);
	    d = phi[c1][c2];
	    s = d / sqrt( SQR(d) + SQR( mag * dx ) );
	    phi1[c1][c2] = phi[c1][c2] - dt / dx * s * ( mag - 1 );
      printf( "mag: %f, %f, (%d,%d) %d/%d\n", mag, phi1[c1][c2], c1, c2, i, garrayNeg->len);
	  }
	  
	  for( i = 0; i < garrayPos->len; ++i)
	  {
	    c2 = aPos[i] / d1;
	    c1 = aPos[i] - d1 * c2;
	    mag = ENO1down( phi, c1, c2);
	    d = phi[c1][c2];
	    s = d / sqrt( SQR(d) + SQR( mag * dx ) );
	    phi1[c1][c2] = phi[c1][c2] - dt / dx * s * ( mag - 1 );
      printf( "mag: %f, %f, (%d,%d) %d/%d\n", mag, phi1[c1][c2], c1, c2, i, garrayPos->len);
	  }
	  swap = phi;
	  phi = phi1;
	  phi1 = swap;
  }
  
  g_array_free(garrayPos, TRUE);
  g_array_free(garrayNeg, TRUE);
  
  PetscLogEventEnd(EVENT_BCFromDist,0,0,0,0);
  PetscFunctionReturn(0);
}
/*
  double maxXm = MAX(f*(phi[i][j]-phi[i][j-1]),0.),
         maxXp = MAX(f*(phi[i][j]-phi[i][j+1]),0.),
         maxX  = MAX(maxXm,maxXp),
         maxYm = MAX(f*(phi[i][j]-phi[i-1][j]),0.),
         maxYp = MAX(f*(phi[i][j]-phi[i+1][j]),0.),
         maxY  = MAX(maxYm,maxYp);
*/
inline double ENO1(double** phi, int i, int j, double f)
{
	double maxXm = phi[i][j-1]==INF ? 0. : (phi[i][j]-phi[i][j-1]),
				 maxXp = phi[i][j+1]==INF ? 0. : (phi[i][j]-phi[i][j+1]),
				 maxX  = MAX(f*maxXm,f*maxXp),
				 maxYm = phi[i-1][j]==INF ? 0. : (phi[i][j]-phi[i-1][j]),
				 maxYp = phi[i+1][j]==INF ? 0. : (phi[i][j]-phi[i+1][j]),
				 maxY  = MAX(f*maxYm,f*maxYp);

  
  if( maxX > 10 || maxY > 10 )
  {
    printf("phi: %f\ti:%d\tj:%d\tf:%f\n",phi[i][j], i, j, f);
    exit(1);
  }
  
//  if( sqrt(maxX * maxX + maxY * maxY) < .0001 )
//    exit(5);
    
//  printf("maxX: %f\tmaxY: %f\n", maxX, maxY);
	
	return sqrt(maxX * maxX + maxY * maxY);
}

inline double ENO1up(double **phi, int i, int j)
{
  return ENO1( phi, i, j, -1.); 
}

inline double ENO1down(double **phi, int i, int j)
{
  return ENO1( phi, i, j, 1.); 
}