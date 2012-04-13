//#include "ImmersedInterfaceMethod.h"
//#include "MyCheck.h"
#include "Main.h"
#include "petscda.h"
int RunCheck(){}
//int BenchmarkBoundaryChecks( int n, PetscLogDouble t_bulk, PetscLogDouble )
int PetscMain()
{
  int i,j,k, n=128;
  PetscLogDouble t1,t2,s1,s2;
  Vec U,V,W,DIV1,DIV2;
  PetscReal ***u,***v,***w,***div1, ***div2;
  DA da;
  DALocalInfo info;
  PetscErrorCode ierr;
  ierr = DACreate3d(PETSC_COMM_SELF,//MPI Communicator   
      DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
      DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
      n,n,n, //Global array dimension
      1,1,1,//Number procs per dim
      1,    //Number of chemical species
      1,           //stencil width
      0,0,0,       //specific array of nodes
      &da); CHKERRQ(ierr);
  DACreateGlobalVector(da,&U);
  DACreateGlobalVector(da,&V);
  DACreateGlobalVector(da,&W);
  DACreateGlobalVector(da,&DIV1);
  DACreateGlobalVector(da,&DIV2);
  VecSet(DIV1,0);
  VecSet(DIV2,0);
  DAVecGetArray(da,U,&u);
  DAVecGetArray(da,V,&v);
  DAVecGetArray(da,W,&w);
  DAVecGetArray(da,DIV1,&div1);
  DAVecGetArray(da,DIV2,&div2);
  DAGetLocalInfo(da,&info);
  PetscBarrier(0);
  for( k = 0; k < n; ++k)
  {
    for( j = 0; j < n; ++j)
    {
      for( i = 0; i < n; ++i)
      {
        u[k][j][i] = i * (i-n) * j * (j-n) * k * (k-n);
        v[k][j][i] = 1 - u[k][j][i];
        w[k][j][i] = u[k][j][i] * v[k][j][i];
      }
    }
  }
  PetscBarrier(0);
  PetscGetTime(&t1);
  for( k = 1; k < n-1; ++k)
  {
    for( j = 1; j < n-1; ++j)
    {
      for( i = 1; i < n-1; ++i)
      {
        div1[k][j][i] = u[k][j][i+1] - u[k][j][i-1] +
                        v[k][j+1][i] - v[k][j-1][i] +
                        w[k+1][j][i] - w[k-1][j][i];
        div1[k][j][i] /= 2;
      }
    }
  }
  PetscGetTime(&t2);
  PetscReal uE,uW,vN,vS,wF,wB;
  PetscReal hx,hy,hz;
  PetscBarrier(0);
  PetscGetTime(&s1);
  for( k = 1; k < n-1; ++k)
  {
    for( j = 1; j < n-1; ++j)
    {
      for( i = 1; i < n-1; ++i)
      {
/*        uE = i == info.mx-1 ? u[k][j][i] : u[k][j][i+1]; 
        uW = i == 0         ? u[k][j][i] : u[k][j][i-1];
        vN = j == info.my-1 ? v[k][j][i] : v[k][j+1][i];
        vS = j == 0         ? v[k][j][i] : v[k][j-1][i];
        wB = k == info.mz-1 ? w[k][j][i] : w[k+1][j][i];
        wF = k == 0         ? w[k][j][i] : w[k-1][j][i];
  */      
if( i == info.mx-1 ) { uE = u[k][j][i]; hx= 1;  }else{ uE = u[k][j][i+1];hx= 2;}
if( i == info.mx-1 ) { uE = u[k][j][i]; hx= 1;  }else{ uE = u[k][j][i+1];hx= 2;}
if( j == info.mx-1 ) { uE = u[k][j][i]; hy= 1;  }else{ uE = u[k][j][i+1];hy= 2;}
if( j == info.mx-1 ) { uE = u[k][j][i]; hy= 1;  }else{ uE = u[k][j][i+1];hy= 2;}
if( k == info.mx-1 ) { uE = u[k][j][i]; hz= 1;  }else{ uE = u[k][j][i+1];hz= 2;}
if( k == info.mx-1 ) { uE = u[k][j][i]; hz= 1;  }else{ uE = u[k][j][i+1];hz= 2;}

        div2[k][j][i] = uE - uW + vN - vS + wB - wF;
        div2[k][j][i] /= 2;
//        printf("%f\t%f\t%f\n",div1[k][j][i], div2[k][j][i], div2[k][j][i] - div1[k][j][i]);        
      }
    }
  }
  PetscGetTime(&s2);
  DAVecRestoreArray(da,DIV1,&div1);
  DAVecRestoreArray(da,DIV2,&div2);
  VecAXPY(DIV1,-1,DIV2);
  PetscReal norm;
  VecNorm(DIV1,NORM_2,&norm);
  
  printf("BULK: %f\nIF's: %f\nDIFF:\t%f\nRATIO:\t%f\nnorm: %f\n", (t2-t1), (s2-s1), (s2-s1)-(t2-t1),(s2-s1)/(t2-t1),norm);
}

/*
Suite* FluidField_Suite(void)
{
  Suite *s = suite_create ("Fluid Field Suite");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core, time_BoundaryChecks );
  //  tcase_add_test( tc_core,  );  
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);

  return s;
}
*/