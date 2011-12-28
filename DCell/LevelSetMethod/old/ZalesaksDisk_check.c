#include "LevelSetMethod.h"
#include "MyCheck.h"

START_TEST( test_LevelSetAdvect2D )
{
  int i,j;
  PetscInt d1 = 65, d2 = 65, dd = d1*d1; 
  PetscReal dt = 1./66.;
  Grid2D gvx, gvy;
  LevelSet2D ls, lsnew, lstemp;
  char outname[25];
  PetscErrorCode ierr;
  
//Read the ZD levelset
  int fd;
  ierr = CreateLevelSet2D( d1, d1, &ls ); CHKERRQ(ierr);
  ierr = CreateLevelSet2D( d1, d1, &lsnew ); CHKERRQ(ierr);
  ierr = PetscBinaryOpen("zd.Real64",FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, ls->g2d->v1, dd, PETSC_REAL); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
    
//Create circular velocity field
  CreateGrid2D(d1, d1, &gvx);
  CreateGrid2D(d1, d1, &gvy);  
  
  for( j = 0; j < d1; j++ )
  {
    for( i = 0; i < d1; i++ )
    {
      gvx->v2[i][j] = -33;
      gvy->v2[i][j] = 33;
      gvx->v2[i][j] = ( j - d2/2 - 0.1);
      gvy->v2[i][j] = (-i + d1/2 - 0.1);
    }
  }
  
  PetscTruth nan;
  sprintf(outname, "%s%d", "result.",0);
  WriteVectorArray(outname, ls->g2d->len, ls->g2d->v1);
  for( i = 1; i < 10; i++)
  {
    printf("%d\n", i);
    UpdateIrregularNodeList( ls );
    IrregularNodeListWrite( ls, i );
printf("LINE:%d\n", __LINE__);
LevelSet2DCheck( ls, &nan );
    ReinitializeLevelSet( ls );
    sprintf(outname, "%s%d", "reinit.",i);
    WriteVectorArray(outname, ls->g2d->len, ls->g2d->v1);
printf("LINE:%d\n", __LINE__);
LevelSet2DCheck( ls, &nan );
    LevelSetAdvect2D( dt, gvx, gvy, ls, lsnew);
printf("LINE:%d\n", __LINE__);
LevelSet2DCheck( lsnew, &nan );
    sprintf(outname, "%s%d", "result.",i);
    WriteVectorArray(outname, lsnew->g2d->len, lsnew->g2d->v1);
    
    lstemp = ls;
    ls = lsnew;
    lsnew = lstemp;
  }
}
END_TEST

START_TEST( test_ReadZD )
{
  PetscErrorCode ierr;
  LevelSet2D ls, lsnew;
  PetscInt d1, d2, dd;
  char* name = "zd.Real64";
  char outname[25];
  int fd;
  
  PetscFunctionBegin;
  
  d1 = d2 = 65;
  dd = d1 * d2;
  ierr = CreateLevelSet2D( d1, d2, &ls ); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, ls->g2d->v1, dd, PETSC_REAL); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
//  VecView( ls->g2d->v, PETSC_VIEWER_STDOUT_SELF);
  CreateLevelSet2D( d1, d2, &lsnew);
  VecCopy( ls->g2d->v, lsnew->g2d->v);
  UpdateIrregularNodeList( ls );
  ReinitializeLevelSet( ls  );
  PetscReal **temp, **phi = ls->g2d->v2, **fi = lsnew->g2d->v2, 
            px, py, vx, vy, dt = 10*0.1 / d1;
  int i, j, t;
  
  vx = vy = 1. / sqrt(2.);
  for( t = 0; t < 1000; t++ )
  {
    for( j = 1; j < d2-1; j++ )
    {
      for( i = 1; i < d1-1; i++ )
      {
        vx = (j - d2/2);
        vy = (-i + d1/2);
        px = vx > 0. ? phi[i+1][j]-phi[i][j] : phi[i][j]-phi[i-1][j];
        py = vy > 0. ? phi[i][j+1]-phi[i][j] : phi[i][j]-phi[i][j-1];
        fi[i][j] = dt * (vx * px + vy * py) + phi[i][j]; 
      }
    }
//    UpdateIrregularNodeList( lsnew );
//    ReinitializeLevelSet( lsnew  );
    VecCopy( lsnew->g2d->v, ls->g2d->v );
    sprintf(outname, "%s%d", "result.",t);
    WriteVectorArray(outname, lsnew->g2d->len, lsnew->g2d->v1);
  }
  ierr = UpdateIrregularNodeList( lsnew ); CHKERRQ(ierr);
//  ierr = IrregularNodeListWrite( lsnew ); CHKERRQ(ierr);

  DestroyLevelSet2D(ls);
  DestroyLevelSet2D(lsnew);
}
END_TEST

Suite* ZD_suite (void)
{
  Suite *s = suite_create ("Zalesak's Disk Check");
  TCase *tc_core = tcase_create("Core");
//  tcase_add_checked_fixture( tc_core, setup, teardown);
//  tcase_add_test( tc_core,  test_ReadZD );
  tcase_add_test( tc_core,  test_LevelSetAdvect2D );  
//  tcase_add_test( tc_core,  );
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);
  return s;
}