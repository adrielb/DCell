#include <stdlib.h>
#include "MyCheck.h"
#include "ActiveContoursWithoutEdges.h"

#undef CHKERRQ
#define CHKERRQ(n) if (n) { \
  const char *text;\
  char *msg;\
  PetscErrorMessage(n,&text,&msg);\
  fail("\n%s(%d): %s (%s)",__FILE__,__LINE__, text, msg);}
/*
  PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,n,0," "); 

#define myfail_unless( test, text ) if( test ) { \
  printf("%s(%d): %s\n", __FILE__, __LINE__, text); \
  fail_unless(0==1, text);} \ 
*/
#undef fail_unless
#define fail_unless(expr, ...) {\
  char buffer[256], msg[256]; \
  sprintf(msg, ##__VA_ARGS__); \
  sprintf(buffer, "\n%s(%d): %s",__FILE__, __LINE__,msg); \
        _fail_unless(expr, __FILE__, __LINE__,\
        "Assertion '"#expr"' failed" , buffer, NULL);}


 
void setup();
void teardown();

START_TEST( DMMG_usage )
{
  PetscErrorCode ierr;
  DMMG *dmmg;
  DA da;
  PC pc;
  
  DMMGCreate(PETSC_COMM_SELF, 3, PETSC_NULL, &dmmg);
  DACreate1d(PETSC_COMM_SELF, DA_XPERIODIC, 5, 1, 1, PETSC_NULL, &da);
  DMMGSetDM( dmmg, (DM) da );
  
  ierr = DMMGSetKSP(dmmg, ComputeRHS, ComputeJacobian); CHKERRQ(ierr);
  ierr = DMMGSetUp(dmmg); CHKERRQ(ierr);
  KSP ksp;
/*  for( int i = 0; i < DMMGGetLevels(dmmg); i++)
  {
    ksp = DMMGGetKSP(&dmmg[i]);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
//    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  }*/
  ksp = DMMGGetKSP(dmmg);
  KSPSetType(ksp, KSPCG);
  KSPGetPC(ksp, &pc);
//  PCMGGetSmoother(pc, i, ksp);
  PCView(pc, PETSC_VIEWER_STDOUT_SELF);
  ierr = DMMGView(dmmg, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = DMMGSolve(dmmg); CHKERRQ(ierr);
  
  DMMGDestroy(dmmg); 
}
END_TEST

START_TEST( DA_usage )
{
  PetscErrorCode ierr;
  DA da, daf;
  Mat M;
  
  ierr = DACreate2d(PETSC_COMM_SELF,//MPI Communicator   
    DA_NONPERIODIC,         //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
    DA_STENCIL_STAR,       //DA_STENCIL_BOX or DA_STENCIL_STAR
    4, 4,     //Global dimension
    1, 1,      //Number procs per dim
    1,        //dof
    1,       //stencil width
    0,0,    //specific array of nodes
    &da); CHKERRQ(ierr);

  
  Vec v, r, scale;
  PetscReal **array;
  PetscInt mx, my, i,j, count;
  ierr = DAGetGlobalVector(da, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(da, v, &array); CHKERRQ(ierr);
  ierr = DAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0); CHKERRQ(ierr);
  count = 1;
  for (i = 0; i < mx; ++i)
	{
		for (j = 0; j < my; ++j)
		{
			array[i][j] = count;
			count++;
		}
	}
  ierr = DAVecRestoreArray(da, v, &array); CHKERRQ(ierr);
  ierr = VecView(v, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
	ierr = DARefine(da, PETSC_COMM_SELF, &daf); CHKERRQ(ierr);
	ierr = DAGetInterpolation(da, daf, &M, &scale); CHKERRQ(ierr);
	ierr = DAGetGlobalVector(daf,&r); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE); CHKERRQ(ierr);
	ierr = MatView(M, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = MatInterpolate(M, v, r); CHKERRQ(ierr);
//	ierr = MatMult(M, v, r); CHKERRQ(ierr);
	ierr = VecView(r, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = VecView(scale, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = MatRestrict(M, r, v); CHKERRQ(ierr);
//	ierr = MatMultTranspose(N, r, v); CHKERRQ(ierr);
	ierr = VecView(v, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = VecPointwiseMult(v, scale, v); CHKERRQ(ierr);
  ierr = VecView(v, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
	ierr = MatDestroy(M); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  ierr = DADestroy(daf); CHKERRQ(ierr);
}
END_TEST

START_TEST( FMM_last_phi )
{
  PetscErrorCode ierr;
  PetscInt len = 512*512;
  int fd;
  PetscReal *dd, **m, **p;
  Vec mask, phi;
  ierr = VecCreate(PETSC_COMM_SELF, &mask); CHKERRQ(ierr);
  ierr = VecSetSizes(mask, PETSC_DECIDE, len); CHKERRQ(ierr);
  ierr = VecSetType(mask, VECSEQ); CHKERRQ(ierr);
  ierr = VecDuplicate(mask, &phi); CHKERRQ(ierr);
  ierr = VecGetArray(mask, &dd); CHKERRQ(ierr);
  ierr = PetscBinaryOpen("phi_final.img",FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, dd, len, PETSC_REAL); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  ierr = VecRestoreArray(mask, &dd); CHKERRQ(ierr);
  
  ierr = VecGetArray2d(mask, 512,512, 0,0,&m); CHKERRQ(ierr);
  ierr = VecGetArray2d(phi, 512, 512, 0,0, &p); CHKERRQ(ierr);
  
  HEAP *q1 = HeapAlloc(100,&MyDoubleCompMAX);
  HEAP *q2 = HeapAlloc(100,&MyDoubleCompMIN);
  
  
  BCFromMask( 512,512, m, p, q1, q2);
  
  SolveEikonal(q1, -1, 512, 512, p);
  SolveEikonal(q2,  1, 512, 512, p);
  
//  CheckEikonalProperty( 512, 512, p);
  
  HeapFree(q1);
  HeapFree(q2);
  
  ierr = VecRestoreArray2d(phi, 512, 512, 0, 0, &p); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(mask,512, 512, 0,0, &m); CHKERRQ(ierr);

	ierr = VecGetArray(phi, &dd); CHKERRQ(ierr);
  ierr = PetscBinaryOpen("phi_final.bin",FILE_MODE_WRITE,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd, dd, len, PETSC_REAL, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  ierr = VecRestoreArray(phi, &dd); CHKERRQ(ierr);
  
  ierr = VecDestroy(mask); CHKERRQ(ierr);
  ierr = VecDestroy(phi); CHKERRQ(ierr);
}
END_TEST

START_TEST( Single_time_step )
{
  PetscErrorCode ierr;
  UserContext *uc;
  Image *img;
  Vec phi_init, phi_old;
  PetscReal *dd, **m, norm, norm_old;
  char *name =
//  "/home/abergman/Images/image.a.1.gray";
  "/home/abergman/Images/Joanne/DIR.IGF1 Adaptation-10.spl/image.a.1.gray";
//	"testimage.img";
  ierr = PetscLogBegin(); CHKERRQ(ierr);
  
  ierr = ReadCascade512File( name, &img); CHKERRQ(ierr);
	ierr = UserContextCreate(img->v, &uc); CHKERRQ(ierr);
  
  VecDuplicate(uc->phi, &phi_old);
  
  VecCopy(img->v, uc->phi);
  VecShift(uc->phi, -5000.5);
  Reinitialize(uc->phi);
  ierr = WriteVector("phi_0.bin", uc->phi); CHKERRQ(ierr);
  
  char buf[256];
  for (int i = 1; i < 100; ++i)
  {
    norm_old = uc->c1;
    ierr = VecCopy(uc->phi, phi_old); CHKERRQ(ierr);
    for (int j = 0; j < 10; ++j)
    {
      ierr = SingleStep( img, uc); CHKERRQ(ierr);
    }
    sprintf(buf, "phi_%d.bin", i);
    ierr = Reinitialize(uc->phi); CHKERRQ(ierr);
//    ierr = WriteVector(buf, uc->phi); CHKERRQ(ierr);
    ierr = VecAXPY(phi_old, -1., uc->phi); CHKERRQ(ierr);
//    norm_old = norm;
//    ierr = VecNorm(phi_old,NORM_2, &norm); CHKERRQ(ierr);
    norm = uc->c1;
    LINE("\nc1: %f\nc2: %f\nNorm:%f\n",uc->c1,uc->c2, norm-norm_old);
  } 
  ierr = WriteVector("phi_final.bin", uc->phi); CHKERRQ(ierr);  
  ierr = PetscLogPrintSummary(PETSC_COMM_WORLD,PETSC_NULL); CHKERRQ(ierr);
  ierr = DestroyImage(img); CHKERRQ(ierr);
	ierr = UserContextDestroy(uc); CHKERRQ(ierr);
}
END_TEST

START_TEST( Curvature_check )
{
	PetscErrorCode ierr;
  int d1 = 50, d2 = 50, len = d1*d2, i, j;
  Vec phi, curv;
  
  PetscReal **data;
  
  ierr = VecCreate(PETSC_COMM_SELF, &phi); CHKERRQ(ierr);
  ierr = VecSetSizes(phi, len, len); CHKERRQ(ierr);
  ierr = VecSetType(phi, VECSEQ ); CHKERRQ(ierr);
  ierr = VecGetArray2d(phi, d1, d2, 0, 0, &data); CHKERRQ(ierr);
  
  for( i = 0; i < d1; ++i)
  {
    for( j = 0; j < d2; ++j)
    {
      data[i][j] = sqrt(PetscSqr(((double)i)-d1/2)+PetscSqr(((double)j)-d2/2));
    }
  }
  ierr = VecRestoreArray2d(phi, d1, d2, 0, 0, &data); CHKERRQ(ierr);
  ierr = VecDuplicate(phi, &curv); CHKERRQ(ierr);
  ierr = Curvature(phi, d1, d2, curv); CHKERRQ(ierr);
  ierr = WriteVector("curvature_check.bin", curv); CHKERRQ(ierr);
}
END_TEST


START_TEST( CalculateC1C2_check )
{
	Image *img;
	short data[9] = {-1,-1,-1,1,1,1,-1,-1,-1};
	PetscReal p[9] = {-1,-1,-1,1,1,1,-1,-1,-1};
	CreateImage(3,3,data,&img);
	Vec phi;
	VecCreate(PETSC_COMM_SELF,&phi);
	VecSetSizes(phi, 9,9);
	VecSetType(phi,VECSEQ);
	VecPlaceArray(phi,p);
	PetscReal c1,c2;
	CalcuateC1C2(phi, img->v, &c1, &c2);
	fail_unless(c1==0.2,"c1 = %f",c1);
	fail_unless(PetscAbsReal(c2- -0.714286)<.001, "c2 = %f", c2);
}
END_TEST

START_TEST( ReadingBinaryFile )
{
  PetscErrorCode ierr;
  
  #if !defined(PETSC_WORDS_BIGENDIAN)
  fail_unless(0,"little endian, cannot read Cascade512 images");
  #endif 
  
  int fd;
  char *name = 
"/home/abergman/Images/image.a.1.gray";
//	"/home/abergman/Images/Joanne/DIR.IGF1 Adaptation-10.spl/image.a.1.gray";
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  short *d;
  PetscMalloc(512*512*sizeof(short),&d);
  int r = read(fd,d,512*512*sizeof(short));
  //PetscByteSwapShort(d,512*512);
//  for(int i = 0; i < 10; i++)
//    PetscPrintf(PETSC_COMM_SELF,"%d\t%d\n",i,d[i]);
  PetscFree(d);
}
END_TEST

START_TEST( ReadCascade512File_check )
{
	PetscErrorCode ierr;
  Image *img;
  
  char *name = "/home/abergman/Images/image.a.1.gray";
//name = "/home/abergman/Images/Joanne/DIR.IGF1 Adaptation-10.spl/image.a.1.gray"; 
  ierr = ReadCascade512File( name, &img); CHKERRQ(ierr);
	
  fail_unless(img->dim_1==512,"Dim 1: %d", img->dim_1 );
  fail_unless(img->dim_2==512,"Dim 2: %d", img->dim_2 );
  fail_unless(img->total_length==512*512,"Total length: %d", img->total_length);
 	
 	PetscReal *a;
	VecGetArray(img->v,&a);
  fail_unless(a!=PETSC_NULL,"array is null");

  PetscReal ans2[10] = {3465, 3288, 3238, 3414, 3344, 3442, 3462, 3313, 3431, 3322};
  for(int i = img->total_length-10; i < img->total_length; i++)
  {
//  	LINE("i:%d\n",i-img->total_length+10);
    fail_unless(ans2[i-img->total_length+10]==a[i],"incorrect endian: %d;%g,%g", i, a[i], ans2[i-img->total_length+10]);
  }
  PetscReal ans[10] = {3252, 3235, 3291, 3307, 3274, 3252, 3313, 3341, 3293, 3269};
  for(int i = 0; i < 10; i++)
    fail_unless(ans[i]==a[i],"incorrect endian: %d;%g,%g", i, a[i], ans[i]);
  VecRestoreArray(img->v, &a);
  DestroyImage(img);
}
END_TEST

START_TEST( CreateImage_check )
{
	PetscErrorCode ierr;
	
	Image *img;
	short *data;
	PetscInt len = 9;
	ierr = PetscMalloc(len*sizeof(PetscReal),&data); CHKERRQ(ierr);
	for( int i = 0; i < len; i++ )
		data[i] = i;
	fail_unless( data[8]==8,"mem allocation failed" );
	
	
	ierr = CreateImage(3,3,data,&img); CHKERRQ(ierr);
	
	fail_unless(img->dim_1==3,"Dim 1: %d", img->dim_1 );
  fail_unless(img->dim_2==3,"Dim 2: %d", img->dim_2 );
  fail_unless(img->total_length==len,"Total length: %d", img->total_length);
	
  PetscInt l;
  VecGetSize(img->v,&l);
  fail_unless(l==len,"lengths");
  PetscReal *arr;
  VecGetArray(img->v,&arr);
  for( int i = 0; i < l; i++)
    fail_unless(arr[i]==((PetscReal)data[i]), "element %d", i);
  VecRestoreArray(img->v,&arr);
  
	PetscFree(data);
	ierr = DestroyImage(img); CHKERRQ(ierr);
  
}
END_TEST

START_TEST( GetVecArrayTwice )
{
  PetscErrorCode ierr;
  PetscInt len=10;
  Vec v;
  VecCreate(PETSC_COMM_WORLD, &v);
  VecSetSizes(v,PETSC_DECIDE,len);
  VecSetType(v, VECSEQ);
  PetscReal *a;
  VecGetArray(v, &a);
  for( int i = 0; i < len; i++)
    a[i] = i/7.;
  VecRestoreArray(v, &a);
  PetscReal *b;
  VecGetArray(v, &b);
  for( int i = 0; i< len; i++)
    fail_unless(i/7.==b[i],"");
  VecRestoreArray(v,&b);
  VecDestroy(v);
}
END_TEST

Suite* suite (void)
{
  Suite *s = suite_create ("Active Contours Without Edges Check");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture( tc_core, setup, teardown);
  tcase_set_timeout(tc_core, 0);
//  tcase_add_test( tc_core, Heaviside_check);
  tcase_add_test( tc_core, CreateImage_check);
//  tcase_add_test( tc_core, ReadCascade512File_check);
//  tcase_add_test( tc_core, ReadingBinaryFile);
  tcase_add_test( tc_core, GetVecArrayTwice );
  tcase_add_test( tc_core, CalculateC1C2_check );
  tcase_add_test( tc_core, Curvature_check );
//  tcase_add_test( tc_core, Single_time_step );
//  tcase_add_test( tc_core, FMM_last_phi );
//  tcase_add_test( tc_core, DA_usage );
  tcase_add_test( tc_core, DMMG_usage );
  
  suite_add_tcase( s, tc_core);

  return s;
}

int main (void)
{ 
  PetscInitializeNoArguments();
  
  int number_failed;
  Suite *s = suite ();
  SRunner *sr = srunner_create (s);
  srunner_set_fork_status(sr, CK_FORK);
//  srunner_run_all (sr, CK_VERBOSE);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
//  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
  
  PetscFinalize();
   return 0;
}


void setup()
{ 
	PetscErrorCode ierr;
}

void teardown()
{
  PetscErrorCode ierr;
}