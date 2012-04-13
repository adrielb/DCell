#include "Main.h"
#include "MyCheck.h"
#include "ImmersedInterfaceMethod.h"
#include "gts.h"
#include "LevelSetMethod.h"
#include "PoissonSolver.h"

START_TEST( viz_LocalStencil )
{
  FILE *file, *afile;
  file = fopen( "/home/abergman/Research/DCell/temp/viz_LocalStencil","w");
  afile = fopen( "/home/abergman/Research/DCell/temp/viz_LocalStencilACoeff","w");
  int i, j;
  IrregularNode *n;
  IIM iim;
  IIMCreate(&iim, 12);
  PetscReal *eta, *xi;
  PetscInt d1=63,d2=d1;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
  LevelSetInitializeToCircle(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  PetscReal **phi = ls->g2d->v2;
  PetscReal a2,a4,a6,a8,a10,a12; 
  PetscErrorCode ierr;

  IIMUpdateIrregularNodes( iim , ls);
  LocalCoor2DGetVecs(iim->lc, &eta, &xi);
  IrregularNodeListWrite(ls,0);
fprintf(file,"{");
fprintf(afile,"{");
  for( i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    //Rotate the stencil into local coor space (eta, xi)
    LocalCoor2DSolveStencil( iim->lc, n );
    
    //calc a2 - a12 from 3.17
    a2=a4=a6=a8=a10=a12=0;
    for( j = 0; j < 5; j++)
    {
      if( phi[Stencil2Dx[j]+n->x][Stencil2Dy[j]+n->y] < 0)
        continue; // determine which points in the stencil are inside the interface
      a2 += DiscretePoisson[j];
      a4 += DiscretePoisson[j] * eta[j];
      a6 += DiscretePoisson[j] *  xi[j];
      a8 += DiscretePoisson[j] * eta[j] * eta[j];
      a10+= DiscretePoisson[j] *  xi[j] *  xi[j];
      a12+= DiscretePoisson[j] * eta[j] *  xi[j];
    }
    a8  /= two;
    a10 /= two;
    
    fprintf(afile, "{%f,%f,%f,%f,%f,%f},",a2,a4,a6,a8,a10,a12);
    
    fprintf(file, "{");
    for( j = 0; j < 5; ++j)
    {
      if( phi[Stencil2Dx[j]+n->x][Stencil2Dy[j]+n->y] < 0)
        continue; // determine which points in the stencil are inside the interface
      fprintf(file, "{{%f,%f},{%d,%d}},",eta[j],xi[j],Stencil2Dx[j],Stencil2Dy[j]);
    }
    fprintf(file, "{}},\n");
  }
fprintf(file,"{}}");
fprintf(afile,"{}}");
}
END_TEST

START_TEST( test_LocalCoorSolve )
{
  PetscReal *eta, *xi;
  IrregularNode n;
  PetscErrorCode ierr;
  n.x=n.ox=n.y=n.oy=0;
  n.nx=0.995004; n.ny=0.0998334;
  LocalCoor2D lc;
  LocalCoor2DCreate(12, &lc);
  LocalCoor2DGetVecs(lc, &eta, &xi);
  for( int i = 0; i < 6; i++ )
  {
    eta[i] = i;
    xi[i]  = 0;
  }
  eta[0]=n.nx; xi[0]=n.ny;
  LocalCoor2DSolve(lc, &n);
  for( int i = 0; i < 6; i++ )
    printf("{%f,%f},",eta[i],xi[i]);

  LocalCoor2DDestroy(lc);
}
END_TEST

START_TEST( viz_SurfaceDerivatives )
{
  FILE *file;
  file = fopen( "/home/abergman/Research/DCell/temp/viz_SurfaceDerivatives","w");
  PetscErrorCode ierr;
  int i, j, len;
  PetscReal *eta, *xi, *s, *g, x[12],y[12], f1[12];
  PetscReal w,dw,ddw,v,dv;
  IrregularNode *n, *node;
  GSList *slist, *iter;
  
  IIM iim;
  ierr = IIMCreate(&iim, 12); CHKERRQ(ierr);
  
  int d1=64, d2=d1;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
//  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  LevelSetInitializeToCircle(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  
  IIMUpdateIrregularNodes( iim , ls);
  LocalCoor2DGetVecs(iim->lc, &eta, &xi);
  LeastSqGetVecs(iim->lsq, &s, &g, PETSC_NULL);
fprintf(file,"{");
  for( i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    KDTreeRange(iim->kdtree, n, &slist );
    iter = slist;
    for( j = 0; j < iim->Np; j++ )
    {
      if( !iter ) break;
      node = KDTreeGetIrregularNode(iter);
      eta[j] = node->x + node->ox;
      xi[j]  = node->y + node->oy;
      f1[j] = node->f1;
      x[j] = eta[j] - (n->x + n->ox);
      y[j] = xi[j]  - (n->y + n->oy);
      iter = g_slist_next(iter);
    }
    len = j;
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);

    // redefine eta and xi to local coordiates
    LocalCoor2DSetLength(iim->lc, len);
    LocalCoor2DSolve( iim->lc, n);
    LocalCoor2DToArcLength(iim->lc, s);
    
    JumpConditionPressure( iim->lsq, n, slist, &w, &dw, &ddw, &v, &dv );
    
    iter = slist;
fprintf(file,"{");
    for( j = 0; j < len; j++ )
    {                                         
      node = KDTreeGetIrregularNode(iter);  //  0     1     2       3      4     5
      fprintf(file,"{%f, %f, %f, %f, %f, %f},", x[j], y[j], eta[j], xi[j], s[j], f1[j]);
      iter = g_slist_next(iter);
    }
//    fprintf(file, "{%f, %f, %f, %f, %f}}",w,dw,ddw,v,dv);
    fprintf(file, "{}}");
    if( i != ls->irregularNodes->len - 1)
      fprintf(file, ",\n");
  }
fprintf(file,"}");
  fclose(file);
  IrregularNodeListWrite( ls, 0);
  DestroyLevelSet2D(ls);
  IIMDestroy(iim);
}
END_TEST

START_TEST( test_IIMContext )
{
  IIM iim;
  int d1=63, d2=d1;
  Grid2D rhsC, p, px, py, u, v;
  LevelSet2D ls, lstemp;
  PetscErrorCode ierr;
  
  Mat m;
  KSP ksp;
  PC pc;
  Generate2DLaplacianPeriodicBC( d1, d2, &m);
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp); CHKERRQ(ierr);
  KSPGetPC(ksp, &pc);
  
//  KSPSetType(ksp, KSPCG);
//  PCSetType(pc, PCICC);
//  PCFactorSetLevels(pc, 0);
//  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  
  KSPSetType(ksp, KSPPREONLY);
  PCSetType(pc, PCCHOLESKY);
  PCFactorSetMatOrderingType(pc, MATORDERING_ND);
  
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp, m, m, SAME_PRECONDITIONER);
  PCSetUp(pc);
  
  ierr = IIMCreate(&iim, 12); CHKERRQ(ierr);
//  IIMSetForceComponents(iim, ForceComponentNormalSimple, ForceComponentTangentialSimple);
  CreateGrid2D(d1, d2, &rhsC);
  CreateGrid2D(d1, d2, &p);
  CreateGrid2D(d1, d2, &px);
  CreateGrid2D(d1, d2, &py);
  CreateGrid2D(d1, d2, &u);
  CreateGrid2D(d1, d2, &v);
  ierr = CreateLevelSet2D(d1,d2,&ls); CHKERRQ(ierr);  
  CreateLevelSet2D(d1, d2, &lstemp);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
//  LevelSetInitializeToCircle(ls, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE);
  
  Vec u2,v2,magVel;
  VecCreateSeq(PETSC_COMM_SELF,d1*d2,&magVel);
  VecDuplicate(magVel, &u2);
  VecDuplicate(magVel, &v2);
  PetscReal maxVel, dt;
  
  WriteVec wv, wvC, wvP, wvPX, wvPY, wvU, wvV;
  WriteVecCreate(lstemp->g2d->v, "ls",&wv);
  WriteVecCreate(rhsC->v, "c", &wvC);
  WriteVecCreate(p->v, "p", &wvP);
  WriteVecCreate(px->v, "px", &wvPY);
  WriteVecCreate(py->v, "py", &wvPX);
  WriteVecCreate(u->v, "u", &wvU);
  WriteVecCreate(v->v, "v", &wvV);
  
  for( int t = 0; t < 100; t++ )
  {
printf("TIME: %d ",t);
    VecSet(rhsC->v, 0.); // TODO: this is important when reusing vecs, need to reset back to zero after time step
    ierr = IIMComputeCorrection2D(iim,JumpConditionPressure,ls,rhsC); CHKERRQ(ierr);
    IIMDiscreteCompatabilityCondition(ls, rhsC);
    ierr = KSPSolve(ksp, rhsC->v, p->v); CHKERRQ(ierr);
    IrregularNodeListWrite(ls,t);
    
    ierr = IIMPressureGradient(ls,p,px,py); CHKERRQ(ierr);
    
    ierr = IIMComputeCorrection2D(iim,JumpConditionXVelocity,ls,px); CHKERRQ(ierr);
    px->v2[d1/2][d2/2] = 0.;
    ierr = KSPSolve(ksp, px->v, u->v); CHKERRQ(ierr);
    
    ierr = IIMComputeCorrection2D(iim,JumpConditionYVelocity,ls,py); CHKERRQ(ierr);
    py->v2[d1/2][d2/2] = 0.;
    ierr = KSPSolve(ksp, py->v, v->v); CHKERRQ(ierr);
     
    VecPointwiseMult(u2, u->v, u->v);
    VecPointwiseMult(v2, u->v, u->v);
    VecWAXPY(magVel, 1, u2, v2);
    VecSqrt( magVel );

    VecNorm(magVel, NORM_INFINITY, &maxVel);
    dt = PetscMin(1/(2*maxVel),.1);
    LevelSetAdvect2D(dt, u, v, ls, lstemp);
printf("\tmaxVel: %f",maxVel);
printf("\tdt: %f\n",dt);

    
    WriteVecToDisk(wv);
    WriteVecToDisk(wvC);
    WriteVecToDisk(wvP);
    WriteVecToDisk(wvPX);
    WriteVecToDisk(wvPY);
    WriteVecToDisk(wvU);
    WriteVecToDisk(wvV);
    VecCopy( lstemp->g2d->v, ls->g2d->v);
    UpdateIrregularNodeList(ls);
//    ReinitializeLevelSet(ls);
  }
  
  ierr = VecDestroy(magVel); CHKERRQ(ierr);
  ierr = VecDestroy(u2); CHKERRQ(ierr);
  ierr = VecDestroy(v2); CHKERRQ(ierr);
  ierr = DestroyLevelSet2D(ls); CHKERRQ(ierr);
  ierr = DestroyGrid2D(rhsC); CHKERRQ(ierr);
  ierr = DestroyGrid2D(p); CHKERRQ(ierr);
  ierr = DestroyGrid2D(px); CHKERRQ(ierr);
  ierr = DestroyGrid2D(py); CHKERRQ(ierr);
  ierr = DestroyGrid2D(u); CHKERRQ(ierr);
  ierr = DestroyGrid2D(v); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
}
END_TEST

START_TEST( test_NormalDirection )
{
  PetscErrorCode ierr;
  int d1 = 210, d2 = 210;
  LevelSet2D ls;
  
  ierr = CreateLevelSet2D(d1, d2, &ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToCircle(ls, 10, 5, 3.2); CHKERRQ(ierr);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  ierr = WriteVector("initls",ls->g2d->v); CHKERRQ(ierr);
  ierr = UpdateIrregularNodeList(ls); CHKERRQ(ierr);
  GArray *ga = ls->irregularNodes;
  LocalCoor2D lc;
  ierr = LocalCoor2DCreate(12, &lc); CHKERRQ(ierr);
  IrregularNode *n;
  for( int i = 0; i < ga->len; i++ )
  {
    n = &g_array_index( ga, IrregularNode, i); 
    PetscReal nx, ny, h;
      
    nx = Bilinear2D(GridFunction2D_DerivX, ls->g2d, n->x+n->ox, n->y+n->oy);
    ny = Bilinear2D(GridFunction2D_DerivY, ls->g2d, n->x+n->ox, n->y+n->oy);
    h = sqrt( nx*nx + ny*ny );
      
    n->nx = nx / h;
    n->ny = ny / h;
  }
  ierr = IrregularNodeListWrite(ls,0); CHKERRQ(ierr);
}
END_TEST

START_TEST( testGTSCast )
{
  struct _MyPoint {
    GtsObject object;

    gdouble x, y, z;
    
    void *data;
  };
  
  struct _MyPoint *p;
  PetscNew(struct _MyPoint, &p );
  
  gts_point_set((GtsPoint*)p,1,2,3);
  printf("%f, %f, %f\n", p->x,p->y,p->z);
  PetscFree(p);
}
END_TEST

START_TEST( viz_KDTree )
{
  PetscErrorCode ierr;
  int d1=64,d2=d1;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  KDTree kdt;
  ierr = KDTreeCreate(ls->irregularNodes, &kdt); CHKERRQ(ierr);
  GSList *slist, *iter;
  IrregularNode *n, *p;
  
  WriteVector("ls",ls->g2d->v);
  IrregularNodeListWrite(ls,0);
  LevelSet2D lskd;
  CreateLevelSet2D(d1,d2,&lskd);
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index(ls->irregularNodes,IrregularNode,i);
    ierr = KDTreeRange(kdt, n, &slist); CHKERRQ(ierr);
    iter = slist;
    while( iter  )
    {
      p = KDTreeGetIrregularNode(iter);
      g_array_append_val(lskd->irregularNodes,*p);
      iter = g_slist_next(iter);
    }
    IrregularNodeListWrite(lskd,i+1);
    g_array_remove_range(lskd->irregularNodes, 0, lskd->irregularNodes->len);
  }
    
  
  ierr = KDTreeDestroy(kdt); CHKERRQ(ierr);
}
END_TEST

START_TEST( test_KDTreeCreate )
{
  PetscErrorCode ierr;
  IrregularNode n;
  n.x=1;n.y=1;n.z=1;
  GArray *g = g_array_new( FALSE, FALSE, sizeof(IrregularNode) );
  g_array_append_val(g, n);
  
  KDTree kdt;
  
  ierr = KDTreeCreate( g, &kdt); CHKERRQ(ierr);
  mark_point();
  ierr = KDTreeDestroy(kdt); CHKERRQ(ierr);
}
END_TEST


START_TEST( Stencil2D_pointers )
{
  PetscErrorCode ierr;
  
  for( int j = 0; j < 5; j++ )
  {
    fail_unless( Stencil2D[0][j] == Stencil2Dx[j], "stencil x pointer failed");
    fail_unless( Stencil2D[1][j] == Stencil2Dy[j], "stencil y pointer failed");
  }
}
END_TEST

Suite* CreateSuite (void)
{
  Suite *s = suite_create ("IIM Check");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core,  Stencil2D_pointers );
  tcase_add_test( tc_core, test_KDTreeCreate );
//  tcase_add_test( tc_core, testGTSCast );
//  tcase_add_test( tc_core, test_NormalDirection );
//  tcase_add_test( tc_core,  viz_KDTree );
  tcase_add_test( tc_core, test_IIMContext );
//  tcase_add_test( tc_core,  viz_SurfaceDerivatives );
//  tcase_add_test( tc_core,  test_LocalCoorSolve );
//  tcase_add_test( tc_core,  viz_LocalStencil );
  suite_add_tcase( s, tc_core);
  return s;
}

int RunCheck()
{ 
  Suite *s = CreateSuite ();
  SRunner *sr = srunner_create (s);
  
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}

PetscErrorCode PetscMain()
{
  PetscFunctionReturn(0);
}