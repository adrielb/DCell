#include "LevelSetMethod.h"
#include "LSM_private.h"

/* Trying static variables for val, coor, idx arrays and
 * velocity grid.
 */
#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvect"
PetscErrorCode LevelSetAdvect( LevelSet ls, int ga, PetscReal dt )
{
  static Grid velgrid;
//static int count = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( velgrid == NULL ) {
    Grid phi = ls->phi;
    int dof = phi->is2D?2:3;
    ierr = GridCreate(phi->d,phi->p,phi->n,dof,&velgrid); CHKERRQ(ierr);
    ierr = GridSetName(velgrid,"velgrid"); CHKERRQ(ierr);
  }
  ierr = LevelSetGetVelocity( ls, ga, velgrid ); CHKERRQ(ierr);
//  ierr = LevelSetGatherVelocity( ls, ga, velgrid ); CHKERRQ(ierr);

  ierr = GridWrite(velgrid, 0); CHKERRQ(ierr);

  ls->Advect( ls, velgrid, dt);

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectAndReinit"
PetscErrorCode LevelSetAdvectAndReinit(LevelSet ls, Grid velgrid, PetscReal dt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridCopy(ls->phi,ls->phi0); CHKERRQ(ierr);
  ierr = LevelSetAdvectSL( ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = LevelSetCFLIncrement( ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  if( ls->AdvectCount > ls->AdvectThres || ls->CFLcount > ls->CFLthres ) {
    // TODO: specifically say which level set object is reinitializing (ls->ID)
    ierr = PetscInfo4(0,"CFLcount: %f:%f  AdvectCount: %d:%d\n", ls->CFLcount, ls->CFLthres, ls->AdvectCount,ls->AdvectThres); CHKERRQ(ierr);
    ls->CFLcount = 0;
    ls->AdvectCount = 0;
    ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectSL"
PetscErrorCode LevelSetAdvectSL(LevelSet ls, Grid velgrid, PetscReal dt)
{
  int b;
  iCoor *band;
  Coor X; // Grid point
  Coor S; // Projected grid point
  Coor V; // Interpolated velocity
  Coor dh = ls->phi->d;
  PetscReal ***vel, **phi, **phi0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->phi0,&phi0); CHKERRQ(ierr);

  for ( b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,b,&band); CHKERRQ(ierr);
    X.x = band->x;
    X.y = band->y;
    // Interpolate velocity at grid point
    ierr = InterpolateVelocity2D( 0, vel, X, &V ); CHKERRQ(ierr);
    // Project grid point back in time
    S.x = X.x - V.x * dt / dh.x;
    S.y = X.y - V.y * dt / dh.y;
    // Interpolate phi at projection
    phi[band->y][band->x] = Bilinear2D(GridFunction2D_Identity,phi0,dh, S.x, S.y );
  } // for b in band
  ierr = PetscLogEventEnd(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectSL_3D"
PetscErrorCode LevelSetAdvectSL_3D(LevelSet ls, Grid velgrid, PetscReal dt)
{
  int b;
  iCoor *band;
  Coor X; // Grid point
  Coor S; // Projected grid point
  Coor V; // Interpolated velocity
  Coor dh = ls->phi->d;
  PetscReal ****vel, ***phi, ***tmp;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->tmp,&tmp); CHKERRQ(ierr);
  ierr = VecCopy(ls->phi->v,ls->tmp->v); CHKERRQ(ierr);
  for ( b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,b,&band); CHKERRQ(ierr);
    X.x = band->x;
    X.y = band->y;
    X.z = band->z;

    ierr = InterpolateVelocity3D( 0, vel, X, &V ); CHKERRQ(ierr);

    S.x = X.x - V.x * dt / dh.x;
    S.y = X.y - V.y * dt / dh.y;
    S.z = X.z - V.z * dt / dh.z;

    phi[band->z][band->y][band->x] = Bilinear3D(GridFunction3D_Identity,tmp,dh, S );
  } // for b in band

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetCFLIncrement"
PetscErrorCode LevelSetCFLIncrement( LevelSet ls, Grid velgrid, PetscReal dt )
{
  int i;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  PetscReal mag,maxVel = 0;
  PetscReal ***vel;
  Coor X,V;
  Coor dh = ls->phi->d;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  for ( i = 0; i < len; ++i) {
    X.x = nodes[i].pos.x;
    X.y = nodes[i].pos.y;
    // Interpolate velocity at grid point
    ierr = InterpolateVelocity2D( 0, vel, X, &V ); CHKERRQ(ierr);

    mag = PetscSqrtScalar( V.x*V.x + V.y*V.y );

    maxVel = mag > maxVel ? mag : maxVel;
  }

  ls->CFLcount += dt * maxVel / dh.x;
  ls->AdvectCount++;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetGetVelocityGA"
PetscErrorCode LevelSetGetVelocityGA(LevelSet ls, int ga, Grid velgrid)
{
  int lo[4],hi[4],ld[3] = {0,0,0};
  iCoor p,q;
  int type,dim,dims[4];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  NGA_Inquire(ga,&type,&dim,dims);
  ierr = GridResize(velgrid,ls->phi->p,ls->phi->n); CHKERRQ(ierr);
  ierr = GridGetBounds(velgrid, &p, &q); CHKERRQ(ierr);
  lo[0] = p.y < 0 ? 0 : p.y;
  lo[1] = p.x < 0 ? 0 : p.x;
  lo[2] = 0;
  hi[0] = q.y-1;// > dims[0] ? dims[0] : q.y-1;
  hi[1] = q.x-1;// > dims[1] ? dims[1] : q.x-1;
  hi[2] = 1;

  printf("%d %d %d\n", lo[0], lo[1], lo[2] );
  printf("%d %d %d\n", hi[0], hi[1], hi[2] );

  ld[0] = q.x-p.x;
  ld[1] = 2;
  NGA_Get(ga,lo,hi,velgrid->v1,ld);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetGetVelocity"
PetscErrorCode LevelSetGetVelocity(LevelSet ls, int ga, Grid velgrid)
{
  int lo[4],hi[4],ld[3] = {0,0,0};
  iCoor p,q;
  int type,dim,dims[4];
  PetscReal ***grid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetGetVelocity,0,0,0,0); CHKERRQ(ierr);
  NGA_Inquire(ga,&type,&dim,dims);
  ierr = GridResize(velgrid,ls->phi->p,ls->phi->n); CHKERRQ(ierr);
  ierr = GridGetBounds(velgrid, &p, &q); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&grid); CHKERRQ(ierr);
  lo[0] = p.y < 0 ? 0 : p.y;
  lo[1] = p.x < 0 ? 0 : p.x;
  lo[2] = U_FACE;
  hi[0] = q.y-1 > dims[0]-1 ? dims[0]-1 : q.y-1;
  hi[1] = q.x-1 > dims[1]-1 ? dims[1]-1 : q.x-1;
  hi[2] = V_FACE;
/*
  printf("%d %d %d\n", lo[0], lo[1], lo[2] );
  printf("%d %d %d\n", hi[0], hi[1], hi[2] );
  printf("%d %d %d\n", dims[0], dims[1], dims[2] );
*/
  ld[0] = q.x - p.x;
  ld[1] = V_FACE;
  NGA_Get(ga,lo,hi,&grid[lo[0]][lo[1]][0], ld);
  ierr = PetscLogEventEnd(EVENT_LevelSetGetVelocity,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//TODO: Try block ga_get instead of ga_gather.
//      Gather seems to have too much overhead since cell volume too small.
#undef __FUNCT__
#define __FUNCT__ "LevelSetGatherVelocity"
PetscErrorCode LevelSetGatherVelocity(LevelSet ls, int ga, Grid velgrid)
{
  static Array val, coor, idx;
  PetscReal     *v;
  int *c, **i;
  int b,d,j;
  const int dof = velgrid->dof;
  const int blen = ArrayLength(ls->band);
  int len = dof * blen;
  iCoor *band;
  int type,dim,dims[4];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( val == NULL ) {
    ierr = ArrayCreate("vel",sizeof(PetscReal),len,&val); CHKERRQ(ierr);
    ierr = ArrayCreate("idx",sizeof(int*)     ,len,&idx); CHKERRQ(ierr);
    ierr = ArrayCreate("coor",sizeof(int)*4    ,len,&coor); CHKERRQ(ierr); // [z,y,x,d]
  }
  ierr = ArraySetSize(val,len); CHKERRQ(ierr);
  ierr = ArraySetSize(idx,len); CHKERRQ(ierr);
  ierr = ArraySetSize(coor,len); CHKERRQ(ierr);
  band = ArrayGetData(ls->band);
  i = ArrayGetData(idx);
  j = 0;
  NGA_Inquire(ga,&type,&dim,dims);

  for (b = 0; b < blen; ++b) {
    if( band[b].x <= 0 || dims[1]-1 < band[b].x ||
        band[b].y <= 0 || dims[0]-1 < band[b].y )
      continue;
    for (d = 0; d < dof; ++d) {
      ierr = ArrayGet(coor,j,&c); CHKERRQ(ierr);
      c[0] = d;
      c[1] = band[b].x;
      c[2] = band[b].y;
      c[3] = band[b].z;
      i[j] = c;
      j++;
    }
  }
  len = j;
  ierr = ArraySetSize(val,len); CHKERRQ(ierr);
  ierr = ArraySetSize(idx,len); CHKERRQ(ierr);
  ierr = ArraySetSize(coor,len); CHKERRQ(ierr);
  ierr = GAGather(ga,idx,val); CHKERRQ(ierr);
  ierr = GridResize(velgrid,ls->phi->p,ls->phi->n); CHKERRQ(ierr);

  v = ArrayGetData(val);
  if( ls->phi->is2D ) {
    PetscReal ***vel;
    ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
    for ( j = 0; j < len; j++) {
      ierr = ArrayGet(coor,j,&c); CHKERRQ(ierr);
      vel[c[2]][c[1]][c[0]] = v[j];
    }
  } else {
    PetscReal ****vel;
    ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
    for ( j = 0; j < len; j++) {
      ierr = ArrayGet(coor,j,&c); CHKERRQ(ierr);
      vel[c[3]][c[2]][c[1]][c[0]] = v[j];
    }
  }
//  static int count;
//  ierr = GridSetName(velgrid,"gavel"); CHKERRQ(ierr);
//  ierr = GridWrite(velgrid,count++); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
