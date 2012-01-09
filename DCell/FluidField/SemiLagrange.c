#include "FluidField.h"

PetscReal InterpolateField2D(PetscReal  ***field, int dof, const PetscReal *shift, const Coor X );
PetscReal InterpolateField3D(PetscReal ****field, int dof, const PetscReal *shift, const Coor X );

#undef __FUNCT__
#define __FUNCT__ "AdvectSL_2D"
PetscErrorCode AdvectSL_2D( DM daVel, Vec vecVel, DM daBuf, Vec vecBuf, DM daPhi, Vec vecPhi, PetscReal *phiShift, int f, PetscReal dh[3], PetscReal dt )
{
  int i,j;
  PetscReal ***vel;
  PetscReal ***phi;
  PetscReal  **buf;
  Coor X; // Index coordinate
  Coor S; // Characteristic curve
  Coor V; // Interpolated velocity
  DMDALocalInfo a;
  PetscLogDouble t1,t2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);

  Vec lvel;
  ierr = DMGetLocalVector(daVel,&lvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daVel,vecVel,INSERT_VALUES,lvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (daVel,vecVel,INSERT_VALUES,lvel); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(daVel,lvel,&vel); CHKERRQ(ierr);

  Vec lphi;
  ierr = DMGetLocalVector(daPhi,&lphi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daPhi,vecPhi,INSERT_VALUES,lphi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (daPhi,vecPhi,INSERT_VALUES,lphi); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(daPhi,lphi,&phi); CHKERRQ(ierr);

  ierr = DMDAVecGetArray( daBuf,vecBuf,&buf); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(daPhi,&a); CHKERRQ(ierr);

  for (  j = a.ys; j < a.ys+a.ym && j < a.my-2; ++j) {
    for (i = a.xs; i < a.xs+a.xm && i < a.mx-2; ++i) {
      // Staggered grid coordinate
      X.x = i - phiShift[0];
      X.y = j - phiShift[1];

      // Velocity interpolation at grid coordinate
      ierr = InterpolateVelocity2D( U_FACE, vel, X, &V ); CHKERRQ(ierr);

      //TODO: use RK2 like particle LS
      // Euler step
      S.x = X.x - V.y * dt / dh[0];
      S.y = X.y - V.y * dt / dh[1];

      // Test if trajectory within local grid and global domain
      if( a.xs-1 <= S.x && S.x < a.xs+a.xm && 0 <= S.x && S.x <= a.mx-2 &&
          a.ys-1 <= S.y && S.y < a.ys+a.ym && 0 <= S.y && S.y <= a.my-2 ){

        // Strain interpolation at take off point
        buf[j][i] = InterpolateField2D( phi, f, phiShift, S );

      } else {
        buf[j][i] = 0;
//        SETERRQ4(PETSC_ERR_USER,"Need to implement BC [%d, %d], (%f, %f)", i, j, S[0], S[1]);
      }
    }  // i
  }   // j

  ierr = DMDAVecRestoreArray(daPhi, lphi, &phi); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daVel, lvel, &vel); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daBuf, vecBuf, &buf); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daVel,&lvel); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daPhi,&lphi); CHKERRQ(ierr);

  // Scatter single component parallel buffer into multi-component parallel vector.
  ierr = VecStrideScatter(vecBuf,f,vecPhi,INSERT_VALUES); CHKERRQ(ierr);

  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Advect SL: %f sec\n", t2-t1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AdvectSL_3D"
PetscErrorCode AdvectSL_3D( DM daVel, Vec vecVel, DM daBuf, Vec vecBuf, DM daPhi, Vec vecPhi, PetscReal *phiShift, int f, PetscReal dh[3], PetscReal dt )
{
  int i,j,k;
  PetscReal ****vel;
  PetscReal ****phi;
  PetscReal  ***buf;
  Coor X; // Index coordinate
  Coor S; // Characteristic curve
  Coor V; // Interpolated velocity
  DMDALocalInfo a;
  PetscLogDouble t1,t2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);

  Vec lvel;
  ierr = DMGetLocalVector(daVel,&lvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daVel,vecVel,INSERT_VALUES,lvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(daVel,vecVel,INSERT_VALUES,lvel); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(daVel,lvel,&vel); CHKERRQ(ierr);

  Vec lphi;
  ierr = DMGetLocalVector(daPhi,&lphi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daPhi,vecPhi,INSERT_VALUES,lphi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (daPhi,vecPhi,INSERT_VALUES,lphi); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(daPhi,lphi,&phi); CHKERRQ(ierr);

  ierr = DMDAVecGetArray( daBuf,vecBuf,&buf); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(daPhi,&a); CHKERRQ(ierr);

  for     (k = a.zs; k < a.zs+a.zm && k < a.mz-2; ++k) {
    for   (j = a.ys; j < a.ys+a.ym && j < a.my-2; ++j) {
      for (i = a.xs; i < a.xs+a.xm && i < a.mx-2; ++i) {
        // Staggered grid coordinate
        X.x = i - phiShift[0];
        X.y = j - phiShift[1];
        X.z = k - phiShift[2];

        // Velocity interpolation at grid coordinate
        ierr = InterpolateVelocity3D( U_FACE, vel, X, &V ); CHKERRQ(ierr);

        // Euler step
        S.x = X.x - V.x * dt / dh[0];
        S.y = X.y - V.y * dt / dh[1];
        S.z = X.z - V.z * dt / dh[2];

        // Test if trajectory within local grid and global domain
        if( a.xs-1 <= S.x && S.x < a.xs+a.xm && 0 <= S.x && S.x <= a.mx-2 &&
            a.ys-1 <= S.y && S.y < a.ys+a.ym && 0 <= S.y && S.y <= a.my-2 &&
            a.zs-1 <= S.z && S.z < a.zs+a.zm && 0 <= S.z && S.z <= a.mz-2 ){

// Strain interpolation at take off point
          buf[k][j][i] = InterpolateField3D( phi, f, phiShift, S );
        } else {
          buf[k][j][i] = 0;
        }
      }  // i
    }   // j
  }    // k
  ierr = DMDAVecRestoreArray(daPhi, lphi, &phi); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daVel, lvel, &vel); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daBuf, vecBuf, &buf); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daVel,&lvel); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daPhi,&lphi); CHKERRQ(ierr);

  // Scatter single component parallel buffer into multi-component parallel vector.
  ierr = VecStrideScatter(vecBuf,f,vecPhi,INSERT_VALUES); CHKERRQ(ierr);

  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Advect SL: %f sec\n", t2-t1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscReal InterpolateField2D(PetscReal ***field, int dof, const PetscReal *shift, const Coor X )
{
  int i,j;
  int xs = (int)floor(X.x + shift[0]),
      ys = (int)floor(X.y + shift[1]);
  PetscReal sx = X.x - xs,
            sy = X.y - ys;
  PetscReal sum = 0;

  for( j = 0; j < 2; ++j)
  {
    for( i = 0; i < 2; ++i)
    {
      sum += ((1 - sx) * (1 - i) + i * sx) *
             ((1 - sy) * (1 - j) + j * sy) * field[ys+j][xs+i][dof];
    }
  }

  return sum;
}

PetscReal InterpolateField3D(PetscReal ****field, int dof, const PetscReal *shift, const Coor X )
{
  int i,j,k;
  int xs = (int)floor(X.x + shift[0]),
      ys = (int)floor(X.y + shift[1]),
      zs = (int)floor(X.z + shift[2]);
  PetscReal sx = X.x - xs,
            sy = X.y - ys,
            sz = X.z - zs;
  PetscReal sum = 0;

  for( k = 0; k < 2; ++k)
  {
    for( j = 0; j < 2; ++j)
    {
      for( i = 0; i < 2; ++i)
      {
        sum += ((1 - sx) * (1 - i) + i * sx) *
               ((1 - sy) * (1 - j) + j * sy) *
               ((1 - sz) * (1 - k) + k * sz) * field[zs+k][ys+j][xs+i][dof];
      }
    }
  }

  return sum;
}
