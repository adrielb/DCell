#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

void MyForce( IrregularNode *n, void *ctx );
PetscErrorCode WriteCorrections( IIM iim, int i, const char du[4], int sign, iCoor *band );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1;
  Coor dh = {dx, dx, 0};
  Jump jump ;

  const char *phinames[3] = {"phi_-1", "", "phi_1"};
  LevelSet levelsets[3];
  ierr = LevelSetCreate(dh,(iCoor){-2,-2,0},(iCoor){6,6,0},&levelsets[0]); CHKERRQ(ierr);
  ierr = LevelSetCreate(dh,(iCoor){-2,-2,0},(iCoor){6,6,0},&levelsets[2]); CHKERRQ(ierr);
  ierr = GridSetName(levelsets[0]->phi,phinames[0]); CHKERRQ(ierr);
  ierr = GridSetName(levelsets[2]->phi,phinames[2]); CHKERRQ(ierr);
  PetscReal **phi[3];
  ierr = GridGet( levelsets[0]->phi, &phi[0]); CHKERRQ(ierr);
  ierr = GridGet( levelsets[2]->phi, &phi[2]); CHKERRQ(ierr);

  iCoor box[4] = { {.x = 0, .y = 0 },
                   {.x = 1, .y = 0 },
                   {.x = 0, .y = 1 },
                   {.x = 1, .y = 1 } };

  IIM iim;
  ierr = IIMCreate(PETSC_TRUE,20,dh,&iim); CHKERRQ(ierr);
  ierr = IIMSetViscosity(iim, 1.0); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim, MyForce ); CHKERRQ(ierr);

  iCoor *band;
  IrregularNode *node;
  int i, k, b, sign, s;
  LevelSet ls;
  const int len = 5;

  for( sign = -1; sign <= 1; sign += 2 ) {
    s = sign + 1;
    ls = levelsets[s];
    ierr = VecSet( ls->phi->v, -1*sign ); CHKERRQ(ierr);
    for (b = 0; b < 4; ++b) {
      ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
      band->x = box[b].x;
      band->y = box[b].y;
      for( i = 0; i <= len; i++ ) {
        phi[s][0][0] = sign * ( (len-i)*(0.2) + i*(5) ) / len;
        ierr = GridWrite(ls->phi, i); CHKERRQ(ierr);
        ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
        ierr = IIMUpdateSurfaceQuantities( iim, ls); CHKERRQ(ierr);
        ierr = LevelSetWriteIrregularNodeList(ls, i); CHKERRQ(ierr);

        // C{ px = p[i] - p[i-1] }
        for (k = 0; k < ArrayLength(ls->irregularNodes); ++k) {
          ierr = ArrayGet(ls->irregularNodes,k,(void*)&node); CHKERRQ(ierr);
          if( node->axis == -1 ) continue;
          JumpPressure( node, &jump );
          ierr = IIMPressureGradientCorrection( iim, node, jump ); CHKERRQ(ierr);
        }
        ierr = WriteCorrections(iim, i, "dp", sign, band); CHKERRQ(ierr);

        // C{ laplace(u) }
        for (k = 0; k < ArrayLength(ls->irregularNodes); ++k) {
          ierr = ArrayGet(ls->irregularNodes,k,(void*)&node); CHKERRQ(ierr);
          if( node->axis == -1 ) continue;
          int vel = node->shift == CELL_CENTER ? node->axis : node->shift - U_FACE;
          JumpVelocity( iim->mu, node, &jump, vel );
          ierr = IIMLaplaceCorrection(iim, node, jump); CHKERRQ(ierr);
        }
        ierr = WriteCorrections(iim, i, "ddu", sign, band); CHKERRQ(ierr);

        // C{ ux = u[i+1] - u[i] }
        for (k = 0; k < ArrayLength(ls->irregularNodes); ++k) {
          ierr = ArrayGet(ls->irregularNodes,k,(void*)&node); CHKERRQ(ierr);
          if( node->axis == -1 ) continue;
          int vel = node->shift == CELL_CENTER ? node->axis : node->shift - U_FACE;
          JumpVelocity( iim->mu, node, &jump, vel );
          ierr = IIMVelocityGradientCorrection( iim, node, jump ); CHKERRQ(ierr);
        }
        ierr = WriteCorrections(iim, i, "du", sign, band); CHKERRQ(ierr);

      } // phi[0][0]
      ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);
    } // band [x,y]
  } // sign

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteCorrections"
PetscErrorCode WriteCorrections( IIM iim, int i, const char du[4], int sign, iCoor *b )
{
  char name[64];
  PetscErrorCode  ierr;
  PetscFunctionBegin;

  sprintf(name,  "%s_%d_%d_%d", du, sign, b->x, b->y );
  ierr = ArraySetName(iim->debug, name); CHKERRQ(ierr);
  ierr = ArrayWrite(  iim->debug, i); CHKERRQ(ierr);
  ierr = ArraySetSize(iim->debug, 0); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

void MyForce( IrregularNode *n, void *ctx ) {
  n->F1 = 100;
}
