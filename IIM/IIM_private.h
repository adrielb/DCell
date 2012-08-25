//#ifndef IIM_PRIVATE_H_
//#define IIM_PRIVATE_H_

struct _LocalCoor {
  Array coor;
};

typedef struct IIMDebug {
  PetscReal h;
  PetscReal dd;
  PetscReal j, jp, jpp;
  PetscReal C;
  int dof;
  int sign;
  int axis;
  iCoor X;
} IIMDebug;

PetscLogEvent EVENT_IIMUpdateIrregularNodeGrid;
PetscLogEvent EVENT_IIMUpdateSurfaceQuantities;
PetscLogEvent EVENT_IIMUpdateSurfaceDerivatives;
PetscLogEvent EVENT_IIMUpdateRHS;
PetscLogEvent EVENT_IIMIrregularNodes;
PetscErrorCode IIMRegisterEvents( void );

void IIMRotation1st_2D( IrregularNode *n, Jump *j );
void IIMRotation2nd_2D( IrregularNode *n, Jump *j );

//#endif /* IIM_PRIVATE_H_ */
