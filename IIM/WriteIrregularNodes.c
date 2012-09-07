#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMWriteIrregularNodeList"
PetscErrorCode IIMWriteIrregularNodeList( IIM iim, char *name, int idx )
{
  int i;
  int rowlen;
  int len = ArrayLength( iim->irregularNodes );
  IIMIrregularNode *nodes = ArrayGetData( iim->irregularNodes );
  IIMIrregularNode *node;
  char filename[PETSC_MAX_PATH_LEN], wd[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  PetscErrorCode ierr;

  PetscFunctionBegin;
//  ierr = PetscLogEventBegin(EVENT_LevelSetWriteIrregularNodeList,0,0,0,0); CHKERRQ(ierr);

  // TODO: label with levelset id tag:
//  ierr = ArrayWrite( irreg, "irregNode", idx ); CHKERRQ(ierr);

  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/%s.irregNode.%d.array",wd,name,idx+FILE_COUNT_START); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);

  for ( i = 0; i < len; ++i) {
    node = &nodes[i];
    if( iim->is2D ) {
      PetscReal row[] = {
          node->X.x,
          node->X.y,
          node->nx,
          node->ny,
          node->f1,
          node->fa1,
          node->k,
          node->f1_n,
          node->f1_nn
      };
      rowlen = sizeof(row)/sizeof(PetscReal);
      ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
    } else {
      PetscReal row[] = {
          node->X.x,    // 0
          node->X.y,
          node->X.z,
          node->nx,     // 6
          node->ny,
          node->nz,
          node->k,      // 9
          node->f1,
          node->f1_n,
          node->f1_t,   //12
          node->f1_nn,
          node->f1_tt,
          node->f1_nt,  //15
          node->sx,     //16
          node->sy,
          node->sz,
          node->rx,     //19
          node->ry,
          node->rz,     //21
          node->numNei  //22
      };
      rowlen = sizeof(row)/sizeof(PetscReal);
      ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
    }
  }

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
//  ierr = PetscLogEventEnd(EVENT_LevelSetWriteIrregularNodeList,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
