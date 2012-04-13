#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "IIMSurfaceNeighbors_2D"
PetscErrorCode IIMSurfaceNeighbors_2D( IIM iim, LevelSet ls, IrregularNode *n, IrregularNode **nodes, int *l )
{
  int i,j,m;
  int I,J;
  IrregularNode *node;
  IrregularNode *nodelist = ArrayGetData(ls->irregularNodes);
  int e = iim->eps+1;
  int len = 0;
  iCoor p, q;
  PetscBool unique = PETSC_TRUE;
  GridPoint *gridpoint;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  for (j = -e; j <= e; ++j) {
    for (i = -e; i <= e; ++i) {
      I = i+n->pos.x;
      J = j+n->pos.y;

      // if node out of bounds, skip it
      if( I < p.x || q.x <= I ||
          J < p.y || q.y <= J ) continue;

      ierr = IIMGetIrregularNodeGridPoint_2D( iim, I, J, &gridpoint); CHKERRQ(ierr);

      // each node has 4 neighbors
      for (m = 0; m < gridpoint->len; ++m) {
        // get (i,j;m)
        node = &nodelist[m + gridpoint->idx];

        // If nei null, skip it
        if( node == 0 ) continue;

        // if node outside ball, skip it
        if(sqrt(PetscSqr(node->X.x - n->X.x) +
                PetscSqr(node->X.y - n->X.y) ) > iim->eps ) continue;

        // if ortho-proj, skip it
        if( node->axis == -1 ) continue;

        /* Filter irregnodes that correspond to the same point
         * when irregnodes are stencil intersections
         *
         * NOW ALL IRREGNODES ARE UNIQUE
        unique = PETSC_TRUE;
        for( c = 0; c < len; c++) {
          if( node->x + node->axis.x == nodes[c]->x && node->shift.x == nodes[c]->shift.x &&
              node->y + node->axis.y == nodes[c]->y && node->shift.y == nodes[c]->shift.y )
          {
            unique = PETSC_FALSE;
            break;
          }
        } // O(n^2) test for uniqueness
        */

        // if unique, add it to nei nodes list
        if( unique ) {
          nodes[len] = node;
          len++;
        } //if unique

        // if buffer full, stop adding nodes
        if( len >= iim->Np ) {
          printf("IIM Buffer Full: Np = %d\n", iim->Np);
          printf("{%f,%f}\n",n->X.x,n->X.y);
          printf("{");
          for(i = 0; i < len; i++)
          {
            printf("{%f,%f},",nodes[i]->X.x,nodes[i]->X.y);
          }
          printf("}\n");
          exit(1);
          *l = len;
          PetscFunctionReturn(0);
        }
      } //m
    } //i
  } //j

  *l = len;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMSurfaceNeighbors_3D"
PetscErrorCode IIMSurfaceNeighbors_3D( IIM iim, LevelSet ls, IrregularNode *n, IrregularNode **nodes, int *l )
{
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
