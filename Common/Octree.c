
typedef struct _Octree {
  Octree parent;
  Octree leafs[8];
  int level;
  iCoor a,b; [a.x---b.x]
  void *ptr;
} *Octree;

typedef struct _OctreeRoot {
  Octree root;
  int maxDepth;
} *OctreeRoot;

#undef __FUNCT__
#define __FUNCT__ "OctreeCreate"
PetscErrorCode OctreeCreate( int maxDepth, iCoor a, iCoor b, Octree *root )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}
