#include "ImmersedInterfaceMethod.h"
#include "gts.h"

typedef struct _KDTree *KDTree;
typedef struct _GtsPointWrapper GtsPointWrapper;
PetscErrorCode KDTreeCreate( KDTree *kdt );
PetscErrorCode KDTreeDestroy( KDTree kdt );
PetscErrorCode KDTreeBuild( KDTree kdt, Array g );
PetscErrorCode KDTreeRange(KDTree kdt, IrregularNode *n,
    IrregularNode ***nodes, int *len );
void KDTreeSetEps( KDTree kdt, PetscReal eps );
//IrregularNode* KDTreeGetIrregularNode( GSList *p );


struct _KDTree
{
  GNode *kdtree;
  GPtrArray *array;
  GtsBBox *bbox;
  PetscReal eps;
  IrregularNode **nodes;
  int Np;  // max number of nodes
  int len; // actual number of nodes < Np 
};

struct _GtsPointWrapper {
  GtsObject object;
  gdouble x, y, z;
  IrregularNode *node;
};

#undef __FUNCT__
#define __FUNCT__ "KDTreeCreate"
PetscInt EVENT_KDTreeCreate;  
PetscErrorCode KDTreeCreate( KDTree *kdt )
{
  PetscErrorCode ierr;
  
  KDTree t;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_KDTreeCreate,0,0,0,0);
  
  ierr = PetscNew(struct _KDTree, &t); CHKERRQ(ierr);
  
  t->eps = PetscSqrtScalar(3) + .1; //In 3D, this will capture the smallest 3x3 planar grid
  ierr = PetscOptionsGetReal(0,"-kdtree_eps",&t->eps,0); CHKERRQ(ierr);

  t->Np = 32;
  ierr = PetscOptionsGetInt(0,"-kdtree_Np",&t->Np,0); CHKERRQ(ierr);
  ierr = PetscMalloc(t->Np*sizeof(IrregularNode*), &t->nodes); CHKERRQ(ierr);
   
  t->bbox = gts_bbox_new( gts_bbox_class(), NULL, 
                            0, 0, 0, 0, 0, 0 );
  *kdt = t;
  
  PetscLogEventEnd(EVENT_KDTreeCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KDTreeDestroy"
PetscInt EVENT_KDTreeDestroy;
PetscErrorCode KDTreeDestroy( KDTree kdt )
{
  PetscErrorCode ierr;
//  GtsPoint *p;
  GtsPointWrapper *p;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_KDTreeDestroy,0,0,0,0);

  if( kdt )
  {
    ierr = PetscFree(kdt->nodes); CHKERRQ(ierr);
    gts_object_destroy((GtsObject*)kdt->bbox);
    gts_kdtree_destroy( kdt->kdtree );
    
    int i;
    for( i = 0; i < kdt->array->len; i++ )
    {
      p = g_ptr_array_index(kdt->array, i);
      ierr = PetscFree(p); CHKERRQ(ierr);
    }
    
    g_ptr_array_free(kdt->array, TRUE);
    
    ierr = PetscFree(kdt); CHKERRQ(ierr);
  }
  PetscLogEventEnd(EVENT_KDTreeDestroy,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KDTreeBuild"
PetscInt EVENT_KDTreeBuild;
PetscErrorCode KDTreeBuild(KDTree kdt, GArray *g)
{
  GtsPointWrapper *p;
  IrregularNode *n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_KDTreeBuild,0,0,0,0);
//  PetscLogEventRegister(&EVENT_KDTreeBuild,"KDTreeBuild", 0);
  if( g->len == 0 )
      SETERRQ(PETSC_ERR_ARG_SIZ,"Irregular nodes list not updated! (len = 0)");
 
  //First destroy the old tree
  if( kdt->kdtree )
  {
    gts_kdtree_destroy( kdt->kdtree );
    int i;
    for( i = 0; i < kdt->array->len; i++ )
    {
      p = g_ptr_array_index(kdt->array, i);
      ierr = PetscFree(p); CHKERRQ(ierr);
    }
    g_ptr_array_free(kdt->array, TRUE);
  }
   
  //Then rebuild the new tree
  kdt->array = g_ptr_array_new ();
  int j;
  for( j = 0; j < g->len; j++ )
  {
    n = &g_array_index( g, IrregularNode, j);
    //TODO: too many calls to malloc here, need to optimize this
    ierr = PetscNew(GtsPointWrapper, &p); CHKERRQ(ierr); 
    p->x = n->x + n->ox;
    p->y = n->y + n->oy;
    p->z = n->z + n->oz;
    p->node = n;
    g_ptr_array_add( kdt->array, p);
  }
  kdt->kdtree = gts_kdtree_new(kdt->array, NULL );
  
  PetscLogEventEnd(EVENT_KDTreeBuild,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KDTreeRange"
PetscInt EVENT_KDTreeRange;
PetscErrorCode KDTreeRange(KDTree kdt, IrregularNode *n, 
    IrregularNode ***nodes, int *len )
{
  PetscTruth unique;
  PetscReal eps = kdt->eps, 
    X = n->x+n->ox, Y = n->y+n->oy, Z = n->z+n->oz;
  int count = 0;
  IrregularNode *node;
  GSList *slist, *iter;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_KDTreeRange,0,0,0,0);
  
  gts_bbox_set(kdt->bbox, NULL,
        X-eps, Y-eps, Z-eps, 
        X+eps, Y+eps, Z+eps );
  
  slist = gts_kdtree_range(kdt->kdtree, kdt->bbox, NULL);
  
  iter = slist;
  while( iter )
  {
    node = KDTreeGetIrregularNode(iter);
    //Spherically filter points from box
    if( sqrt(PetscSqr(node->x+node->ox - X) +
             PetscSqr(node->y+node->oy - Y) + 
             PetscSqr(node->z+node->oz - Z)) < eps && 
        count < kdt->Np )
    {
      unique = PETSC_TRUE;
      /* Filter nodes that correspond to the same point when 
       * irreg nodes are stencil intersections (not orthogonal projections)
       */
      int i;
      for( i = 0; i < count; i++) 
      {
        if( node->x + node->axis.x == kdt->nodes[i]->x &&
            node->y + node->axis.y == kdt->nodes[i]->y &&
            node->z + node->axis.z == kdt->nodes[i]->z )
        {
          unique = PETSC_FALSE;
          break;
        }
      }
      if( unique )
      {
        kdt->nodes[count] = node;
        count++;
      }
    }
    iter = g_slist_next(iter);
  }
  
  g_slist_free(slist);
  
  *nodes = kdt->nodes;
  *len = count;
  
  PetscLogEventEnd(EVENT_KDTreeRange,0,0,0,0);
  PetscFunctionReturn(0);
}
//TODO: eps in KDTree is reset since tree is destroyed and created after every iteration!!!
void KDTreeSetEps( KDTree kdt, PetscReal eps )
{
  kdt->eps = eps;
}

IrregularNode* KDTreeGetIrregularNode( GSList *p ) 
{
  return ((GtsPointWrapper*)(p->data))->node;
}
