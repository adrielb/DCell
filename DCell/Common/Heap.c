#include "Common.h"

// Data Structures and Algorithms in Java, Goodrich.

struct _Heap {
  Comparator difference;
  Array tree; // Binary tree using array implementation
};

#undef __FUNCT__
#define __FUNCT__ "HeapCreate"
PetscErrorCode HeapCreate(Comparator cmp, Heap *heap)
{
  Heap h;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _Heap, &h); CHKERRQ(ierr);
  h->difference = cmp;
  //TODO: tune initial array size
  ierr = ArrayCreate("heap",sizeof(void*),1e5,&h->tree); CHKERRQ(ierr);
  // First index for array will start at one for an empty tree
  ierr = ArraySetSize(h->tree,1); CHKERRQ(ierr);
  *heap = h;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapDestroy"
PetscErrorCode HeapDestroy(Heap h)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayDestroy(h->tree); CHKERRQ(ierr);
  ierr = PetscFree(h); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapClear"
PetscErrorCode HeapClear( Heap h, MemCache mc )
{
  int i;
  int len = ArrayLength(h->tree);
  void **elem = ArrayGetData(h->tree);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 1; i < len; ++i) {
    ierr = MemCacheFree(mc,elem[i]); CHKERRQ(ierr);
  }
  ierr = ArraySetSize(h->tree,1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapInsert"
static PetscLogEvent EVENT_HeapInsert;
PetscErrorCode HeapInsert( Heap h, void *item )
{
  void* parent;
  int p; // parent idx
  int c; // child idx
  void** tree;
  PetscErrorCode ierr;

  PetscFunctionBegin;
//  ierr = PetscLogEventBegin(EVENT_HeapInsert,0,0,0,0); CHKERRQ(ierr);
  //Add item to end of list (position z)
  c = ArrayLength(h->tree);
  ierr = ArraySetSize(h->tree,c+1); CHKERRQ(ierr);
  tree = ArrayGetData(h->tree);
  tree[c] = item;

  //Bubble item up to correct position in heap
  while( c != 1 ) { // while not at root node
    p = c / 2;
    parent = tree[p];
    // compare parent with new item
    // if heap property satisfied, break;
    if( h->difference( parent, item ) < 0 ) break;
    // if not, swap parent with new item
    tree[c] = parent;
    tree[p] = item;
    c = p;
  }
//  ierr = PetscLogEventEnd(EVENT_HeapInsert,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "HeapPop"
static PetscLogEvent EVENT_HeapPop;
PetscErrorCode HeapPop( Heap h, void *topItem )
{
  int p, l, r, s;
  void **tree;
  void *bubble;
  int len = ArrayLength(h->tree);
  PetscErrorCode ierr;

  PetscFunctionBegin;
//  ierr = PetscLogEventBegin(EVENT_HeapPop,0,0,0,0); CHKERRQ(ierr);
  // (1-base indexing)
  ierr = ArrayGetP(h->tree,1,topItem); CHKERRQ(ierr);
  // if tree only contains one item
  if( len == 2 ) {
    // pop it off and return
    ierr = ArraySetSize(h->tree, 1); CHKERRQ(ierr);
    return 0;
  } // else restore the heap property
  tree = ArrayGetData(h->tree);
  // move last item to top of heap
  p = 1;
  len--;
  bubble = tree[p] = tree[len];
  ierr = ArraySetSize(h->tree, len); CHKERRQ(ierr);
  l = 2*p;  // left child
  while( l < len )
  {
    r = 2*p+1; // right child
    s = l;
    // if the right child exists and is smaller than the left child, use the right child
    if( r < len && h->difference(tree[r],tree[l]) < 0 ) s = r;
    // if the smallest child is greater than the parent,
    // heap property satisfied, return;
    if( h->difference(tree[p],tree[s]) <= 0 ) break;
    // heap property not satisfied,
    // swap parent with smallest child
    tree[p] = tree[s];
    tree[s] = bubble;
    p = s;
    l = 2*p;
  }
//  ierr = PetscLogEventEnd(EVENT_HeapPop,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapPeek"
PetscErrorCode HeapPeek( Heap h, void *topItem )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayGetP(h->tree,1,topItem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapCheck"
PetscErrorCode HeapCheck( Heap h )
{
  int p;
  int l, r;
  int len = ArrayLength(h->tree);
  void **tree = ArrayGetData(h->tree);

  PetscFunctionBegin;
  for ( p = 1; p < len; ++p) {
    l = 2*p;
    r = 2*p+1;
    if( l < len && h->difference( tree[p], tree[l] ) > 0 ) {
      printf("Parent [%d] > Left child [%d]", p, l);
    }
    if( r < len && h->difference( tree[p], tree[r] ) > 0 ) {
      printf("Parent [%d] > Right child [%d]", p, l);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeapPrint"
PetscErrorCode HeapPrint( Heap h, HeapPrintNode printNode )
{
  int i;
  void **tree = ArrayGetData(h->tree);
  int len = ArrayLength(h->tree);
//  int d, depth = ceil(log2(len+1))+1;

  PetscFunctionBegin;
  for ( i = 1; i < len; ++i) {
/*    d = ceil(log2(i+1));
    for (int j = 0; j < depth-d; ++j) {
      printf("   ");
    }
//    printf("%d:", i);
*/
    printNode(tree[i]);
    if ( (i & (i+1)) == 0 ) printf("\n");
  }
  printf("\n");
  PetscFunctionReturn(0);
}

PetscBool HeapIsEmpty( Heap h )
{
  return ArrayLength(h->tree) == 1;
}

int HeapSize( Heap h )
{
  return ArrayLength(h->tree) - 1;
}

#undef __FUNCT__
#define __FUNCT__ "HeapRegisterEvents"
PetscErrorCode  HeapRegisterEvents()
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;

  ierr = PetscLogEventRegister("HeapInsert", 0, &EVENT_HeapInsert); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("HeapPop", 0, &EVENT_HeapPop); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
