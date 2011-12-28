typedef struct _RLEBlock *RLEBlock;
struct _RLEBlock {
  int extents[2];  // [min, max+1) along axis
  Array start;     // 1D pointer to runcode array
  Array runcodes;  // either +, -, or pointer to defined values
  Array runbreaks; // absolute coordinates where each run starts
};

typedef struct _HRLELevelSet *HRLELevelSet;
struct _HRLELevelSet {
  RLEBlock xblock, yblock, zblock;
  Array phi;
};

PetscErrorCode RLEBlockCreate( RLEBlock *block ) {
  int initsize = 32; // initial size of arrays
  RLEBlock rle;
  ierr = PetscNew(sizeof(struct _RLEBlock), &rle); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int),initsize,&rle->start); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int),initsize,&rle->runcodes); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int),initsize,&rle->runbreaks); CHKERRQ(ierr);
}

PetscErrorCode HRLELevelSetCreate( HRLELevelSet *hrlels ) {
  HRLELevelSet ls;
  ierr = PetscNew(struct _HRLELevelSet, &ls); CHKERRQ(ierr);
  ierr = RLEBlockCreate(ls->xblock); CHKERRQ(ierr);
  ierr = RLEBlockCreate(ls->yblock); CHKERRQ(ierr);
}

PetscErrorCode EncodeLevelSet( LevelSet ls, HRLELevelSet hls ) {
  HRLELevelSet hls;
  ierr = HRLELevelSetCreate( &hls ); CHKERRQ(ierr);

  hls->xblock->extents[0] = p.x;
  hls->xblock->extents[1] = q.x;
  hls->yblock->extents[0] = p.y;
  hls->yblock->extents[1] = q.y;



  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {

    }
  }
}


int binarySearch(int[] a, int key) {
  int low = 0;
  int high = a.length - 1;
  while (low <= high) {
    int mid = low + ((high - low) / 2)
    int midVal = a[mid];
    if (midVal < key)
      low = mid + 1
    else if (midVal > key)
      high = mid - 1;
    else
      return mid; // key found
  }
  return -(low + 1);  // key not found.
}
