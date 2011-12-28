#include "DCell.h"

#undef __FUNCT__
#define __FUNCT__ "DCellCreate"
PetscInt EVENT_DCellCreate;
PetscErrorCode DCellCreate(Reaction rxn, DCell *cell)
{
  DCell c;
  char *prefix = "dcell"; // for command line dcell options
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DCellCreate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_DCellCreate,"DCellCreate", 0);
  ierr = PetscNew( struct _DCell, &c); CHKERRQ(ierr);
  c->rxn = rxn;
  
  //TODO: sphere or circle default state?
  iCoor s = {16, 16, 0};
  ierr = GridCreate(s, &c->u); CHKERRQ(ierr);
  ierr = GridCreate(s, &c->v); CHKERRQ(ierr);
  ierr = LevelSetCreate( s, &c->lsPlasmaMembrane ); CHKERRQ(ierr);
//  ierr = ; CHKERRQ(ierr);
  
  *cell = c;
  PetscLogEventEnd(EVENT_DCellCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellDestroy"
PetscInt EVENT_DCellDestroy;
PetscErrorCode DCellDestroy( DCell c )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DCellDestroy,0,0,0,0);
//  PetscLogEventRegister(&EVENT_DCellDestroy,"DCellDestroy", 0);
  ierr = LevelSetDestroy(c->lsPlasmaMembrane); CHKERRQ(ierr);
  ierr = GridDestroy(c->u); CHKERRQ(ierr);
  ierr = GridDestroy(c->v); CHKERRQ(ierr);
  ierr = GridDestroy(c->w); CHKERRQ(ierr);
  ierr = ReactionDestroy(c->rxn); CHKERRQ(ierr);
  ierr = PetscFree(c); CHKERRQ(ierr);
  PetscLogEventEnd(EVENT_DCellDestroy,0,0,0,0);
  PetscFunctionReturn(0);
}

void DCellResizeMemory()
{
  //Make a new DCell at the
  //Find the minimally containing box
  //Recopy only the internal concentrations
  //Destroys old DCell
}

void DCellRecenterMemory()
{
  //Create buffer for only one layer of chemical species
  //Calculate center of mass using orthoprojections
  //Calc displacement between CoM and CoB (center of box)
  //Copy shifted chem_1 to buffer
  //Copy buffer to conc vec
  //Zero out buffer
}

void DCellAssembleMassTransport(DCell cell)
{
  PetscInt s, S = 5;
  MatStencil row, col[S];
  iCoor n = cell->lsPlasmaMembrane->g->n;
  PetscInt i, j;
  PetscReal **phi = cell->lsPlasmaMembrane->g->v2;
  PetscReal **u = cell->u->v2,
            **v = cell->v->v2;
  PetscReal val[5];
  PetscReal xx = 2, yy = 2; // TODO: fix length of discretization (2 h)
  
  for (row.j = 0; row.j < n.y; ++row.j)
  {
    for (row.i = 0; row.i < n.x; ++row.i)
    {
      i = row.i; j = row.j;
      
      if( phi[i][j] < 0 ) //inside cell
      {
        /* second order flux form
         * w'[t] = ( a1 (w0 + w1) - a2 (w1 + w2) ) / (2 h) 
         */
        val[0] =    u[i][j] / xx;
        val[1] = -u[i+1][j] / xx;
        val[2] =    v[i][j] / yy;
        val[3] = -v[i][j+1] / yy;
        val[4] = ( u[i][j] - u[i+1][j] ) / xx + 
                 ( v[i][j] - v[i][j+1] ) / yy;
      } else {
        for( s = 0; s < S; s++ )
          val[s] = 0;
        val[S-1] = 1;
      }
      for( row.c = 0; row.c < cell->rxn->dof; ++row.c)
      {
        for( s = 0; s < S; s++ )
          col[s].c = row.c;
//        MatSetValuesStencil(cell->mat, 1, &row, S, col, val, ADD_VALUES); 
      }
    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "DCellGrow"
PetscInt EVENT_DCellGrow;
PetscErrorCode DCellGrow(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DCellGrow,0,0,0,0);
  
  
   
  
  PetscLogEventEnd(EVENT_DCellGrow,0,0,0,0);
  PetscFunctionReturn(0);
}


DCellGetInterfaceVelocity()
{
  ierr = IIMUpdateSurfaceQuantities(iim,CELL_CENTER,ls[i]); CHKERRQ(ierr);
  ArrayCreate( 4*sizeof(int), 16, &a);
  ArraySetSize( a, 4 * ls->irreg->len ); // indexes 
  ArraySetSize( val, 4 * ls->irreg->len ); //values
  iCoor = idx = ArrayGetData(a);
  int d = 0;
  for( int i = 0; i <  ) //For each node, get its four neighbors for bilinear interpolation
  {
    n = IrregGet(i);
    X = n->x + n->ox;
    Y = n->y + n->oy; Z
    
    //TODO: have shifts linked with X_FACE iCoor constants in IIM and irreg nodes
    int xU[2] = { (int)floor(X+0.5), (int)ceil(X+0.5) }; 
    int yU[2] = { (int)floor(Y), (int)ceil(Y) };
    int xV[2] = { (int)floor(X), (int)ceil(X) }; 
    int yV[2] = { (int)floor(Y+0.5), (int)ceil(Y+0.5) };
    
    for( y = 0; y < 2; ++y)
    {
      for( x = 0; x < 2; ++x)
      {
        idxU[d][0] = yU[y];
        idxU[d][1] = xU[x];
        idxV[d][0] = yV[y];
        idxV[d][1] = xV[x];
        d++;
      }
    }
  }
  
  NGA_Gather(gaU, valU, idxU, len);
  NGA_Gather(gaV, valV, idxV, len);
  
  for( irreg )
  {
    
    vel.x = IIMInterpolation( idxU[d], valU, i);
    vel.y = IIMInterpolation( idxV[d], valV, i);
    d += 4; //four corners of the bilinear interpolation
    n->Vn = vel.x * n->nx + vel.y * n->ny; //dot product of the velocity vector with the surface normal
  }
}

IIMInterpolation()
{
  xl = X - idx[d][0];
  yl = Y - idx[d][1];
  
  vel.x = (1-xl)*(1-yl)*val[i+0] +
            xl  *(1-yl)*val[i+1] +
            xl  *  yl  *val[i+2] +
          (1-xl)*  yl  *val[i+3];
}