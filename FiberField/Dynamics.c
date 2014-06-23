#include "FiberField.h"
#include "FiberField_private.h"

/*#define DEBUG_DYNAMICS*/

/* 1) AO: renumber obj so that contiguous entries are adjacent
 * 2) determine needed neighbors
 * 3) generate local numbering
 * 4) generate local and global vectors
 * 5) Create VecScatter
 */
/* 
 * 1) spatially partition verticies and edges
 *    global to global mapping
 *    adding and removing verticies
 * 2) local ID's and ghosted ID's
 * 5) copy 
 * local # ->  vID  -> global #
 * [0,m-1]    [a,b]    [0, N-1]
 */

#undef __FUNCT__
#define __FUNCT__ "FiberFieldSolve"
PetscErrorCode FiberFieldSolve( FiberField f )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = FiberField_SpatiallyBalance( f ); CHKERRQ(ierr);
  ierr = FiberField_Init( f ); CHKERRQ(ierr);
  ierr = FiberField_ToVec( f ); CHKERRQ(ierr);
  ierr = FiberField_AssembleDisplacementMatrix( f ); CHKERRQ(ierr);
  ierr = FiberField_Step( f ); CHKERRQ(ierr);
  ierr = FiberField_FromVec( f ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Take local vert/edge list data structures and create global petsc vectors
#undef __FUNCT__
#define __FUNCT__ "FiberField_Init"
PetscErrorCode FiberField_Init( FiberField field )
{
  int local_len_verts = ArrayLength( field->verts ); 
  int local_len_edges = ArrayLength( field->edges );
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInfo1(0, "local_len_verts = %d\n", local_len_verts); CHKERRQ(ierr);
  ierr = PetscInfo1(0, "local_len_edges = %d\n", local_len_edges); CHKERRQ(ierr);
  ierr = ArraySetSize( field->vIDs, local_len_verts ); CHKERRQ(ierr);
  ierr = ArraySetSize( field->eIDs, local_len_edges ); CHKERRQ(ierr);

  // Create global vertex vectors
  // {x1,x2,x3, x1,x2,x3,...}
  ierr = VecDestroy(&field->x0); CHKERRQ(ierr);
  ierr = VecDestroy(&field->xi); CHKERRQ(ierr);
  ierr = VecDestroy(&field->x1); CHKERRQ(ierr);
  ierr = VecDestroy(&field->v0); CHKERRQ(ierr);
  ierr = VecDestroy(&field->vi); CHKERRQ(ierr);
  ierr = VecDestroy(&field->v1); CHKERRQ(ierr);
  ierr = VecDestroy(&field->vf); CHKERRQ(ierr);
  ierr = VecDestroy(&field->Fe_v); CHKERRQ(ierr);
  ierr = VecCreateMPI(field->comm, 3*local_len_verts, PETSC_DETERMINE, &field->Fe_v); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->x0); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->xi); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->x1); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->v0); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->vi); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->v1); CHKERRQ(ierr);
  ierr = VecDuplicate(field->Fe_v, &field->vf); CHKERRQ(ierr);
  
  // Create global edge vectors
  ierr = VecDestroy(&field->l0); CHKERRQ(ierr);
  ierr = VecDestroy(&field->u); CHKERRQ(ierr);
  ierr = VecDestroy(&field->w); CHKERRQ(ierr);
  ierr = VecDestroy(&field->Fe_e); CHKERRQ(ierr);
  ierr = VecCreateMPI(field->comm,   local_len_edges, PETSC_DETERMINE, &field->l0); CHKERRQ(ierr);
  ierr = VecCreateMPI(field->comm, 3*local_len_edges, PETSC_DETERMINE, &field->u); CHKERRQ(ierr);
  ierr = VecDuplicate(field->u, &field->w); CHKERRQ(ierr);
  ierr = VecDuplicate(field->l0, &field->Fe_e); CHKERRQ(ierr);

  // Create displacment and force matrix
  ierr = MatDestroy(&field->matXtoU); CHKERRQ(ierr);
  ierr = MatDestroy(&field->matFe2Fv); CHKERRQ(ierr);
  ierr = MatCreate( field->comm, &field->matXtoU); CHKERRQ(ierr);
  ierr = MatCreate( field->comm, &field->matFe2Fv ); CHKERRQ(ierr);
  ierr = MatSetType( field->matXtoU, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetType( field->matFe2Fv, MATMPIAIJ); CHKERRQ(ierr);
  //                                   # of local rows,  # of local cols,   # of global rows,# of global cols
  ierr = MatSetSizes( field->matXtoU,  3*local_len_edges,3*local_len_verts, PETSC_DETERMINE, PETSC_DETERMINE ); CHKERRQ(ierr);
  ierr = MatSetSizes( field->matFe2Fv, 3*local_len_verts,  local_len_edges, PETSC_DETERMINE, PETSC_DETERMINE ); CHKERRQ(ierr);
  //                                          # of nonzeros in diagonal, # of nonzeros in off-diagonal
  //                                                    d_nz,     d_nnz,    o_nz,     o_nnz 
  ierr = MatMPIAIJSetPreallocation( field->matXtoU,        2,PETSC_NULL,       2,PETSC_NULL); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation( field->matFe2Fv,MAXEDGES,PETSC_NULL,MAXEDGES,PETSC_NULL); CHKERRQ(ierr);

  ierr = MatSetOption( field->matFe2Fv, MAT_ROW_ORIENTED, PETSC_FALSE ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_ToVec"
PetscErrorCode FiberField_ToVec( FiberField f )
{
  int i;
  struct _Vertex *verts = ArrayGetData( f->verts );
  struct _Edge *edges   = ArrayGetData( f->edges );
  int *vIDs             = ArrayGetData( f->vIDs );
  int *eIDs             = ArrayGetData( f->eIDs );
  const int vlen = ArrayLength( f->verts );
  const int elen = ArrayLength( f->edges );
  PetscReal *x0;
  PetscReal *v0;
  PetscReal *l0;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = VecGetArray(f->x0,&x0);CHKERRQ(ierr);
  ierr = VecGetArray(f->v0,&v0);CHKERRQ(ierr);
  for (i = 0; i < vlen; i++) {
    vIDs[i] = verts[i].vID;
    x0[3*i+0] = verts[i].X.x;
    x0[3*i+1] = verts[i].X.y;
    x0[3*i+2] = verts[i].X.z;
    v0[3*i+0] = verts[i].V.x;
    v0[3*i+1] = verts[i].V.y;
    v0[3*i+2] = verts[i].V.z;
  }
  ierr = VecRestoreArray(f->x0,&x0);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->v0,&v0);CHKERRQ(ierr);

  ierr = VecGetArray(f->l0,&l0);CHKERRQ(ierr);
  for (i = 0; i < elen; i++) {
    eIDs[i] = edges[i].eID;
    l0[i] = edges[i].l0;
  }
  ierr = VecRestoreArray(f->l0,&l0);CHKERRQ(ierr);

  // vIDs to global petsc ordering
  ierr = AODestroy(&f->aoVerts); CHKERRQ(ierr);
  ierr = AOCreateMapping(PETSC_COMM_WORLD, vlen, vIDs, PETSC_NULL, &f->aoVerts); CHKERRQ(ierr); 
  // eIDs to global petsc ordering
  ierr = AODestroy(&f->aoEdges); CHKERRQ(ierr);
  ierr = AOCreateMapping(PETSC_COMM_WORLD, elen, eIDs, PETSC_NULL, &f->aoEdges); CHKERRQ(ierr); 

  int v; // vertex indexing
  int e; // edge indexing
  // for each vertex, convert eID into petsc ordering
  for (v = 0; v < vlen; v++) {
    for (e = 0; e < MAXEDGES; e++) {
      verts[v].ePO[e] = verts[v].eID[e];
    }
    ierr = AOApplicationToPetsc(f->aoEdges, MAXEDGES, verts[v].ePO); CHKERRQ(ierr);
    verts[v].vPO = verts[v].vID;
    ierr = AOApplicationToPetsc(f->aoVerts, 1, &verts[v].vPO); CHKERRQ(ierr);
  }
   
  // for each edge, convert vID into petsc ordering
  for (e = 0; e < elen; e++) {
    edges[e].vPO[0] = edges[e].vID[0];
    edges[e].vPO[1] = edges[e].vID[1];
    ierr = AOApplicationToPetsc(f->aoVerts, 2, edges[e].vPO); CHKERRQ(ierr);
    edges[e].ePO = edges[e].eID;
    ierr = AOApplicationToPetsc(f->aoEdges, 1, &edges[e].ePO); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_FromVec"
PetscErrorCode FiberField_FromVec( FiberField f )
{
  int i;
  PetscReal *x1;
  PetscReal *v1;
  const int vlen = ArrayLength( f->verts );
  struct _Vertex *verts = ArrayGetData( f->verts );
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = VecGetArray(f->x1,&x1);CHKERRQ(ierr);
  ierr = VecGetArray(f->v1,&v1);CHKERRQ(ierr);
  for (i = 0; i < vlen; i++) {
    verts[i].X.x = x1[3*i+0];
    verts[i].X.y = x1[3*i+1];
    verts[i].X.z = x1[3*i+2];
    verts[i].V.x = v1[3*i+0];
    verts[i].V.y = v1[3*i+1];
    verts[i].V.z = v1[3*i+2];
  }
  ierr = VecRestoreArray(f->x1,&x1);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->v1,&v1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_Step"
PetscErrorCode FiberField_Step( FiberField field )
{ 
  // Unpack FiberField struct
  PetscReal dt = field->dt;
  PetscReal ti = field->t + 0.5*dt; // Intermediate time
  Vec x0 = field->x0;
  Vec x1 = field->x1;
  Vec xi = field->xi;
  Vec v0 = field->v0;
  Vec vi = field->vi;
  Vec v1 = field->v1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // vi = v0 + dt/2 * a( ti, x0, vi ) 
  ierr = FiberField_Solve( field, ti, x0 ); CHKERRQ(ierr);
  
  // x1 = x0 + dt * vi
  ierr = VecWAXPY(x1, dt, vi, x0); CHKERRQ(ierr);

  // body collisions modify x1 and v0
  ierr = FiberField_BodyCollisions( field ); CHKERRQ(ierr);

  // xi = (x0 + x1) / 2.0
  ierr = VecAXPBYPCZ( xi, 0.5, 0.5, 0.0, x0, x1); CHKERRQ(ierr);

  // vi = v0 + dt/2 * a( ti, xi, vi )
  ierr = FiberField_Solve( field, ti, xi ); CHKERRQ(ierr);

  // v1 = 2 * vi - v0
  ierr = VecAXPBYPCZ( v1, 2.0, -1.0, 0.0, vi, v0 ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FiberField_Solve"
PetscErrorCode FiberField_Solve( FiberField f, PetscReal ti, Vec x )
{
  int i; // 3-component edge vector
  int j; // single component edge vector, j = i / 3; 
  int nlocal; // local array size of edge list
  PetscReal norm_u;
  PetscReal strain;
  PetscReal *u; 
  PetscReal *w;
  PetscReal *Fe_e;
  PetscReal *Fe_v;
  PetscReal *l0;
  PetscReal *v0;
  PetscReal *vi;
  PetscReal *vf; // Fluid velocity at vertex
  const PetscReal dt2m = 0.5 * f->dt / f->mass;
  const PetscReal E = 1; // Young's Modulus
  const PetscReal drag = f->fluidDrag; 
  PetscErrorCode ierr;

  PetscFunctionBegin;
#ifdef DEBUG_DYNAMICS
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nx\n"); CHKERRQ(ierr);
  ierr = VecView( x, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
#endif


  // u: displacement
  // u = A . x
  // u_k = x_j - x_i
  ierr = MatMult( f->matXtoU, x, f->u); CHKERRQ(ierr); 

#ifdef DEBUG_DYNAMICS
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nu\n"); CHKERRQ(ierr);
  ierr = VecView( f->u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
#endif

  // norm_u = ( u1^2 + u2^2 + u3^2 ) / sqrt( norm_u )
  // normalized: w = u / norm_u
  // l = norm_u
  // strain: e = (l - l0) / l0 = l / l0 - 1
  // hook's law: stress = E * e

  ierr = VecGetLocalSize(f->u,&nlocal);CHKERRQ(ierr);
  ierr = VecGetArray(f->u,&u);CHKERRQ(ierr);
  ierr = VecGetArray(f->w,&w);CHKERRQ(ierr);
  ierr = VecGetArray(f->l0,&l0);CHKERRQ(ierr);
  ierr = VecGetArray(f->Fe_e,&Fe_e);CHKERRQ(ierr);
  for (i=0; i<nlocal; i+=3) {
    norm_u = u[i+0]*u[i+0] +
             u[i+1]*u[i+1] +
             u[i+2]*u[i+2]; 
    norm_u = PetscSqrtReal(norm_u);
    w[i+0] = u[i+0] / norm_u;
    w[i+1] = u[i+1] / norm_u;
    w[i+2] = u[i+2] / norm_u;
    // if spring has collapsed to a point
    if ( PetscIsInfOrNanReal(w[i+0]) ) {
      // arbitrary offset
      // TODO: use a random direction
      w[i+0] = 0.333;
      w[i+1] = 0.333;
      w[i+2] = 0.333;
    }
    j = i / 3;
    strain = norm_u / l0[j] - 1;
    Fe_e[j] = E * strain; 
  }
  ierr = VecRestoreArray(f->u,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->w,&w);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->l0,&l0);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->Fe_e,&Fe_e);CHKERRQ(ierr);

#ifdef DEBUG_DYNAMICS
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nw\n"); CHKERRQ(ierr);
  ierr = VecView( f->w, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nFe_e\n"); CHKERRQ(ierr);
  ierr = VecView( f->Fe_e, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
#endif

  // f = B . e
  // f_i = Sum( Fe_i )
  ierr = FiberField_AssembleForceMatrix( f ); CHKERRQ(ierr);
  ierr = MatMult( f->matFe2Fv, f->Fe_e, f->Fe_v); CHKERRQ(ierr); 

#ifdef DEBUG_DYNAMICS
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nFe_v\n"); CHKERRQ(ierr);
  ierr = VecView( f->Fe_v, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
#endif

  // Evaluate fluid velocity to compute drag force
  //   returns f->vf
  ierr = FiberFieldEvaluateFluidVelocity( f, ti, x ); CHKERRQ(ierr);

  // Solve for vi
  // vi = v0 + dt/2 * a( ti, x, vi )
  // vi = v0 + dt/2 * ( Fe(x) + Fd(vi) )
  ierr = VecGetLocalSize(f->vi,&nlocal);CHKERRQ(ierr);
  ierr = VecGetArray(f->vi,&vi);CHKERRQ(ierr);
  ierr = VecGetArray(f->v0,&v0);CHKERRQ(ierr);
  ierr = VecGetArray(f->vf,&vf);CHKERRQ(ierr);
  ierr = VecGetArray(f->Fe_v,&Fe_v);CHKERRQ(ierr);
  for (i = 0; i < nlocal; i++) {
    vi[i] = v0[i] + dt2m * Fe_v[i] + dt2m * drag * vf[i];
    vi[i] = vi[i] / (1 + dt2m * drag);
  }
  ierr = VecRestoreArray(f->vi,&vi);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->v0,&v0);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->vf,&vf);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->Fe_v,&Fe_v);CHKERRQ(ierr);

#ifdef DEBUG_DYNAMICS
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nvi\n"); CHKERRQ(ierr);
  ierr = VecView( f->vi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
#endif

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_BodyCollisions"
PetscErrorCode FiberField_BodyCollisions( FiberField f )
{
  int i;
  int len;
  PetscReal *x1;
  PetscReal *v0;
  Coor gmin = f->globalBounds.min;
  Coor gmax = f->globalBounds.max;
  const PetscReal delta = PETSC_SMALL;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // add small delta
  gmin.x += delta;
  gmin.y += delta;
  gmin.z += delta;
  gmax.x -= delta;
  gmax.y -= delta;
  gmax.z -= delta;

  // body collisions modify x1 and v0
  
  ierr = VecGetLocalSize( f->x1, &len); CHKERRQ(ierr);
  ierr = VecGetArray(f->x1,&x1);CHKERRQ(ierr);
  ierr = VecGetArray(f->v0,&v0);CHKERRQ(ierr);
  for (i = 0; i < len; i+=3) {
    // upper boundary 
    x1[i+0] = x1[i+0] > gmax.x ? gmax.x : x1[i+0];
    x1[i+1] = x1[i+1] > gmax.y ? gmax.y : x1[i+1];
    x1[i+2] = x1[i+2] > gmax.z ? gmax.z : x1[i+2];

    // lower boundary
    x1[i+0] = x1[i+0] < gmin.x ? gmin.x : x1[i+0];
    x1[i+1] = x1[i+1] < gmin.y ? gmin.y : x1[i+1];
    x1[i+2] = x1[i+2] < gmin.z ? gmin.z : x1[i+2];

    //TODO: modify velocity if collision with global bbox occurs 
    //      friction? no-slip?
    /*v1[i+0] = */
    /*v1[i+1] = */
    /*v1[i+2] = */
  }

  ierr = VecRestoreArray(f->x1,&x1);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->v0,&v0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FiberField_AssembleDisplacementMatrix"
PetscErrorCode FiberField_AssembleDisplacementMatrix( FiberField field )
{
  int i; // index over components {x,y,z}
  int e; // index over edges
  PetscInt  row[2];
  PetscInt  col[2];
  PetscReal val[2];
  const Mat XtoU = field->matXtoU;
  const int elen = ArrayLength( field->edges );
  const struct _Edge *edges = ArrayGetData( field->edges );
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // Assemble: u = A.x
  // for each edge e 
  //   get v0 and v1 
  //   convert PO -> vector index = i + 3*PO
 
  for (e = 0; e < elen; e++) {
    /*sign = edges[e].vID[1] - edges[e].vID[0] > 0 ? 1 : -1; */
    for (i = 0; i < 3; i++) {
      row[0]    = i + 3*edges[e].ePO;
      col[0] = i + 3*edges[e].vPO[0]; 
      col[1] = i + 3*edges[e].vPO[1];
      val[0] = -1.0;
      val[1] =  1.0; 
      /*u[e] = x[1] - x[0];*/
      ierr = MatSetValues(XtoU,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(XtoU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(XtoU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_AssembleForceMatrix"
PetscErrorCode FiberField_AssembleForceMatrix( FiberField f )
{
  int i; // index over components {x,y,z}
  int e; // index over edges
  PetscReal *w;
  PetscInt  row[2];
  PetscInt  col[2];
  PetscReal val[2];
  const Mat Fe2Fv = f->matFe2Fv;
  const int elen = ArrayLength( f->edges );
  const struct _Edge *edges = ArrayGetData( f->edges );
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // f = B.e 
  // for each edge e
  //   for each dim i
  //     Fv = Sum_i( w_i * Fe_i )
  ierr = VecGetArray(f->w,&w);CHKERRQ(ierr);
  for (e = 0; e < elen; e++) {
    for (i = 0; i < 3; i++) {
      col[0] =       edges[e].ePO;    
      row[0] = i + 3*edges[e].vPO[0];
      row[1] = i + 3*edges[e].vPO[1];
      val[0] =  w[i+3*e];
      val[1] = -w[i+3*e];
      ierr = MatSetValues(Fe2Fv,2,row,1,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = VecRestoreArray(f->w,&w);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Fe2Fv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Fe2Fv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef DEBUG_DYNAMICS

#define __BUGGY__
#ifndef __BUGGY__ 

  // f = B.e 
  // for each local vertex v (rows of B)
  //   for each edge e (cols of B)
  //     f_i = (+/-)w_i
  EdgeID eID = verts[v].eID[e];
  for (v = 0; v < vlen; v++) {
    for (i = 0; i < 3; i++) {
      for (e = 0; e < MAXEDGES; e++) {
        if( verts[v].eID[e] == FIBERFIELD_NO_EDGE ) 
          col[e] = -1;
        else
          col[e] = 3*verts[v].ePO[e] + i;
        vID = edge->vID;
        sign = vID[1] - vID[0] > 0 ? 1 : 0; 
        val[e] = sign * w[i];
      }
      row = 3*verts[v].vPO + i
      /* Fe_v[vIDp] +=  Fe_e[eIDp] * w[eIDp] */
      ierr   = MatSetValues(B,1,&row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

#endif
