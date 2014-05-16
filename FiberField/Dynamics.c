#include "FiberField.h"

PetscErrorCode FiberField_Solve( FiberField f, PetscReal ti, Vec x );
PetscErrorCode FiberField_Assemble( FiberField field );

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

// Take local vert/edge list data structures and create global petsc vectors
#undef __FUNCT__
#define __FUNCT__ "FiberField_Setup"
PetscErrorCode FiberField_Setup( FiberField field )
{
  int local_len_verts = ArrayLength( field->verts ); 
  int local_len_edges = ArrayLength( field->edges );
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // Create global vertex vectors
  // {x1,x2,x3, x1,x2,x3,...}
  ierr = VecCreateMPI(field->comm, 3*local_len_verts, PETSC_DETERMINE, &field->x0); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->xi); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->x1); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->v0); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->vi); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->v1); CHKERRQ(ierr);
  ierr = VecDuplicate(field->x0, &field->Fe_v); CHKERRQ(ierr);
  
  // Create global edge vectors
  ierr = VecCreateMPI(field->comm,   local_len_edges, PETSC_DETERMINE, &field->l0); CHKERRQ(ierr);
  ierr = VecCreateMPI(field->comm, 3*local_len_edges, PETSC_DETERMINE, &field->u); CHKERRQ(ierr);
  ierr = VecDuplicate(field->u, &field->w); CHKERRQ(ierr);

  Mat A, B;
  ierr = MatCreate( field->comm, &A); CHKERRQ(ierr);
  ierr = MatCreate( field->comm, &B); CHKERRQ(ierr);
  ierr = MatSetType( A, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetType( B, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetSizes( A, 3*local_len_verts,3*local_len_edges, PETSC_DETERMINE, PETSC_DETERMINE ); CHKERRQ(ierr);
  ierr = MatSetSizes( B, 3*local_len_edges,3*local_len_verts, PETSC_DETERMINE, PETSC_DETERMINE ); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation( A,       2,PETSC_NULL,       2,PETSC_NULL); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation( B,MAXEDGES,PETSC_NULL,MAXEDGES,PETSC_NULL); CHKERRQ(ierr);
  field->matXtoU  = A;
  field->matFe2Fv = B;


  int *vIDs;
  int *eIDs;
  // vIDs to global petsc ordering
  ierr = AOCreateMapping(PETSC_COMM_WORLD, local_len_verts, vIDs, 0, &field->aoVerts); CHKERRQ(ierr); 
  // eIDs to global petsc ordering
  ierr = AOCreateMapping(PETSC_COMM_WORLD, local_len_verts, eIDs, 0, &field->aoEdges); CHKERRQ(ierr); 

  int v; // vertex indexing
  int e; // edge indexing
  const int vlen = ArrayLength( field->verts );
  const int elen = ArrayLength( field->edges );
  struct _Vertex *verts = ArrayGetData( field->verts );
  struct _Edge *edges = ArrayGetData( field->edges );
  // for each vertex, convert eID into petsc ordering
  for (v = 0; v < vlen; v++) {
    for (e = 0; e < MAXEDGES; e++) {
      verts[v].ePO[e] = verts[v].eID[e];
    }
    ierr = AOApplicationToPetsc(field->aoEdges, MAXEDGES, verts[v].ePO); CHKERRQ(ierr);
  }
   
  // for each edge, convert vID into petsc ordering
  for (e = 0; e < elen; e++) {
    edges[e].vPO[0] = edges[e].vID[0];
    edges[e].vPO[1] = edges[e].vID[1];
    ierr = AOApplicationToPetsc(field->aoVerts, 2, edges[e].vPO); CHKERRQ(ierr);
  }

  ierr = FiberField_Assemble( field ); CHKERRQ(ierr);



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
  const PetscReal dt2m = 0.5 * f->dt / f->mass;
  const PetscReal E = 1; // Young's Modulus
  const PetscReal vf = 0; // Fluid velocity
  const PetscReal drag = f->drag; 
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // u = A . x
  // u_i = x_j - x_i
  ierr = MatMult( f->matXtoU, x, f->u); CHKERRQ(ierr); 
// u: displacement
// l = norm( u )
// strain: e = (l - l0) / l0 = l / l0 - 1
// hook's law: stress = E * e
// normalized: w = u / norm_u

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
    j = i / 3;
    strain = norm_u / l0[j] - 1;
    Fe_e[j] = E * strain; 
  }
  ierr = VecRestoreArray(f->u,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->w,&w);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->l0,&l0);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->Fe_e,&Fe_e);CHKERRQ(ierr);

  // f = B . e
  // f_i = Sum( Fe_i )
  ierr = MatMult( f->matFe2Fv, f->Fe_e, f->Fe_v); CHKERRQ(ierr); 

  // Solve for vi
  // vi = v0 + dt/2 * a( ti, x, vi )
  // vi = v0 + dt/2 * ( Fe(x) + Fd(vi) )
  ierr = VecGetLocalSize(f->vi,&nlocal);CHKERRQ(ierr);
  ierr = VecGetArray(f->vi,&vi);CHKERRQ(ierr);
  ierr = VecGetArray(f->v0,&v0);CHKERRQ(ierr);
  ierr = VecGetArray(f->Fe_v,&Fe_v);CHKERRQ(ierr);
  for (i = 0; i < nlocal; i++) {
    vi[i] = v0[i] + dt2m * Fe_v[i] + dt2m * drag * vf; // vf = vf[i]
    vi[i] = vi[i] / (1 + dt2m * drag);
  }
  ierr = VecRestoreArray(f->vi,&vi);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->v0,&v0);CHKERRQ(ierr);
  ierr = VecRestoreArray(f->Fe_v,&Fe_v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_Assemble"
PetscErrorCode FiberField_Assemble( FiberField field )
{
  int i; // index over components {x,y,z}
  int e; // index over edges
  int v; // index over verts
  int row;
  PetscReal sign;
  PetscInt  col[MAXEDGES];
  PetscReal val[MAXEDGES];
  VertexID vIDp[2]; // petsc index of vertex
  EdgeID   eIDp; // petsc index of edge
  const Mat A = field->matXtoU;
  const Mat B = field->matFe2Fv;
  const int vlen = ArrayLength( field->verts );
  const int elen = ArrayLength( field->edges );
  const struct _Vertex *verts = ArrayGetData( field->verts );
  const struct _Edge *edges = ArrayGetData( field->edges );
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // Assemble: u = A.x
  // for each edge i 
  //   get v0 and v1 vertexIDs
  //   convert vID -> petsc ordering

  for (e = 0; e < elen; e++) {
    eIDp = edges[e].eID;
    ierr = AOApplicationToPetsc(field->aoEdges, 1, &eIDp); CHKERRQ(ierr);
    vIDp[0] = edges[e].vID[0];
    vIDp[1] = edges[e].vID[1];
    sign = vIDp[1] - vIDp[0] > 0 ? 1 : 0; 
    ierr = AOApplicationToPetsc(field->aoVerts, 2, vIDp); CHKERRQ(ierr);
    for (i = 0; i < 3; i++) {
      row    = i + 3*eIDp;
      col[0] = i + 3*vIDp[0]; 
      col[1] = i + 3*vIDp[1];
      val[0] = -1.0 * sign;
      val[1] =  1.0 * sign; 
      /*u[e] = x[a] - x[b];*/
      ierr = MatSetValues(A,1,&row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // global vIDs -> global petsc -> local + ghosted 
  
  int low;   // start index in petsc ordering
  int high;
  ierr = VecGetOwnershipRange( x, &low, &high ); CHKERRQ(ierr);

  // f = B.e 
  // for each local vertex v
  //   vID -> pID
  //   for each edge e
  //     eID -> pID
  //     f_i = (+/-)w_i
  const Edge *edge;
  for (v = 0; v < vlen; v++) {
    // petscIndex = v + low
    row = verts[v].petscIndex;
    for (i = 0; i < 3; i++) {
      for (e = 0; e < MAXEDGES; e++) {
        edge = &verts[v].edges[e];
        eIDp = edges.petscIndex;
        if( edge == FIBERFIELD_NO_EDGE ) 
          col[ ] = -1;
        else
          col[ ] = eIDp
        vID = edge->vID;
        sign = vID[1] - vID[0] > 0 ? 1 : 0; 
        col[i] = ePetsc + k;
        val[i] = sign * w;
      }
      row = vIDp + i;
      /* Fe_v[vIDp] +=  Fe_e[eIDp] * w[eIDp] */
      ierr   = MatSetValues(B,1,&row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#define __BUGGY__
#ifndef __BUGGY__ 

#endif
