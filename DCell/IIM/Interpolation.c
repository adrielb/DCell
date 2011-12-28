#include "ImmersedInterfaceMethod.h"

PetscReal IIMInterpolation_2D(PetscReal h, Jump *j, Grid phi, const PetscReal *shift, iCoor p, Coor s, int **idx, PetscReal *val);
PetscReal IIMInterpolation3D(Coor s, int *idx, PetscReal *val, int d);

/* Interpolates the velocity component of an irregular node at its orthogonal projection
 * and takes the 
 */
#undef __FUNCT__
#define __FUNCT__ "IIMInterfaceVelocity"
PetscInt EVENT_IIMInterfaceVelocity;
PetscErrorCode IIMInterfaceVelocity( IIM iim, DA da, Vec vecVel, int *ga, LevelSet ls )
{
  iCoor p = ls->g->p;
  PetscReal X,Y,Z; //orthogonal projection of irregular node
  int i, d;
  int xs,ys,zs; // bottom corner of interpolation box
  int x, y, z; // loop index for coordinate
  int len; // number of irreg nodes
  IrregularNode *n;
  int **idx, *coor;
  PetscReal *val;
  PetscReal vel; //interpolated velocity
  PetscReal *normal; //normal vector at irregular node
  DALocalInfo info;
  iCoor CELL_CENTER = {0,0,0};
  const PetscReal shift[3][3] = {{0.5,  0,  0},
                                 {0  ,0.5,  0},   // TODO: use common.h Tensor[][]
                                 {0  ,  0,0.5}};


  PetscErrorCode ierr;
//  PetscLogEventRegister(&EVENT_IIMInterfaceVelocity,"IIMInterfaceVelocity", 0);
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMInterfaceVelocity,0,0,0,0); CHKERRQ(ierr);
  ierr = IIMUpdateSurfaceQuantities(iim,CELL_CENTER,ls); CHKERRQ(ierr);
  len = ArrayLength(ls->irregularNodes);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  const int size = info.dim * (info.dim == 2 ? 4 : 8 ) * len;
  ierr = ArraySetSize( iim->aIdx, size ); CHKERRQ(ierr);
  ierr = ArraySetSize( iim->aCoor, size ); CHKERRQ(ierr);
  ierr = ArraySetSize( iim->aVal,  size ); CHKERRQ(ierr);
  idx = ArrayGetData( iim->aIdx );
  coor = ArrayGetData( iim->aCoor );
  val = ArrayGetData( iim->aVal );

  //Assemble idx array
  int c = 0; // counter for idx[] arrays, every irreg node has 4 (2D) or 8 (3D) neighbors
  for( i = 0; i < len; i++ ) //For each node, get its neighbors
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
    X = n->x + n->op.x;
    Y = n->y + n->op.y;
    Z = n->z + n->op.z;

    for( d = 0; d < info.dim; d++) // For each dimension, 2 or 3
    {
      //TODO: have shifts linked with X_FACE iCoor constants in IIM and irreg nodes
      xs = (int)floor(X+shift[d][0]);
      ys = (int)floor(Y+shift[d][1]);

      // u                    v
      // [x y  x y  x y  x y] [x y  x y  x y  x y]
      for( y = 0; y < 2; ++y)
      {
        for( x = 0; x < 2; ++x)
        {
          idx[c] = &coor[4*c]; // [x y z d]
          idx[c][0] = ys+y;
          idx[c][1] = xs+x;
          idx[c][2] = d;
          c++;
        }
      }

      // u                                        v     w
      // [xyz  xyz  xyz  xyz  xyz  xyz  xyz  xyz] [...] [...]
    }
  }

//  Submit idx array to obtain values
//  ierr = VecGetValuesGA( da, vec, ga, val, idx, 4*len); CHKERRQ(ierr);

  //Single process
  {
  PetscReal ***vel;
  ierr = DAVecGetArrayDOF(da,vecVel,&vel); CHKERRQ(ierr);
  for (int m = 0; m < c; ++m) {
    val[m] = vel[idx[m][0]][idx[m][1]][idx[m][2]];
  }
  ierr = DAVecRestoreArrayDOF(da,vecVel,&vel); CHKERRQ(ierr);
  }

//  DEBUG: write coor and val to disk
//  ierr = ArrayWrite(iim->aCoor,"coor",0); CHKERRQ(ierr);
//  ierr = ArrayWrite(iim->aVal,"val",0); CHKERRQ(ierr);


  Jump j;
  Coor s, S;
  c = 0; //Reset the counter
  for( i = 0; i < len; i++ ) // For each node, perform the bilinear interpolation
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
    n->Vn = 0;
    s.x = n->x + n->op.x + p.x;
    s.y = n->y + n->op.y + p.y;
    s.z = n->z + n->op.z + p.z;

    for (d = 0; d < info.dim; ++d) {
      S.x = s.x + shift[d][0];
      S.y = s.y + shift[d][1];
      S.z = s.z + shift[d][2];

      JumpVelocity( *iim->mu, n, &j, d);
      IIMLocalToGlobal_1st( n, &j );

      if( ls->g->is2D )
      {
        vel = IIMInterpolation_2D( iim->dh.x, &j, ls->g, shift[d], p, S, &idx[c], &val[c]);
        c += 4; //next four corners of the bilinear interpolation
      } else {
        vel = IIMInterpolation3D( S, idx[c], val, c);
        c += 8; //next eight corners of the bilinear interpolation
      }
      //dot product of the velocity vector with the surface normal
      normal = &n->nx;
      n->Vn += vel * normal[d];

//      if( d == 0 )  n->Vn = vel;
    }
//    n->Vn = -1;
//printf("{%d,%d,%f},",n->x,n->y,n->Vn);
  }
//printf("\n");
  ierr = PetscLogEventEnd(EVENT_IIMInterfaceVelocity,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscReal IIMInterpolation_2D(PetscReal h, Jump *j, Grid phi, const PetscReal *shift, iCoor p, Coor s, int **idx, PetscReal *val)
{
  PetscReal C1,C2,C3,C4;

  s.y -= idx[0][0];
  s.x -= idx[0][1];

  C1 =  h * (1-s.x) * (1-s.y) * (   s.x *j->x +    s.y *j->y);
  C2 = -h *    s.x  * (1-s.y) * ((1-s.x)*j->x -    s.y *j->y);
  C3 = -h *    s.x  *    s.y  * ((1-s.x)*j->x + (1-s.y)*j->y);
  C4 =  h * (1-s.x) *    s.y  * (   s.x *j->x - (1-s.y)*j->y);

  C1 = Bilinear2D( GridFunction2D_Identity, phi, idx[0][1]-shift[0], idx[0][0]-shift[1]) > 0 ? C1 : 0;
  C2 = Bilinear2D( GridFunction2D_Identity, phi, idx[1][1]-shift[0], idx[1][0]-shift[1]) > 0 ? C2 : 0;
  C3 = Bilinear2D( GridFunction2D_Identity, phi, idx[3][1]-shift[0], idx[3][0]-shift[1]) > 0 ? C3 : 0;
  C2 = Bilinear2D( GridFunction2D_Identity, phi, idx[2][1]-shift[0], idx[2][0]-shift[1]) > 0 ? C4 : 0;

//  C1=0;C2=0;C3=0;C4=0;
  return (1-s.x)*(1-s.y)*val[0] + C1 +
           s.x  *(1-s.y)*val[1] + C2 +
           s.x  *  s.y  *val[3] + C3 +
         (1-s.x)*  s.y  *val[2] + C4;
}


/*
#undef __FUNCT__
#define __FUNCT__ "IIMInterfaceVelocity"
PetscInt EVENT_IIMInterfaceVelocity;
PetscErrorCode IIMInterfaceVelocity( IIM iim, DA da, Vec *vec, int *ga, LevelSet ls )
{
  PetscReal X,Y,Z; //orthogonal projection of irregular node
  int x, y, z; // loop index for coordinate
  int len = ArrayLength(ls->irregularNodes); // number of irreg nodes
  IrregularNode *n;
  int **idx, *coor;
  PetscReal *val;
  PetscReal vel; //interpolated velocity
  PetscReal *normal; //normal vector at irregular node
  DALocalInfo info;
  iCoor CELL_CENTER = {0,0,0};
  const PetscReal shift[3][3] = {{0.5,  0,  0},
                                 {0,  0.5,  0},
                                 {0,    0,0.5}};
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMInterfaceVelocity,0,0,0,0);; CHKERRQ(ierr);
//  PetscLogEventRegister(&EVENT_IIMInterfaceVelocity,"IIMInterfaceVelocity", 0);

  ierr = IIMUpdateSurfaceQuantities(iim,CELL_CENTER,ls); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = ArraySetSize( iim->aIdx,  4 * len ); CHKERRQ(ierr);
  ierr = ArraySetSize( iim->aCoor, 4 * len ); CHKERRQ(ierr);
  ierr = ArraySetSize( iim->aVal,  4 * len ); CHKERRQ(ierr);
  
  for( int d = 0; d < info.dim; d++) // For each dimension, 2 or 3
  {
    //Assemble idx array
    int i, c = 0; // counter for idx[] arrays, every irreg node has 4 (2D) or 8 (3D) neighbors
    for( i = 0; i < len; i++ ) //For each node, get its neighbors
    {
      ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
      X = n->x + n->op.x;
      Y = n->y + n->op.y;
      Z = n->z + n->op.z; 
      
      //TODO: have shifts linked with X_FACE iCoor constants in IIM and irreg nodes
      int xs[2] = { (int)floor(X+shift[d][0]), (int)ceil(X+shift[d][0]) };
      int ys[2] = { (int)floor(Y+shift[d][1]), (int)ceil(Y+shift[d][1]) };
      int zs[2] = { (int)floor(Z+shift[d][2]), (int)ceil(Z+shift[d][2]) };
  
      if( ls->g->is2D )
      {
        for( y = 0; y < 2; ++y)
        {
          for( x = 0; x < 2; ++x)
          {
            ierr = ArrayGet(iim->aCoor,c,(void*)&coor); CHKERRQ(ierr);
            ierr = ArrayGet(iim->aIdx,c,(void*)&idx); CHKERRQ(ierr);
            idx = &coor;
            idx[0] = ys[y];
            idx[1] = xs[x];
            c++;
          }
        }  
      } else {
        for( z = 0; z < 2; ++z)
        {
          for( y = 0; y < 2; ++y)
          {
            for( x = 0; x < 2; ++x)
            {
              idx[c] = &coor[c];
              idx[c][0] = zs[z];
              idx[c][1] = ys[y];
              idx[c][2] = xs[x];
              c++;
            }
          }
        }
      }
    }

    //DEBUG: write coor and val to disk
    int fd;
    char *file = "/tmp/coor-%i.%i.Int32",filename[128];
    sprintf(filename, file, 4*len,d);
    ierr = PetscBinaryOpen(filename, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
    ierr = PetscBinaryWrite(fd,coor,4*len,PETSC_INT,PETSC_FALSE); CHKERRQ(ierr);
    ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
    

    //Submit idx array to obtain values
    ierr = VecGetValuesGA( da, vec[d], ga[d], val, idx, 4*len); CHKERRQ(ierr);
    
    Coor s;
    c = 0; //Reset the counter
    for( i = 0; i < len; i++ ) // For each node, perform the bilinear interpolation
    {
      ierr = ArrayGet(ls->irregularNodes,i,&n); CHKERRQ(ierr);
      s.x = n->x + n->op.x;
      s.y = n->y + n->op.y;
      s.z = n->z + n->op.z;
      
      if( ls->g->is2D )
      {
        vel = IIMInterpolation2D( s, idx[c], val, c);
        c += 4; //next four corners of the bilinear interpolation
      } else {
        vel = IIMInterpolation3D( s, idx[c], val, c);
        c += 8; //next eight corners of the bilinear interpolation
      }
      //dot product of the velocity vector with the surface normal
      normal = &n->nx;
      n->Vn += vel * normal[d];
    }
  }
  ierr = PetscLogEventEnd(EVENT_IIMInterfaceVelocity,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
*/
//TODO: add correction terms to IIM iterpolation

PetscReal IIMInterpolation3D(Coor s, int *idx, PetscReal *val, int d)
{
  int i, j, k;
  PetscReal sum = 0;
  
  s.z -= idx[0];
  s.y -= idx[1];
  s.x -= idx[2];
  
  for( k = 0; k < 2; ++k)
  {
    for( j = 0; j < 2; ++j)
    {
      for( i = 0; i < 2; ++i)
      {
        sum += ((1 - s.x) * (1 - i) + i * s.x) *
               ((1 - s.y) * (1 - j) + j * s.y) *
               ((1 - s.z) * (1 - k) + k * s.z) * val[d++];
      }
    }
  }
  return sum;
}
