#ifndef DWORLD_H_
#define DWORLD_H_

typedef struct _DWorld *DWorld;

void DWorldCreate(double gxs, double gxe, double gys, double gye, DWorld *world);
void DWorldDestroy( DWorld world );

void DWorldBalance(DWorld world);
void DWorldTimeStep(DWorld world, double dt );

#endif /* DWORLD_H_ */
