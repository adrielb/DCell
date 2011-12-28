#ifndef INTELFFT_H_
#define INTELFFT_H_

#include "Utilities.h"
#include "petscda.h"

typedef struct _IntelFFT *IntelFFT;

PetscErrorCode PCIntelFFTSetPC( PC pc, DALocalInfo g, IS is, iCoor s, iCoor e );

#endif /*INTELFFT_H_*/
