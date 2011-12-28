#ifndef MYCHECK_H_
#define MYCHECK_H_

#include "check.h"

#undef CHKERRQ
#define CHKERRQ(n) if (n) { \
  const char *text;\
  char *msg;\
  PetscErrorMessage(n,&text,&msg);\
  fail("\n%s(%d): %s (%s)",__FILE__,__LINE__, text, msg);}

#undef fail_unless
#define fail_unless(expr, ...) {\
  char buffer[256], msg[256]; \
  sprintf(msg, ##__VA_ARGS__); \
  sprintf(buffer, "\n%s(%d): %s",__FILE__, __LINE__,msg); \
        _fail_unless(expr, __FILE__, __LINE__,\
        "Assertion '"#expr"' failed" , buffer, NULL);}

#define LINE(...) { \
  char buffer[256]; \
  sprintf(buffer, ##__VA_ARGS__); \
  PetscPrintf(PETSC_COMM_SELF,"%s(%d):%s\n",__FILE__,__LINE__,buffer);}
  
void setup();
void teardown();

#endif /*MYCHECK_H_*/
