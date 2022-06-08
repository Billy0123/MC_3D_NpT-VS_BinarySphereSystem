#ifndef MTGENERATOR_H_INCLUDED
#define MTGENERATOR_H_INCLUDED

void InitMT (unsigned int);

int InitRandomMT ();

unsigned int MTGenerate (double*);

double MTRandom0to1 (double*);

#endif
