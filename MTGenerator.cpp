#include "MTGenerator.h"
#include <time.h>

// Generator Mersenne Twister
// Data   : 16.04.2008
// (C)2012 mgr Jerzy Wałaszek
// Period: 2^19937-1
// http://edu.i-lo.tarnow.pl/inf/alg/001_search/0022.php
//---------------------------

typedef unsigned long long ulong;
unsigned int MT[624];
int mti;

void InitMT (unsigned int x0) {
    ulong x;
    mti = 0;

    MT[0] = x0;
    for(int i = 1; i < 623; i++) {
        x = MT[i-1];
        x = (23023 * x) & 0xffffffffull;
        x = (    3 * x) & 0xffffffffull;
        MT[i] = x;
    }
}

int InitRandomMT () {
    int timeStart = time(0);
    InitMT((unsigned int)timeStart);
    return timeStart;
}

unsigned int MTGenerate (double *randomStartStep) {randomStartStep[1]++;
    const unsigned int MA[] = {0,0x9908b0df};
    long int y;
    int i1,i397;

    i1      = mti +   1; if(  i1 > 623) i1 = 0;
    i397    = mti + 397; if(i397 > 623) i397 -= 624;
    y       = (MT[mti] & 0x80000000) | (MT[i1] & 0x7fffffff);
    MT[mti] = MT[i397] ^ (y >> 1) ^ MA[y & 1];
    y       = MT[mti];
    y       ^=  y >> 11;
    y       ^= (y <<  7) & 0x9d2c5680;
    y       ^= (y << 15) & 0xefc60000;
    y       ^=  y >> 18;
    mti      = i1;
    return y;
}

double MTRandom0to1 (double *randomStartStep) {
    return MTGenerate(randomStartStep)/4294967295.0;  //32-bitowe slowo: maxValue=2^32-1=4294967295
}
