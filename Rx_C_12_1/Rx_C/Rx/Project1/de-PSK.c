#include <stdio.h>
#include <string.h>
#include <math.h>
#include "complex.h"

void decBPSK(const dcomplex* pIn, double* pOutNRZ, int outBitLen)
{
    for (int i = 0; i < outBitLen; ++i)
    {
        pOutNRZ[i] = (pIn[i].x + pIn[i].y) * 0.5;
    }
}

void decQPSK(dcomplex* pIn, double* pOutNRZ, int outBitLen)
{
    memcpy(pOutNRZ, pIn, sizeof(dcomplex) * outBitLen / 2);
}

void dec8PSK(const dcomplex* pIn, double* pOutNRZ, int outBitLen)
{
    const int inSampLen = outBitLen / 3;

    for (int i = 0; i < inSampLen; ++i)
    {
        double* pRef = &pOutNRZ[3 * i];

        pRef[0] = pIn[i].y;// * SQRT2;
        pRef[1] = pIn[i].x;// * SQRT2;
        pRef[2] = fabs(pIn[i].x) - fabs(pIn[i].y);
    }
}