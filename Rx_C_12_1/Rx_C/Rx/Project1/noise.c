#include <stdlib.h>
#include <math.h>
#include "complex.h"

double* MLS_gaussian2(double d)
{
    static double  res[2];

    double U1;
    double U2;
    double r;

    static const double C_INV_RAND_MAX = 1.0 / RAND_MAX;

    do
    {
        U1 = (double)(2 * rand()) * C_INV_RAND_MAX - 1.0;
        U2 = (double)(2 * rand()) * C_INV_RAND_MAX - 1.0;
        r = U1 * U1 + U2 * U2;

    } while (r >= 1.0);

    r = sqrt((-2.0 * log(r)) / r);

    res[0] = U1 * r * d;
    res[1] = U2 * r * d;

    return res;
}

void Complex_Batch_Guass(dcomplex* Out, double noisePower, int dataLen)
{

    const double noiseVolt = sqrt(noisePower * 0.5);
    static const double C_INV_RAND_MAX = 1.0 / RAND_MAX;

    for (int i = 0; i < dataLen; i++)
    {
        const double* pAWGN = MLS_gaussian2(noiseVolt);
        Out[i].x = pAWGN[0];
        Out[i].y = pAWGN[1];
    }
}

void batchAdd(dcomplex* output, const dcomplex* in0, const dcomplex* in1, int size)
{

    for (int i = 0; i < size; ++i)
    {
        output[i] = cAdd(&in0[i], &in1[i]);
        output[i].x = (int)output[i].x;
        output[i].y = (int)output[i].y;
    }
}

void doFading(dcomplex* pRxOut, dcomplex* pTxIn, double SNR, int numOfSamp)
{
    double noisePower = 2 * pow(2, 9) * pow(2, 9) * pow(10, -14 * 0.1);
    double sigPower = 2 * pow(2, 9) * pow(2, 9) * pow(10, (-14 + SNR) * 0.1);

    dcomplex* m_awgn_per_sinr;
    m_awgn_per_sinr = (dcomplex*)malloc(numOfSamp * sizeof(dcomplex));
    for (int i = 0; i < numOfSamp; i++)
    {
        m_awgn_per_sinr[i].x = 0;
        m_awgn_per_sinr[i].y = 0;
    }

    dcomplex* tmpTxIn;
    tmpTxIn = (dcomplex*)malloc(numOfSamp * sizeof(dcomplex));
    for (int i = 0; i < numOfSamp; i++)
    {
        tmpTxIn[i].x = sqrt(sigPower) * pTxIn[i].x;
        tmpTxIn[i].y = sqrt(sigPower) * pTxIn[i].y;
    }

    Complex_Batch_Guass(m_awgn_per_sinr, noisePower, numOfSamp);
    batchAdd(pRxOut, tmpTxIn, m_awgn_per_sinr, numOfSamp);
    for (int i = 0; i < numOfSamp; i++)
    {
        if (pRxOut[i].x > 511)
        {
            pRxOut[i].x = 511;
        }
        if (pRxOut[i].x < -511)
        {
            pRxOut[i].x = -511;
        }

        if (pRxOut[i].y > 511)
        {
            pRxOut[i].y = 511;
        }
        if (pRxOut[i].y < -511)
        {
            pRxOut[i].y = -511;
        }
    }

    free(m_awgn_per_sinr);
    free(tmpTxIn);

}