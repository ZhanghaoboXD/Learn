#include <stdio.h>
#include <math.h>
#include "complex.h"

void getPseudoRandomSeq(int* pOut_NRZ, int Cinit, int len)
{
	unsigned int X1 = 1;
	unsigned int X2 = Cinit;

	int tmpBit = 0;

	//循环移位 1600次
	for (int i = 0; i < 1600; ++i)
	{
		X1 |= ((((X1 >> 3) + (X1 >> 0)) & 0x01) << 31);    // 0x09
		//printf("0x%x ", X1);
		X1 >>= 1;

		X2 |= ((((X2 >> 3) + (X2 >> 2) + (X2 >> 1) + (X2 >> 0)) & 0x01) << 31);         // 0x0F
		X2 >>= 1;
	}
	//printf("0x%x, 0x%x\n", X1, X2);
	//到了后期，开始输出
	for (int i = 0; i < len; ++i)
	{
		tmpBit = ((X1 & 0x01) ^ (X2 & 0x01)) & 0x01;

		X1 |= ((((X1 >> 3) + (X1 >> 0)) & 0x01) << 31);    // 0x09
		X1 >>= 1;

		X2 |= ((((X2 >> 3) + (X2 >> 2) + (X2 >> 1) + (X2 >> 0)) & 0x01) << 31);         // 0x0F
		X2 >>= 1;

		pOut_NRZ[i] = 1 - 2 * tmpBit;
		//printf("%d\n", pOut_NRZ[i]);
	}
}

void nrzDeSfProcess_Msk(double* outNrzData, double* inSfedData, int outLen, int* sfNrzSeq, int sfFactor)
{
	int inSfed = 0;
	for (int bit = 0; bit < outLen; bit++)
	{
		outNrzData[bit] = 0;
		for (int sf = 0; sf < sfFactor; sf++)
		{
			outNrzData[bit] += inSfedData[inSfed] * sfNrzSeq[inSfed] * (-1);
			inSfed++;
		}
		outNrzData[bit] /= sfFactor;
		//取整
		outNrzData[bit] = (int)(outNrzData[bit]);
		//printf("%f\n", outNrzData[bit]);
	}
}

void nrzDeSfProcess_QPSK(dcomplex* outNrzData, dcomplex* inSfedData, int outLen, int* sfNrzSeq, int sfFactor)
{
	int inSfed = 0;
	for (int bit = 0; bit < outLen; bit++)
	{
		outNrzData[bit].x = 0;
		outNrzData[bit].y = 0;
		for (int sf = 0; sf < sfFactor; sf++)
		{
			outNrzData[bit].x += inSfedData[inSfed].x * sfNrzSeq[inSfed] * (-1);
			outNrzData[bit].y += inSfedData[inSfed].y * sfNrzSeq[inSfed] * (-1);
			inSfed++;
		}
		outNrzData[bit].x /= sfFactor;
		outNrzData[bit].y /= sfFactor;
		//取整,右移3bit 8->5
		outNrzData[bit].x = (int)outNrzData[bit].x ;
		outNrzData[bit].y = (int)outNrzData[bit].y ;
		//if(fabs(outNrzData[bit].x) > 31 | fabs(outNrzData[bit].y) > 31)
		//printf("%f %f\n", outNrzData[bit].x, outNrzData[bit].y);

	}
}