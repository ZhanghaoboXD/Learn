#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#include <math.h>

const double PI_2 = 6.283185307179586476925286766559;

dcomplex m_phase[32] = { {8191,0}, {8033,1597}, {7567,3134}, {6810,4550}, {5791,5791}, {4550,6810}, {3134,7567}, {1597,8033},
	{0,8191}, {-1597,8033}, {-3134,7567}, {-4550,6810}, {-5791,5791}, {-6810,4550}, {-7567,3134}, {-8033,1597},
	{-8191,0}, {-8033,-1597}, {-7567,-3134}, {-6810,-4550}, {-5791,-5791}, {-4550,-6810}, {-3134,-7567}, {-1597,-8033},
	{0,-8191}, {1597,-8033}, {3134,-7567}, {4550,-6810}, {5791,-5791}, {6810,-4550}, {7567,-3134}, {8033,-1597} };
int m_sysBps;
int m_sampRate;
int m_numOfCycleSamp;


//dcomplex POLAR_INTEFER(double sita, int i)
//{
//	dcomplex res;
//
//#ifdef USE_F87_ASSEMBLY
//	__asm fld  sita;  // sita - > ST0
//	__asm fild i;     // sita - > ST1  i -> ST0
//	__asm fmul ST(0), ST(1)
//	__asm fsincos;
//	__asm fstp res.x;
//	__asm fstp res.y;
//	__asm ffree ST(0)
//#else
//	const double phi = sita * i;
//	res.x = cos(phi);
//	res.y = sin(phi);
//#endif
//
//	return res;
//}

void initPhase(int sysBps, int divation, int sampRate)
{
	m_sysBps = sysBps;
	m_sampRate = sampRate;
	m_numOfCycleSamp = sampRate * divation;

	//m_phase = (dcomplex*)malloc(sizeof(dcomplex) * (m_numOfCycleSamp + 1));

	//double phaseData = PI_2 / (m_numOfCycleSamp);

	//for (int i = 0; i < m_numOfCycleSamp + 1; i++)
	//{
	//	m_phase[i] = POLAR_INTEFER(phaseData, i);

	//	//量化
	//	m_phase[i].x = (int)(m_phase[i].x * (pow(2, 5) - 1));
	//	m_phase[i].y = (int)(m_phase[i].y * (pow(2, 5) - 1));
	//	//printf("%d  %f  %f\n", i, m_phase[i].x, m_phase[i].y);
	//}

}

void deMskFsk(dcomplex* pIn, double* pOutNRZ, int outBitLen, int type)
{
	const int modFactor = (type == 0) ? 1 : 2;

	dcomplex tempSamp[2];
	tempSamp[0] = set(0.0, 0.0);
	tempSamp[1] = set(0.0, 0.0);
	int glbSamp = 0;
	for (int j = 0; j < outBitLen; j++)
	{
		tempSamp[0] = set(0.0, 0.0);  //bit 0
		tempSamp[1] = set(0.0, 0.0);  //bit 1

		for (int samp = 0; samp < m_sampRate; samp++)
		{
			dcomplex temp2, temp3, temp4;
			temp2 = cMult(&pIn[j * m_sampRate + samp], &m_phase[glbSamp]);
			//右移7bit 19->12
			temp2.x = (int)temp2.x >> 7;
			temp2.y = (int)temp2.y >> 7;
			//if (fabs(temp2.x) > 2047 | fabs(temp2.y) > 2047)
			//	printf("%f  %f\n", temp2.x, temp2.y);
			tempSamp[0] = cAdd(&tempSamp[0], &temp2);  //bit 0
			temp3 = getConj(&m_phase[glbSamp]);
			temp4 = cMult(&pIn[j * m_sampRate + samp], &temp3);
			//右移7bit 19->12
			temp4.x = (int)temp4.x >> 7;
			temp4.y = (int)temp4.y >> 7;
			//if (fabs(temp4.x) > 2047 | fabs(temp4.y) > 2047)
			//	printf("%f  %f\n", temp4.x, temp4.y);
			tempSamp[1] = cAdd(&tempSamp[1], &temp4);
			glbSamp = (glbSamp + modFactor) % m_numOfCycleSamp;
		}

		tempSamp[0].x = (int)tempSamp[0].x >> 2;
		tempSamp[0].y = (int)tempSamp[0].y >> 2;
		tempSamp[1].x = (int)tempSamp[1].x >> 2;
		tempSamp[1].y = (int)tempSamp[1].y >> 2;
		//if (fabs(tempSamp[1].x) > 2047 | fabs(tempSamp[1].y) > 2047)
		//	printf("%f  %f\n", tempSamp[1].x, tempSamp[1].y);


		tempSamp[0] = cDiv(&tempSamp[0], m_sampRate);
		tempSamp[1] = cDiv(&tempSamp[1], m_sampRate);
		//取整操作
		tempSamp[0].x = (int)tempSamp[0].x;
		tempSamp[0].y = (int)tempSamp[0].y;
		tempSamp[1].x = (int)tempSamp[1].x;
		tempSamp[1].y = (int)tempSamp[1].y;
		// 0/1 相差软判决
		pOutNRZ[j] = getPower(&tempSamp[0]) - getPower(&tempSamp[1]);
		//移位11bit 16->5
		pOutNRZ[j] = (int)pOutNRZ[j] >> 11;//理论应该11bit，但实际貌似10bit就够
		//if (fabs(pOutNRZ[j]) > 63)
		//	printf("%d %f\n", j, pOutNRZ[j]);
		// 0/1 直接软判决
		// pOutNRZ[j] = (getPower(&tempSamp[0]) > getPower(&tempSamp[1])) ? getPower(&tempSamp[0] : (-1) * getPower(&tempSamp[1];
		// 硬判决
		// pOutNRZ[j] = (getPower(&tempSamp[0] > getPower(&tempSamp[1]) ? 1 : -1;
	}

}
