#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"

//double m_filterCoef[65] = {
//   -0.00021993, -0.00096935, -0.00075026, 0.00033277, 0.00125512, 0.00098512, -0.00045171, -0.00180972,
//   -0.00167292, 0.00013745, 0.00215351, 0.00243849, 0.00045473, -0.00213200, -0.00273638, -0.00035831,
//   0.00297467, 0.00356461, -0.00045181, -0.00642028, -0.00827894, -0.00157099, 0.01125395, 0.02022974,
//   0.01414922, -0.00952653, -0.03840097, -0.04872477, -0.01875039, 0.05567169, 0.15383624, 0.23759494,
//   0.27049297,
//   0.23759494, 0.15383624, 0.05567169, -0.01875039, -0.04872477, -0.03840097, -0.00952653, 0.01414922,
//   0.02022974, 0.01125395, -0.00157099, -0.00827894, -0.00642028, -0.00045181, 0.00356461, 0.00297467,
//   -0.00035831, -0.00273638, -0.00213200, 0.00045473, 0.00243849, 0.00215351, 0.00013745, -0.00167292,
//   -0.00180972, -0.00045171, 0.00098512, 0.00125512, 0.00033277, -0.00075026, -0.00096935, -0.00021993 };

int m_filterCoef[65] = {
    -2, -8, -6, 3, 10, 8, -4, -15,
    -14, 1, 18, 20, 4, -17, -22, -3,
    24, 29, -4, -53, -68, -13, 92, 166,
    116, -78, -315, -399, -154, 456, 1260, 1946,
    2216,
    1946, 1260, 456, -154, -399, -315, -78, 116,
    166, 92, -13, -68, -53, -4, 29, 24,
    -3, -22, -17, 4, 20, 18, 1, -14,
    -15, -4, 8, 10, 3, -6, -8, -2 };

int m_isFixed = 0;      // default is float
int m_fixFac = 0;
static const int C_FILTER_SAMP = 8;
static const int C_INTERPOLA_OFF = 4;


void decFilter(dcomplex* pIn, dcomplex* pOut, int inLen)
{
    const double C_OutFactor = sqrt(0.5);  // 无噪音, 纯信号时保证 pIn/pOut 功率一致

    dcomplex* tmpSampIn;
    tmpSampIn = (dcomplex*)malloc((inLen + 64) * sizeof(dcomplex));

    for (int i = 0; i < inLen + 64; i++)
    {
        tmpSampIn[i].x = 0;
        tmpSampIn[i].y = 0;
    }
    memcpy(&(tmpSampIn[32]), pIn, inLen * sizeof(dcomplex));

    for (int j = 0; j < inLen / C_FILTER_SAMP; j++)
    {
        pOut[j].x = 0;
        pOut[j].y = 0;

        for (int slid = 0; slid < 65; slid++)
        {
            double x = tmpSampIn[j * C_FILTER_SAMP + slid + C_INTERPOLA_OFF].x * m_filterCoef[slid];
            double y = tmpSampIn[j * C_FILTER_SAMP + slid + C_INTERPOLA_OFF].y * m_filterCoef[slid];
            //右移5bit 17->12
            x = (int)x >> 4;
            y = (int)y >> 4;
            //if (fabs(y) > 2047)
            //    printf("%f\n", y);

            pOut[j].x = pOut[j].x + x;
            pOut[j].y = pOut[j].y + y;           
        }  
        pOut[j].x = (int)pOut[j].x >> 2;
        pOut[j].y = (int)pOut[j].y >> 2;

        pOut[j].x = pOut[j].x * C_OutFactor;
        pOut[j].x = pOut[j].y * C_OutFactor;

        pOut[j].x = (int)pOut[j].x >> 4;
        pOut[j].y = (int)pOut[j].y >> 4;
        //if(fabs(pOut[j].x) > 127 | fabs(pOut[j].y) > 127)
        //printf("%f  %f\n", pOut[j].x, pOut[j].y);

        if (1 == m_isFixed)
        {
            pOut[j].x /= (1 << m_fixFac);
            pOut[j].y /= (1 << m_fixFac);
            pOut[j].x = (double)((int)pOut[j].x);
            pOut[j].y = (double)((int)pOut[j].y);
        }
    }
    free(tmpSampIn);
}