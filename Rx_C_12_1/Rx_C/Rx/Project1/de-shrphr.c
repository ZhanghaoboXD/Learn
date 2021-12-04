#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "tp680.h"

int getCrc8Buf(const int* inputBit, int rawByteLen);
int decCrc8_UNB(const int* inputBit, int rawByteLen);

dcomplex g_ss1_zcSeq[128 * 8] = { 0 };
dcomplex g_ss2_zcSeq[128 * 8] = { 0 };
const double PI = 3.1415926535897932384626433832795;
const double PI2 = 6.283185307179586476925286766559;

dcomplex POLAR(double sita)
{
    dcomplex res;

#ifdef USE_F87_ASSEMBLY
    __asm fld sita;
    __asm fsincos;
    __asm fstp res.x;
    __asm fstp res.y;
#else
    res.x = cos(sita);
    res.y = sin(sita);
#endif

    return res;
}

void getZadoffChuSeq_upSampling(dcomplex* seq, int u, int Nzc, int Kos)
{

    const double phaseUnit = -PI * u / (Kos * Kos * Nzc);

    for (int i = 0; i < Kos * Nzc; ++i)
    {
        const double sita = phaseUnit * i * (i + Kos);

        seq[i] = POLAR(sita);
    }
}

void getZadoffChuSeq_upSampling_v082(dcomplex* seq, int u, int Nzc, int Kos)
{

    const double phaseUnit = -PI * u / (Kos * Kos * Nzc);

    for (int i = 0; i < Kos * Nzc; ++i)
    {
        const double sita = phaseUnit * i * (i + Kos - Kos * Nzc);

        seq[i] = POLAR(sita);
    }
}

void batchConj(dcomplex* output, dcomplex* input, int size)
{

    for (int i = 0; i < size; ++i)
    {
        output[i].x = input[i].x;
        output[i].y = -input[i].y;
    }
}

dcomplex batchMuxAndSum(dcomplex* in0, dcomplex* in1, int size)
{
    dcomplex sum;

    sum = set(0.0, 0.0);

    for (int i = 0; i < size; ++i)
    {
        dcomplex temp;
        temp = cMult(&in0[i], &in1[i]);
        temp.x = (int)temp.x >> 7;
        temp.y = (int)temp.y >> 7;
        //if (fabs(temp.x) > 2047 | fabs(temp.y) > 2047)
        //    printf("%f  %f\n", temp.x, temp.y);
        sum = cAdd(&sum, &temp);
    }

    return sum;
}

void findMax(double* maxVal, int* maxOffset, const double* input, int size)
{

    int tmpOffset = 0;
    double  tmpVal = input[0];

    for (int i = 1; i < size; ++i)
    {
        if (tmpVal < input[i])
        {
            tmpVal = input[i];
            tmpOffset = i;
        }
    }
    if (maxVal) *maxVal = tmpVal;
    if (maxOffset) *maxOffset = tmpOffset;
}

void inrzToByte(int* inputNrzBit, int* outByte, int bitLen)
{
    const int BYTELEN = bitLen / 8;
    for (int byteIdx = 0; byteIdx < BYTELEN; byteIdx++)
    {
        outByte[byteIdx] = 0;
        for (int bitIdx = 0; bitIdx < 8; bitIdx++)
        {
            outByte[byteIdx] ^= ((((1 - inputNrzBit[byteIdx * 8 + bitIdx]) >> 1) & 0x01) << (7 - bitIdx));
        }
    }
}



dcomplex ss1Seq[128 * 8] = { 0 };
dcomplex ss2Seq[128 * 8] = { 0 };
int decShrPhrSeq(dcomplex* rxAirSamp, int rxSampLen)
{

    const int u_ss1 = 4;
    const int u_ss2 = 2;

    const int SS_ZC_Samp_LEN = g_tp680Para.shr_N_zc;
    const int K_os_UpSampling = g_tp680Para.sys_sampRate;

    int estSSPos[2] = { 0 };       // for ss1 pos / ss2 pos
    double estPeak[2] = { 0 };    // for ss1 peak / ss2 peak


    getZadoffChuSeq_upSampling(g_ss1_zcSeq, u_ss1, SS_ZC_Samp_LEN, K_os_UpSampling);
    getZadoffChuSeq_upSampling(g_ss2_zcSeq, u_ss2, SS_ZC_Samp_LEN, K_os_UpSampling);
    //getZadoffChuSeq_upSampling_v082(g_ss1_zcSeq, u_ss1, SS_ZC_Samp_LEN, K_os_UpSampling);
    //getZadoffChuSeq_upSampling_v082(g_ss2_zcSeq, u_ss2, SS_ZC_Samp_LEN, K_os_UpSampling);

    batchConj(ss1Seq, g_ss1_zcSeq, SS_ZC_Samp_LEN * K_os_UpSampling);     //共轭
    batchConj(ss2Seq, g_ss2_zcSeq, SS_ZC_Samp_LEN * K_os_UpSampling);

    //对ss1Seq和ss2Seq同输入量化
    for (int i = 0; i < SS_ZC_Samp_LEN * K_os_UpSampling; i++)
    {
        ss1Seq[i].x = (int)(ss1Seq[i].x * (pow(2, 9) - 1));
        ss1Seq[i].y = (int)(ss1Seq[i].y * (pow(2, 9) - 1));
        ss2Seq[i].x = (int)(ss2Seq[i].x * (pow(2, 9) - 1));
        ss2Seq[i].y = (int)(ss2Seq[i].y * (pow(2, 9) - 1));
    }

    // 4. find the Max
   // 4.1.1 get ss1 max
    int searchStartPos = 0;
    dcomplex corrRes = batchMuxAndSum(&rxAirSamp[searchStartPos], ss1Seq, SS_ZC_Samp_LEN * K_os_UpSampling);   // 函数内部对乘法进行量化，累加和在外面量化
    //移位10bit 22->12
    //printf("%f  %f\n", corrRes.x, corrRes.y);
    corrRes.x = (int)corrRes.x >> 7;
    corrRes.y = (int)corrRes.y >> 7;
    //if (fabs(corrRes.x) > 2047 | fabs(corrRes.y) > 2047)
    //    printf("%f  %f\n", corrRes.x, corrRes.y);

    double tmpCorrPwr = getPower(&corrRes);
    //移位11bit 23->12
    tmpCorrPwr = (int)tmpCorrPwr >> 11;
    //printf("%f\n", tmpCorrPwr);

    estSSPos[0] = searchStartPos;
    estPeak[0] = tmpCorrPwr;

    // 1s blind dectect, find ss1 the max
    for (int i = searchStartPos + 1; i < 76800 * g_tp680Para.sys_sampRate; i++)
    {
        corrRes = batchMuxAndSum(&rxAirSamp[i], ss1Seq, SS_ZC_Samp_LEN * K_os_UpSampling);
        //移位10bit 22->12
        corrRes.x = (int)corrRes.x >> 7;
        corrRes.y = (int)corrRes.y >> 7;
        //if (fabs(corrRes.x) > 2047 | fabs(corrRes.y) > 2047)
//    printf("%f  %f\n", corrRes.x, corrRes.y);

        tmpCorrPwr = getPower(&corrRes);
        //移位11bit 23->12
        tmpCorrPwr = (int)tmpCorrPwr >> 11;

        if (estPeak[0] < tmpCorrPwr)
        {
            estPeak[0] = tmpCorrPwr;
            estSSPos[0] = i;
        }

        if (i - estSSPos[0] >= SS_ZC_Samp_LEN * K_os_UpSampling - 1)
        {
            // TODO: judge estPeak[0] is valid??
            break;
        }
    }

    // 2 get ss2 max, blind detect at ss2's ideal position from -MD_OFFSET samp to 3*MD_OFFSET samp
    const int BD_OFFSET = SS_ZC_Samp_LEN * K_os_UpSampling / 8;
    const int BD_WIN_LEN = 4 * BD_OFFSET;

    double* corrPwr;
    corrPwr = (double*)malloc(sizeof(double) * BD_WIN_LEN);
    for (int i = 0; i < BD_WIN_LEN; i++)
    {
        corrPwr[i] = 0;
    }

    // ahead 10 zc sym to detect ss2 max
    const int ssOffset = estSSPos[0] + SS_ZC_Samp_LEN * K_os_UpSampling - BD_OFFSET;
    // find from ZC1 start position
    for (int i = 0; i < BD_WIN_LEN; i++)
    {
        corrRes = batchMuxAndSum(&rxAirSamp[ssOffset + i], ss2Seq, SS_ZC_Samp_LEN * K_os_UpSampling);
        //移位10bit 22->12
        corrRes.x = (int)corrRes.x >> 7;
        corrRes.y = (int)corrRes.y >> 7;

        corrPwr[i] = getPower(&corrRes);
        //移位11bit 23->12
        corrPwr[i] = (int)corrPwr[i] >> 11;        
    }

    double maxValue = 0.0;
    findMax(&estPeak[1], &estSSPos[1], corrPwr, BD_WIN_LEN);
    estSSPos[1] += ssOffset;

    // 3. estimate the freq offset
    double delta_f = ((double)u_ss1 * u_ss2 / (u_ss1 - u_ss2)) * (estSSPos[1] - estSSPos[0] - SS_ZC_Samp_LEN * K_os_UpSampling) / (K_os_UpSampling * K_os_UpSampling * SS_ZC_Samp_LEN);
    int startSS1 = (u_ss2 * estSSPos[1] - u_ss1 * estSSPos[0] - u_ss2 * SS_ZC_Samp_LEN * K_os_UpSampling) / (u_ss2 - u_ss1);
    startSS1 = (startSS1 < 0) ? 0 : startSS1;

    if (0 != g_tp680Para.shr_freq_bias)
    {
        // 4. correct the freq bias
        const double phaseOffset = -PI2 * delta_f;
        for (int i = 0; i < rxSampLen; i++)
        {
            dcomplex temp = POLAR(phaseOffset * i);
            temp.x = (int)(temp.x * (pow(2, 9) - 1));
            temp.y = (int)(temp.y * (pow(2, 9) - 1));

            rxAirSamp[startSS1 + i] = cMult(&rxAirSamp[startSS1 + i], &temp);
            rxAirSamp[startSS1 + i].x = (int)rxAirSamp[startSS1 + i].x >> 9;
            rxAirSamp[startSS1 + i].y = (int)rxAirSamp[startSS1 + i].y >> 9;
        }
    }


    // dec phr
    const int N_bitOfPhr = 64;    
    int phrOff = startSS1 + 2 * SS_ZC_Samp_LEN * K_os_UpSampling;
    int estZCPos = 0;       // for phr zc pos each grp
    double estZCPeak = 0;   // for phr zc peak each grp

    int* rxBytePhr;
    rxBytePhr = (int*)malloc(sizeof(int) * N_bitOfPhr / 8);
    for (int i = 0; i < N_bitOfPhr / 8; i++)
    {
        rxBytePhr[i] = 0;
    }

    dcomplex estZCVal;
    double* corrMaxAngle;
    corrMaxAngle = (double*)malloc(g_tp680Para.rx_L_numOfGroup * sizeof(double));
    for (int i = 0; i < g_tp680Para.rx_L_numOfGroup; i++)
    {
        corrMaxAngle[i] = 0;
    }

    int* rxNrzPhr;
    rxNrzPhr = (int*)malloc(sizeof(int) * g_tp680Para.rx_I_zc * g_tp680Para.rx_L_numOfGroup); // g_tp680Para.rx_I_zc,g_tp680Para.rx_L_numOfGroup 从发送端得知
    for (int i = 0; i < g_tp680Para.rx_I_zc * g_tp680Para.rx_L_numOfGroup; i++)
    {
        rxNrzPhr[i] = 0;
    }

    for (int grp = 0; grp < g_tp680Para.rx_L_numOfGroup; grp++)
    {
        // 4.1.2 get zc1 max, notice that must cycle shift process
        corrRes = batchMuxAndSum(&rxAirSamp[phrOff], ss2Seq, SS_ZC_Samp_LEN * K_os_UpSampling);
        //移位10bit 22->12
        corrRes.x = (int)corrRes.x >> 7;
        corrRes.y = (int)corrRes.y >> 7;

        estZCPos = 0;
        estZCPeak = getPower(&corrRes);
        //移位11bit 23->12
        estZCPeak = (int)estZCPeak >> 11;

        estZCVal = corrRes;

        dcomplex* ptrRxSamp = &rxAirSamp[phrOff];

        //for (int i = 1; i < SS_ZC_Samp_LEN * K_os_UpSampling; i++)
        for (int i = K_os_UpSampling; i < SS_ZC_Samp_LEN * K_os_UpSampling; i += K_os_UpSampling)
        {
            corrRes = set(0.0, 0.0);
            for (int j = 0; j < SS_ZC_Samp_LEN * K_os_UpSampling; ++j)
            {
                dcomplex temp;
                temp = cMult(&ptrRxSamp[(i + j) % (SS_ZC_Samp_LEN * K_os_UpSampling)], &ss2Seq[j]);
                temp.x = (int)temp.x >> 7;
                temp.y = (int)temp.y >> 7;
                //if (fabs(temp.x) > 2047 | fabs(temp.y) > 2047)
                //    printf("%f  %f\n", temp.x, temp.y);

                corrRes = cAdd(&corrRes, &temp);
            }
            //移位10bit 22->12
            //printf("%f  %f\n", corrRes.x, corrRes.y);
            corrRes.x = (int)corrRes.x >> 7;
            corrRes.y = (int)corrRes.y >> 7;
            //if (fabs(corrRes.x) > 2047 | fabs(corrRes.y) > 2047)
            //    printf("%f  %f\n", corrRes.x, corrRes.y);

            tmpCorrPwr = getPower(&corrRes);
            //移位11bit 23->12
            tmpCorrPwr = (int)tmpCorrPwr >> 11;
            //if(fabs(tmpCorrPwr) > 2047)
            //    printf("%f\n", tmpCorrPwr);

            if (estZCPeak < tmpCorrPwr)
            {
                estZCPeak = tmpCorrPwr;
                estZCPos = i;
                estZCVal = corrRes;
            }
        }

        corrMaxAngle[grp] = atan2(estZCVal.y, estZCVal.x);

        int tmpPos = (SS_ZC_Samp_LEN * K_os_UpSampling - estZCPos) % (SS_ZC_Samp_LEN * K_os_UpSampling);   // as tx's (Cv * K_os_UpSampling)
        int v_bit = (int)((g_tp680Para.phr_N_cs * K_os_UpSampling / 2 + tmpPos) / (g_tp680Para.phr_N_cs * K_os_UpSampling));

        for (int bit = 0; bit < g_tp680Para.rx_I_zc; bit++)
        {
            rxNrzPhr[grp * g_tp680Para.rx_I_zc + bit] = 1 - 2 * (int)((v_bit >> (g_tp680Para.rx_I_zc - 1 - bit)) & 0x01);
        }

        phrOff += SS_ZC_Samp_LEN * K_os_UpSampling;
    }

    double avgDiff = 0;
    int numOfDiff = 0;
    double tmpDiff = 0;
    for (int grp = 1; grp < g_tp680Para.rx_L_numOfGroup; grp++)
    {
        tmpDiff = corrMaxAngle[grp] - corrMaxAngle[grp - 1];

        if ((tmpDiff > -3) && (tmpDiff < 3))
        {
            avgDiff += tmpDiff;
            numOfDiff++;
        }
    }

    if (numOfDiff > 0)
    {
        avgDiff /= numOfDiff;
    }

    if (0 != g_tp680Para.shr_freq_bias)
    {
        // 5. correct the freq bias according PHR est
        const double phrPhaseOffset = -PI2 * avgDiff / (PI2 * SS_ZC_Samp_LEN * K_os_UpSampling);
        for (int i = 0; i < rxSampLen; i++)
        {
            dcomplex temp = POLAR(phrPhaseOffset * i);
            temp.x = (int)(temp.x * (pow(2, 9) - 1));
            temp.y = (int)(temp.y * (pow(2, 9) - 1));
            rxAirSamp[startSS1 + i] = cMult(&rxAirSamp[startSS1 + i], &temp);
            //printf("%f  %f\n", rxAirSamp[startSS1 + i].x, rxAirSamp[startSS1 + i].y);
            //取整
            rxAirSamp[startSS1 + i].x = (int)rxAirSamp[startSS1 + i].x >> 9;
            rxAirSamp[startSS1 + i].y = (int)rxAirSamp[startSS1 + i].y >> 9;
        }
    }

    inrzToByte(rxNrzPhr, rxBytePhr, N_bitOfPhr);
    if (0 == decCrc8_UNB(rxBytePhr, 7))
    {
        //printf(" *************  PHR Dec Success **************** payload off[%d]\r\n", phrOff);
        return phrOff;
    }
    else
    {
        return -1;
    }

    free(corrPwr);
    free(rxBytePhr);
    free(corrMaxAngle);
    free(rxNrzPhr);

}