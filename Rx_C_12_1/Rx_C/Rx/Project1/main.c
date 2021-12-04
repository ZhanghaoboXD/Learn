#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"
#include "tp680.h"

int decShrPhrSeq(dcomplex* rxAirSamp, int rxSampLen);
void getPseudoRandomSeq(int* pOut_NRZ, int Cinit, int len);
void initPhase(int sysBps, int divation, int sampRate);
void deMskFsk(dcomplex* pIn, double* pOutNRZ, int outBitLen, int type);
void nrzDeSfProcess_Msk(double* outNrzData, double* inSfedData, int outLen, int* sfNrzSeq, int sfFactor);
void decFilter(dcomplex* pIn, dcomplex* pOut, int inLen);
void nrzDeSfProcess_QPSK(dcomplex* outNrzData, dcomplex* inSfedData, int outLen, int* sfNrzSeq, int sfFactor);
void decBPSK(dcomplex* pIn, double* pOutNRZ, int outBitLen);
void decQPSK(dcomplex* pIn, double* pOutNRZ, int outBitLen);
void dec8PSK(dcomplex* pIn, double* pOutNRZ, int outBitLen);
void nrzGetPN9SeqUNB(double* ioBuff, int bitLen);
void interleaveProc(double* ccOutNrzByte, unsigned int ccOutBitLen);
void Quantify(double* input, int* output, double delta, int length);
void decodeCC(int rate, int* din1, int len, int* output);
void inrzToByte(const int* inputNrzBit, int* outByte, int bitLen);
int decCrc16_UNB(const int* inputBit, int rawByteLen);
int decCrc24_UNB(const int* inputBit, int rawByteLen);
dcomplex POLAR(double sita);
void doFading(dcomplex* pRxOut, dcomplex* pTxIn, double SNR, int numOfSamp);

dcomplex* g_txDataSamp = NULL;

void initGlobalPara(void)
{
	g_tp680Para.sys_symBps = 76800;         // 1 sym = 307200 / 76800 = 4 Ts
	g_tp680Para.sys_sampRate = 8;
	g_tp680Para.shr_freq_bias = 0;       // 4k freq bias

	g_tp680Para.shr_N_zc = 128;             // 16, 64 or 128
	g_tp680Para.phr_N_cs = 1;               // 1, 2, 4 or 8

	g_tp680Para.n_sf = 4;
	g_tp680Para.modType = 5;             // 0: MSK          5: QPSK
	g_tp680Para.codeType = 0;
	g_tp680Para.ccPreLen = 20;

	g_tp680Para.rx_L_numOfGroup = 10;
	g_tp680Para.rx_I_zc = 7;
	g_tp680Para.rx_NumBit = 2816;     // MSK:5632  QPSK:2816
}

double g_glbDelta = 0;
double g_glbPhi = 0;

void decFreqTracking(dcomplex* rxSampSeq, dcomplex* txRawSamp)
{
	const int numOfDataInTTI = 2944 * g_tp680Para.sys_symBps / 307200;
	const int numOfPilotInTTI = 128 * g_tp680Para.sys_symBps / 307200;

	const int lpSize = 8;
	const int N_numOfEstSamp = lpSize * g_tp680Para.sys_sampRate;
	const int numIter = numOfPilotInTTI / lpSize;

	// 1st pilot check and adjust freq bias
	double prePhase = 0;
	double cfo_est = 0;
	for (int pilotIdx = 0; pilotIdx < numIter; pilotIdx++)
	{
		const int pilotOff = pilotIdx * N_numOfEstSamp;

		double z_n = 0;
		for (int i = 0; i < N_numOfEstSamp; i++)
		{
			dcomplex temp1 = POLAR(prePhase);
			//量化
			temp1.x = (int)(temp1.x * (pow(2, 4) - 1));
			temp1.y = (int)(temp1.y * (pow(2, 4) - 1));

			dcomplex temp2 = cMult(&txRawSamp[pilotOff + i], &temp1);
			dcomplex temp3 = getConj(&temp2);
			dcomplex cfo_part = cMult(&rxSampSeq[pilotOff + i], &temp3);

			//右移1bit 13->12
			cfo_part.x = (int)cfo_part.x >> 2;
			cfo_part.y = (int)cfo_part.y >> 2;
			//if (fabs(cfo_part.x) > 2047 | fabs(cfo_part.y) > 2047)
			//	printf("%f  %f\n", cfo_part.x, cfo_part.y);
			double curPhase = atan2(cfo_part.y, cfo_part.x);
			z_n += curPhase;
		}
		// per each bit, phi_k / N
		z_n /= N_numOfEstSamp;

		prePhase += z_n;

		//printf("phi :%.5f ", z_n / (N_numOfEstSamp * PI2));
		if (0 != pilotIdx)
		{
			cfo_est += z_n / N_numOfEstSamp;
		}
	}
	cfo_est /= (numIter - 1);
	//printf("cfo_est :%.5f \n", cfo_est / PI2);

	const double delta_f = -cfo_est;

	// sigma(delta_f_0 .. delta_f_tti) pilot est result, for cur tti's data
	g_glbDelta += delta_f;

	// jump the cur tti's pilot freq phi
	g_glbPhi += delta_f * numOfPilotInTTI * g_tp680Para.sys_sampRate;

	// samp freq adjust for cur tti's data + next tti's pilot
	for (int samp = 0; samp < (numOfPilotInTTI + numOfDataInTTI) * g_tp680Para.sys_sampRate; samp++)
	{
		dcomplex temp = POLAR(g_glbPhi);
		//量化
		temp.x = (int)(temp.x * (pow(2, 4) - 1));
		temp.y = (int)(temp.y * (pow(2, 4) - 1));
		rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate] = cMult(&rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate], &temp);
		//右移3bit 8->5
		rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].x = (int)rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].x >> 3;
		rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].y = (int)rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].y >> 3;
		//if (fabs(rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].x) > 31 | fabs(rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].y) > 31)
		//	printf("%f  %f\n", rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].x, rxSampSeq[samp + numOfPilotInTTI * g_tp680Para.sys_sampRate].y);
		// store the g_glbPhi = sigma(phi) at end of the next's data
		g_glbPhi += g_glbDelta;
	}
}

void decPhaseAdjust(dcomplex* rxSampSeq, dcomplex* txRawSamp, int ttiIdx)
{
	const int numOfDataInTTI = 2944 * g_tp680Para.sys_symBps / 307200;
	const int numOfPilotInTTI = 128 * g_tp680Para.sys_symBps / 307200;
	const int numOfTTI = (g_tp680Para.rx_NumBit + numOfDataInTTI - 1) / numOfDataInTTI; // rx_NumBit 由发送端得知

	const int lpSize = numOfPilotInTTI;
	const int N_numOfEstSamp = lpSize * g_tp680Para.sys_sampRate;
	const int numIter = numOfPilotInTTI / lpSize;

	// 1st pilot check and adjust freq bias
	double prePhase = 0;
	for (int pilotIdx = 0; pilotIdx < numIter; pilotIdx++)
	{
		const int pilotOff = pilotIdx * N_numOfEstSamp;

		double z_n = 0;
		for (int i = 0; i < N_numOfEstSamp; i++)
		{
			dcomplex temp1 = POLAR(prePhase);
			//量化
			temp1.x = (int)(temp1.x * (pow(2, 4) - 1));
			temp1.y = (int)(temp1.y * (pow(2, 4) - 1));

			dcomplex temp2 = cMult(&txRawSamp[pilotOff + i], &temp1);
			dcomplex temp3 = getConj(&temp2);
			dcomplex cfo_part = cMult(&rxSampSeq[pilotOff + i], &temp3);
			//右移3bit 15->12
			cfo_part.x = (int)cfo_part.x >> 2;
			cfo_part.y = (int)cfo_part.y >> 2;
	//		if (fabs(cfo_part.x) > 2047 | fabs(cfo_part.y) > 2047)
	//printf("%f  %f\n", cfo_part.x, cfo_part.y);
			double curPhase = atan2(cfo_part.y, cfo_part.x);
			
			z_n += curPhase;
		}
		// per each bit, phi_k / N
		z_n /= N_numOfEstSamp;

		prePhase += z_n;
	}
	prePhase /= numIter;

	for (int samp = 0; samp < (numOfPilotInTTI + numOfDataInTTI) * (numOfTTI - ttiIdx) * g_tp680Para.sys_sampRate; samp++)
	{
		dcomplex temp = POLAR(-prePhase);
		temp.x = (int)(temp.x * (pow(2, 4) - 1));
		temp.y = (int)(temp.y * (pow(2, 4) - 1));

		rxSampSeq[samp] = cMult(&rxSampSeq[samp], &temp);
		//取整
		rxSampSeq[samp].x = (int)rxSampSeq[samp].x >> 4;
		rxSampSeq[samp].y = (int)rxSampSeq[samp].y >> 4;
	//	if (fabs(rxSampSeq[samp].x) > 31 | fabs(rxSampSeq[samp].y) > 31)
	//printf("%f  %f\n", rxSampSeq[samp].x, rxSampSeq[samp].y);
	}
}


int g_spreadFactor[76800] = { 0 };
int decPayload(dcomplex* rxSampSeq, int* output)
{
	const int numOfDataInTTI = 2944 * g_tp680Para.sys_symBps / 307200;
	const int numOfPilotInTTI = 128 * g_tp680Para.sys_symBps / 307200;
	const int numOfTTI = (g_tp680Para.rx_NumBit + numOfDataInTTI - 1) / numOfDataInTTI; // rx_NumBit 由发送端得知
	int postModOff = 0;
	int preDataOff = 0;

	double rxPilotNrzBit[128] = { 0 };
	const int C_SF_FACTOR = 1 << g_tp680Para.n_sf;

	getPseudoRandomSeq(g_spreadFactor, 2, g_tp680Para.sys_symBps);

	int postBitLen = 0;
	double* deSfBitData = NULL;

	if (0 == g_tp680Para.modType || 2 == g_tp680Para.modType) // 0:MSK  2:FSK
	{

		double* rxNrzBit;
		rxNrzBit = (double*)malloc(numOfTTI * numOfDataInTTI * sizeof(double));
		for (int i = 0; i < numOfTTI * numOfDataInTTI; i++)
		{
			rxNrzBit[i] = 0;
		}

		g_glbDelta = 0;
		g_glbPhi = 0;

		initPhase(g_tp680Para.sys_symBps, 4, g_tp680Para.sys_sampRate);
		for (int tti = 0; tti < numOfTTI; tti++)
		{
			// 1st pilot
			decFreqTracking(&rxSampSeq[postModOff], &g_txDataSamp[postModOff]);
			postModOff += numOfPilotInTTI * g_tp680Para.sys_sampRate;

			// 2nd data, decMskOnly or decFskMsk
			deMskFsk(&rxSampSeq[postModOff], &rxNrzBit[preDataOff], numOfDataInTTI, g_tp680Para.modType);
			preDataOff += numOfDataInTTI;
			postModOff += numOfDataInTTI * g_tp680Para.sys_sampRate;
		}

		g_glbDelta = 0;
		g_glbPhi = 0;

		postBitLen = g_tp680Para.rx_NumBit / C_SF_FACTOR;
		deSfBitData = (double*)malloc(postBitLen * sizeof(double));
		for (int i = 0; i < postBitLen; i++)
		{
			deSfBitData[i] = 0;
		}

		nrzDeSfProcess_Msk(deSfBitData, rxNrzBit, postBitLen, g_spreadFactor, C_SF_FACTOR);

		free(rxNrzBit);
	}
	else if (4 == g_tp680Para.modType || 5 == g_tp680Para.modType || 6 == g_tp680Para.modType)
	{
		const int C_PSK_MODRATE = g_tp680Para.modType - 3;

		for (int tti = 0; tti < numOfTTI; tti++)
		{
			// 1st pilot tracking
			decPhaseAdjust(&rxSampSeq[postModOff], &g_txDataSamp[postModOff], tti);
			postModOff += numOfPilotInTTI * g_tp680Para.sys_sampRate;
			// 2nd data
			postModOff += numOfDataInTTI * g_tp680Para.sys_sampRate;
		}

		const int C_RE_SYM = (numOfDataInTTI + numOfPilotInTTI) * numOfTTI;

		dcomplex* reMapSym;
		reMapSym = (dcomplex*)malloc(C_RE_SYM * sizeof(dcomplex));
		for (int i = 0; i < C_RE_SYM; i++)
		{
			reMapSym[i].x = 0;
			reMapSym[i].y = 0;
		}

		decFilter(rxSampSeq, reMapSym, C_RE_SYM * g_tp680Para.sys_sampRate);

		dcomplex* rxSfdeMod;
		rxSfdeMod = (dcomplex*)malloc(numOfDataInTTI * numOfTTI * sizeof(dcomplex));
		for (int i = 0; i < numOfDataInTTI * numOfTTI; i++)
		{
			rxSfdeMod[i].x = 0;
			rxSfdeMod[i].y = 0;
		}

		dcomplex* rxPilotSym;
		rxPilotSym = (dcomplex*)malloc(numOfPilotInTTI * sizeof(dcomplex));
		for (int i = 0; i < numOfPilotInTTI; i++)
		{
			rxPilotSym[i].x = 0;
			rxPilotSym[i].y = 0;
		}

		postModOff = 0;
		preDataOff = 0;
		for (int tti = 0; tti < numOfTTI; tti++)
		{
			// 1st pilot
			memcpy(rxPilotSym, &reMapSym[postModOff], numOfPilotInTTI * sizeof(dcomplex));
			postModOff += numOfPilotInTTI;

			// 2nd data
			memcpy(&rxSfdeMod[preDataOff], &reMapSym[postModOff], numOfDataInTTI * sizeof(dcomplex));
			postModOff += numOfDataInTTI;

			preDataOff += numOfDataInTTI;
		}

		dcomplex* rxDeSfMod;
		rxDeSfMod = (dcomplex*)malloc(g_tp680Para.rx_NumBit / C_SF_FACTOR * sizeof(dcomplex));
		for (int i = 0; i < g_tp680Para.rx_NumBit / C_SF_FACTOR; i++)
		{
			rxDeSfMod[i].x = 0;
			rxDeSfMod[i].y = 0;
		}

		nrzDeSfProcess_QPSK(rxDeSfMod, rxSfdeMod, g_tp680Para.rx_NumBit / C_SF_FACTOR, g_spreadFactor, C_SF_FACTOR);

		postBitLen = g_tp680Para.rx_NumBit * C_PSK_MODRATE / C_SF_FACTOR;
		deSfBitData = (double*)malloc(postBitLen * sizeof(double));
		for (int i = 0; i < postBitLen; i++)
		{
			deSfBitData[i] = 0;
		}

		if (C_PSK_MODRATE == 1)
		{
			decBPSK(rxDeSfMod, deSfBitData, postBitLen);
		}
		else if (C_PSK_MODRATE == 2)
		{
			decQPSK(rxDeSfMod, deSfBitData, postBitLen);
		}
		else if (C_PSK_MODRATE == 3)
		{
			dec8PSK(rxDeSfMod, deSfBitData, postBitLen);
		}

		free(reMapSym);
		free(rxSfdeMod);
		free(rxPilotSym);
		free(rxDeSfMod);
	}

	nrzGetPN9SeqUNB(deSfBitData, postBitLen);

	interleaveProc(deSfBitData, postBitLen);

	//for (int i = 0; i < postBitLen; i++)
	//{
	//	printf("%f ",deSfBitData[i]);
	//	if ((i + 1) % 8 == 0)printf("\n");
	//}

	int Len = (((postBitLen / 2) - 1) / 64 + 1) * 64;
	int* cbInput;
	cbInput = (int*)malloc(Len * sizeof(int));
	for (int i = 0; i < Len; i++)
	{
		cbInput[i] = 0;
	}

	int* rxRawByte;
	int rawDataLen_B = postBitLen / 8 / 2;
	rxRawByte = (int*)malloc(rawDataLen_B * sizeof(int));
	for (int i = 0; i < rawDataLen_B; i++)
	{
		rxRawByte[i] = 0;
	}

	int* Q;
	Q = (int*)malloc(postBitLen * sizeof(int));
	Quantify(deSfBitData, Q, 1, postBitLen);  //delta适配暂定1

	int* decodein;
	decodein = (int*)malloc(sizeof(int) * postBitLen / 2);
	for (int i = 0; i < postBitLen / 2; i++)
	{
		decodein[i] = (Q[2 * i] << 4) + Q[2 * i + 1];
	}

	decodeCC(0, decodein, postBitLen / 2, cbInput);   //0: 1/2  1: 2/3  2: 3/4
	for (int i = 0; i < Len; i++)
	{
		cbInput[i] = 2 * cbInput[i] - 1;
	}

	inrzToByte(cbInput, rxRawByte, postBitLen / 2);

	for (int i = 0; i < 160; i++)
	{
		output[i] = cbInput[i];
	}


	int result = -1;
	if (rawDataLen_B <= 130)
	{
		result = decCrc16_UNB(rxRawByte, rawDataLen_B - 2);
	}
	else
	{
		result = decCrc24_UNB(rxRawByte, rawDataLen_B - 3);
	}

	if (0 == result)
	{
		//printf("  ****************  Payload Dec Success  ****************\r\n");
		return 0;
	}
	else
	{
		return -1;
	}

	free(deSfBitData);
	free(cbInput);
	free(rxRawByte);
	free(Q);
	free(decodein);
}

double tmp[36864][2] = { 0 };//MSK:61440 QPSK:36864
int main()
{
	initGlobalPara(); // 选MSK还是QPSK需要更改参数

	FILE* fp = NULL;
	if (g_tp680Para.modType == 0)
	{
		fp = fopen("1_2_MSK.txt", "r");
	}
	else if (g_tp680Para.modType == 5)
	{
		fp = fopen("1_2_QPSK.txt", "r");
	}
	

	dcomplex* simTxSeq;
	simTxSeq = (dcomplex*)malloc((g_tp680Para.sys_symBps * g_tp680Para.sys_sampRate + 50) * sizeof(dcomplex));
	for (int i = 0; i < g_tp680Para.sys_symBps * g_tp680Para.sys_sampRate + 50; i++)
	{
		simTxSeq[i].x = 0;
		simTxSeq[i].y = 0;
	}

	if (fp == NULL)
	{
		printf("READ ERROR");
		return -1;
	}

	for (int i = 0; i < 36864; i++)  //MSK:61440 QPSK:36864
	{
		for (int j = 0; j < 2; j++)
		{
			fscanf(fp, "%lf", &tmp[i][j]);
		}
	}
	fclose(fp);

	for (int i = 0; i < 36864; i++)  //MSK:61440 QPSK:36864
	{
		simTxSeq[i].x = tmp[i][0];
		simTxSeq[i].y = tmp[i][1];
	}

	int SHR_SAMP_LEN = 2048;
	int PAYLOAD_SAMP_OFF = 12288;
	int txSampOff = 36864;   //接收数据长度, 由发送端得知     测试1/2码率 MSK:61440 QPSK:36864

	g_txDataSamp = (dcomplex*)malloc((txSampOff - PAYLOAD_SAMP_OFF) * sizeof(dcomplex));
	for (int i = 0; i < txSampOff - PAYLOAD_SAMP_OFF; i++)
	{
		g_txDataSamp[i].x = simTxSeq[PAYLOAD_SAMP_OFF + i].x;
		g_txDataSamp[i].y = simTxSeq[PAYLOAD_SAMP_OFF + i].y;
		//printf("%f %f\n", g_txDataSamp[i].x, g_txDataSamp[i].y);
	}

	//量化
	for (int i = 0; i < txSampOff - PAYLOAD_SAMP_OFF; i++)
	{
		g_txDataSamp[i].x = (int)(g_txDataSamp[i].x * (pow(2, 4) - 1));
		g_txDataSamp[i].y = (int)(g_txDataSamp[i].y * (pow(2, 4) - 1));
	}


	dcomplex* rxAirSamp;
	rxAirSamp = (dcomplex*)malloc((txSampOff + SHR_SAMP_LEN) * sizeof(dcomplex));
	rxAirSamp = (dcomplex*)malloc((txSampOff + SHR_SAMP_LEN) * sizeof(dcomplex));
	for (int i = 0; i < txSampOff + SHR_SAMP_LEN; i++)
	{
		rxAirSamp[i].x = 0;
		rxAirSamp[i].y = 0;
	}

	int input[160] = { 1,1,-1,1,-1,1,1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,1,1,1,-1,-1,1,-1,-1,1,1,-1,1,1,-1,1,1 };
	int output[160] = { 0 };
	for (double SNR = -20; SNR < 0; SNR++)
	{
		int error = 0;
		double BER = 0;
		for (double frame_num = 1; frame_num <= 1000; frame_num++)
		{
			doFading(rxAirSamp, simTxSeq, SNR, txSampOff);
			//printf("SNR = %f\n", SNR);

			int payloadPos = decShrPhrSeq(&rxAirSamp[0], txSampOff);

			if (payloadPos > 0)
			{
				for (int i = 0; i < txSampOff; i++)  // 暂定移5位
				{
					rxAirSamp[i].x = (int)rxAirSamp[i].x >> 5;
					rxAirSamp[i].y = (int)rxAirSamp[i].y >> 5;
				}
				decPayload(&rxAirSamp[payloadPos], output);
			}

			for (int i = 0; i < 160; i++)
			{
				if (output[i] != input[i])
				{
					error++;
				}
			}
			
			BER = (double)error / (double)(160 * frame_num);
		}

		printf("SNR = %f  BER = %e\n", SNR, BER);
	}
	free(simTxSeq);
	free(g_txDataSamp);
	free(rxAirSamp);
}