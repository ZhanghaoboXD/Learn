#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Quantify(double* input, int* output, double delta, int length)
{
	for (int i = 0; i < length; i++)
	{
		int j = floor(input[i] / delta);
		if (j >= 7)
		{
			output[i] = 15;
		}
		else if (j < -7)
		{
			output[i] = 0;
		}
		else
		{
			output[i] = j + 8;
		}
	}
}

int viterbi(int len, int* din, int swit, int** mem, int times, int* si, int** so)
{
	int final_state = 0;
	int d[64] = { 0 };
	int so_tmp[64] = { 0 };
	int min_s;

	for (int i = 0; i < len; i++)
	{
		int r0 = 0;
		int r1 = 0;
		int bm[4] = { 0 };
		//每次输入(按r1低位r0高位来)			
		r0 = (din[i] >> 4) & 0xf;
		r1 = din[i] & 0xf;
		//每次分支度量
		bm[0] = r0 + r1;
		bm[1] = 15 + r0 - r1;
		bm[2] = 15 - r0 + r1;
		bm[3] = 30 - r0 - r1;
		//确定当前时刻路径度量值
		for (int j = 0; j < 64; j++)
		{
			so[j][0] = si[j] + bm[((j & 0x01) + ((j & 0x02) >> 1) + ((j & 0x04) >> 2) + ((j & 0x20) >> 5)) % 2 + (((j & 0x02) >> 1) + ((j & 0x04) >> 2) + ((j & 0x10) >> 4) + ((j & 0x20) >> 5)) % 2 * 2];
			so[j][1] = si[j] + bm[3 - (((j & 0x01) + ((j & 0x02) >> 1) + ((j & 0x04) >> 2) + ((j & 0x20) >> 5)) % 2 + (((j & 0x02) >> 1) + ((j & 0x04) >> 2) + ((j & 0x10) >> 4) + ((j & 0x20) >> 5)) % 2 * 2)];
		}
		//比较出d，将d存放在mem中
		for (int k = 0; k < 32; k++)
		{
			so_tmp[2 * k] = (so[k][0] <= so[k + 32][0]) ? so[k][0] : so[32 + k][0];
			so_tmp[2 * k + 1] = (so[k][1] <= so[k + 32][1]) ? so[k][1] : so[32 + k][1];
			d[2 * k] = (so[k][0] > so[k + 32][0]) ? 1 : 0;
			d[2 * k + 1] = (so[k][1] > so[k + 32][1]) ? 1 : 0;
			if (swit == 1)
			{
				if (times == 0)
				{
					mem[i][2 * k] = d[2 * k];
					mem[i][2 * k + 1] = d[2 * k + 1];
				}
				else
				{
					mem[(i + 64 * times - 64) % 128][2 * k] = d[2 * k];
					mem[(i + 64 * times - 64) % 128][2 * k + 1] = d[2 * k + 1];
				}
			}
		}
		//比较最小值min_s，最终路径final_state
		int m[64] = { 0 };
		int final_idx = 0;
		for (int k = 0; k < 32; k++)
		{
			m[2 * k] = (d[2 * k] == 1) ? so[k + 32][0] : so[k][0];
			m[2 * k + 1] = (d[2 * k + 1] == 1) ? so[k + 32][1] : so[k][1];
		}
		int mm[32] = { 0 };
		for (int k = 0; k < 32; k++)
		{
			mm[k] = (m[2 * k] < m[2 * k + 1]) ? m[2 * k] : m[2 * k + 1];
		}
		int mmm[16] = { 0 };
		for (int k = 0; k < 16; k++)
		{
			mmm[k] = (mm[2 * k] < mm[2 * k + 1]) ? mm[2 * k] : mm[2 * k + 1];
		}
		int mmmm[8] = { 0 };
		for (int k = 0; k < 8; k++)
		{
			mmmm[k] = (mmm[2 * k] < mmm[2 * k + 1]) ? mmm[2 * k] : mmm[2 * k + 1];
		}
		int mmmmm[4] = { 0 };
		for (int k = 0; k < 4; k++)
		{
			mmmmm[k] = (mmmm[2 * k] < mmmm[2 * k + 1]) ? mmmm[2 * k] : mmmm[2 * k + 1];
		}
		int z[2] = { 0 };
		for (int k = 0; k < 2; k++)
		{
			z[k] = (mmmmm[2 * k] < mmmmm[2 * k + 1]) ? mmmmm[2 * k] : mmmmm[2 * k + 1];
		}
		min_s = (z[0] < z[1]) ? z[0] : z[1];
		final_state = (min_s == z[1]) ? 0x20 : 0x00;
		final_state += (min_s == mmmmm[1] || min_s == mmmmm[3]) ? 0x10 : 0x00;
		final_state += (min_s == mmmm[1] || min_s == mmmm[3] || min_s == mmmm[5] || min_s == mmmm[7]) ? 0x08 : 0x00;
		final_state += (min_s == mmm[1] || min_s == mmm[3] || min_s == mmm[5] || min_s == mmm[7] || min_s == mmm[9] || min_s == mmm[11] || min_s == mmm[13] || min_s == mmm[15]) ? 0x04 : 0x00;
		final_state += (min_s == mm[1] || min_s == mm[3] || min_s == mm[5] || min_s == mm[7] || min_s == mm[9] || min_s == mm[11] || min_s == mm[13] || min_s == mm[15] ||
			min_s == mm[17] || min_s == mm[19] || min_s == mm[21] || min_s == mm[23] || min_s == mm[25] || min_s == mm[27] || min_s == mm[29] || min_s == mm[31]) ? 0x02 : 0x00;
		final_state += (min_s == m[1] || min_s == m[3] || min_s == m[5] || min_s == m[7] || min_s == m[9] || min_s == m[11] || min_s == m[13] || min_s == m[15] || min_s == m[17] || min_s == m[19] ||
			min_s == m[21] || min_s == m[23] || min_s == m[25] || min_s == m[27] || min_s == m[29] || min_s == m[31] || min_s == m[33] || min_s == m[35] || min_s == m[37] || min_s == m[39] || min_s == m[41] ||
			min_s == m[43] || min_s == m[45] || min_s == m[47] || min_s == m[49] || min_s == m[51] || min_s == m[53] || min_s == m[55] || min_s == m[57] || min_s == m[59] || min_s == m[61] || min_s == m[63]) ? 0x01 : 0x00;

		//防止溢出，减去min_s,给下一次循环赋值
		for (int n = 0; n < 64; n++)
		{
			si[n] = (so_tmp[n] - min_s) & 0x7ff;
		}
	}

	return final_state;
}

void decodeCC(int rate, int* din1, int len, int* output)
{
	int len0 = 0;
	if (rate == 0)
	{
		len0 = len;
	}
	else if (rate == 1)
	{
		if ((len * 4 / 3) % 4 == 0)
		{
			len0 = len * 4 / 3;
		}
		if ((len * 4 / 3) % 4 == 1)
		{
			len0 = len * 4 / 3 + 3;
		}
		if ((len * 4 / 3) % 4 == 2)
		{
			len0 = len * 4 / 3 + 2;
		}
	}
	else if (rate == 2)
	{
		len0 = len * 3 / 2;
	}

	int* din;
	din = (int*)malloc((len0) * sizeof(int));
	for (int i = 0; i < len0; i++)
	{
		din[i] = 0;
	}
	//获取不同码率的输入
	if (rate == 0)
	{
		for (int i = 0; i < len0; i++)
		{
			din[i] = din1[i];
		}
	}
	else if (rate == 1)
	{
		if ((len * 4 / 3) % 4 == 0)
		{
			int j = 0;
			for (int i = 0; i < len0; i++)
			{
				if (i % 4 == 0)
				{
					din[i] = din1[j % len];
					j++;
				}
				else if (i % 4 == 1)
				{
					din[i] = (din1[j % len] & 0xf0) + 0x07;
				}
				else if (i % 4 == 2)
				{
					din[i] = ((din1[j % len] & 0x0f) << 4) + ((din1[(j + 1) % len] & 0xf0) >> 4);
					j++;
				}
				else if (i % 4 == 3)
				{
					din[i] = ((din1[j % len] & 0x0f) << 4) + 0x08;
					j++;
				}
			}
		}
		if ((len * 4 / 3) % 4 == 1)
		{
			int j = 0;
			for (int i = 0; i < len0; i++)
			{
				if (i < len0 - 3)
				{
					if (i % 4 == 0)
					{
						din[i] = din1[j % len];
						j++;
					}
					else if (i % 4 == 1)
					{
						din[i] = (din1[j % len] & 0xf0) + 0x07;
					}
					else if (i % 4 == 2)
					{
						din[i] = ((din1[j % len] & 0x0f) << 4) + ((din1[(j + 1) % len] & 0xf0) >> 4);
						j++;
					}
					else if (i % 4 == 3)
					{
						din[i] = ((din1[j % len] & 0x0f) << 4) + 0x08;
						j++;
					}
				}
				else
				{
					if (i % 4 == 1)
					{
						din[i] = (din1[len] & 0xf0) + 0x07;
					}
					else if (i % 4 == 2)
					{
						din[i] = ((din1[len - 2] & 0x0f) << 4) + ((din1[len - 1] & 0xf0) >> 4);
					}
					else if (i % 4 == 3)
					{
						din[i] = ((din1[len - 1] & 0x0f) << 4) + 0x08;
					}
				}
			}
		}
		if ((len * 4 / 3) % 4 == 2)
		{
			int j = 0;
			for (int i = 0; i < len0; i++)
			{
				if (i < len0 - 2)
				{
					if (i % 4 == 0)
					{
						din[i] = din1[j % len];
						j++;
					}
					else if (i % 4 == 1)
					{
						din[i] = (din1[j % len] & 0xf0) + 0x07;
					}
					else if (i % 4 == 2)
					{
						din[i] = ((din1[j % len] & 0x0f) << 4) + ((din1[(j + 1) % len] & 0xf0) >> 4);
						j++;
					}
					else if (i % 4 == 3)
					{
						din[i] = ((din1[j % len] & 0x0f) << 4) + 0x08;
						j++;
					}
				}
				else
				{
					if (i % 4 == 2)
					{
						din[len0 - 2] = ((din1[len - 1] & 0x0f) << 4) + ((din1[len - 3] & 0xf0) >> 4);
					}
					else if (i % 4 == 3)
					{
						din[len0 - 1] = ((din1[len - 3] & 0x0f) << 4) + 0x08;
					}
				}
			}
		}
	}
	else if (rate == 2)
	{
		if (len % 2 == 1)
		{
			int j = 0;
			for (int i = 0; i < len0 - 1; i++)
			{
				if (i % 3 == 0)
				{
					din[i] = din1[j % len];
					j++;
				}
				else if (i % 3 == 1)
				{
					din[i] = (((din1[j % len] >> 4) & 0x0f) << 4) + 0x07;
				}
				else if (i % 3 == 2)
				{
					din[i] = (din1[j % len] & 0x0f) + (0x08 << 4);
					j++;
				}
				din[len0 - 1] = (din1[len - 2] & 0x0f) + (0x08 << 4);
			}
		}
		else
		{
			int j = 0;
			for (int i = 0; i < len0; i++)
			{

				if (i % 3 == 0)
				{
					din[i] = din1[j % len];
					j++;
				}
				else if (i % 3 == 1)
				{
					din[i] = (((din1[j % len] >> 4) & 0x0f) << 4) + 0x07;
				}
				else if (i % 3 == 2)
				{
					din[i] = (din1[j % len] & 0x0f) + (0x08 << 4);
					j++;
				}
			}

		}
	}
	int times = 0;
	int gogo;
	int tmp_din[128] = { 0 };
	int tmp_din1[64] = { 0 };
	int swit = 0;
	int addr = 0;
	int len1;
	int len2;
	int final_state = 0;
	int state = 0;
	int* tmp_d;

	tmp_d = (int*)malloc((len0 + 6) * sizeof(int));
	//状态传输si,so
	int* si;
	si = (int*)malloc(64 * sizeof(int));
	////对si初始化
	for (int i = 0; i < 64; i++)
	{
		si[i] = 0;
	}
	//**********************************//
	int** so;
	so = (int**)malloc(64 * sizeof(int*));
	for (int i = 0; i < 64; i++)
	{
		so[i] = (int*)malloc(2 * sizeof(int));
	}
	////对so初始化
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			so[i][j] = 0;
		}
	}
	//***********************************//
	//存储器
	int** mem;
	mem = (int**)malloc(128 * sizeof(int*));
	for (int i = 0; i < 128; i++)
	{
		if (mem == NULL)
		{
			printf("Error\n");
		}
		else
		{
			mem[i] = (int*)malloc(64 * sizeof(int));
		}
	}
	////对存储器初始化
	for (int i = 0; i < 128; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			mem[i][j] = 0;
		}

	}
	//*********************************************************//
	//第一次viterbi	
	for (int i = 0; i < len0 + 6; i++)
	{
		addr = i % len0;
		tmp_d[i] = din[addr];
	}
	addr++;

	viterbi(len0 + 6, tmp_d, swit, mem, times, si, so);
	//*****************************************************//
	//开始回溯的viterbi
	if (len0 <= 64) //一共只有两次viterbi
	{
		swit = 1;
		//第二次viterbi，回溯前64位
		for (int i = 1; i < 129; i++)
		{
			len1 = (addr - 1 + i) % len0;
			tmp_din[i - 1] = din[len1];
		}
		final_state = viterbi(128, tmp_din, swit, mem, times, si, so);
		//取mem中前len位回溯
		state = final_state;
		int state1 = 0;
		for (int i = 127; i >= 0; i--)
		{
			state1 = (mem[i % 128][state] == 0) ? (state >> 1) : (state >> 1) + 0x20;
			if (i < 64)
			{
				output[i] = mem[i % 128][state];
			}
			state = state1;
		}
	}
	//*******************************************************//
	else if (len0 > 64)//一共有times+2次
	{
		times = (len0 - 1) / 64;
		swit = 1;
		//第二次viterbi，整体存RAM128位长，回溯前64位
		for (int i = 1; i < 129; i++)
		{
			len1 = (addr - 1 + i) % len0;
			tmp_din[i - 1] = din[len1];
		}
		final_state = viterbi(128, tmp_din, swit, mem, 0, si, so);
		//取mem中前64位回溯
		state = final_state;
		int state1 = 0;
		for (int i = 127; i >= 0; i--)
		{
			state1 = (mem[i % 128][state] == 0) ? (state >> 0x01) : (state >> 0x01) + 0x20;
			if (i < 64)
			{
				output[i] = mem[i % 128][state];
			}
			state = state1;
		}
		//每次往后延64bit长	
		int go = 1;
		int type = 0;
		for (gogo = 1; gogo <= times; gogo++)
		{
			for (int i = 1; i < 65; i++)
			{
				len2 = (len1 + i) % len0;
				tmp_din1[i - 1] = din[len2];
			}
			len1 = len2;

			final_state = viterbi(64, tmp_din1, swit, mem, gogo, si, so);
			//对更新的64位回溯
			state = final_state;
			int state1 = 0;
			for (int i = 127; i >= 0; i--)
			{
				state1 = (mem[(i + 64 * gogo) % 128][state] == 0) ? (state >> 0x01) : (state >> 1) + 0x20;
				if (i < 64)
				{
					output[i + 64 * gogo] = mem[(i + 64 * gogo) % 128][state];
				}
				state = state1;
			}
		}
	}
	free(tmp_d);
	free(din);
	free(mem);
	free(so);
	free(si);
}