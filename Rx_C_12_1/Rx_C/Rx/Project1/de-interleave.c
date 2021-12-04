#include <stdio.h>

#define UL_INTLV_INNER_SIZE 64

void interleaveProc(double* ccOutNrzByte, unsigned int ccOutBitLen)  //uint32 -> unsigned int  ����ʱ��int
{
    static const int UL_INTLV_TRANS[UL_INTLV_INNER_SIZE] =
    {
        0, 8, 16, 24, 32, 40, 48, 56,
        1, 9, 17, 25, 33, 41, 49, 57,
        2, 10, 18, 26, 34, 42, 50, 58,
        3, 11, 19, 27, 35, 43, 51, 59,
        4, 12, 20, 28, 36, 44, 52, 60,
        5, 13, 21, 29, 37, 45, 53, 61,
        6, 14, 22, 30, 38, 46, 54, 62,
        7, 15, 23, 31, 39, 47, 55, 63
    };
    double tmpIntlv[UL_INTLV_INNER_SIZE];
    const int intlvFloorLen = ccOutBitLen / UL_INTLV_INNER_SIZE;   //һ�ν�֯��Ҫ64bit���ж�Ҫ������֯��

    for (int idx = 0; idx < intlvFloorLen; idx++)   //ÿ��
    {
        for (int bit = 0; bit < UL_INTLV_INNER_SIZE; bit++)   //��ÿ�齻֯ǰ����ŵ� tmpIntlv��
        {
            tmpIntlv[bit] = ccOutNrzByte[idx * UL_INTLV_INNER_SIZE + bit];
        }

        for (int bit = 0; bit < UL_INTLV_INNER_SIZE; bit++)   //��ÿ��tmpIntlv����֯˳��Ż�������У��齻֯���
        {
            ccOutNrzByte[idx * UL_INTLV_INNER_SIZE + bit] = tmpIntlv[UL_INTLV_TRANS[bit]];
            //printf("%f\n", ccOutNrzByte[idx * UL_INTLV_INNER_SIZE + bit]);
        }
    }
}