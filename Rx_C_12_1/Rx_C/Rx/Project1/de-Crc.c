#include <stdio.h>
int getCrc8Buf(const int* inputBit, int rawByteLen)
{
    const int CRC8_MASK_UNB = 0x07;
    int buf = 0x0;
    for (int i = 0; i < rawByteLen; i++)
    {
        buf ^= inputBit[i];
        for (int j = 0; j < 8; j++)
        {
            if ((buf & 0x80u) > 0)
            {
                buf = (buf << 1u) ^ CRC8_MASK_UNB;
            }
            else
            {
                buf <<= 1u;
            }
        }
    }
    // buf ^= 0x55;

    return buf;
}

int decCrc8_UNB(const int* inputBit, int rawByteLen)
{
    int buf = getCrc8Buf(inputBit, rawByteLen);

    if (inputBit[rawByteLen] == (buf & 0xFF))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int getCrc16Buf(const int* inputBit, int rawByteLen)
{
    const int CRC16_MASK_UNB = 0x1021;
    int buf = 0xFFFF;
    for (int i = 0; i < rawByteLen; i++)
    {
        buf ^= (inputBit[i] << 8);
        for (int j = 0; j < 8; j++)
        {
            if ((buf & 0x8000u) > 0)
            {
                buf = (buf << 1u) ^ CRC16_MASK_UNB;
            }
            else
            {
                buf <<= 1u;
            }
        }
    }

    return buf;
}

int decCrc16_UNB(const int* inputBit, int rawByteLen)
{
    int buf = getCrc16Buf(inputBit, rawByteLen);

    if (inputBit[rawByteLen] == (buf & 0xFF)
        && inputBit[rawByteLen + 1] == ((buf >> 8) & 0xFF))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int getCrc24Buf(const int* inputBit, int rawByteLen)
{
    const int CRC24_MASK_UNB = 0x800063;
    int buf = 0;
    for (int i = 0; i < rawByteLen; i++)
    {
        buf ^= (inputBit[i] << 16);
        for (int j = 0; j < 8; j++)
        {
            if ((buf & 0x800000u) > 0)
            {
                buf = (buf << 1u) ^ CRC24_MASK_UNB;
            }
            else
            {
                buf <<= 1u;
            }
        }
    }
    return buf;
}

int decCrc24_UNB(const int* inputBit, int rawByteLen)
{
    int buf = getCrc24Buf(inputBit, rawByteLen);

    if (inputBit[rawByteLen] == (buf & 0xFF)
        && inputBit[rawByteLen + 1] == ((buf >> 8) & 0xFF)
        && inputBit[rawByteLen + 2] == ((buf >> 16) & 0xFF))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}