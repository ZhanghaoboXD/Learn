#include <stdio.h>


int updateWhiting_data(int whiting_pn9)  //  uint16_t -> int
{
    int i, c; //uint8_t -> int

    for (i = 0; i < 8; i++)
    {
        c = whiting_pn9 & 0x21;

        whiting_pn9 >>= 1;    //ÓÒÒÆÒ»Î»  

        if ((c == 0x21) || (c == 0))
        {
            whiting_pn9 &= 0x0FF;
        }
        else
        {
            whiting_pn9 |= 0x100;
        }

    }
    return whiting_pn9;
}

void nrzGetPN9SeqUNB(double* ioBuff, int bitLen) 
{
    int whiting_pn9 = 0x1FF;  //uint16 -> int
    int index = 0;

    while (index < bitLen)
    {
        // every times process 8 nrz bit
        for (int i = 0; i < 8; i++)
        {
            ioBuff[index] = ioBuff[index] * (1 - 2 * ((whiting_pn9 >> (7 - i)) & 0x01));
            //printf("%f\n", ioBuff[index]);
            index++;
        }

        whiting_pn9 = updateWhiting_data(whiting_pn9);
    }
}