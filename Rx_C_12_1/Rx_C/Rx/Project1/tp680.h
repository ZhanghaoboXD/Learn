
typedef struct
{
    int sys_symBps;
    int sys_sampRate;

    int shr_N_zc;    // 16, 64, 128
    int shr_freq_bias;

    int phr_N_cs;

    int n_sf;
    int modType;     // 0: MSK,  1: GMSK,   2: FSK,  3: GFSK,   4:BPSK,    5: QPSK,  6: 8PSK
    int codeType;    // 0: CC 1/2,       1: CC 2/3         2: CC 3/4          3: Turbo
    int ccPreLen;

    int turboNtti;
    int turboSizeIdx;

    // rx Para
    int rx_L_numOfGroup;
    int rx_I_zc;     // bit len of ZC in each group

    int rx_NumBit;
}Tp680xPara;
Tp680xPara g_tp680Para;
