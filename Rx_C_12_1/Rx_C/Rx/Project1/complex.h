#ifndef COMPLEX_H
#define COMPLEX_H

#ifdef _MSC_VER
    #ifndef _M_X64
        #define USE_F87_ASSEMBLY
    #endif
#endif

typedef struct DCOMPLEX
{
    double x;
    double y;

} dcomplex;

dcomplex cAdd(dcomplex* a, dcomplex* b); //�ӷ�
dcomplex cSub(dcomplex* a, dcomplex* b);//����
dcomplex cMult(dcomplex* a, dcomplex* b);//�˷�
dcomplex cDiv(dcomplex* a, int b);//����
dcomplex set(double x, double y);
dcomplex getConj(dcomplex* a);
double getPower(dcomplex* a);

#endif
