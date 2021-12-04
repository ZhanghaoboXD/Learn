#include "complex.h"

dcomplex cAdd(dcomplex* a, dcomplex* b)
{
	dcomplex c;
	c.x = a->x + b->x;
	c.y = a->y + b->y;
	return c;
}

dcomplex cSub(dcomplex* a, dcomplex* b)
{
	dcomplex c;
	c.x = a->x - b->x;
	c.y = a->y - b->y;
	return c;
}

dcomplex cMult(dcomplex* a, dcomplex* b)
{
	dcomplex c;
	c.x = a->x * b->x - a->y * b->y;
	c.y = a->x * b->y + a->y * b->x;
	return c;
}

dcomplex cDiv(dcomplex* a, int b)
{
	dcomplex c;
	c.x = a->x / b;
	c.y = a->y / b;
	return c;
}

dcomplex getConj(dcomplex* a)
{
	dcomplex c;
	c.x = a->x;
	c.y = -(a->y);
	return c;
}

dcomplex set(double x, double y)
{
	dcomplex c;
	c.x = x;
	c.y = y;
	return c;
}
double getPower(dcomplex* a)
{
	return a->x * a->x + a->y * a->y;
}