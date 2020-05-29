/*
 * correction.c
 *
 *  Created on: 13 июня 2019 г.
 *      Author: Олег
 */
#include <math.h>
#include "correction.h"
#include "const.h"

float xdata Lin_A1[N-1], Lin_B1[N-1], Lin_C1[N-1], Lin_D1[N-1];
float xdata K1_low, B1_low, K1_high, B1_high;
int xdata aa[N], bb[N], cc[N];
float xdata B_[N],  res_normal[N];

int _round(float x) {
	float r;
	int y;
	y = (int) x;
	r = x - y;
	r = fabs(r);
	if (r > 0.5)
		if (x > 0)
			y += 1;
		else
			y -= 1;
	return y;
}
//**************************************************************
//--------------------------------------------------------------
//**************************************************************
void progon(void)
{
    //straight run
    float xdata alph[N], bett[N];

    int i;
    alph[0] = -cc[0]/bb[0];
    bett[0] = B_[0]/bb[0];
    for (i = 1; i < N-1; i++){
        alph[i]=-cc[i]/(bb[i]+aa[i]*alph[i-1]);
        bett[i]=(B_[i]-aa[i]*bett[i-1])/(bb[i]+aa[i]*alph[i-1]);
    }
    alph[N-1]=0;
    bett[N-1]=(B_[N-1]-aa[N-1]*bett[N-2])/(bb[N-1]+aa[N-1]*alph[N-2]);
    //back run
    res_normal[N-1]=bett[N-1];
    for (i = N-2; i >= 0; i--)
    	res_normal[i]=alph[i]*res_normal[i+1]+bett[i];
    return;
}

//**************************************************************
//-------------------linear coeffs init-------------------------
//**************************************************************
void linear_init(void)
{
	int xdata h[N];
	int i;
	//---------------------------------------
	K1_low = (float)(i_grad[1] - i_grad[0])/(grad_1[1] - grad_1[0]);
	B1_low = i_grad[0] - K1_low*grad_1[0];

	K1_high = (float)(i_grad[N-1] - i_grad[N-2])/(grad_1[N-1] - grad_1[N-2]);
	B1_high = i_grad[N-1] - K1_high*grad_1[N-1];
	//---------------------------------------
	h[0] = 0;
	for(i = 1; i < N; i++)
        h[i]=grad_1[i]-grad_1[i-1];

	B_[0] = 0;
	B_[N-1] = 0;
	for(i = 1; i < N-1; i++)
			B_[i]=3*((float)(i_grad[i+1]-i_grad[i])/h[i+1]-(float)(i_grad[i]-i_grad[i-1])/h[i]);
	//main diagonaly
	bb[0] = 1;
	bb[N-1] = 1;
	for(i = 1; i < N-1; i++)
			bb[i]=2*(h[i]+h[i+1]);
	//upper and down diagonaly
	for(i = 0; i < N-1; i++)
			aa[i]=h[i];
	aa[N-1] = 0;
	cc[0] = 0;
	for(i = 0; i < N-2; i++)
			cc[i+1]=h[i+2];
	cc[N-1] = 0;
	progon();
	//getting spline coeefs
	for(i = 0; i < N-1; i++)
			Lin_C1[i]=res_normal[i];
	for(i = 0; i < N-1; i++){
			Lin_D1[i]=(Lin_C1[i+1]-Lin_C1[i])/(3*h[i+1]);
			Lin_A1[i]=i_grad[i];
	}
	Lin_D1[N-2]=-Lin_C1[N-2]/(3*h[N-1]);
	Lin_A1[N-2]=i_grad[N-2];
	for(i = 1; i < N-1; i++)
			Lin_B1[i-1]=(float)(i_grad[i]-i_grad[i-1])/h[i]-h[i]*(Lin_C1[i]+2*Lin_C1[i-1])/3;
	Lin_B1[N-2]=(float)(i_grad[N-1]-i_grad[N-2])/h[N-1]-2*h[N-1]*Lin_C1[N-2]/3;

	return;
}
//**************************************************************
//-------------------linear get result--------------------------
//**************************************************************
uint16_t linear(float XU)
{
//рассчётная часть
	float T;
	int res;
	uint16_t pos = 0;
	uint16_t i;
	//----------------------------------
    int up = N - 1,
		down = 0;
	//----------------------------------

	if (XU < grad_1[0]){
		res = _round(K1_low*XU+B1_low);
		goto end;
	}
	else
    if (XU > grad_1[N-1]){
			res = _round(K1_high*XU+B1_high);
			goto end;
		}

/*	pos = N-2;
	for (i = 1; i < N-1; i++)
		if (XU < grad_1[i]){
			pos = i-1;
			break;
		}*/

    while (down != (up-1)){
		pos = (up - down)/2 + down;

		if (XU >= grad_1[pos])
			down = pos;
		else //if (point <= grad[pos])
			up = pos;
	}

    T = (XU - grad_1[down]);
		res = _round((T*(T*(Lin_D1[down]*T+Lin_C1[down])+Lin_B1[down])+Lin_A1[down]));

end:if (res < 180) return 180;
		else if (res > i_grad[N-1]) return i_grad[N-1];
				 else return (uint16_t)res;
}
