/*
 * const.h
 *
 *  Created on: 13 ���� 2019 �.
 *      Author: ����
 */

#ifndef INC_CONST_H_
#define INC_CONST_H_

#define SETTING = 1

#define NUMB_POINT  	13
#define N 				NUMB_POINT

//						   6    5     4	    3	  2     1	  0	   -1	 -2	   -3	  -4     -5	    -6
int code i_grad[N]		= {200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600};//{0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000};

int code grad_1[N]	_at_ 0x3C00;//	= {1381, 1585, 1826, 2095, 2386, 2690, 3000, 3308, 3608, 3895, 4162, 4403, 4606};//{1578, 2411, 3464, 4614, 5863, 7164,  8500, 9835, 11130, 12364, 13519, 14542, 15402};

#endif /* INC_CONST_H_ */
