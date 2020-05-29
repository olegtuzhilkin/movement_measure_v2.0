//=========================================================
// src/movement_measure_main.c: generated by Hardware Configurator
//
// This file will be updated when saving a document.
// leave the sections inside the "$[...]" comment tags alone
// or they will be overwritten!!
//=========================================================

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <math.h>
#include <SI_EFM8LB1_Register_Enums.h>                  // SFR declarations
#include "InitDevice.h"
#include "movement_measure_main.h"
#include "correction.h"

#define SETTING = 1
#define DELAY 		40		//12us@48MHz
#define DELAY_2		10		//2us@72MHz
#define DELAY_3		1		//1us@48MHz
#define K			1.092
#define	B			-164
SI_SFR16 (DAC2_16,    0x8B);
// $[Generated Includes]
// [Generated Includes]$

//-----------------------------------------------------------------------------
// SiLabs_Startup() Routine
// ----------------------------------------------------------------------------
// This function is called immediately after reset, before the initialization
// code is run in SILABS_STARTUP.A51 (which runs before main() ). This is a
// useful place to disable the watchdog timer, which is enable by default
// and may trigger before main() in some instances.
//-----------------------------------------------------------------------------
void SiLabs_Startup(void) {
	// $[SiLabs Startup]
	// [SiLabs Startup]$
}

int __round(float x) {
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
//-----------------------------------------------------------------------------
// main() Routine
// ----------------------------------------------------------------------------
int main(void) {
	uint16_t U1, U2;
	int delta, sigma, txRel, result, i;
	float U1f, U2f;
	float rel;
	float K1, K2;
	int dac, shft;

	// Call hardware initialization routine
	enter_DefaultMode_from_RESET();

	SYNC = 1;

	linear_init();

	KEY = COIL1;
	U1 = 0;
	U2 = 0;
	U1f = 0.0;
	U2f = 0.0;

	K1 = (float) (FIL - 1) / FIL;
	K2 = (float) 1 / FIL;

	while (!SYNC)
		;
	while (SYNC)
		;
	while (1) {
// $[Generated Run-time code]
// [Generated Run-time code]$

		while (!SYNC)
			;

		ADC0CN0_ADINT = 0;
		ADC0CN0_ADBUSY = 1;
		KT_ADC = 1;
		TMR2CN0_TR2 = 1;

		while (!ADC0CN0_ADINT)
			;

		KT_ADC = 0;
		U1 = ADC0;
		U1f = U1f * K1 + U1 * K2;
		U1 = (uint16_t) __round(U1f);

		result = linear(rel);

		dac = K*result;//+B;
//		shft = dac>>8;

		while (SYNC)
			;

		while (!SYNC)
			;

		ADC0CN0_ADINT = 0;
		ADC0CN0_ADBUSY = 1;
		KT_ADC = 1;
		TMR2CN0_TR2 = 1;

		while (!ADC0CN0_ADINT)
			;

		KT_ADC = 0;
		U2 = ADC0;
		U2f = U2f * K1 + U2 * K2;
		U2 = (uint16_t) __round(U2f);

		delta = U1 - U2;
		sigma = U1 + U2;
		rel = (float) (delta) / (sigma) * K_REL + 3000;//8500;

		txRel = __round(rel);

#ifdef SETTING
		/*		U1 = 0;
		 U2 = 0x1555;
		 txRel = 0x2AAA;
		 result = 0x3FFF;*/
		setAddr(0x00);
		setData(U1);
		for (i = 0; i < DELAY_3; ++i)
			;
		READY = 1;
		for (i = 0; i < DELAY; ++i)
			;
		READY = 0;

		for (i = 0; i < DELAY_3; ++i)
			;

		setAddr(0x01);
		setData(U2);
		for (i = 0; i < DELAY_3; ++i)
			;
		READY = 1;
		for (i = 0; i < DELAY; ++i)
			;
		READY = 0;

		for (i = 0; i < DELAY_3; ++i)
			;

		setAddr(0x02);
		setData(txRel);
		for (i = 0; i < DELAY_3; ++i)
			;
		READY = 1;
		for (i = 0; i < DELAY; ++i)
			;
		READY = 0;

		for (i = 0; i < DELAY_3; ++i)
			;
#endif

		setAddr(0x03);
		setData(result);
		for (i = 0; i < DELAY_3; ++i)
			;
		READY = 1;
		for (i = 0; i < DELAY; ++i)
			;
		READY = 0;

		SFRPAGE = 0x30;
		dac = dac+B;
		DAC2L = dac;
		DAC2H = dac>>8;
		SFRPAGE = 0x00;

		while (SYNC)
			;

	}
}