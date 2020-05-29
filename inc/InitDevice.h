//=========================================================
// inc/InitDevice.h: generated by Hardware Configurator
//
// This file will be regenerated when saving a document.
// leave the sections inside the "$[...]" comment tags alone
// or they will be overwritten!!
//=========================================================
#ifndef __INIT_DEVICE_H__
#define __INIT_DEVICE_H__

// USER CONSTANTS
sbit SYNC = P0 ^ 4;
sbit KEY = P0 ^ 2;
sbit READY = P3 ^ 1;
sbit KT_ADC = P3 ^ 0;

sbit ADDR0 = P0 ^ 5;
sbit ADDR1 = P0 ^ 6;
sbit ADDR2 = P0 ^ 7;

#define setAddr(addr) P0 &= 0x1F; P0 |= ((addr&0x07) << 5);
#define setData(data) P1 = (data&0xFF); P2 &= 0x40; P2 |= ((data&0x3F00) >> 8);
// USER PROTOTYPES

// $[Mode Transition Prototypes]
extern void enter_DefaultMode_from_RESET(void);
// [Mode Transition Prototypes]$

// $[Config(Per-Module Mode)Transition Prototypes]
extern void WDT_0_enter_DefaultMode_from_RESET(void);
extern void PORTS_0_enter_DefaultMode_from_RESET(void);
extern void PORTS_1_enter_DefaultMode_from_RESET(void);
extern void PORTS_2_enter_DefaultMode_from_RESET(void);
extern void PORTS_3_enter_DefaultMode_from_RESET(void);
extern void PBCFG_0_enter_DefaultMode_from_RESET(void);
extern void ADC_0_enter_DefaultMode_from_RESET(void);
extern void DAC_2_enter_DefaultMode_from_RESET(void);
extern void DACGCF_0_enter_DefaultMode_from_RESET(void);
extern void VREF_0_enter_DefaultMode_from_RESET(void);
extern void EXTOSC_0_enter_DefaultMode_from_RESET(void);
extern void CIP51_0_enter_DefaultMode_from_RESET(void);
extern void CLOCK_0_enter_DefaultMode_from_RESET(void);
extern void TIMER16_2_enter_DefaultMode_from_RESET(void);
extern void TIMER_SETUP_0_enter_DefaultMode_from_RESET(void);
extern void INTERRUPT_0_enter_DefaultMode_from_RESET(void);
// [Config(Per-Module Mode)Transition Prototypes]$

#endif

