#ifndef __XOSHIRO256PLUSPLUS_H

#define __XOSHIRO256PLUSPLUS_H

#define MAX32 0xFFFFFFFF
#define MAX64 0xFFFFFFFFFFFFFFFF

#include <stdint.h>

static uint32_t s32[8];
uint32_t next32(void);
void jump32(void);
void long_jump32(void);

//static uint64_t s64[4];
//uint64_t next64(void);
//void jump64(void);
//void long_jump64(void);

#endif