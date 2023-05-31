/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "xoshiro256plusplus.h"

/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
   It has excellent (sub-ns) speed, a state (256 bits) that is large
   enough for any parallel application, and it passes all tests we are
   aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static inline uint32_t rotl32(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}


uint32_t next32(void) {

	const uint32_t result1 = rotl32(s32[0] + s32[6], 11) + s32[0];
	const uint32_t result2 = rotl32(s32[1] + s32[7], 13) + s32[1];

	const uint32_t t1 = s32[2] << 7;
	const uint32_t t2 = s32[3] << 9;

	s32[4] ^= s32[0];
	s32[5] ^= s32[1];
	s32[6] ^= s32[2];
	s32[7] ^= s32[3];
	s32[3] ^= s32[4];
	s32[2] ^= s32[5];
	s32[1] ^= s32[6];
	s32[0] ^= s32[7];

	s32[4] ^= t1;
	s32[5] ^= t2;

	s32[6] = rotl32(s32[6], 19);
	s32[7] = rotl32(s32[7], 23);

	return result1+result2;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void jump32(void) {
	static const uint32_t JUMP32[] = { 0x180ec6d3, 0x3cfd0aba, 0xd5a61266, 0xf0c9392c, 0xa9582618, 0xe03fc9aa, 0x39abdc45, 0x29b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	uint64_t s4 = 0;
	uint64_t s5 = 0;
	uint64_t s6 = 0;
	uint64_t s7 = 0;

	for(uint32_t i = 0; i < sizeof JUMP32 / sizeof *JUMP32; i++) {
		for(uint32_t b = 0; b < 32; b++) {
			if (JUMP32[i] & UINT32_C(1) << b) {
				s0 ^= s32[0];
				s1 ^= s32[1];
				s2 ^= s32[2];
				s3 ^= s32[3];
				s4 ^= s32[4];
				s5 ^= s32[5];
				s6 ^= s32[6];
				s7 ^= s32[7];
			}
			next32();	
		}
	}
		
	s32[0] = s0;
	s32[1] = s1;
	s32[2] = s2;
	s32[3] = s3;
	s32[4] = s4;
	s32[5] = s5;
	s32[6] = s6;
	s32[7] = s7;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump32(void) {
	static const uint32_t LONG_JUMP32[] = { 0x76e15d3e, 0xfefdcbbf, 0xc5004e44, 0x1c522fb3, 0x77710069, 0x854ee241, 0x39109bb0, 0x2acbe635 };

	uint32_t s0 = 0;
	uint32_t s1 = 0;
	uint32_t s2 = 0;
	uint32_t s3 = 0;
	uint32_t s4 = 0;
	uint32_t s5 = 0;
	uint32_t s6 = 0;
	uint32_t s7 = 0;
	
	for(uint32_t i = 0; i < sizeof LONG_JUMP32 / sizeof *LONG_JUMP32; i++)
		for(uint32_t b = 0; b < 32; b++) {
			if (LONG_JUMP32[i] & UINT32_C(1) << b) {
				s0 ^= s32[0];
				s1 ^= s32[1];
				s2 ^= s32[2];
				s3 ^= s32[3];
				s4 ^= s32[4];
				s5 ^= s32[5];
				s6 ^= s32[6];
				s7 ^= s32[7];
			}
			next32();	
		}
		
	s32[0] = s0;
	s32[1] = s1;
	s32[2] = s2;
	s32[3] = s3;
	s32[4] = s4;
	s32[5] = s5;
	s32[6] = s6;
	s32[7] = s7;
}

/*
int main() {

	printf("RNG Xoshiro 256 ++ 32 bits TEST\n");
    
	// Semilla inicial (puedes cambiarla para obtener diferentes secuencias de números aleatorios)

	s32[0] = (uint32_t)time(NULL); // 0x01234567;
	s32[4] = (uint32_t)time(NULL); // 0x89ABCDEF;
	s32[1] = (uint32_t)time(NULL); // 0xABCDEF01;
	s32[5] = (uint32_t)time(NULL); // 0x23456789;
	s32[2] = (uint32_t)time(NULL); // 0x13579BDF;
	s32[6] = (uint32_t)time(NULL); // 0x2468ACE0;
	s32[3] = (uint32_t)time(NULL); // 0xE0AC4682;
	s32[7] = (uint32_t)time(NULL); // 0xDF9B7531;

	printf("MAX32: %u\n", MAX32);

	// Generar número aleatorio

	for(unsigned int i = 0; i < 10; ++i) {

		uint32_t num = next32();
		jump32();
		// long_jump32();
		printf("RNG num: %d\tnum/MAX32: %5.5f\n", num, num/(float)MAX32);
		// printf("s32[0]: %x\ts32[1]: %x\ts32[2]: %x\ts32[3]: %x\ts32[4]: %x\ts32[5]: %x\ts32[6]: %x\ts32[7]: %x\n", s32[0], s32[1], s32[2], s32[3], s32[4], s32[5], s32[6], s32[7]);

		// printf("s32[0]: %x\tMAX32: %x\ts32[0]/MAX32: %.5f\n", s32[0], MAX32, (float)s32[0]/MAX32);
		// printf("s32[1]: %x\tMAX32: %x\ts32[1]/MAX32: %.5f\n", s32[1], MAX32, (float)s32[1]/MAX32);
		// printf("s32[2]: %x\tMAX32: %x\ts32[2]/MAX32: %.5f\n", s32[2], MAX32, (float)s32[2]/MAX32);
		// printf("s32[3]: %x\tMAX32: %x\ts32[3]/MAX32: %.5f\n", s32[3], MAX32, (float)s32[3]/MAX32);
		// printf("s32[4]: %x\tMAX32: %x\ts32[4]/MAX32: %.5f\n", s32[4], MAX32, (float)s32[4]/MAX32);
		// printf("s32[5]: %x\tMAX32: %x\ts32[5]/MAX32: %.5f\n", s32[5], MAX32, (float)s32[5]/MAX32);
		// printf("s32[6]: %x\tMAX32: %x\ts32[6]/MAX32: %.5f\n", s32[6], MAX32, (float)s32[6]/MAX32);
		// printf("s32[7]: %x\tMAX32: %x\ts32[7]/MAX32: %.5f\n", s32[7], MAX32, (float)s32[7]/MAX32);

	};

	return 0;
}
*/