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

static inline uint64_t rotl64(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

uint64_t next64(void) {
	const uint64_t result = rotl64(s64[0] + s64[3], 23) + s64[0];

	const uint64_t t = s64[1] << 17;

	s64[2] ^= s64[0];
	s64[3] ^= s64[1];
	s64[1] ^= s64[2];
	s64[0] ^= s64[3];

	s64[2] ^= t;

	s64[3] = rotl64(s64[3], 45);

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void jump64(void) {
	static const uint64_t JUMP64[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(uint64_t i = 0; i < sizeof JUMP64 / sizeof *JUMP64; i++)
		for(uint64_t b = 0; b < 64; b++) {
			if (JUMP64[i] & UINT64_C(1) << b) {
				s0 ^= s64[0];
				s1 ^= s64[1];
				s2 ^= s64[2];
				s3 ^= s64[3];
			}
			next64();
		}
		
	s64[0] = s0;
	s64[1] = s1;
	s64[2] = s2;
	s64[3] = s3;
}



/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */


void long_jump64(void) {
	static const uint64_t LONG_JUMP64[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;

	for(uint64_t i = 0; i < sizeof LONG_JUMP64 / sizeof *LONG_JUMP64; i++)
		for(uint64_t b = 0; b < 64; b++) {
			if (LONG_JUMP64[i] & UINT64_C(1) << b) {
				s0 ^= s64[0];
				s1 ^= s64[1];
				s2 ^= s64[2];
				s3 ^= s64[3];
			}
			next64();	
		}
		
	s64[0] = s0;
	s64[1] = s1;
	s64[2] = s2;
	s64[3] = s3;
}

int main() {

	printf("RNG Xoshiro 256 ++ 64 bits TEST\n");
    
	// Semilla inicial (puedes cambiarla para obtener diferentes secuencias de números aleatorios)

	s64[0] = (uint64_t)time(NULL); // 0x0123456789ABCDEF;
	s64[1] = (uint64_t)time(NULL); // 0xABCDEF0123456789;
	s64[2] = (uint64_t)time(NULL); // 0x13579BDF2468ACE0;
	s64[3] = (uint64_t)time(NULL); // 0xE0AC4682DF9B7531;

	printf("MAX64: %lu\n", MAX64);

	// Generar número aleatorio

	// uint64_t random_number64 = next64();
	// Usa el número aleatorio aquí como desees
	// printf("Number: %lx\t(size: %ld bytes)\n", random_number64, sizeof(random_number64));

	for(unsigned int i = 0; i < 10; ++i) {

		uint64_t num = next64();
		jump64();
		// long_jump64();
		printf("RNG num: %d\tnum/MAX64: %5.5f\n", num, num/(double)MAX64);

		// printf("s64[0]: %lx\ts64[1]: %lx\ts64[2]: %lx\ts64[3]: %lx\n", s64[0], s64[1], s64[2], s64[3]);
		
		// printf("s64[0]: %lx\tMAX64: %lx\ts64[0]/MAX64: %.5f\n", s64[0], MAX64, (double)s64[0]/MAX64, s64[0]);
		// printf("s64[1]: %lx\tMAX64: %lx\ts64[1]/MAX64: %.5f\n", s64[1], MAX64, (double)s64[1]/MAX64, s64[0]);
		// printf("s64[2]: %lx\tMAX64: %lx\ts64[2]/MAX64: %.5f\n", s64[2], MAX64, (double)s64[2]/MAX64, s64[0]);
		// printf("s64[3]: %lx\tMAX64: %lx\ts64[3]/MAX64: %.5f\n", s64[3], MAX64, (double)s64[3]/MAX64, s64[0]);

	};

	return 0;
}