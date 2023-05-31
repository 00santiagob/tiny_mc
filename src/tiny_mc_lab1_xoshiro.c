/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#define _XOPEN_SOURCE 500  // M_PI

#define MAX32 0xFFFFFFFF

#include "params.h"
#include "wtime.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
   It has excellent (sub-ns) speed, a state (256 bits) that is large
   enough for any parallel application, and it passes all tests we are
   aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static uint32_t s32[8];

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


char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";


// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];


/***
 * Photon
 ***/

static void photon()
{
    const float albedo = MU_S / (MU_S + MU_A);
    const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

    /* STEP 1: Launching a photon packet */
    // Initial position
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    // Initial direction of propagation
    float dir_x = 0.0f;
    float dir_y = 0.0f;
    float dir_z = 1.0f;
    // Initial weight of photon
    float weight = 1.0f;

    for (;;) {
        /* Step 2: Step size selection and photon packet movement */
        // Distance the photon packet travels between interaction sites
        float t = -logf(next32()/(float)MAX32); /* move */
       
        /* Roulette */
		if (weight < 0.005f) { /* roulette */ 
            if (next32()/(float)MAX32 > 0.1f)
                break;
            weight /= 0.1f;
        }
	
		x += t * dir_x;
        y += t * dir_y;
        z += t * dir_z;

        /* Step 3: Absorption and scattering */
        unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;

        float xi1, xi2;
		      
       	if (shell > SHELLS - 1) {
            shell = SHELLS - 1;
        }

        /* New direction, rejection method */
        /* float xi1, xi2; // La declaración de las variables se movio hacia arriba*/
        do {
            xi1 = 2.0f * next32()/(float)MAX32 - 1.0f;
            xi2 = 2.0f * next32()/(float)MAX32 - 1.0f;
            t = xi1 * xi1 + xi2 * xi2;
        } while (1.0f < t);
        

			     
	
        heat[shell] += (1.0f - albedo) * weight;
        heat2[shell] += (1.0f - albedo) * (1 - albedo) * weight * weight; /* add up squares */
       
       
        dir_x = 2.0f * t - 1.0f;
        weight *= albedo; 
		dir_y = xi1 * sqrtf((1.0f - dir_x * dir_x) / t);
		dir_z = xi2 * sqrtf((1.0f - dir_x * dir_x) / t);

        /* roulette: Se agrando el valor en la condicional de 0.001 a 0.005 */
        // if (weight < 0.005f) {
        //     if (genRandInt(&rand) / (float)RAND_MAX > 0.1f) {
        //         break;
        //     };
        //     weight /= 0.1f;
        // }
	}
}


/***
 * Main matter
 ***/

int main(void)
{
    // heading
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n#\n", PHOTONS);

    // configure RNG Xoshiro 256 ++ para 32 bits
    s32[0] = (uint32_t)SEED;
    s32[4] = (uint32_t)SEED;
    s32[1] = (uint32_t)SEED;
    s32[5] = (uint32_t)SEED;
    s32[2] = (uint32_t)SEED;
    s32[6] = (uint32_t)SEED;
    s32[3] = (uint32_t)SEED;
    s32[7] = (uint32_t)SEED;

    // start timer
    double start = wtime();
    
    // simulation
    for (unsigned int i = 0; i < PHOTONS; ++i) {
        photon();
    }
    // stop timer
    double end = wtime();
    assert(start <= end);
    double elapsed = end - start;

    int len = sizeof(heat) / sizeof(heat[0]); 
    FILE *fp;
    fp = fopen("dati_mod_1.bin", "wb");
    fwrite(heat, sizeof(float), len, fp);
    fclose(fp);

    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);

    printf("# Radius\tHeat\n");
    printf("# [microns]\t[W/cm^3]\tError\n");
    float t = 4.0f * M_PI * powf(MICRONS_PER_SHELL, 3.0f) * PHOTONS / 1e12;
    for (unsigned int i = 0; i < SHELLS - 1; ++i) {
        printf("%6.0f\t%12.5f\t%12.5f\n", i * (float)MICRONS_PER_SHELL,
               heat[i] / t / (i * i + i + 1.0 / 3.0),
               sqrt(heat2[i] - heat[i] * heat[i] / PHOTONS) / t / (i * i + i + 1.0f / 3.0f));
    }
    printf("# extra\t%12.5f\n", heat[SHELLS - 1] / PHOTONS);

    return 0;
}