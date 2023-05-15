/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

// #include "params.h"
#include "wtime.h"
// #include "mtwister.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h> // time

/* params.h */
#define SHELLS 101              // discretization level
#define PHOTONS 32768           // 32K photons
#define MU_A 2.0f               // Absorption Coefficient in 1/cm !!non-zero!!
#define MU_S 20.0f              // Reduced Scattering Coefficient in 1/cm
#define MICRONS_PER_SHELL 50    // Thickness of spherical shells in microns
#define SEED (time(NULL))       // random seed

/* mtwister.h */
#define STATE_VECTOR_LENGTH 624
#define STATE_VECTOR_M      397 /* changes to STATE_VECTOR_LENGTH also require changes to this */

#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7fffffff

#define TEMPERING_MASK_B    0x9d2c5680
#define TEMPERING_MASK_C    0xefc60000

#define _XOPEN_SOURCE 500       // M_PI

#define albedo (MU_S / (MU_S + MU_A))
#define shells_per_mfp (1e4 / MICRONS_PER_SHELL / (MU_A + MU_S))

#define N_MIN_FOR 56
#define N_MAX_FOR 225

#define N_MIN_WHILE 2
#define N_MAX_WHILE 16

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";

// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];

typedef struct tagMTRand {
  unsigned long mt[STATE_VECTOR_LENGTH];
  int index;
} MTRand;

/* Photon */
// static void photon(MTRand * restrict rand) {}

/* Main matter */
int main(void) {
    // heading
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n#\n", PHOTONS);

    // configure RNG
    unsigned long rng;
    static unsigned long mag[2] = {0x0, 0x9908b0df};
    MTRand rand;
    rand.mt[0] = SEED & 0xffffffff;
    for(rand.index=1; rand.index<STATE_VECTOR_LENGTH; rand.index++) {
        rand.mt[rand.index] = (6069 * rand.mt[rand.index-1]) & 0xffffffff;
    }

    // start timer
    double start = wtime();

    // simulation
    for (unsigned int k = 0; k < PHOTONS; ++k) {
        // photon(&rand);

        /* Photon */

        // const float albedo = MU_S / (MU_S + MU_A);
        // const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

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

        // float xi1, xi2;
        float xi[2] = {0, 0};

        // for (;;) {
        for (int i = 0; i<N_MAX_FOR; ++i) {
    
            /* Step 2: Step size selection and photon packet movement */

            // Distance the photon packet travels between interaction sites

            // RNG M. Twister
            // mag[2] = {0x0, 0x9908b0df}; /* mag[x] = x * 0x9908b0df for x = 0,1 */
            mag[0] = 0x0;
            mag[1] = 0x9908b0df;
            if(rand.index >= STATE_VECTOR_LENGTH || rand.index < 0) {
                /* generate STATE_VECTOR_LENGTH words at a time */
                int kk;
                if(rand.index >= STATE_VECTOR_LENGTH+1 || rand.index < 0) {
                    // rand.mt[0] = 4357 & 0xffffffff;
                    rand.mt[0] = SEED & 0xffffffff;
                    for(rand.index=1; rand.index<STATE_VECTOR_LENGTH; rand.index++) {
                        rand.mt[rand.index] = (6069 * rand.mt[rand.index-1]) & 0xffffffff;
                    }
                }
                for(kk=0; kk<STATE_VECTOR_LENGTH-STATE_VECTOR_M; kk++) {
                    rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                    rand.mt[kk] = rand.mt[kk+STATE_VECTOR_M] ^ (rng >> 1) ^ mag[rng & 0x1];
                }
                for(; kk<STATE_VECTOR_LENGTH-1; kk++) {
                    rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                    rand.mt[kk] = rand.mt[kk+(STATE_VECTOR_M-STATE_VECTOR_LENGTH)] ^ (rng >> 1) ^ mag[rng & 0x1];
                }
                rng = (rand.mt[STATE_VECTOR_LENGTH-1] & UPPER_MASK) | (rand.mt[0] & LOWER_MASK);
                rand.mt[STATE_VECTOR_LENGTH-1] = rand.mt[STATE_VECTOR_M-1] ^ (rng >> 1) ^ mag[rng & 0x1];
                rand.index = 0;
            }
            rng = rand.mt[rand.index++];
            rng ^= (rng >> 11);
            rng ^= (rng << 7) & TEMPERING_MASK_B;
            rng ^= (rng << 15) & TEMPERING_MASK_C;
            rng ^= (rng >> 18);
            // (unsigned int)(rng % (RAND_MAX + 1LL))
            
            float t = -logf((unsigned int)(rng % (RAND_MAX + 1LL)) / (float)RAND_MAX);
           
            x += t * dir_x;
            y += t * dir_y;
            z += t * dir_z;

            /* Step 3: Absorption and scattering */

            unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;
            
            /*
                if (shell > SHELLS - 1) {
                    shell = SHELLS - 1;
                }
            */
            shell = ((SHELLS - 1) * (shell > SHELLS - 1)) + (shell * (shell <= SHELLS - 1));

            heat[shell] += (1.0f - albedo) * weight;
            heat2[shell] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */
           
            weight *= albedo;
            
            /* Step 4: Photon termination */

            // RNG M. Twister
            // mag[2] = {0x0, 0x9908b0df}; /* mag[x] = x * 0x9908b0df for x = 0,1 */
            mag[0] = 0x0;
            mag[1] = 0x9908b0df;
            if(rand.index >= STATE_VECTOR_LENGTH || rand.index < 0) {
                /* generate STATE_VECTOR_LENGTH words at a time */
                int kk;
                if(rand.index >= STATE_VECTOR_LENGTH+1 || rand.index < 0) {
                    // rand.mt[0] = 4357 & 0xffffffff;
                    rand.mt[0] = SEED & 0xffffffff;
                    for(rand.index=1; rand.index<STATE_VECTOR_LENGTH; rand.index++) {
                        rand.mt[rand.index] = (6069 * rand.mt[rand.index-1]) & 0xffffffff;
                    }
                }
                for(kk=0; kk<STATE_VECTOR_LENGTH-STATE_VECTOR_M; kk++) {
                    rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                    rand.mt[kk] = rand.mt[kk+STATE_VECTOR_M] ^ (rng >> 1) ^ mag[rng & 0x1];
                }
                for(; kk<STATE_VECTOR_LENGTH-1; kk++) {
                    rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                    rand.mt[kk] = rand.mt[kk+(STATE_VECTOR_M-STATE_VECTOR_LENGTH)] ^ (rng >> 1) ^ mag[rng & 0x1];
                }
                rng = (rand.mt[STATE_VECTOR_LENGTH-1] & UPPER_MASK) | (rand.mt[0] & LOWER_MASK);
                rand.mt[STATE_VECTOR_LENGTH-1] = rand.mt[STATE_VECTOR_M-1] ^ (rng >> 1) ^ mag[rng & 0x1];
                rand.index = 0;
            }
            rng = rand.mt[rand.index++];
            rng ^= (rng >> 11);
            rng ^= (rng << 7) & TEMPERING_MASK_B;
            rng ^= (rng << 15) & TEMPERING_MASK_C;
            rng ^= (rng >> 18);
            // (unsigned int)(rng % (RAND_MAX + 1LL))

            // Roulette
            // if (weight < 0.005f) {
            //     if (genRandInt(&rand) / (float)RAND_MAX > 0.1f) {
            //         break;
            //     };
            //     weight /= 0.1f;
            // }
            weight = (weight < 0.005f) ? ((((unsigned int)(rng % (RAND_MAX + 1LL)) / (float)RAND_MAX) > 0.1f) ? i = N_MAX_FOR : weight / 0.1f) : weight;


            /* New direction, rejection method */
            // do {
            for(int j = 0; (j < N_MAX_WHILE) && (1.0f < t); ++j) {
                // RNG M. Twister
                // mag[2] = {0x0, 0x9908b0df}; /* mag[x] = x * 0x9908b0df for x = 0,1 */
                mag[0] = 0x0;
                mag[1] = 0x9908b0df;
                if(rand.index >= STATE_VECTOR_LENGTH || rand.index < 0) {
                    /* generate STATE_VECTOR_LENGTH words at a time */
                    int kk;
                    if(rand.index >= STATE_VECTOR_LENGTH+1 || rand.index < 0) {
                        // rand.mt[0] = 4357 & 0xffffffff;
                        rand.mt[0] = SEED & 0xffffffff;
                        for(rand.index=1; rand.index<STATE_VECTOR_LENGTH; rand.index++) {
                            rand.mt[rand.index] = (6069 * rand.mt[rand.index-1]) & 0xffffffff;
                        }
                    }
                    for(kk=0; kk<STATE_VECTOR_LENGTH-STATE_VECTOR_M; kk++) {
                        rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                        rand.mt[kk] = rand.mt[kk+STATE_VECTOR_M] ^ (rng >> 1) ^ mag[rng & 0x1];
                    }
                    for(; kk<STATE_VECTOR_LENGTH-1; kk++) {
                        rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                        rand.mt[kk] = rand.mt[kk+(STATE_VECTOR_M-STATE_VECTOR_LENGTH)] ^ (rng >> 1) ^ mag[rng & 0x1];
                    }
                    rng = (rand.mt[STATE_VECTOR_LENGTH-1] & UPPER_MASK) | (rand.mt[0] & LOWER_MASK);
                    rand.mt[STATE_VECTOR_LENGTH-1] = rand.mt[STATE_VECTOR_M-1] ^ (rng >> 1) ^ mag[rng & 0x1];
                    rand.index = 0;
                }   
                rng = rand.mt[rand.index++];
                rng ^= (rng >> 11);
                rng ^= (rng << 7) & TEMPERING_MASK_B;
                rng ^= (rng << 15) & TEMPERING_MASK_C;
                rng ^= (rng >> 18);
                // (unsigned int)(rng % (RAND_MAX + 1LL))
                xi[0] = 2.0f * (unsigned int)(rng % (RAND_MAX + 1LL)) / (float)RAND_MAX - 1.0f;
                // RNG M. Twister
                // mag[2] = {0x0, 0x9908b0df}; /* mag[x] = x * 0x9908b0df for x = 0,1 */
                mag[0] = 0x0;
                mag[1] = 0x9908b0df;
                if(rand.index >= STATE_VECTOR_LENGTH || rand.index < 0) {
                    /* generate STATE_VECTOR_LENGTH words at a time */
                    int kk;
                    if(rand.index >= STATE_VECTOR_LENGTH+1 || rand.index < 0) {
                        // rand.mt[0] = 4357 & 0xffffffff;
                        rand.mt[0] = SEED & 0xffffffff;
                        for(rand.index=1; rand.index<STATE_VECTOR_LENGTH; rand.index++) {
                            rand.mt[rand.index] = (6069 * rand.mt[rand.index-1]) & 0xffffffff;
                        }
                    }
                    for(kk=0; kk<STATE_VECTOR_LENGTH-STATE_VECTOR_M; kk++) {
                        rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                        rand.mt[kk] = rand.mt[kk+STATE_VECTOR_M] ^ (rng >> 1) ^ mag[rng & 0x1];
                    }
                    for(; kk<STATE_VECTOR_LENGTH-1; kk++) {
                        rng = (rand.mt[kk] & UPPER_MASK) | (rand.mt[kk+1] & LOWER_MASK);
                        rand.mt[kk] = rand.mt[kk+(STATE_VECTOR_M-STATE_VECTOR_LENGTH)] ^ (rng >> 1) ^ mag[rng & 0x1];
                    }
                    rng = (rand.mt[STATE_VECTOR_LENGTH-1] & UPPER_MASK) | (rand.mt[0] & LOWER_MASK);
                    rand.mt[STATE_VECTOR_LENGTH-1] = rand.mt[STATE_VECTOR_M-1] ^ (rng >> 1) ^ mag[rng & 0x1];
                    rand.index = 0;
                }
                rng = rand.mt[rand.index++];
                rng ^= (rng >> 11);
                rng ^= (rng << 7) & TEMPERING_MASK_B;
                rng ^= (rng << 15) & TEMPERING_MASK_C;
                rng ^= (rng >> 18);
                // (unsigned int)(rng % (RAND_MAX + 1LL))
                xi[1] = 2.0f * (unsigned int)(rng % (RAND_MAX + 1LL)) / (float)RAND_MAX - 1.0f;
                t = xi[0] * xi[0] + xi[1] * xi[1];
            }; // while (1.0f < t);   

            dir_x = 2.0f * t - 1.0f;
            dir_y = xi[0] * sqrtf((1.0f - dir_x * dir_x) / t);
            dir_z = xi[0] * sqrtf((1.0f - dir_x * dir_x) / t);
        };
    };

    // stop timer
    double end = wtime();

    assert(start <= end);
    double elapsed = end - start;

    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);

    printf("# Radius\tHeat\n");
    printf("# [microns]\t[W/cm^3]\tError\n");

    float t = 4.0f * M_PI * powf(MICRONS_PER_SHELL, 3.0f) * PHOTONS / 1e12;

    for (unsigned int i = 0; i < SHELLS - 1; ++i) {
        printf("%6.0f\t%12.5f\t%12.5f\n",
            i * (float)MICRONS_PER_SHELL,
            heat[i] / t / (i * i + i + 1.0 / 3.0),
            sqrt(heat2[i] - heat[i] * heat[i] / PHOTONS) / t / (i * i + i + 1.0f / 3.0f)
        );
    }

    printf("# extra\t%12.5f\n", heat[SHELLS - 1] / PHOTONS);

    return 0;
}
    