/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#define _XOPEN_SOURCE 500  // M_PI

#include "params.h"
#include "wtime.h"
#include "mtwister.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


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

static void photon(MTRand *rand)
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
        float t = -logf(genRngMTInt(rand) / (float)RAND_MAX); /* move */
       
        /* Roulette */
		if (weight < 0.005f) { /* roulette */ 
            if (genRngMTInt(rand) / (float)RAND_MAX > 0.1f)
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
            xi1 = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
            xi2 = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
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

    // configure RNG
    MTRand rand = seedRand(SEED);

    // start timer
    double start = wtime();
    
    // simulation
    for (unsigned int i = 0; i < PHOTONS; ++i) {
        photon(&rand);
    }
    // stop timer
    double end = wtime();
    assert(start <= end);
    double elapsed = end - start;

    /*  */
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
