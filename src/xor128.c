/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#define _XOPEN_SOURCE 500  // M_PI

#include "params.h"
#include "wtime.h"
#include "simdxorshift128plus.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <immintrin.h>

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";


// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];
//static float merged_heat[SHELLS*2];
int random_n[4];

/***
 * Photon
 ***/

static void photon(avx_xorshift128plus_key_t *my_key)
{
    const float albedo = MU_S / (MU_S + MU_A);
    const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

    /* launch */
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float w = 1.0f;
    float weight = 1.0f;
    //bool found = false;
    for (;;) {
    __m256i vect = avx_xorshift128plus(my_key);
   
     random_n[0] = _mm256_extract_epi32(vect, 0);
     random_n[1] = _mm256_extract_epi32(vect, 1);
     random_n[2] = _mm256_extract_epi32(vect, 2);
     random_n[3] = _mm256_extract_epi32(vect, 3);

     //Scale each integer value to the range [0, RAND_MAX]
    random_n[0] = (int)(((unsigned int)random_n[0] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
    random_n[1] = (int)(((unsigned int)random_n[1] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
    random_n[2] = (int)(((unsigned int)random_n[2] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
    random_n[3] = (int)(((unsigned int)random_n[3] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
   

        float t = -logf(random_n[0] / (float)RAND_MAX); /* move */
        
	if (weight < 0.005f) { /* roulette */ 
            if (random_n[1] / (float)RAND_MAX > 0.1f)
                break;
            weight /= 0.1f;
        }
	
	x += t * u;
        y += t * v;
        z += t * w;
        
        unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp; /* absorb */
        float xi1, xi2;
		      
       	if (shell > SHELLS - 1) {
            shell = SHELLS - 1;
        }

        /* New direction, rejection method */   
               do {
            xi1 = 2.0f * random_n[2] / (float)RAND_MAX - 1.0f;
            xi2 = 2.0f * random_n[3] / (float)RAND_MAX - 1.0f;
            t = xi1 * xi1 + xi2 * xi2;
				      
        } while (1.0f < t);   
			     
	
        heat[shell] += (1.0f - albedo) * weight;
        heat2[shell] += (1.0f - albedo) * (1 - albedo) * weight * weight; /* add up squares */
       
       
        u = 2.0f * t - 1.0f;
        weight *= albedo; 
	v = xi1 * sqrtf((1.0f - u * u) / t);

	w = xi2 * sqrtf((1.0f - u * u) / t);

        
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
    
    // init RNG
    avx_xorshift128plus_key_t my_key;


    avx_xorshift128plus_init(time(NULL), time(NULL)+1, &my_key);
    // start timer
    double start = wtime();
    // simulation
  
    for (unsigned int i = 0; i < PHOTONS; ++i) {
        photon(&my_key);
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
