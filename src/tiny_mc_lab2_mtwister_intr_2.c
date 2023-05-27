/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#include "params.h"
#include "wtime.h"
#include "mtwister.h"
#include "simdxorshift128plus.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h> // time

#define _XOPEN_SOURCE 500       // M_PI

#define albedo (MU_S / (MU_S + MU_A))
#define shells_per_mfp (1e4 / MICRONS_PER_SHELL / (MU_A + MU_S))

#define N_MIN_FOR 56
#define N_MAX_FOR 220//225//225/8   153

#define N_MIN_WHILE 2
#define N_MAX_WHILE 50   //9-10

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";
            
            
void transformVector(__m256i x, __m256 *result);                           
bool areAllFalse(bool array[], int size);                        
// global state, heat and heat square in each shell              
static float heat[SHELLS];
static float heat2[SHELLS];
//int count_f = 0;
//int count_w = 0;

/* Photon */
static void photon(avx_xorshift128plus_key_t * restrict my_key) {
    // const float albedo = MU_S / (MU_S + MU_A);
    // const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);
    
    /* STEP 1: Launching a photon packet */

    // Initial position
    
    float random_n[1];
    __m256 xin = _mm256_set1_ps(0.0f);
    
    __m256 yin = _mm256_set1_ps(0.0f);
    
    __m256 zin = _mm256_set1_ps(0.0f);
    
    __m256 dir_xin = _mm256_set1_ps(0.0f);
    
    __m256 dir_yin = _mm256_set1_ps(0.0f);
   
    __m256 dir_zin = _mm256_set1_ps(1.0f);
    
    float weight[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    
    bool flags[8] = {true, true, true, true, true, true, true, true};

     float dir_x[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    //float dir_y = 0.0f;
    float dir_y[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    //float dir_z = 1.0f;
    float dir_z[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

         bool stop = false; 
           //int c = 0;
           __m256i vect;
     for(;stop == false;) {
    
       __m256i t_g = avx_xorshift128plus(my_key);
       __m256 tin;
       transformVector(t_g, &tin);
       __m256 r_m = _mm256_set1_ps(RAND_MAX);
       tin = _mm256_div_ps (tin, r_m);
       
        xin = _mm256_fmadd_ps(dir_xin, tin, xin);
        
        yin = _mm256_fmadd_ps(dir_yin, tin, yin);
        
        zin = _mm256_fmadd_ps(dir_zin, tin, zin);     
                              
      float x[8];
     _mm256_storeu_ps(x, xin);

     float y[8];
    _mm256_storeu_ps(y, yin);

    float z[8];
    _mm256_storeu_ps(z, zin);
    
    float t[8];
    _mm256_storeu_ps(t, tin);    
        /* Step 3: Absorption and scattering */
        
        
        
        unsigned int shell[8];
        for (int k=0; k < 8 ; ++k) {
             shell[k] = sqrtf(x[k] * x[k] + y[k] * y[k] + z[k] * z[k]) * shells_per_mfp;
        //unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;
        			}
        			
        float xi1[8];
        float xi2[8];

        /*
            if (shell > SHELLS - 1) {
                shell = SHELLS - 1;
            }
        */
        for(int k=0; k < 8; ++k){
            shell[k] = ((SHELLS - 1) * (shell[k] > SHELLS - 1)) + (shell[k] * (shell[k] <= SHELLS - 1));
                                }
        for(int k=0; k < 8; ++k){
            //variabile = condizione ? valore_se_vero : valore_se_falso;
            heat[shell[k]] = flags[k] ? heat[shell[k]]+(1.0f - albedo) * weight[k] : heat[shell[k]] ;
            heat2[shell[k]] = flags[k] ? heat2[shell[k]]+(1.0f - albedo) * (1.0f - albedo) * weight[k] * weight[k] : heat2[shell[k]] ;
             
            //heat[shell[k]] += (1.0f - albedo) * weight ;
            //heat2[shell[k]] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */
                                }
        
        for(int k=0; k<8; ++k){                        
        weight[k] *= albedo ;
                              }
        /* Step 4: Photon termination */

        /* roulette: Se agrando el valor en la condicional de 0.001 a 0.005 */
        // if (weight < 0.005f) {
        //     if (genRngMTInt(&rand) / (float)RAND_MAX > 0.1f) {
        //         break;
        //     };
        //     weight /= 0.1f;
        // }
        // weight = (weight < 0.005f) ? (((genRngMTInt(rand) / (float)RAND_MAX) > 0.1f) ? i = N_MAX_FOR : weight / 0.1f) : weight;


        /* New direction, rejection method */
       //int iterazioni = 0;

        for(int k=0; k<8; ++k){ 
        do {
           
                                            // && (1.0f < t[k])
       //for(int j = 0; (j < N_MAX_WHILE) ; ++j) {
            //iterazioni+=1;
                vect = avx_xorshift128plus(my_key);
                random_n[0] = _mm256_extract_epi32(vect, 0);
                random_n[0] = (int)(((unsigned int)random_n[0] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
                         
            xi1[k] = 2.0f * random_n[0] / (float)RAND_MAX - 1.0f;
                 vect = avx_xorshift128plus(my_key);
                random_n[0] = _mm256_extract_epi32(vect, 0);
                random_n[0] = (int)(((unsigned int)random_n[0] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
            
            
            xi2[k] = 2.0f * random_n[0] / (float)RAND_MAX - 1.0f;
            
            t[k] = xi1[k] * xi1[k]  + xi2[k] * xi2[k] ;
            }  while (1.0f < t[k]);
                              }
       // count_w+= iterazioni/8;
        
        
        
        for(int k=0; k<8; ++k){
        dir_x[k] = 2.0f * t[k] - 1.0f;
        dir_y[k] = xi1[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
        dir_z[k] = xi2[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
                              }
                              
        dir_xin = _mm256_loadu_ps(dir_x);

        // Sovrascrittura di vec2 con i valori di array2
        dir_yin = _mm256_loadu_ps(dir_y);

        // Sovrascrittura di vec3 con i valori di array3
        dir_zin = _mm256_loadu_ps(dir_z);                      
                              
        /* Roulette */
        for(int k=0; k<8; ++k){
        if (weight[k] < 0.005f) {
             vect = avx_xorshift128plus(my_key);
                random_n[0] = _mm256_extract_epi32(vect, 0);
                random_n[0] = (int)(((unsigned int)random_n[0] & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
            if (random_n[0] / (float)RAND_MAX > 0.1f) {
                flags[k] = false;
            };
            weight[k] /= 0.1f;
        }
                               }
       stop = areAllFalse(flags, 8);                        
                               
    };
}

bool areAllFalse(bool array[], int size) {
    for (int i = 0; i < size; i++) {
        if (array[i]) {
            return false;
        }
    }
    return true;
}
void transformVector(__m256i x, __m256 *result) {
    __m256i maxVal = _mm256_set1_epi32(RAND_MAX);
    __m256 floatX = _mm256_cvtepi32_ps(x);
    __m256 normalized = _mm256_div_ps(floatX, _mm256_cvtepi32_ps(maxVal));
    *result = _mm256_mul_ps(normalized, _mm256_cvtepi32_ps(maxVal));
}
/* Main matter */
int main(void) {
    // heading
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n#\n", PHOTONS);

    // configure RNG
 avx_xorshift128plus_key_t my_key;
    avx_xorshift128plus_init(SEED, SEED+1, &my_key);


    // start timer
    double start = wtime();

    // simulation
    for (unsigned int k = 0; k < PHOTONS/8; ++k) {
        photon(&my_key);
    };

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
    printf("#################\n");
 //   printf("%d\n", count_f);
   // printf("%d\n", count_w);
    printf("##################\n");
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
    
