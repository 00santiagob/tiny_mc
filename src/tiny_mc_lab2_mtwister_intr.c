/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#include "params.h"
#include "wtime.h"
#include "mtwister.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <immintrin.h>
#include <limits.h>
#include <time.h> // time

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
bool areAllFalse(bool array[], int size);
void transformVector(__m256i x, __m256 *result);
void sumOfSquares3(__m256 x, __m256 y, __m256 z, __m256 *result);
void sumOfSquares2(__m256 x, __m256 y, __m256 *result);
applyCondition(__m256 shell, float SHELLS, __m256 *result);
// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];

/* Photon */
static void photon(avx_xorshift128plus_key_t * restrict my_key) {
    // const float albedo = MU_S / (MU_S + MU_A);
    // const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

    /* STEP 1: Launching a photon packet */

    // Initial position
    //float x = 0.0f;
    //float x[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    __m256 x = _mm256_set1_ps(0.0f);
    //float y = 0.0f;
    //float y[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    __m256 y = _mm256_set1_ps(0.0f);
    //float z = 0.0f;
    //float z[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    __m256 z = _mm256_set1_ps(0.0f);
    // Initial direction of propagation
    //float dir_x = 0.0f;
    //float dir_x[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    __m256 dir_x = _mm256_set1_ps(0.0f);
    //float dir_y = 0.0f;
    //float dir_y[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    __m256 dir_y = _mm256_set1_ps(0.0f);
    //float dir_z = 1.0f;
    //float dir_z[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    _m256 dir_z = _mm256_set1_ps(1.0f);
    
    // Initial weight of photon
    //float weight = 1.0f;
    //float weight[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    __m256 weight = _mm256_set1_ps(1.0f);
    
    //bool flags[8] = {true, true, true, true, true, true, true, true};
    __m256i flags = _mm256_set1_epi32(0xFFFFFFFF);
        
    //bool stop = false;
    for (;;) {
    // for (int i = 0; i<N_MAX_FOR; ++i) {

        /* Step 2: Step size selection and photon packet movement */

        // Distance the photon packet travels between interaction sites
        //float t = -logf(genRngMTInt(rand) / (float)RAND_MAX);
         //_m256 t = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
        //for (int k=0; k < 8 ; ++k) {
            //t[k] = -logf(genRngMTInt(rand) / (float)RAND_MAX)*flags[k];
         //};
             __m256i t_g = avx_xorshift128plus(my_key);
             __m256 t;
             transformVector(t_g, &t);
             __m256 r_m = _mm256_set1_ps(RAND_MAX);
             t = _mm256_div_ps (t, r_m);
        
       
        //x += t * dir_x;
        //for (int k=0; k < 8 ; ++k) {
          //   x[k] += t[k] * dir_x[k] * flags[k];
        //};
        x = _mm256_fmadd_ps(dir_x, t, x);
        //y += t * dir_y;
        //for (int k=0; k < 8 ; ++k) {
            // y[k] += t[k] * dir_y[k]* flags[k];
        //};
        y = _mm256_fmadd_ps(dir_y, t, y);
        //z += t * dir_z;
        //for (int k=0; k < 8 ; ++k) {
             //z[k] += t[k] * dir_z[k]*flags[k];
               
                                     //}
        z = _mm256_fmadd_ps(dir_z, t, z);                             
        //};                       
        /* Step 3: Absorption and scattering */

        //unsigned int shell[8];
        //for (int k=0; k < 8 ; ++k) {
        //     shell[k] = sqrtf(x[k] * x[k] + y[k] * y[k] + z[k] * z[k]) * shells_per_mfp;
        //unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;
        			//}
        _m256 rad;			
        sumOfSquares3(x, y, z, &rad);
        _m256 shell = _mm256_mul_ps(rad, _m256_set1_ps(shells_per_mfp);
        
        
        			
        //float xi1[8];
        //float xi2[8];
        __m256 xi1;
        __m256 xi2;

        /*
            if (shell > SHELLS - 1) {
                shell = SHELLS - 1;
            }
        */
        //for(int k=0; k < 8; ++k){
         //   shell[k] = ((SHELLS - 1) * (shell[k] > SHELLS - 1)) + (shell[k] * (shell[k] <= SHELLS - 1));
           //                     }
       __m256 shell_c;   
       applyCondition(shell, (float)SHELLS, &shell_c);                   
      
      unsigned int shell_u[8];                          
      _mm256_storeu_ps(shell_u, shell_c);                         
                                
                                
        for(int k=0; k < 8; ++k){
            //variabile = condizione ? valore_se_vero : valore_se_falso;
            heat[shell_c[k]] = _mm256_extract_epi32(mask, k)==1 ? heat[shell_c[k]]+(1.0f - albedo) * _mm256_extract_epi32(weight, k) : heat[shell_c[k]] ;
            heat2[shell_c[k]] = _mm256_extract_epi32(mask, k)==1 ? heat2[shell_c[k]]+(1.0f - albedo) * (1.0f - albedo) * _mm256_extract_epi32(weight, k) * _mm256_extract_epi32(weight, k) : heat2[shell_c[k]];
             
            //heat[shell[k]] += (1.0f - albedo) * weight ;
            //heat2[shell[k]] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */
                                }
        
        //for(int k=0; k<8; ++k){                        
        //weight[k] *= albedo ;
        //                      }
        weight = _mm256_mul_ps(weight, _mm256_set1_ps(albedo));
        
        /* Step 4: Photon termination */

        /* roulette: Se agrando el valor en la condicional de 0.001 a 0.005 */
        // if (weight < 0.005f) {
        //     if (genRngMTInt(&rand) / (float)RAND_MAX > 0.1f) {
        //         break;
        //     };
        //     weight /= 0.1f;
        // }
        // weight = (weight < 0.005f) ? (((genRngMTInt(rand) / (float)RAND_MAX) > 0.1f) ? i = N_MAX_FOR : weight / 0.1f) : weight;

        //__m256 xi1;
        //__m256 xi2;
        
        /* New direction, rejection method */
        
        //for(int k=0; k<8; ++k){
        do {
        // for(int j = 0; (j < N_MAX_WHILE) && (1.0f < t); ++j) {
             __m256i t_g = avx_xorshift128plus(my_key);
             __m256 t;
             transformVector(t_g, &t);
             __m256 a = _mm256_set1_ps(2.0f/(float)RAND_MAX);
             __m256 b = _mm256_set1_ps(-1.0f);
             
            //xi1[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
            xi1 = __mm256_add_ps(b, _mm256_mul_ps(a,t))
            __m256i t_g = avx_xorshift128plus(my_key);
             __m256 t;
             transformVector(t_g, &t);
            //xi2[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
            xi2 = __mm256_add_ps(b, _mm256_mul_ps(a,t))
        
            t = sumOfSquares2(xi1, xi2) ;
        
        
        } while (1.0f < t[k]);
                              }
        // };
        
        for(int k=0; k<8; ++k){
        dir_x[k] = 2.0f * t[k] - 1.0f;
        dir_y[k] = xi1[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
        dir_z[k] = xi2[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
                              }
        /* Roulette */
        for(int k=0; k<8; ++k){
        if (weight[k] < 0.005f) {
            if (genRngMTInt(rand) / (float)RAND_MAX > 0.1f) {
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

void sumOfSquares3(__m256 x, __m256 y, __m256 z, __m256 *result) {
    // Calcola il quadrato di ciascun elemento dei tre vettori
    __m256 xSquared = _mm256_mul_ps(x, x);
    __m256 ySquared = _mm256_mul_ps(y, y);
    __m256 zSquared = _mm256_mul_ps(z, z);

    // Somma gli elementi dei tre vettori di quadrati
    __m256 sum = _mm256_add_ps(xSquared, ySquared);
    sum = _mm256_add_ps(sum, zSquared);

    // Calcola la radice quadrata di ciascun elemento del vettore risultante
    *result = _mm256_sqrt_ps(sum);
}
void sumOfSquares2(__m256 x, __m256 y, __m256 *result) {
    // Calcola il quadrato di ciascun elemento dei tre vettori
    __m256 xSquared = _mm256_mul_ps(x, x);
    __m256 ySquared = _mm256_mul_ps(y, y);
 

    // Somma gli elementi dei tre vettori di quadrati
    __m256 sum = _mm256_add_ps(xSquared, ySquared);
    

    // Calcola la radice quadrata di ciascun elemento del vettore risultante
    *result = _mm256_sqrt_ps(sum);
}

void applyCondition(__m256 shell, float SHELLS, __m256 *result) {
    __m256 shellsMinusOne = _mm256_set1_ps(SHELLS - 1);

    // Esegui la verifica su ogni elemento di shell
    __m256 cmpResult = _mm256_cmp_ps(shell, shellsMinusOne, _CMP_GT_OS);

    // Seleziona il valore corretto per ciascun elemento
    *result = _mm256_blendv_ps(shell, shellsMinusOne, cmpResult);
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
    
