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
#include "simdxorshift128plus.h"
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h> // time

#define _XOPEN_SOURCE 500       // M_PI

#define albedo (MU_S / (MU_S + MU_A))
//#define shells_per_mfp (1e4 / MICRONS_PER_SHELL / (MU_A + MU_S))


#define N_MIN_FOR 56
#define N_MAX_FOR 220//225//225/8   153

#define N_MIN_WHILE 2
#define N_MAX_WHILE 50   //9-10

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";
            
            
__m256 transformVector(__m256i x);
bool areAllFalse(float array[], int size);
__m256 raiz3(__m256 x, __m256 y, __m256 z);
// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];
//int count_f = 0;
//int count_w = 0;

float const shells_per_mfp = 1e4 / 50 / (2.0f + 20.0f);
//const __m256 shells_per_mfp_in = _mm256_set1_ps(shells_per_mfp);

/* Photon */
static void photon(avx_xorshift128plus_key_t * restrict my_key) {
    // const float albedo = MU_S / (MU_S + MU_A);
    // const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);
    
    /* STEP 1: Launching a photon packet */

    // Initial position
    //float x = 0.0f;
    //float x[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    __m256 xin = _mm256_set1_ps(0.0f);
    //float y = 0.0f;
    //float y[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    __m256 yin = _mm256_set1_ps(0.0f);
    //float z = 0.0f;
    //float z[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    __m256 zin = _mm256_set1_ps(0.0f);
    
    // Initial direction of propagation
    //float dir_x = 0.0f;
    float dir_x[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    //__m256 dir_xin = _mm256_set1_ps(0.0f);
    //float dir_y = 0.0f;
    float dir_y[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    //__m256 dir_yin = _mm256_set1_ps(0.0f);
    //float dir_z = 1.0f;
    float dir_z[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    //__m256 dir_zin = _mm256_set1_ps(1.0f);
    // Initial weight of photon
    //float weight = 1.0f;
    float weight[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    
    float flags[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    bool stop = false;
     float random_n;
     __m256i vect;
    for (int j = 0;j < 220 && stop==false; ++j) {
       
        /* Step 2: Step size selection and photon packet movement */

        // Distance the photon packet travels between interaction sites
        //float t = -logf(genRngMTInt(rand) / (float)RAND_MAX);
      
        float t[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; 
        
        for (int k=0; k < 8 ; ++k) {
            vect = avx_xorshift128plus(my_key);
            random_n = _mm256_extract_epi32(vect, 0);
            random_n = (int)(((unsigned int)random_n & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
            t[k] = -logf(random_n / (float)RAND_MAX)*flags[k];
         };
        
        
        
        __m256 flags_in = _mm256_load_ps(flags);
        
        __m256 tin = _mm256_load_ps(t);
        
        //__m256 xin = _mm256_load_ps(x);
        
        //__m256 yin = _mm256_load_ps(y);
        
        //__m256 zin = _mm256_load_ps(z);
        
        __m256 dir_xin = _mm256_load_ps(dir_x);
        
        __m256 dir_yin = _mm256_load_ps(dir_y);
        
        __m256 dir_zin = _mm256_load_ps(dir_z);
        
           dir_xin = _mm256_mul_ps(dir_xin, flags_in);
           xin = _mm256_fmadd_ps(dir_xin, tin, xin);
           dir_yin = _mm256_mul_ps(dir_yin, flags_in);
           yin = _mm256_fmadd_ps(dir_yin, tin, yin);
           dir_zin = _mm256_mul_ps(dir_zin, flags_in);
           zin = _mm256_fmadd_ps(dir_zin, tin, zin);   
      
          // _mm256_store_ps (x, xin);
          // _mm256_store_ps (y, yin);
          // _mm256_store_ps (z, zin);
          
          
        
   //     for (int k=0; k < 8 ; ++k) {
   //          x[k] *= flags[k];
   //          y[k] *= flags[k];
   //          z[k] *= flags[k];
   //     };                       
        
        
        
        /* Step 3: Absorption and scattering */



        unsigned int shell[8];
        //__mm256i shell_x;
        
        //for (int k=0; k < 8 ; ++k) {
        //     shell[k] = sqrtf(x[k] * x[k] + y[k] * y[k] + z[k] * z[k]) * shells_per_mfp;
        //unsigned int shell = sqrtf(x * x + y * y + z * z) * shells_per_mfp;
        //			}
        			
       __m256 rad = raiz3(xin, yin, zin);
       
       __m256 shells_per_mfp_in = _mm256_set1_ps(shells_per_mfp);
        
       
       
       //__m256 shell_in = _mm256_mul_ps(raiz3(xin, yin, zin), shells_per_mfp_in);
       __m256 shell_in = _mm256_mul_ps(rad, shells_per_mfp_in);
       float s[8];
       _mm256_store_ps (s, shell_in); 
       for (int i = 0; i < 8; i++) {
            shell[i] = (unsigned int)s[i];
                                         }
       
       
       
      //  float s[8];
     //_mm256_store_ps (s, shell_in);
     // for (int i = 0; i < 8; i++) {
       // printf("shell in %f\n", s[i]);
       //                                  }
       
       
       // _mm256_storeu_ps((__m256i*) shell, (__m256i)shell_in);
       //_mm256_storeu_ps((float*)shell, shell_in); 
     //  __m256i g = _mm256_castps_si256(shell_in);
    //   int o[8];
    // _mm256_store_si256 ((__m256i*)o, g);
   //   for (int i = 0; i < 8; i++) {
   //     printf("cast1 %d\n", o[i]);
   //                                      }
       
       // _mm256_storeu_si256((__m256i*)shell, _mm256_castps_si256(shell_in));		
      
      
      //  _mm256_store_ps (shell, shell_in);		
        			
        			
        float xi1[8];
        float xi2[8];

        /*
            if (shell > SHELLS - 1) {
                shell = SHELLS - 1;
            }
        */
        for(int k=0; k < 8; ++k){
          //  printf("%u\n", shell[k]);
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
       // int count = 0;  
        for(int k=0; k<8; ++k){ 
        //printf("Inizio a contare\n");
        //count=0;
        
       for(int i = 0;i < 6;i++) {
         
        
        //do {
        //for(int i=0;i<3;++i){
         // count+=1;
                                            // && (1.0f < t[k])
       //for(int j = 0; (j < N_MAX_WHILE) ; ++j) {
            //iterazioni+=1;
            vect = avx_xorshift128plus(my_key);
            random_n = _mm256_extract_epi32(vect, 0);
            random_n = (int)(((unsigned int)random_n & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
            xi1[k] = 2.0f * random_n / (float)RAND_MAX - 1.0f;
            vect = avx_xorshift128plus(my_key);
            random_n = _mm256_extract_epi32(vect, 0);
            random_n = (int)(((unsigned int)random_n & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
            xi2[k] = 2.0f * random_n / (float)RAND_MAX - 1.0f;
            t[k] = xi1[k] * xi1[k]  + xi2[k] * xi2[k] ;
          if ((1.0f >= t[k])==true) {break;}
          }
          //  }  while (1.0f < t[k]);
      //                        if ((1.0f < t[k])==true) {break;}
      
       //printf("%d\n", count);                        
       }
      
                              
       // count_w+= iterazioni/8;
        
        
        
        for(int k=0; k<8; ++k){
        dir_x[k] = 2.0f * t[k] - 1.0f;
        dir_y[k] = xi1[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
        dir_z[k] = xi2[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
                              }
                              
                              
        /* Roulette */
        for(int k=0; k<8; ++k){
        if (weight[k] < 0.005f) {
            vect = avx_xorshift128plus(my_key);
            random_n = _mm256_extract_epi32(vect, 0);
            random_n = (int)(((unsigned int)random_n & 0x7FFFFFFF) / ((double)RAND_MAX + 1.0) * (double)RAND_MAX + 0.5);
          
            if (random_n / (float)RAND_MAX > 0.1f) {
                flags[k] = 0;
            };
            weight[k] /= 0.1f;
        }
                               }
      stop = areAllFalse(flags, 8);                        
                               
    };
}

bool areAllFalse(float array[], int size) {
    for (int i = 0; i < size; i++) {
        if (array[i]==1.0f) {
            return false;
        }
    }
    return true;
}
__m256 raiz3(__m256 x, __m256 y, __m256 z) {
    // Calcola il quadrato di ciascun elemento dei tre vettori
     
    
    __m256 xSquared = _mm256_mul_ps(x, x);
    __m256 ySquared = _mm256_mul_ps(y, y);
    __m256 zSquared = _mm256_mul_ps(z, z);
     
    // Somma gli elementi dei tre vettori di quadrati
    
   __m256 sum = _mm256_add_ps(xSquared, ySquared);
     
   sum = _mm256_add_ps(sum, zSquared);
   // Calcola la radice quadrata di ciascun elemento del vettore risultante
    return _mm256_sqrt_ps(sum);
}

__m256 transformVector(__m256i x) {
    __m256i maxVal = _mm256_set1_epi32(RAND_MAX);
    __m256 floatX = _mm256_cvtepi32_ps(x);
    __m256 normalized = _mm256_div_ps(floatX, _mm256_cvtepi32_ps(maxVal));
    return _mm256_mul_ps(normalized, _mm256_cvtepi32_ps(maxVal));
}

/* Main matter */
int main(void) {
    // heading
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n#\n", PHOTONS);

    // configure RNG
    //MTRand rand = seedRand(SEED);
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
    
