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

//#define albedo (MU_S / (MU_S + MU_A))
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
float const albedo = (MU_S / (MU_S + MU_A));
//const __m256 shells_per_mfp_in = _mm256_set1_ps(shells_per_mfp);

/* Photon */
static void photon(MTRand * restrict rand) {
    
    
    /* STEP 1: Launching a photon packet */

    // Initial position
   
    __m256 xin = _mm256_set1_ps(0.0f);
    
    __m256 yin = _mm256_set1_ps(0.0f);
    
    __m256 zin = _mm256_set1_ps(0.0f);
    
    // Initial direction of propagation
    
    __m256 dir_xin = _mm256_set1_ps(0.0f);
    
    __m256 dir_yin = _mm256_set1_ps(0.0f);
   
    __m256 dir_zin = _mm256_set1_ps(1.0f);
    // Initial weight of photon

    float weight[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

    
    float flags[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    bool stop = false;

    for (int j = 0;j < 220 && stop==false; ++j) {
       
        /* Step 2: Step size selection and photon packet movement */

        // Distance the photon packet travels between interaction sites
        //float t = -logf(genRngMTInt(rand) / (float)RAND_MAX);
      
        float t[8]; 
        
           
        for (int k=0; k < 8 ; ++k) {
            
            t[k] = -logf(genRngMTInt(rand) / (float)RAND_MAX)*flags[k];
         };
      

        
        __m256 flags_in = _mm256_load_ps(flags);
        
        __m256 tin = _mm256_load_ps(t);
        
        
           dir_xin = _mm256_mul_ps(dir_xin, flags_in);
           xin = _mm256_fmadd_ps(dir_xin, tin, xin);
           dir_yin = _mm256_mul_ps(dir_yin, flags_in);
           yin = _mm256_fmadd_ps(dir_yin, tin, yin);
           dir_zin = _mm256_mul_ps(dir_zin, flags_in);
           zin = _mm256_fmadd_ps(dir_zin, tin, zin);   
                           
                
        
        /* Step 3: Absorption and scattering */



        unsigned int shell[8];

        			
       __m256 rad = raiz3(xin, yin, zin);
       
       __m256 shells_per_mfp_in = _mm256_set1_ps(shells_per_mfp);
        
       
       

       __m256 shell_in = _mm256_mul_ps(rad, shells_per_mfp_in);
       float s[8];
       _mm256_store_ps (s, shell_in); 
       for (int i = 0; i < 8; i++) {
            shell[i] = (unsigned int)s[i];
                                         }
       
       
      	
        			
        			
        float xi1[8];
        float xi2[8];

    
        for(int k=0; k < 8; ++k){

            shell[k] = ((SHELLS - 1) * (shell[k] > SHELLS - 1)) + (shell[k] * (shell[k] <= SHELLS - 1));
                                }
        for(int k=0; k < 8; ++k){

            heat[shell[k]] = flags[k] ? heat[shell[k]]+(1.0f - albedo) * weight[k] : heat[shell[k]] ;
            heat2[shell[k]] = flags[k] ? heat2[shell[k]]+(1.0f - albedo) * (1.0f - albedo) * weight[k] * weight[k] : heat2[shell[k]] ;
             
           
                                }
                                
                                
        __m256 weightin = _mm256_load_ps(weight);
        
        weightin = _mm256_mul_ps(weightin, _mm256_set1_ps(albedo));
        

        /* New direction, rejection method */

       
       _mm256_store_ps(t, tin);
         
        for(int k=0; k<8; ++k){ 
       
        
       for(int i = 0;i < 6;i++) {
         
  
            
            xi1[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
            
            xi2[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
            t[k] = xi1[k] * xi1[k]  + xi2[k] * xi2[k] ;
          if ((1.0f >= t[k])==true) {break;}
          }
                                 
       }
      
                              
       
        
        tin = _mm256_load_ps(t);
        __m256 xin1 = _mm256_load_ps(xi1);
        __m256 xin2 = _mm256_load_ps(xi2); 
        dir_xin = _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), tin),_mm256_set1_ps(-1.0f));
        
        __m256 sqt = _mm256_sqrt_ps(_mm256_div_ps(_mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(dir_xin, dir_xin)), tin));
        dir_yin = _mm256_mul_ps(xin1, sqt);          
        dir_zin = _mm256_mul_ps(xin2, sqt);   
                              
                              
        /* Roulette */
      
      _mm256_store_ps(weight, weightin);    
      
       
        for(int k=0; k<8; ++k){
        if (weight[k] < 0.005f) {                      
            if (genRngMTInt(rand) / (float)RAND_MAX > 0.1f) {
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
    float ji[8];
    _mm256_store_ps(ji, _mm256_castsi256_ps(x));
  
    for(int i = 0; i < 8; ++i){    
        printf("RNG %f\n", ji[i]);
    }
    
    __m256 maxVal1 = _mm256_set1_ps(RAND_MAX+1.0f);
    __m256 maxVal2 = _mm256_set1_ps(RAND_MAX+0.5f);
    __m256i unsignedX = _mm256_abs_epi32(x);
    __m256 floatX = _mm256_castsi256_ps(unsignedX);
    __m256  ts = _mm256_set1_ps(0x7FFFFFFF);
    __m256 scaledX = _mm256_mul_ps(floatX, ts);
    __m256 normalized = _mm256_div_ps(scaledX, _mm256_mul_ps(maxVal1, maxVal2));
    return normalized;
}

/* Main matter */
int main(void) {
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
    for (unsigned int k = 0; k < PHOTONS/8; ++k) {
        photon(&rand);
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
    
