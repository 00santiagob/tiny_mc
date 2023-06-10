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
#include <limits.h>
#include <time.h> // time

#define _XOPEN_SOURCE 500       // M_PI

// #define albedo (MU_S / (MU_S + MU_A))
// #define shells_per_mfp (1e4 / MICRONS_PER_SHELL / (MU_A + MU_S))

#define K 8 // K photons per cycle

#define N_MIN_FOR 56
#define N_MAX_FOR 225

#define N_MIN_WHILE 2
#define N_MAX_WHILE 16 // 6 ?

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";

// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];
static const float albedo = MU_S / (MU_S + MU_A);
static const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

// static bool areAllFalse(int array[], int size) {
//     for (int i = 0; i < size; i++) {
//         if (array[i]) {
//             return false;
//         };
//     };
//     return true;
// };

static int photons_lives(int array[], int size) {
    int ret = 0;
    for (int i = 0; i < size; i++) {
        ret += array[i];
    };
    return ret;
}

/* Photon */
static void photon(MTRand * restrict rand) {

    /* STEP 1: Launching a photon packet */

    // Initial position
    // float y = 0.0f;
    // float y = 0.0f;
    // float z = 0.0f;
    // float x[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    // float y[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    // float z[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    float x[K] = { [0 ... K-1] = 0.0f };
    float y[K] = { [0 ... K-1] = 0.0f };
    float z[K] = { [0 ... K-1] = 0.0f };
    
    // Initial direction of propagation
    // float dir_x = 0.0f;
    // float dir_y = 0.0f;
    // float dir_z = 1.0f;
    // float dir_x[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    // float dir_y[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    // float dir_z[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    float dir_x[K] = { [0 ... K-1] = 0.0f };
    float dir_y[K] = { [0 ... K-1] = 0.0f };
    float dir_z[K] = { [0 ... K-1] = 1.0f };

    // Initial weight of photon
    // float weight = 1.0f;
    // float weight[8] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    float weight[K] = { [0 ... K-1] = 1.0f };
    
    // True if the photon is live, False if is die.
    // bool flags[K] = {true, true, true, true, true, true, true, true};
    int flags[K] = { [0 ... K-1] = 1 };
    
//    for (int k = 0; k < K; printf("flags: %d\n", flags[k]), k++);

    // bool stop = false;

    // for (;stop == false;) {
    // for (int i = 0; i < N_MAX_FOR; ++i) {
    for (;photons_lives(flags, K);) {
    // for (int i = 0; (i < N_MAX_FOR) && photons_lives(flags, K); ++i) {

        // Distance the photon packet travels between interaction sites
        float t[K]; // = { [0 ... K-1] * 0.0f };

        unsigned int shell[K];

        float xi1[K];
        float xi2[K];

        /* Step 2: Step size selection and photon packet movement */

        for (int k=0; k < K ; ++k) {
            // t[k] = -logf(genRngMTInt(rand) / (float)RAND_MAX) * flags[k];
            t[k] = -logf(genRngMTInt(rand) / (float)RAND_MAX);    // * flags[k];

            // x[k] += t[k] * dir_x[k] * flags[k];
            // y[k] += t[k] * dir_y[k] * flags[k];
            // z[k] += t[k] * dir_z[k] * flags[k];
            x[k] += t[k] * dir_x[k];    // * flags[k];
            y[k] += t[k] * dir_y[k];    // * flags[k];
            z[k] += t[k] * dir_z[k];    // * flags[k];

            /* Step 3: Absorption and scattering */

            shell[k] = sqrtf(x[k] * x[k] + y[k] * y[k] + z[k] * z[k]) * shells_per_mfp;

            // if (shell[k] > SHELLS - 1) {
            //     shell[k] = SHELLS - 1;
            // }
            shell[k] = ((SHELLS - 1) * (shell[k] > SHELLS - 1)) + (shell[k] * (shell[k] <= SHELLS - 1));

//            printf("shell: %d\n", shell[k]);
//            printf("flag: %d\n", flags[k]);
//            printf("heat old: %f\n", heat[shell[k]]);

            // heat[shell[k]] = flags[k] ? heat[shell[k]]+(1.0f - albedo) * weight[k] : heat[shell[k]] ;
            heat[shell[k]] += (1.0f - albedo) * weight[k] * flags[k];

//             printf("heat new: %f\n", heat[shell[k]]);
            
            // heat2[shell[k]] = flags[k] ? heat2[shell[k]]+(1.0f - albedo) * (1.0f - albedo) * weight[k] * weight[k] : heat2[shell[k]];
            heat2[shell[k]] += (1.0f - albedo) * (1.0f - albedo) * weight[k] * weight[k] * flags[k];

            weight[k] *= albedo;

            /* Step 4: Photon termination */

            // if (weight[k] < 0.005f) {
            //     if (genRngMTInt(rand) / (float)RAND_MAX > 0.1f) {
            //         flags[k] = 0;
            //     };
            //     weight[k] /= 0.1f;
            // };

            // weight[k] = (weight[k] < 0.005f) ? (((genRngMTInt(rand) / (float)RAND_MAX) > 0.1f) ? i = N_MAX_FOR : weight[k] / 0.1f) : weight[k];
            // weight[k] = (weight[k] < 0.005f) ? (((genRngMTInt(rand) / (float)RAND_MAX) > 0.1f) ? flags[k] = 0 : weight[k] / 0.1f) : weight[k];
            
            flags[k] = flags[k] * (weight[k] >= 0.005f) + (weight[k] < 0.005f) * ( ( genRngMTInt(rand) / (float)RAND_MAX ) <= 0.1f );
            weight[k] = weight[k] * (weight[k] >= 0.005f) + (weight[k] < 0.005f) * weight[k] / 0.1f;
            
            // stop = areAllFalse(flags, 8);


            /* New direction, rejection method */
        
            // do {
            for(int j = 0; j < N_MAX_WHILE; ++j) {
                
                xi1[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
                xi2[k] = 2.0f * genRngMTInt(rand) / (float)RAND_MAX - 1.0f;
                t[k] = xi1[k] * xi1[k]  + xi2[k] * xi2[k] ;

                // if (1.0f >= t[k]) {
                //     break;
                // };

                j = N_MAX_WHILE * (1.0f >= t[k]) + j * (1.0f < t[k]);

            // } while (1.0f < t[k]);
            };

            dir_x[k] = 2.0f * t[k] - 1.0f;
            dir_y[k] = xi1[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
            dir_z[k] = xi2[k] * sqrtf((1.0f - dir_x[k] * dir_x[k]) / t[k]);
        };
    };
};

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

    // Binary file with data of heats for benchmark
    int len = sizeof(heat) / sizeof(heat[0]);
    FILE *fp;
    fp = fopen("dati_mod_1.bin", "wb");
    fwrite(heat, sizeof(float), len, fp);
    fclose(fp);

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

    printf("# extra\t%12.5f\n\n", heat[SHELLS - 1] / PHOTONS);

    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n\n", 1e-3 * PHOTONS / elapsed);

    return 0;
}
    
