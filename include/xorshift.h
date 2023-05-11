#ifndef XORSHIFT_H
#define XORSHIFT_H

#include <stdint.h>

struct xorshift16_state {
    uint16_t a;
};

/* The state must be initialized to non-zero */
uint16_t xorshift16(struct xorshift16_state *state);

struct xorshift32_state {
    uint32_t a;
};

/* The state must be initialized to non-zero */
uint32_t xorshift32(struct xorshift32_state *state);

struct xorshift64_state {
    uint64_t a;
};

uint64_t xorshift64(struct xorshift64_state *state);

/* struct xorshift128_state can alternatively be defined as a pair
   of uint64_t or a uint128_t where supported */
struct xorshift128_state {
    uint32_t x[4];
};

/* The state must be initialized to non-zero */
uint32_t xorshift128(struct xorshift128_state *state);

#endif /* XORSHIFT_H */

