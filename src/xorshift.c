#include <stdint.h>
#include "xorshift.h"

/* The state must be initialized to non-zero */

uint16_t xorshift16(struct xorshift16_state *state)
{
	/* Algorithm "xor" by ME */
	uint16_t x = state->a;
	x ^= x << 5;
	x ^= x >> 7;
	x ^= x << 3;
	return state->a = x;
}

uint32_t xorshift32(struct xorshift32_state *state)
{
	/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
	uint32_t x = state->a;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return state->a = x;
}

uint64_t xorshift64(struct xorshift64_state *state)
{
	uint64_t x = state->a;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return state->a = x;
}

/* struct xorshift128_state can alternatively be defined as a pair
   of uint64_t or a uint128_t where supported */

/* The state must be initialized to non-zero */
uint32_t xorshift128(struct xorshift128_state *state)
{
	/* Algorithm "xor128" from p. 5 of Marsaglia, "Xorshift RNGs" */
	uint32_t t  = state->x[3];
    
    uint32_t s  = state->x[0];  /* Perform a contrived 32-bit shift. */
	state->x[3] = state->x[2];
	state->x[2] = state->x[1];
	state->x[1] = s;

	t ^= t << 11;
	t ^= t >> 8;
	return state->x[0] = t ^ s ^ (s >> 19);
}
