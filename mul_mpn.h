#ifndef MUL_MPN_H
#define MUL_MPN_H

#include <gmp.h>
#include <omp.h>

void toom3_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, mp_limb_t* scratch, int v);

void lohi22(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab);

void lomidhi32(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab);


#endif