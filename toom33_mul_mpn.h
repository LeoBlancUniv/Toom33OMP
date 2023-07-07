#ifndef TOOM33_MUL_MPN_H
#define TOOM33_MUL_MPN_H

#include <gmp.h>
#include <omp.h>

void toom3_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, mp_limb_t* scratch, int v);

#endif