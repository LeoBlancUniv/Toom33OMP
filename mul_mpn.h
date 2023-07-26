#ifndef MUL_MPN_H
#define MUL_MPN_H

#include <gmp.h>
#include "gmp-impl.h"
#include <omp.h>

void check_mul_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, int nb_limbs);

void toom3_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, mp_limb_t* scratch, int para);

void lohi22(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para);

void lomidhi32(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para);

void mullo_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para);

void mul_redc_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, mp_limb_t* p_inv, int nb_limbs, mp_limb_t* ab, int para);

void mul_mont_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, mp_limb_t* p_inv, int nb_limbs, mp_limb_t* ab, int para);


#endif