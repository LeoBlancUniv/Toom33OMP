/* mpn_toom33_mul -- Multiply {input_a,size_a} and {p,size_b} where size_a and size_b are close in
   size.  Or more accurately, size_b <= size_a < (3/2)size_b.

   Contributed to the GNU project by Torbjorn Granlund.
   Additional improvements by Marco Bodrato.

   THE FUNCTION IN THIS FILE IS INTERNAL WITH A MUTABLE INTERFACE.  IT IS ONLY
   SAFE TO REACH IT THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT IT WILL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 2006-2008, 2010, 2012, 2015 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */


#include "gmp-impl.h"
#include <omp.h>

/* Evaluate in: -1, 0, +1, +2, +inf

  <-s--><--n--><--n-->
   ____ ______ ______
  |_a2_|___a1_|___a0_|
   |b2_|___b1_|___b0_|
   <-t-><--n--><--n-->

  v0  =  a0         * b0          #   A(0)*B(0)
  v1  = (a0+ a1+ a2)*(b0+ b1+ b2) #   A(1)*B(1)      ah  <= 2  bh <= 2
  vm1 = (a0- a1+ a2)*(b0- b1+ b2) #  A(-1)*B(-1)    |ah| <= 1  bh <= 1
  v2  = (a0+2a1+4a2)*(b0+2b1+4b2) #   A(2)*B(2)      ah  <= 6  bh <= 6
  vinf=          a2 *         b2  # A(inf)*B(inf)
*/

#if TUNE_PROGRAM_BUILD || WANT_FAT_BINARY
#define MAYBE_mul_basecase 1
#define MAYBE_mul_toom33   1
#else
#define MAYBE_mul_basecase						\
  (MUL_TOOM33_THRESHOLD < 3 * MUL_TOOM22_THRESHOLD)
#define MAYBE_mul_toom33						\
  (MUL_TOOM44_THRESHOLD >= 3 * MUL_TOOM33_THRESHOLD)
#endif

/* FIXME: TOOM33_MUL_N_REC is not quite right for a balanced
   multiplication at the infinity point. We may have
   MAYBE_mul_basecase == 0, and still get s just below
   MUL_TOOM22_THRESHOLD. If MUL_TOOM33_THRESHOLD == 7, we can even get
   s == 1 and mpn_toom22_mul will crash.
*/

/*
  this should be what dictate what alg is going to be used when doing the mult
  toom22 is equivalent to kara
  toom33 is our entry point, which should call toom22 which should call basecase
*/
#define TOOM33_MUL_N_REC(p, a, b, n, ws)				\
  do {									\
    if (/*MAYBE_mul_basecase						\
	&& BELOW_THRESHOLD (n, MUL_TOOM22_THRESHOLD )*/1)			\
      mpn_mul_basecase (p, a, n, b, n);					\
    else if (! MAYBE_mul_toom33						\
	     || BELOW_THRESHOLD (n, MUL_TOOM33_THRESHOLD))		\
      mpn_toom22_mul (p, a, n, b, n, ws);				\
    else								\
      mpn_toom33_mul (p, a, n, b, n, ws);				\
  } while (0)

void
mpn_toom33_mul (mp_ptr pp,
		mp_srcptr input_a, mp_size_t size_a,
		mp_srcptr input_b, mp_size_t size_b,
		mp_ptr scratch)
{
  const int __gmpn_cpuvec_initialized = 1;
  mp_size_t nb_block_coeff, s, t;
  int vm1_neg;
  mp_limb_t cy, vinf0;
  mp_ptr gp;
  mp_ptr as1, asm1, as2;
  mp_ptr bs1, bsm1, bs2;

  #define a0  input_a
  #define a1  (input_a + nb_block_coeff)
  #define a2  (input_a + 2*nb_block_coeff)
  #define b0  input_b
  #define b1  (input_b + nb_block_coeff)
  #define b2  (input_b + 2*nb_block_coeff)

  nb_block_coeff = (size_a + 2) / (size_t) 3;

  s = size_a - 2 * nb_block_coeff;
  t = size_b - 2 * nb_block_coeff;

  ASSERT (size_a >= size_b);

  ASSERT (0 < s && s <= nb_block_coeff);
  ASSERT (0 < t && t <= nb_block_coeff);

  as1  = scratch + 4 * nb_block_coeff + 4;
  asm1 = scratch + 2 * nb_block_coeff + 2;
  as2 = pp + nb_block_coeff + 1;

  bs1 = pp;
  bsm1 = scratch + 3 * nb_block_coeff + 3; /* we need 4n+4 <= 4n+s+t */
  bs2 = pp + 2 * nb_block_coeff + 2;

  gp = scratch;

  vm1_neg = 0;

  /* Compute as1 and asm1.  */
  cy = mpn_add (gp, a0, nb_block_coeff, a2, s);
  #if HAVE_NATIVE_mpn_add_n_sub_n
    if (cy == 0 && mpn_cmp (gp, a1, nb_block_coeff) < 0)
      {
        cy = mpn_add_n_sub_n (as1, asm1, a1, gp, nb_block_coeff);
        as1[nb_block_coeff] = cy >> 1;
        asm1[nb_block_coeff] = 0;
        vm1_neg = 1;
      }
    else
      {
        mp_limb_t cy2;
        cy2 = mpn_add_n_sub_n (as1, asm1, gp, a1, nb_block_coeff);
        as1[nb_block_coeff] = cy + (cy2 >> 1);
        asm1[nb_block_coeff] = cy - (cy2 & 1);
      }
  #else
    as1[nb_block_coeff] = cy + mpn_add_n (as1, gp, a1, nb_block_coeff);
    if (cy == 0 && mpn_cmp (gp, a1, nb_block_coeff) < 0)
      {
        mpn_sub_n (asm1, a1, gp, nb_block_coeff);
        asm1[nb_block_coeff] = 0;
        vm1_neg = 1;
      }
    else
      {
        cy -= mpn_sub_n (asm1, gp, a1, nb_block_coeff);
        asm1[nb_block_coeff] = cy;
      }
  #endif

  /* Compute as2.  */
  #if HAVE_NATIVE_mpn_rsblsh1_n
    cy = mpn_add_n (as2, a2, as1, s);
    if (s != nb_block_coeff)
      cy = mpn_add_1 (as2 + s, as1 + s, nb_block_coeff - s, cy);
    
    cy += as1[nb_block_coeff];
    cy = 2 * cy + mpn_rsblsh1_n (as2, a0, as2, nb_block_coeff);
  #else
    #if HAVE_NATIVE_mpn_addlsh1_n
      cy  = mpn_addlsh1_n (as2, a1, a2, s);
      if (s != nb_block_coeff)
        cy = mpn_add_1 (as2 + s, a1 + s, nb_block_coeff - s, cy);
      cy = 2 * cy + mpn_addlsh1_n (as2, a0, as2, nb_block_coeff);
    #else
      cy = mpn_add_n (as2, a2, as1, s);
      if (s != nb_block_coeff)
        cy = mpn_add_1 (as2 + s, as1 + s, nb_block_coeff - s, cy);

      cy += as1[nb_block_coeff];
      cy = 2 * cy + mpn_lshift (as2, as2, nb_block_coeff, 1);
      cy -= mpn_sub_n (as2, as2, a0, nb_block_coeff);
    #endif
  #endif
  as2[nb_block_coeff] = cy;

  /* Compute bs1 and bsm1.  */
  cy = mpn_add (gp, b0, nb_block_coeff, b2, t);
  #if HAVE_NATIVE_mpn_add_n_sub_n
    if (cy == 0 && mpn_cmp (gp, b1, nb_block_coeff) < 0)
      {
        cy = mpn_add_n_sub_n (bs1, bsm1, b1, gp, nb_block_coeff);
        bs1[nb_block_coeff] = cy >> 1;
        bsm1[nb_block_coeff] = 0;
        vm1_neg ^= 1;
      }
    else
      {
        mp_limb_t cy2;
        cy2 = mpn_add_n_sub_n (bs1, bsm1, gp, b1, nb_block_coeff);
        bs1[nb_block_coeff] = cy + (cy2 >> 1);
        bsm1[nb_block_coeff] = cy - (cy2 & 1);
      }
  #else
    bs1[nb_block_coeff] = cy + mpn_add_n (bs1, gp, b1, nb_block_coeff);
    if (cy == 0 && mpn_cmp (gp, b1, nb_block_coeff) < 0)
      {
        mpn_sub_n (bsm1, b1, gp, nb_block_coeff);
        bsm1[nb_block_coeff] = 0;
        vm1_neg ^= 1;
      }
    else
      {
        cy -= mpn_sub_n (bsm1, gp, b1, nb_block_coeff);
        bsm1[nb_block_coeff] = cy;
      }
  #endif

  /* Compute bs2.  */
  #if HAVE_NATIVE_mpn_rsblsh1_n
    cy = mpn_add_n (bs2, b2, bs1, t);
    if (t != nb_block_coeff)
      cy = mpn_add_1 (bs2 + t, bs1 + t, nb_block_coeff - t, cy);

    cy += bs1[nb_block_coeff];
    cy = 2 * cy + mpn_rsblsh1_n (bs2, b0, bs2, nb_block_coeff);
  #else
    #if HAVE_NATIVE_mpn_addlsh1_n
      cy  = mpn_addlsh1_n (bs2, b1, b2, t);
      if (t != nb_block_coeff)
        cy = mpn_add_1 (bs2 + t, b1 + t, nb_block_coeff - t, cy);

      cy = 2 * cy + mpn_addlsh1_n (bs2, b0, bs2, nb_block_coeff);
    #else
      cy  = mpn_add_n (bs2, bs1, b2, t);
      if (t != nb_block_coeff)
        cy = mpn_add_1 (bs2 + t, bs1 + t, nb_block_coeff - t, cy);

      cy += bs1[nb_block_coeff];
      cy = 2 * cy + mpn_lshift (bs2, bs2, nb_block_coeff, 1);
      cy -= mpn_sub_n (bs2, bs2, b0, nb_block_coeff);
    #endif
  #endif
  bs2[nb_block_coeff] = cy;

  ASSERT (as1[nb_block_coeff] <= 2);
  ASSERT (bs1[nb_block_coeff] <= 2);
  ASSERT (asm1[nb_block_coeff] <= 1);
  ASSERT (bsm1[nb_block_coeff] <= 1);
  ASSERT (as2[nb_block_coeff] <= 6);
  ASSERT (bs2[nb_block_coeff] <= 6);

  #define v0    pp				/* 2n */
  #define v1    (pp + 2 * nb_block_coeff)			/* 2n+1 */
  #define vinf  (pp + 4 * nb_block_coeff)			/* s+t */
  #define vm1   scratch				/* 2n+1 */
  #define v2    (scratch + 2 * nb_block_coeff + 1)		/* 2n+2 */
  #define scratch_out  (scratch + 5 * nb_block_coeff + 5)

  //looks like SMALLER_REDUCTION isn't defined

  /*
    forced basecase : 
      {1 2 3 4 5} ok

      {1 2 3 4}{5} no
      {1} {2 3 4 5} no

      {1 2 3}{4}{5} no
      {1 2 3}{4 5} no
      {1}{2 3 4}{5} no
      {1}{2}{3 4 5} no
      {1 2}{3 4 5} no

      {1 2}{3}{4}{5} no
      {1}{2 3}{4}{5} no
      {1}{2}{3 4}{5} no
      {1}{2}{3}{4 5} no

      


      {1}{2}{3}{4}{5} no
      
  */
  omp_set_num_threads(2);
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      /* vm1, 2n+1 limbs */
      #ifdef SMALLER_RECURSION
        
        /*TOOM33_MUL_N_REC (vm1, asm1, bsm1, nb_block_coeff, scratch_out);
        cy = 0;
        if (asm1[nb_block_coeff] != 0)
          cy = bsm1[nb_block_coeff] + mpn_add_n (vm1 + nb_block_coeff, vm1 + nb_block_coeff, bsm1, nb_block_coeff);

        if (bsm1[nb_block_coeff] != 0)
          cy += mpn_add_n (vm1 + nb_block_coeff, vm1 + nb_block_coeff, asm1, nb_block_coeff);

        vm1[2 * nb_block_coeff] = cy;*/
      #else
        
        TOOM33_MUL_N_REC (vm1, asm1, bsm1, nb_block_coeff + 1, scratch_out);
      #endif

      //-------------------------
    //}

    /*#pragma omp section
    {*/
      TOOM33_MUL_N_REC (v2, as2, bs2, nb_block_coeff + 1, scratch_out); /* v2, 2n+1 limbs */

      //-------------------------
    //}

    /*#pragma omp section
    {*/
      /* vinf, s+t limbs */
      if (s > t)  mpn_mul (vinf, a2, s, b2, t);
      else        TOOM33_MUL_N_REC (vinf, a2, b2, s, scratch_out);

      vinf0 = vinf[0];        /* v1 overlaps with this */
    
      //-------------------------
    //}
    /*#pragma omp section
    {*/
      #ifdef SMALLER_RECURSION
        /* v1, 2n+1 limbs */
        /*TOOM33_MUL_N_REC (v1, as1, bs1, nb_block_coeff, scratch_out);
        if (as1[nb_block_coeff] == 1){
            cy = bs1[nb_block_coeff] + mpn_add_n (v1 + nb_block_coeff, v1 + nb_block_coeff, bs1, nb_block_coeff);
        }
        else if (as1[nb_block_coeff] != 0){
          #if HAVE_NATIVE_mpn_addlsh1_n_ip1
                cy = 2 * bs1[nb_block_coeff] + mpn_addlsh1_n_ip1 (v1 + nb_block_coeff, bs1, nb_block_coeff);
          #else
                cy = 2 * bs1[nb_block_coeff] + mpn_addmul_1 (v1 + nb_block_coeff, bs1, nb_block_coeff, CNST_LIMB(2));
          #endif
        }
        else
          cy = 0;

        if (bs1[nb_block_coeff] == 1){
          cy += mpn_add_n (v1 + nb_block_coeff, v1 + nb_block_coeff, as1, nb_block_coeff);
        }
        else if (bs1[nb_block_coeff] != 0){
          #if HAVE_NATIVE_mpn_addlsh1_n_ip1
                cy += mpn_addlsh1_n_ip1 (v1 + nb_block_coeff, as1, nb_block_coeff);
          #else
                cy += mpn_addmul_1 (v1 + nb_block_coeff, as1, nb_block_coeff, CNST_LIMB(2));
          #endif
        }

        v1[2 * nb_block_coeff] = cy;*/

      #else

        

        cy = vinf[1];

        mp_ptr bs1_ = calloc (nb_block_coeff+2, sizeof(mp_limb_t));
        mpn_copyi(bs1_, bs1, nb_block_coeff+2);

        /*gmp_printf("%Nx \n", bs1, nb_block_coeff+2);
        gmp_printf("%Nx \n\n", bs1_, nb_block_coeff+2);*/


        TOOM33_MUL_N_REC (v1, as1, bs1, nb_block_coeff + 1, scratch_out);
        /*if (vinf[1] != cy){
          printf("%ld , ", vinf[1]);
          printf("%ld , \n", cy);
        }*/
        vinf[1] = cy;

        free(bs1_);

         
      #endif
      //-------------------------
    //}

    /*#pragma omp section
    {*/
      TOOM33_MUL_N_REC (v0, input_a, input_b, nb_block_coeff, scratch_out); /* v0, 2n limbs */

      //-------------------------
    }
  }
  

  mpn_toom_interpolate_5pts (pp, v2, vm1, nb_block_coeff, s + t, vm1_neg, vinf0);
}