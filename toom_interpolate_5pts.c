/* mpn_toom_interpolate_5pts -- Interpolate for toom3, 33, 42.

   Contributed to the GNU project by Robert Harley.
   Improvements by Paul Zimmermann and Marco Bodrato.

   THE FUNCTION IN THIS FILE IS INTERNAL WITH A MUTABLE INTERFACE.  IT IS ONLY
   SAFE TO REACH IT THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT IT WILL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 2000-2003, 2005-2007, 2009 Free Software Foundation, Inc.

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

void
mpn_toom_interpolate_5pts (mp_ptr ret, mp_ptr v2, mp_ptr vm1,
			   mp_size_t nb_block_coeff, mp_size_t nb_block_last_coeff, int vm1_sign,
			   mp_limb_t vinf0)
{
    //vm1_sign : (0, pos), (1, neg)

  mp_limb_t cy, saved; //carry stuff ?
  mp_size_t twok;  // = nb_block_coeff * 2
  mp_size_t kk1;   // = nb_block_coeff * 2 + 1
  mp_ptr ab1, ab2, ab3, ab4;  //acces into return value

  twok = nb_block_coeff + nb_block_coeff;
  kk1 = twok + 1;

  ab1 = ret  + nb_block_coeff; //ab1   = ret + nb_block_coeff
  ab2 = ab1 + nb_block_coeff;   //ab2   = ret + nb_block_coeff * 2
  ab3 = ab2 + nb_block_coeff;   //ab3   = ret + nb_block_coeff * 3
  ab4 = ab3 + nb_block_coeff; //ab4 = ret + nb_block_coeff * 4

#define ab0 (ret)              //ab0   = ret


  /* (1) v2 <- v2-vm1 < v2+|vm1|,       (16 8 4 2 1) - (1 -1 1 -1  1) =
     thus 0 <= v2 < 50*B^(2k) < 2^6*B^(2k)             (15 9 3  3  0)
  */
  if (vm1_sign) //we add vm1 if it's neg, we sub it instead if it's pos
    ASSERT_NOCARRY (mpn_add_n (v2, v2, vm1, kk1));
  else
    ASSERT_NOCARRY (mpn_sub_n (v2, v2, vm1, kk1));

  /* {ret,2k} {ret+2k,2k+1} {ret+4k+1,2r-1} {scratch,2k+1} {scratch+2k+1,2k+1} {scratch+4k+2,2r}
       ab0       ab2       hi(ab4)       |vm1|     v2-vm1      EMPTY */

  ASSERT_NOCARRY (mpn_divexact_by3 (v2, v2, kk1));    /* v2 <- v2 / 3 */
						      /* (5 3 1 1 0)*/

  /* {ret,2k} {ret+2k,2k+1} {ret+4k+1,2r-1} {scratch,2k+1} {scratch+2k+1,2k+1} {scratch+4k+2,2r}
       ab0       ab2      hi(ab4)       |vm1|     (v2-vm1)/3    EMPTY */

  /* (2) vm1 <- tm1 := (ab2 - vm1) / 2  [(1 1 1 1 1) - (1 -1 1 -1 1)] / 2 =
     tm1 >= 0                                         (0  1 0  1 0)
     No carry comes out from {ab2, kk1} +/- {vm1, kk1},

    //

     and the division by two is exact.
     If (vm1_sign!=0) the sign of vm1 is negative */
  if (vm1_sign)
    {

      ASSERT_NOCARRY (mpn_add_n (vm1, ab2, vm1, kk1));
      ASSERT_NOCARRY (mpn_rshift (vm1, vm1, kk1, 1));
    }
  else
    {
      ASSERT_NOCARRY (mpn_sub_n (vm1, ab2, vm1, kk1));
      ASSERT_NOCARRY (mpn_rshift (vm1, vm1, kk1, 1));
    }

  /* {ret,2k} {ret+2k,2k+1} {ret+4k+1,2r-1} {scratch,2k+1} {scratch+2k+1,2k+1} {scratch+4k+2,2r}
       ab0       ab2        hi(ab4)       tm1     (v2-vm1)/3    EMPTY */

  /* (3) ab2 <- t1 := ab2 - ab0    (1 1 1 1 1) - (0 0 0 0 1) = (1 1 1 1 0)
     t1 >= 0
  */
  ab4[0] -= mpn_sub_n (ab2, ab2, ret, twok);

  /* {ret,2k} {ret+2k,2k+1} {ret+4k+1,2r-1} {scratch,2k+1} {scratch+2k+1,2k+1} {scratch+4k+2,2r}
       ab0     ab2-ab0        hi(ab4)       tm1     (v2-vm1)/3    EMPTY */

  /* (4) v2 <- t2 := ((v2-vm1)/3-t1)/2 = (v2-vm1-3*t1)/6
     t2 >= 0                  [(5 3 1 1 0) - (1 1 1 1 0)]/2 = (2 1 0 0 0)
  */

  ASSERT_NOCARRY (mpn_sub_n (v2, v2, ab2, kk1));
  ASSERT_NOCARRY (mpn_rshift (v2, v2, kk1, 1));


  /* {ret,2k} {ret+2k,2k+1} {ret+4k+1,2r-1} {scratch,2k+1} {scratch+2k+1,2k+1} {scratch+4k+2,2r}
       ab0     ab2-ab0        hi(ab4)     tm1    (v2-vm1-3t1)/6    EMPTY */

  /* (5) ab2 <- t1-tm1           (1 1 1 1 0) - (0 1 0 1 0) = (1 0 1 0 0)
     result is ab2 >= 0
  */
  ASSERT_NOCARRY (mpn_sub_n (ab2, ab2, vm1, kk1));

  /* We do not need to read the value in vm1, so we add it in {ret+nb_block_coeff, ...} */
  cy = mpn_add_n (ab1, ab1, vm1, kk1);
  MPN_INCR_U (ab3 + 1, nb_block_last_coeff + nb_block_coeff - 1, cy); /* 2n-(3k+1) = 2r+nb_block_coeff-1 */
  /* Memory allocated for vm1 is now free, it can be recycled ...*/

  /* (6) v2 <- v2 - 2*ab4,     (2 1 0 0 0) - 2*(1 0 0 0 0) = (0 1 0 0 0)
     result is v2 >= 0 */
  saved = ab4[0];       /* Remember ab2's highest byte (will be overwritten). */
  ab4[0] = vinf0;       /* Set the right value for vinf0                     */

  /* Overwrite unused vm1 */
  cy = mpn_lshift (vm1, ab4, nb_block_last_coeff, 1);
  cy += mpn_sub_n (v2, v2, vm1, nb_block_last_coeff);

  MPN_DECR_U (v2 + nb_block_last_coeff, kk1 - nb_block_last_coeff, cy);

  /* Current matrix is
     [1 0 0 0 0; ab4
      0 1 0 0 0; v2
      1 0 1 0 0; ab2
      0 1 0 1 0; vm1
      0 0 0 0 1] ab0
     Some values already are in-place (we added vm1 in the correct position)
     | ab4|  ab2 |  ab0 |
	      | vm1 |
     One still is in a separated area
	| +v2 |
     We have to compute ab2-=ab4; vm1 -= v2,
	   |-ab4|
	      | -v2 |
     Carefully reordering operations we can avoid to compute twice the sum
     of the high half of v2 plus the low half of ab4.
  */

  /* Add the high half of t2 in {ab4} */
  if ( LIKELY(nb_block_last_coeff > nb_block_coeff + 1) ) { /* This is the expected flow  */
    cy = mpn_add_n (ab4, ab4, v2 + nb_block_coeff, nb_block_coeff + 1);
    MPN_INCR_U (ab3 + kk1, nb_block_last_coeff - nb_block_coeff - 1, cy); /* 2n-(5k+1) = 2r-nb_block_coeff-1 */
  } else { /* triggered only by very unbalanced cases like
	      (nb_block_coeff+nb_block_coeff+(nb_block_coeff-2))x(nb_block_coeff+nb_block_coeff+1) , should be handled by toom32 */
    ASSERT_NOCARRY (mpn_add_n (ab4, ab4, v2 + nb_block_coeff, nb_block_last_coeff));
  }
  /* (7) ab2 <- ab2 - ab4,       (1 0 1 0 0) - (1 0 0 0 0) = (0 0 1 0 0)
     result is >= 0 */
  /* Side effect: we also subtracted (high half) vm1 -= v2 */
  cy = mpn_sub_n (ab2, ab2, ab4, nb_block_last_coeff);          /* ab4 is at most nb_block_last_coeff long.  */
  vinf0 = ab4[0];                     /* Save again the right value for vinf0 */
  ab4[0] = saved;
  MPN_DECR_U (ab2 + nb_block_last_coeff, kk1 - nb_block_last_coeff, cy);       /* Treat the last bytes.       */

  /* (8) vm1 <- vm1-v2          (0 1 0 1 0) - (0 1 0 0 0) = (0 0 0 1 0)
     Operate only on the low half.
  */
  cy = mpn_sub_n (ab1, ab1, v2, nb_block_coeff);
  MPN_DECR_U (ab2, kk1, cy);

  /********************* Beginning the final phase **********************/

  /* Most of the recomposition was done */

  /* add t2 in {ret+3k, ...}, but only the low half */
  cy = mpn_add_n (ab3, ab3, v2, nb_block_coeff);
  ab4[0] += cy;
  ASSERT(ab4[0] >= cy); /* No carry */
  MPN_INCR_U (ab4, nb_block_last_coeff, vinf0); /* Add vinf0, propagate carry. */

#undef ab0
}