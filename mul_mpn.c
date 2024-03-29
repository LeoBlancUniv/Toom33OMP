#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <gmp.h>

#include "mul_mpn.h"

void toom33_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, mp_limb_t* scratch, int para){
	/*
		gmp inspired toom cook 33 alg for multiplication
	

		the scratch variable will hold all the intermediate values requiered 
	*/
	
	//printf("mul start\n");

	int nb_block_coeff = (nb_limbs + 2) /  3;	//nb_block_coeff is the base unit of size mesurement
	int nb_block_last_coeff = nb_limbs - 2 * nb_block_coeff; //

	//for object of size 8192bits (nb_limbs = 128), they should be respectively 43 and 42

	//splitting of a without any big calculation
	#define a0 a
	#define a1 a + nb_block_coeff
	#define a2 a + 2*nb_block_coeff

	//splitting of b without any big calculation
	#define b0 b
	#define b1 b + nb_block_coeff
	#define b2 b + 2*nb_block_coeff

	//offset used in scratch 
	#define bpts_offset  				nb_block_coeff * 3 + 3
	#define abpts_offset bpts_offset  + nb_block_coeff * 3 + 3
	#define ab_offset 	 abpts_offset + 10*nb_block_coeff + 6
	#define aux_offset   ab_offset 	  + 10*nb_block_coeff + 6



	/*
		data pointer to all of the points calculation of A
	*/

	//A(inf) = last coeff of a, A(0) = const term of a 
	
	#define Apts_inf  a2 									 //	nb_block_last_coeff
	#define Apts_zero a0 									 // nb_block_coeff
							 
	mp_limb_t* Apts_one =  scratch; 						 // nb_block_coeff + 1
	mp_limb_t* Apts_mone = scratch + nb_block_coeff + 1; 	 // nb_block_coeff + 1
	mp_limb_t* Apts_two =  scratch + nb_block_coeff * 2 + 2; // nb_block_coeff + 1
	

	//version of above using distinct memory spaces

	/*
	mp_limb_t* Apts_inf =  calloc(nb_block_last_coeff, sizeof(mp_limb_t));// nb_block_last_coeff
	mp_limb_t* Apts_zero = calloc(nb_block_coeff, sizeof(mp_limb_t));	  // nb_block_coeff

	mpn_copyd(Apts_inf, a2, nb_block_last_coeff);
	mpn_copyd(Apts_zero, a0, nb_block_coeff);

	mp_limb_t* Apts_one =  calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	mp_limb_t* Apts_mone = calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	mp_limb_t* Apts_two =  calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	*/

	/*
		data pointer to all of the points calculation of B
	*/

	//B(inf) = last coeff of b, B(0) = const term of b 
	
	#define Bpts_inf  b2 									 				// nb_block_last_coeff
	#define Bpts_zero b0 													// nb_block_coeff
	
	mp_limb_t* Bpts_one =  scratch + bpts_offset;							// nb_block_coeff + 1
	mp_limb_t* Bpts_mone = scratch + bpts_offset + nb_block_coeff + 1; 		// nb_block_coeff + 1
	mp_limb_t* Bpts_two =  scratch + bpts_offset + nb_block_coeff * 2 + 2;	// nb_block_coeff + 1
	

	//version of above using distinct memory spaces
	/*
	mp_limb_t* Bpts_inf =  calloc(nb_block_last_coeff, sizeof(mp_limb_t));// nb_block_last_coeff
	mp_limb_t* Bpts_zero = calloc(nb_block_coeff, sizeof(mp_limb_t));	  // nb_block_coeff

	mpn_copyd(Bpts_inf, b2, nb_block_last_coeff);
	mpn_copyd(Bpts_zero, b0, nb_block_coeff);

	mp_limb_t* Bpts_one =  calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	mp_limb_t* Bpts_mone = calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	mp_limb_t* Bpts_two =  calloc(nb_block_coeff + 1, sizeof(mp_limb_t)); // nb_block_coeff + 1
	*/

	//the sign of Apts_mone and Bpts_mone are tracked on the side to only have to hold absolute value
	int Apts_mone_sign = 1;
	int Bpts_mone_sign = 1;

	//flag used to skip comparaison in cases where they cannot be done
	int Apts_cmp_skip = 0;
	int Bpts_cmp_skip = 0;


	/*
		data pointer to the points calculations of AB
	*/

	
	mp_limb_t* ABpts_inf =  scratch + abpts_offset; 							// nb_block_coeff * 2
	mp_limb_t* ABpts_zero = scratch + abpts_offset + nb_block_coeff * 2; 		// nb_block_coeff * 2
	mp_limb_t* ABpts_one =  scratch + abpts_offset + nb_block_coeff * 4; 		// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_mone = scratch + abpts_offset + nb_block_coeff * 6 + 2; 	// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_two =  scratch + abpts_offset + nb_block_coeff * 8 + 4; 	// nb_block_coeff * 2 + 2
	

	//version of above using distinct memory spaces
	/*
	mp_limb_t* ABpts_inf  = calloc(nb_block_coeff * 2, sizeof(mp_limb_t));		// nb_block_coeff * 2
	mp_limb_t* ABpts_zero = calloc(nb_block_coeff * 2, sizeof(mp_limb_t));		// nb_block_coeff * 2
	mp_limb_t* ABpts_one  = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_mone = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_two  = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2 + 2?
	*/

	//the sign of ABpts_mone and Bpts_mone are tracked on the side to only have to hold absolute value
	int ABpts_mone_sign;

	/*
		data pointer to the value of the coefficient of AB
	*/

	
	mp_limb_t* ab0 = scratch + ab_offset;							// nb_block_coeff * 2
	mp_limb_t* ab1 = scratch + ab_offset + nb_block_coeff * 2;		// nb_block_coeff * 2 + 2
	mp_limb_t* ab2 = scratch + ab_offset + nb_block_coeff * 4 + 2; 	// nb_block_coeff * 2 + 2
	mp_limb_t* ab3 = scratch + ab_offset + nb_block_coeff * 6 + 4;	// nb_block_coeff * 2 + 2
	mp_limb_t* ab4 = scratch + ab_offset + nb_block_coeff * 8 + 6; 	// nb_block_coeff * 2
	

	//version of above using distinct memory spaces
	/*
	mp_limb_t* ab0  = calloc(nb_block_coeff * 2, sizeof(mp_limb_t));		// nb_block_coeff * 2
	mp_limb_t* ab1  = calloc(nb_block_coeff * 2, sizeof(mp_limb_t));		// nb_block_coeff * 2 + 2
	mp_limb_t* ab2  = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2 + 2
	mp_limb_t* ab3  = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2 + 2
	mp_limb_t* ab4  = calloc(nb_block_coeff * 2 + 2, sizeof(mp_limb_t));	// nb_block_coeff * 2
	*/


	/*
		auxiliary variable used in some calculation
	*/
	
	mp_limb_t* aux_inter_6 = scratch + aux_offset;		// nb_block_coeff * 2
	

	//version of above using distinct memory spaces
	/*
	mp_limb_t* aux_inter_6 = calloc(nb_block_coeff * 2, sizeof(mp_limb_t));		// nb_block_coeff * 2
	*/

	//printf("done pointer setting\n");

	//printf("mpn version\n");

	/*if (!v)
		printf("%d %d %d\n",nb_limbs, nb_block_coeff, nb_block_last_coeff);*/

	//print a and b
	/*gmp_printf("A : %Nd \n\n", a, nb_limbs);
	gmp_printf("B : %Nd \n\n", b, nb_limbs);*/

	//print A0, A1, A2 
	/*gmp_printf("A0 : %Nd \n\n", a0, nb_block_coeff);
	gmp_printf("A1 : %Nd \n\n", a1, nb_block_coeff);
	gmp_printf("A2 : %Nd \n\n", a2, nb_block_last_coeff);

	//print B0, B1, B2
	gmp_printf("B0 : %Nd \n\n", b0, nb_block_coeff);
	gmp_printf("B1 : %Nd \n\n", b1, nb_block_coeff);
	gmp_printf("B2 : %Nd \n\n", b2, nb_block_last_coeff);*/


	/*
		Apts calc, A is always positive so A0, A1, A2 are also always positive
		
		Apts_mone is the only one that can be negative, Apts_mone_sign will track it's sign

		1 : pos, 0 : neg
	*/

	//Apts_inf <- a2 :


	mpn_copyd(Apts_inf, a2, nb_block_last_coeff);
	//should only be usefull if using the non scratch version, or you are just copying a2 into itself

	//Apts_two <- 4*a2 + 2*a1 + a0 : 

	//Apt_two <- a2
	mpn_copyd(Apts_two, a2, nb_block_last_coeff); //nb_block_last_coeff

	//Apts_two <- 2*a2
	if (mpn_lshift(Apts_two, Apts_two, nb_block_last_coeff, 1)){ //nb_block_last_coeff +1
		Apts_two[nb_block_last_coeff] = 1;
	}
	
	//Apts_two <- 2*a2 + a1
	if (mpn_add(Apts_two, a1, nb_block_coeff, Apts_two, nb_block_last_coeff + 1)){ //nb_block_coeff+1
		Apts_two[nb_block_coeff] = 1;
	}

	//Apts_two <- 4*a2 + 2*a1
	mpn_lshift(Apts_two, Apts_two, nb_block_coeff + 1, 1); // nb_block_coeff + 1
		
	//Apts_two <- 4*a2 + 2*a1 + a0
	mpn_add(Apts_two, Apts_two, nb_block_coeff + 1, a0, nb_block_coeff); // nb_block_coeff + 1



	//Apts_mone <- a2 + a0 - a1 :

	//Apts_mone <- a2 + a0, always positive
	if (mpn_add(Apts_mone, a0, nb_block_coeff, a2, nb_block_last_coeff)){ // nb_block_coeff + 1
		Apts_mone[nb_block_coeff] = 1;
		Apts_cmp_skip = 1; //if we have to overflow, (a2 + a0) is always bigger than a1
	}

	//Apts_mone <- (a2 + a0) - a1, can be negative if a1 is bigger than (a2 + a0)

	if (Apts_cmp_skip){ 
		//we know that (a2 + a0) is bigger than a1
		mpn_sub(Apts_mone, Apts_mone, nb_block_coeff + 1, a1, nb_block_coeff);
	}
	else{
		//we know that (a2 + a0) and a1 are the same size
		if (mpn_cmp(Apts_mone, a1, nb_block_coeff) >= 0){
			//(a2 + a0) is bigger than a1, so Apts_mone stays positive
			mpn_sub_n(Apts_mone, Apts_mone, a1, nb_block_coeff);
		}
		else{
			//(a2 + a0) is smaller than a1, so Apts_mone becomes negative
			mpn_sub_n(Apts_mone, a1, Apts_mone, nb_block_coeff);
			Apts_mone_sign = 0;
		}
	}



	//Apts_one <- a2 + a1 + a0 :

	//Apts_one <- a1 + a0
	if (mpn_add_n(Apts_one, a0, a1, nb_block_coeff)){
		Apts_one[nb_block_coeff] = 1;
	}

	//Apts_one <- a2 + a1 + a0
	mpn_add(Apts_one, Apts_one, nb_block_coeff + 1, a2, nb_block_last_coeff);


	
	//Apts_zero <- a0 :

	mpn_copyd(Apts_zero, a0, nb_block_coeff);
	//should only be usefull if using the non scratch version, or you are just copying a0 into itself



	/*gmp_printf("Apts0 : %Nx \n\n", Apts_inf,  nb_block_last_coeff);
	gmp_printf("Apts1 : %Nx \n\n", Apts_two,  nb_block_coeff + 1);

	if (Apts_mone_sign){
		gmp_printf("Apts2 : %Nx \n\n", Apts_mone, nb_block_coeff + 1);
	}
	else{
		gmp_printf("Apts2 : -%Nx \n\n", Apts_mone, nb_block_coeff + 1);
	}
	
	gmp_printf("Apts3 : %Nx \n\n", Apts_one,  nb_block_coeff + 1);
	gmp_printf("Apts4 : %Nx \n\n", Apts_zero, nb_block_coeff);*/


	//printf("done Apts\n");


	/*
		Bpts calc, B is always positive so B0, B1, B2 are also always positive
		
		Bpts_mone is the only one that can be negative, Bpts_mone_sign will track it's sign

		1 : pos, 0 : neg
	*/

	//Bpts_inf <- b2 :

	mpn_copyd(Bpts_inf, b2, nb_block_last_coeff);
	//should only be usefull if using the non scratch version, or you are just copying b2 into itself




	//Bpts_two <- 4*b2 + 2*b1 + b0 : 

	//Bpt_two <- b2
	mpn_copyd(Bpts_two, b2, nb_block_last_coeff); //nb_block_last_coeff

	//Bpts_two <- 2*b2
	if (mpn_lshift(Bpts_two, Bpts_two, nb_block_last_coeff, 1)){ //nb_block_last_coeff +1
		Bpts_two[nb_block_last_coeff] = 1;
	}
	
	//Bpts_two <- 2*b2 + b1
	if (mpn_add(Bpts_two, b1, nb_block_coeff, Bpts_two, nb_block_last_coeff + 1)){ //nb_block_coeff+1
		Bpts_two[nb_block_coeff] = 1;
	}

	//Bpts_two <- 4*b2 + 2*b1
	mpn_lshift(Bpts_two, Bpts_two, nb_block_coeff + 1, 1); // nb_block_coeff + 1
		
	//Bpts_two <- 4*b2 + 2*b1 + b0
	mpn_add(Bpts_two, Bpts_two, nb_block_coeff + 1, b0, nb_block_coeff); // nb_block_coeff + 1



	//Bpts_mone <- b2 + b0 - b1 :

	//Bpts_mone <- b2 + b0, always positive
	if (mpn_add(Bpts_mone, b0, nb_block_coeff, b2, nb_block_last_coeff)){ // nb_block_coeff + 1
		Bpts_mone[nb_block_coeff] = 1;
		Bpts_cmp_skip = 1; //if we have to overflow, (b2 + b0) is always bigger than b1
	}

	//Bpts_mone <- (b2 + b0) - b1, can be negative if b1 is bigger than (b2 + b0)

	if (Bpts_cmp_skip){ 
		//we know that (b2 + b0) is bigger than b1
		mpn_sub(Bpts_mone, Bpts_mone, nb_block_coeff + 1, b1, nb_block_coeff);
	}
	else{
		//we know that (b2 + b0) and b1 are the same size
		if (mpn_cmp(Bpts_mone, b1, nb_block_coeff) >= 0){
			//(b2 + b0) is bigger than b1, so Bpts_mone stays positive
			mpn_sub_n(Bpts_mone, Bpts_mone, b1, nb_block_coeff);
		}
		else{
			//(b2 + b0) is smaller than b1, so Bpts_mone becomes negative
			mpn_sub_n(Bpts_mone, b1, Bpts_mone, nb_block_coeff);
			Bpts_mone_sign = 0;
		}
	}



	//Bpts_one <- b2 + b1 + b0 :

	//Bpts_one <- b1 + b0
	if (mpn_add_n(Bpts_one, b0, b1, nb_block_coeff)){
		Bpts_one[nb_block_coeff] = 1;
	}

	//Bpts_one <- b2 + b1 + b0
	mpn_add(Bpts_one, Bpts_one, nb_block_coeff + 1, b2, nb_block_last_coeff);


	
	//Bpts_zero <- b0 :

	mpn_copyd(Bpts_zero, b0, nb_block_coeff);
	//should only be usefull if using the non scratch version, or you are just copying b0 into itself




	/*gmp_printf("Bpts0 : %Nx \n\n", Bpts_inf,  nb_block_last_coeff);
	gmp_printf("Bpts1 : %Nx \n\n", Bpts_two,  nb_block_coeff + 1);
	
	if (Bpts_mone_sign){
		gmp_printf("Bpts2 : %Nx \n\n", Bpts_mone, nb_block_coeff + 1);
	}
	else{
		gmp_printf("Bpts2 : -%Nx \n\n", Bpts_mone, nb_block_coeff + 1);
	}
	
	gmp_printf("Bpts3 : %Nx \n\n", Bpts_one,  nb_block_coeff + 1);
	gmp_printf("Bpts4 : %Nx \n\n", Bpts_zero, nb_block_coeff);*/


	//printf("done Bpts\n");


	//ABpts calc

	//parallel version using OpenMP with parallel for

	if (para){
		mp_limb_t* Apts_table[5] = {Apts_inf, Apts_two, Apts_mone, Apts_one, Apts_zero};
		mp_limb_t* Bpts_table[5] = {Bpts_inf, Bpts_two, Bpts_mone, Bpts_one, Bpts_zero};
		mp_limb_t* ABpts_table[5] = {ABpts_inf, ABpts_two, ABpts_mone, ABpts_one, ABpts_zero};

		int size_table[5] = {nb_block_last_coeff, nb_block_coeff + 1, nb_block_coeff +1, nb_block_coeff + 1, nb_block_coeff};
		//#pragma omp parallel sections firstprivate(a, b, ABpts_inf, ABpts_two, ABpts_mone, ABpts_one, ABpts_zero, Apts_one, Apts_mone, Apts_two, Apts_inf, Apts_zero, Bpts_one, Bpts_mone, Bpts_two, Bpts_inf, Bpts_zero, nb_block_coeff, nb_block_last_coeff)

		#pragma omp parallel for firstprivate(scratch, nb_block_coeff, nb_block_last_coeff, Apts_table, Bpts_table, ABpts_table, size_table)
		for (int i = 0; i < 5; i++){
			mpn_mul_n(ABpts_table[i], Apts_table[i], Bpts_table[i], size_table[i]);
		}

		ABpts_mone_sign = Apts_mone_sign == Bpts_mone_sign;
	}
	else{
		//sequential version, the calls to mpn_mul should call karatsuba (toom22) which should itelf call basecase
	
		mpn_mul_n(ABpts_inf, Apts_inf, Bpts_inf, nb_block_last_coeff);

		mpn_mul_n(ABpts_two, Apts_two, Bpts_two, nb_block_coeff + 1);

		mpn_mul_n(ABpts_mone, Apts_mone, Bpts_mone, nb_block_coeff + 1); //nb_block_coeff * 2 + 1
		ABpts_mone_sign = Apts_mone_sign == Bpts_mone_sign;

		mpn_mul_n(ABpts_one, Apts_one, Bpts_one, nb_block_coeff + 1); //nb_block_coeff * 2 + 1


		mpn_mul_n(ABpts_zero, Apts_zero, Bpts_zero, nb_block_coeff);

	}

	


	/*gmp_printf("ABpts0 : %Nx \n\n", ABpts_inf,  2 * nb_block_coeff);
	gmp_printf("ABpts1 : %Nx \n\n", ABpts_two,  2 * nb_block_coeff + 2); // ?

	if (ABpts_mone_sign){
		gmp_printf("ABpts2 : %Nx \n\n", ABpts_mone, 2 * nb_block_coeff + 1);//+1
	}
	else{
		gmp_printf("ABpts2 : -%Nx \n\n", ABpts_mone, 2 * nb_block_coeff + 1);//+1
	}
	
	gmp_printf("ABpts3 : %Nx \n\n", ABpts_one,  2 * nb_block_coeff + 1);//+1
	gmp_printf("ABpts4 : %Nx \n\n", ABpts_zero, 2 * nb_block_coeff);*/





	// trying to recreate steps from interpolate5 (mpn version):

	/* 
		(1) v2 <- v2-vm1 < v2+|vm1| : (16 8 4 2 1) - (1 -1 1 -1  1) = (15 9 3  3  0)
		we need can do v2 <- v2 - |vm1|

		then

		v2 <- v2 / 3 : (15 9 3  3  0)/3 =  (5 3 1 1 0)
	*/

	//should double check if this is right with the handling of sign being separate

	if (ABpts_mone_sign){ //pos
		mpn_sub_n(ABpts_two, ABpts_two, ABpts_mone, nb_block_coeff * 2 + 1);
	}
	else{ //neg
		mpn_add_n(ABpts_two, ABpts_two, ABpts_mone, nb_block_coeff * 2 + 1);
	}

	mpn_divexact_by3 (ABpts_two, ABpts_two, nb_block_coeff * 2 + 1);
	//since 15, 9, and 3 are all multiple of 3, this operation doesn't loose any information

	//test (1)

	//gmp_printf("(1) v2 : %Nd\n\n", ABpts_two, nb_block_coeff * 2 + 1);



	/*
		(2) vm1 <- tm1 := (v1 - vm1) / 2 : [(1 1 1 1 1) - (1 -1 1 -1 1)] / 2 = (0  1 0  1 0)
		we can do v1 - |vm1|
	*/

	//should double check if this is right with the handling of sign being separate

	if (ABpts_mone_sign){ //pos
		mpn_sub_n(ABpts_mone, ABpts_one, ABpts_mone, nb_block_coeff * 2 + 1);
	}
	else //neg
	{
		mpn_add_n(ABpts_mone, ABpts_one, ABpts_mone, nb_block_coeff * 2 + 1);
	}

	mpn_rshift(ABpts_mone, ABpts_mone, nb_block_coeff * 2 + 1, 1);

	//gmp_printf("(2) vm1 : %Nd\n\n", ABpts_mone, nb_block_coeff * 2 + 2);




	/*
		(3) v1 <- t1 := v1 - v0 : (1 1 1 1 1) - (0 0 0 0 1) = (1 1 1 1 0)
	*/

	mpn_sub(ABpts_one, ABpts_one, nb_block_coeff * 2 + 1, ABpts_zero, nb_block_coeff * 2);

	//gmp_printf("(3) v1 : %Nd\n\n", ABpts_one, nb_block_coeff * 2 + 2);




	/*
		(4) v2 <- t2 := ((v2-vm1)/3-t1)/2 = (v2-vm1-3*t1)/6 : [(5 3 1 1 0) - (1 1 1 1 0)]/2 = (2 1 0 0 0)
		(v2 - vm1)/3 is v2 currently
		t1 is v1 currently
		so we do v2 <- (v2 - v1) / 2
	*/

	mpn_sub_n(ABpts_two, ABpts_two, ABpts_one, nb_block_coeff * 2 + 1);
	mpn_rshift(ABpts_two, ABpts_two, nb_block_coeff * 2 + 1, 1);

	//gmp_printf("(4) v2 : %Nd\n\n", ABpts_two, nb_block_coeff * 2 + 2);




	/*
		(5) v1 <- t1-tm1 : (1 1 1 1 0) - (0 1 0 1 0) = (1 0 1 0 0)
		t1 is v1 currently
		tm1 is vm1 currently
		so we do v1 <- v1 - vm1
	*/

	mpn_sub_n(ABpts_one, ABpts_one, ABpts_mone, nb_block_coeff * 2 + 1);

	//gmp_printf("(5) v1 : %Nd\n\n", ABpts_one, nb_block_coeff * 2 + 2);



	
	/*
		(6) v2 <- v2 - 2*vinf : (2 1 0 0 0) - 2*(1 0 0 0 0) = (0 1 0 0 0)
	*/

	mpn_copyd(aux_inter_6, ABpts_inf, nb_block_coeff * 2);

	mpn_lshift(aux_inter_6, aux_inter_6, nb_block_coeff * 2, 1);

	mpn_sub(ab3, ABpts_two, nb_block_coeff * 2 + 1, aux_inter_6, nb_block_coeff * 2);
	//result is (0 1 0 0 0) so we can directly put it in ab3

	//gmp_printf("(6) v2 : %Nd\n\n", ABpts_two, nb_block_coeff * 2 + 2);




	/*
		(7) v1 <- v1 - vinf : (1 0 1 0 0) - (1 0 0 0 0) = (0 0 1 0 0)
	*/

	mpn_sub(ab2, ABpts_one, nb_block_coeff * 2 + 1, ABpts_inf, nb_block_coeff * 2);
	//result is (0 0 1 0 0) so we can directly put it in ab2


	//gmp_printf("(7) v1 : %Nd\n\n", ABpts_one, nb_block_coeff * 2 + 2);




	/*
		(8) vm1 <- vm1-v2 : (0 1 0 1 0) - (0 1 0 0 0) = (0 0 0 1 0)
	*/

	mpn_sub_n(ab1, ABpts_mone, ab3, nb_block_coeff * 2 + 1);
	//result is (0 0 0 1 0) so we can directly put it in ab1


	//gmp_printf("(8) vm1 : %Nd\n\n", ABpts_mone, nb_block_coeff * 2 + 2);


	//AB calc

	//AB0 <- P(0)
	mpn_copyd(ab0, ABpts_zero, nb_block_coeff * 2);


	//AB4 <- P(inf)
	mpn_copyd(ab4, ABpts_inf, nb_block_coeff * 2);



	/*gmp_printf("AB0 : %Nx \n\n", ab0, 2 * nb_block_coeff);
	gmp_printf("AB1 : %Nx \n\n", ab1, 2 * nb_block_coeff + 1);
	gmp_printf("AB2 : %Nx \n\n", ab2, 2 * nb_block_coeff + 1);
	gmp_printf("AB3 : %Nx \n\n", ab3, 2 * nb_block_coeff + 1);
	gmp_printf("AB4 : %Nx \n\n", ab4, 2 * nb_block_coeff);

	gmp_printf("C    : %Nx \n\n", ab0, nb_block_coeff * 10 + 4);*/

	/*
		mpz version (eval5)
	*/
	
	mp_limb_t* acc = calloc(nb_block_coeff * 10 + 4, sizeof(mp_limb_t));

	

	int size_table_[5] = {nb_block_coeff * 2, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2};
	mp_limb_t* ab_table[5] = {ab0, ab1, ab2, ab3, ab4};

	#define current ab_table[i]
	#define current_size size_table_[i]

	for (int i = 0; i < 5; i++){

		mpn_zero(acc, nb_block_coeff * 10 + 4);
		
		mpn_copyd(acc + 43 * i, current, current_size);

		

		/*for (int j = 0; j < i; j++){
			for (int k = 0; k < 43; k++){
				mpn_lshift(acc, acc, nb_block_coeff * 10 + 4, 63);
				mpn_lshift(acc, acc, nb_block_coeff * 10 + 4, 1);

			}
		}*/


		mpn_add_n(ab, ab, acc, nb_block_coeff * 10 + 4);

	}


	/*if (!v){
		gmp_printf("r : %Nx \n\n", 	ab, nb_limbs * 2);
	}*/
	

	free(acc);

	

	//free for non scratch version only
	/*free(Apts_one);
	free(Apts_mone);
	free(Apts_two);

	free(Bpts_one);
	free(Bpts_mone);
	free(Bpts_two);

	free(ABpts_inf);
	free(ABpts_zero);
	free(ABpts_one);
	free(ABpts_mone);
	free(ABpts_two);

	free(ab0);
	free(ab1);
	free(ab2);
	free(ab3);
	free(ab4);

	free(aux_inter_6);*/

	

}

void schoolbook22(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para){



	mp_limb_t accs[4][128]; // will hold intermediate values 

	int nb_limbs_half = nb_limbs / 2;

	if (para){
		#pragma omp parallel for
	    for (int i = 0; i < 4; i++){
	    	mpn_mul_n(accs[i], a + nb_limbs_half * (i&1), b + nb_limbs_half * ((i&2)>0), nb_limbs_half);
	    }
	}
	else{
		for (int i = 0; i < 4; i++){
	    	mpn_mul_n(accs[i], a + nb_limbs_half * (i&1), b + nb_limbs_half * ((i&2)>0), nb_limbs_half);
	    }
	}
	
    

  	
	mpn_copyd(ab, accs[0], nb_limbs);
	ab[nb_limbs + nb_limbs_half] += mpn_add_n(ab + nb_limbs_half, ab + nb_limbs_half, accs[1], nb_limbs);
	ab[nb_limbs + nb_limbs_half] += mpn_add_n(ab + nb_limbs_half, ab + nb_limbs_half, accs[2], nb_limbs);
	mpn_add_n(ab + nb_limbs, ab + nb_limbs, accs[3], nb_limbs);


}

void schoolbook32(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para){

	mp_limb_t accs[6][128]; // will hold intermediate values 

	int nb_limbs_half = nb_limbs / 2;

	int nb_block_coeff = (nb_limbs + 2) /  3;	
	int nb_block_last_coeff = nb_limbs - 2 * nb_block_coeff; 

	mp_limb_t* a_lo  = a;						// nb_block_coeff	   
	mp_limb_t* a_mid = a + nb_block_coeff;		// nb_block_coeff
	mp_limb_t* a_hi  = a + nb_limbs_half * 2;   // nb_block_last_coeff

	mp_limb_t* b_lo = b; 				    // nb_limbs_half
	mp_limb_t* b_hi = b + nb_limbs_half; 	// nb_limbs_half

	/*mpn_mul(accs[0], a_lo, nb_block_coeff, b_lo, nb_limbs_half);

	mpn_mul(accs[1], a_lo, nb_block_coeff, b_hi, nb_limbs_half);

	mpn_mul(accs[2], a_mid, nb_block_coeff, b_lo, nb_limbs_half);

	mpn_mul(accs[3], a_mid, nb_block_coeff, b_hi, nb_limbs_half);

	mpn_mul(accs[4], a_hi, nb_block_coeff, b_lo, nb_limbs_half);

	mpn_mul(accs[5], a_hi, nb_block_coeff, b_hi, nb_limbs_half);*/


	/*#pragma omp parallel for
	for (int i = 0; i < 6; i++)
	{
		mpn_mul_n
	}*/
	
	omp_set_num_threads(6);
	#pragma omp parallel sections
  	{
  		
		#pragma omp section
    	{
    		mpn_mul(accs[0], a_lo, nb_block_coeff, b_lo, nb_limbs_half);
    	}
    	#pragma omp section
    	{
			mpn_mul(accs[1], a_lo, nb_block_coeff, b_hi, nb_limbs_half);
    	}
    	#pragma omp section
    	{
    		mpn_mul(accs[2], a_mid, nb_block_coeff, b_lo, nb_limbs_half);
    	}
    	#pragma omp section
    	{
    		mpn_mul(accs[3], a_mid, nb_block_coeff, b_hi, nb_limbs_half);
    	}
    	#pragma omp section
    	{
    		mpn_mul(accs[4], a_hi, nb_block_last_coeff, b_lo, nb_limbs_half);
    	}
    	#pragma omp section
    	{
    		mpn_mul(accs[5], a_hi, nb_block_last_coeff, b_hi, nb_limbs_half);
    	}

    }

	mpn_copyd(ab, accs[0], nb_limbs_half + nb_block_coeff);
	
	ab[nb_limbs] += mpn_add_n(ab + nb_block_coeff, ab + nb_block_coeff, accs[1], nb_limbs);
	ab[nb_limbs] += mpn_add_n(ab + nb_block_coeff, ab + nb_block_coeff, accs[2], nb_limbs);
	
	ab[nb_limbs + 2 * nb_block_coeff] += mpn_add_n(ab + nb_limbs, ab + nb_limbs, accs[3], nb_limbs);
	ab[nb_limbs + 2 * nb_block_coeff] += mpn_add_n(ab + nb_limbs, ab + nb_limbs, accs[4], nb_limbs);
	mpn_add_n(ab + nb_limbs + 2 * nb_block_coeff, ab + nb_limbs + 2 * nb_block_coeff, accs[5], nb_limbs);

}

void schoolbooklow_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, int para){

	mp_limb_t accs[3][128];

	if (para){
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				mpn_mul_n(accs[0], a, b, nb_limbs / 2);
			}

			#pragma omp section
			{
				mpn_mullo_n(accs[1], a, b + nb_limbs / 2, nb_limbs / 2);
			}

			#pragma omp section
			{
				mpn_mullo_n(accs[2], a + nb_limbs / 2, b, nb_limbs / 2);
			}
		}
	}
	else{
		mpn_mul_n(accs[0], a, b, nb_limbs / 2);
		mpn_mullo_n(accs[1], a, b + nb_limbs / 2, nb_limbs / 2);
		mpn_mullo_n(accs[2], a + nb_limbs / 2, b, nb_limbs / 2);
	}


	mpn_copyd(ab, accs[0], nb_limbs);
	mpn_add_n(ab + nb_limbs / 2, ab + nb_limbs / 2, accs[1], nb_limbs / 2);
	mpn_add_n(ab + nb_limbs / 2, ab + nb_limbs / 2, accs[2], nb_limbs / 2);

}

void mul_mont_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, mp_limb_t* p_inv, int nb_limbs, mp_limb_t* ab, int para)
{
	mp_limb_t c[2 * nb_limbs];
	mp_limb_t q[nb_limbs];
	mp_limb_t q2[2 * nb_limbs];
	mp_limb_t r[2 * nb_limbs];
	mp_limb_t cy;
	
	schoolbook22(a, b, nb_limbs, c, para); //c = ab

	schoolbooklow_mpn(c, p_inv, nb_limbs, q, para); //q = c * p_inv mod 2^8192

	schoolbook22(q, p, nb_limbs, q2, para); //q2 = qp

	/*cy = mpn_sub_n (q2 + nb_limbs, q2, c, nb_limbs);
	mpn_sub(q2 + nb_limbs, q2 + nb_limbs, nb_limbs, &cy, 1);*/

	cy = mpn_sub_n(r, c, q2, nb_limbs * 2); //r = c - q2

	if (cy != 0){
    	mpn_add_n (r, r, p, nb_limbs);
	}

	mpn_copyd(ab, r + nb_limbs, nb_limbs);

	

}
	

void mul_redc_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, mp_limb_t* p_inv, int nb_limbs, mp_limb_t* ab, int para){
	mp_limb_t c[2 * nb_limbs];
	

	//mpn_mul_n(c, a, b, nb_limbs); //to change to other mul later
	schoolbook22(a, b, nb_limbs, c, para);
	
	mp_ptr x, y;
	mp_ptr scratch;
  	mp_limb_t cy; //carry
 	mp_size_t rn, rn_4; //should be nb_limbs and nb_limbs/2
  	TMP_DECL;
  	TMP_MARK;
  	
  	rn = mpn_mulmod_bnm1_next_size (nb_limbs); //128
  	rn_4 = mpn_mulmod_bnm1_next_size (nb_limbs / 2); //64

  	//printf("%d %d\n\n", rn, rn_4);

  	// x(nb_limbs) + y(rn) + y0(rn_4) + y1(rn_4) + y2(rn_4) + 3(rn_4) + 4 * scratch for bnm1(rn, nb_limbs/2)
  	scratch = TMP_ALLOC_LIMBS (nb_limbs + rn + rn_4 * 4 + 
  		4 * mpn_mulmod_bnm1_itch (rn, nb_limbs / 2, nb_limbs / 2));

  	int scratch_s = mpn_mulmod_bnm1_itch (rn, nb_limbs / 2, nb_limbs / 2); //size of one scratch space

  	x = scratch; //scratch
  	y = scratch + nb_limbs; //scratch + nb_limbs
  	mp_limb_t* y_ = y + rn; //scratch + nb_limbs + rn
  	mp_limb_t* scratch_ = y_ + rn_4 * 4; //scratch + nb_limbs + rn + 4 * rn_4

  	schoolbooklow_mpn(c, p_inv, nb_limbs, x, para);


  	//split of mpn_mulmod_bnm1
  	
  	if (para){
  		#pragma omp parallel for
  		for (int i = 0; i < 4; i++){
	    	mpn_mulmod_bnm1(y_ + rn_4 * i, rn_4, 
    		x + (nb_limbs / 2) * (i&1), nb_limbs / 2,
    		p + (nb_limbs / 2) * ((i&2)>0), nb_limbs / 2, 
    		scratch_ + (scratch_s) * i);
	    }
  	}
  	else{
  		for (int i = 0; i < 4; i++){
	    	mpn_mulmod_bnm1(y_ + rn_4 * i, rn_4, 
    		x + (nb_limbs / 2) * (i&1), nb_limbs / 2,
    		p + (nb_limbs / 2) * ((i&2)>0), nb_limbs / 2, 
    		scratch_ + (scratch_s) * i);
	    }
  	}

    
    

    //y reconstruction
    mpn_copyd(y, y_, rn_4);
	y[rn + rn_4] += mpn_add_n(y + rn_4, y + rn_4, y_ + rn_4, rn);
	y[rn + rn_4] += mpn_add_n(y + rn_4, y + rn_4, y_ + (rn_4 * 2), rn);
	mpn_add_n(y + rn, y + rn, y_ + (rn_4 * 3), rn);

	//final step
	cy = mpn_sub_n (y + rn, y, c, 2*nb_limbs - rn);
	//MPN_DECR_U (y + 2*nb_limbs - rn, rn, cy);
	mpn_sub(y + 2*nb_limbs - rn, y + 2*nb_limbs - rn, rn, &cy, 1);

	cy = mpn_sub_n (ab, c + nb_limbs, y + nb_limbs, nb_limbs);
	if (cy != 0)
		mpn_add_n (ab, ab, p, nb_limbs);


  	TMP_FREE;




}

void check_mul_mpn(mp_limb_t* a, mp_limb_t* b, mp_limb_t* p, int nb_limbs){

	mp_limb_t* c = calloc(nb_limbs * 2, sizeof(mp_limb_t));
	mp_limb_t* q_limbs = (mp_limb_t*) calloc ((nb_limbs+1), sizeof(mp_limb_t));
	mp_limb_t* r_limbs = (mp_limb_t*) calloc ((nb_limbs), sizeof(mp_limb_t));

	mp_limb_t* p_inv = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	mp_limb_t* scratch_p_inv = calloc(nb_limbs * 2, sizeof(mp_limb_t));

	mp_limb_t* c_redc = calloc(nb_limbs, sizeof(mp_limb_t));

	mpn_binvert(p_inv, p, nb_limbs, scratch_p_inv);
	free(scratch_p_inv);

	mpn_mul_n(c, a, b, nb_limbs);
	mpn_tdiv_qr(q_limbs, r_limbs, 0, c, nb_limbs * 2, p, nb_limbs);

	mpn_redc_n(c_redc, c, p, nb_limbs, p_inv);




	gmp_printf("c : %Nx\n\n", c, nb_limbs * 2);
	gmp_printf("c_lo : %Nx\n\n", c, nb_limbs);
	gmp_printf("r : %Nx\n\n", r_limbs, nb_limbs);
	gmp_printf("redc : %Nx\n\n", c_redc, nb_limbs);

	mp_limb_t* c_mullo = calloc(nb_limbs * 2, sizeof(mp_limb_t));
	schoolbooklow_mpn(a, b, nb_limbs, c_mullo, 0);
	gmp_printf("clo : %Nx\n\n", c_mullo, nb_limbs);
	free(c_mullo);

	mp_limb_t* c_mul_redc_mpn = calloc(nb_limbs * 2, sizeof(mp_limb_t));

	mul_redc_mpn(a, b, p, p_inv, nb_limbs, c_mul_redc_mpn, 0);
	//gmp_printf("c : %Nx\n\n", c_mul_redc_mpn, 2 * nb_limbs);
	free(c_mul_redc_mpn);

	mp_limb_t* c_mul_mont_mpn = calloc(nb_limbs * 2, sizeof(mp_limb_t));

	mul_mont_mpn(a, b, p, p_inv, nb_limbs, c_mul_mont_mpn, 0);
	gmp_printf("redc : %Nx\n\n", c_mul_mont_mpn, nb_limbs);
	free(c_mul_mont_mpn);



	free(c_redc);
	free(r_limbs);
	free(q_limbs);
	free(p_inv);
	free(c);

}