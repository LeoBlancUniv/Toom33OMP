#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <gmp.h>
//#include "toom33_mul.c"

#define NTEST 501
#define NSAMPLES 1001

#define UNUSED(X) (void)(X)

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

const int size = 8192;
int toom3_block_size;
mpz_t powl;

void quicksort(uint64_t* t, int n)
{
	if (n > 0)
	{
	/* partitionning */
	int i, j, temp;
	
	j=0;
	for (i=0; i<n-1; i++)
	{
		/* at least as big as the pivot */
		if (t[i] < t[n-1])
		{
		temp = t[j];
		t[j] = t[i];
		t[i] = temp;
		j++;
		}
	}
	/*replacing the pivot */
	temp = t[j];
	t[j] = t[n-1];
	t[n-1] = temp;
	
	quicksort(t, j);
	quicksort(&t[j+1], n-j-1);
	}
}

uint64_t *quartiles(uint64_t *tab, uint64_t size)
{
	uint64_t *result ;
	uint64_t aux ;
	
	result = malloc(3*sizeof(uint64_t));
	quicksort(tab, size);
	aux = size >> 2;
	if (size % 4) aux++;
	// Q1
	result[0] = tab[aux-1];
	// Mediane
	// size is odd hence it's easy
	result[1]	= tab[(size+1)/2 - 1];
	// Q3
	aux = (3*size) >> 2;
	if ((3*size) % 4) aux++;
	result[2]	= tab[aux - 1];
	
	return result;
}

static inline uint64_t cpucyclesStart(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n    "
			"RDTSC\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

static inline uint64_t cpucyclesStop(void) {
	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			"CPUID\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

void toPoly(mpz_t x, mpz_t p[3]){

	

	for (int i = 0; i < 3; i++){
		mpz_init(p[i]);

		mpz_tdiv_r_2exp(p[i], x, toom3_block_size);

		mpz_tdiv_q_2exp(x, x, toom3_block_size);

	}

}

void clearPoly3(mpz_t* p){
	for (int i = 0; i < 3; i++){
		mpz_clear(p[i]);
	}
}

void clearPoly5(mpz_t* p){
	for (int i = 0; i < 5; i++){
		mpz_clear(p[i]);
	}
}

void eval(mpz_t p[3], mpz_t x, mpz_t ret){
	mpz_init(ret);

	for (int i = 0; i < 3; i++){
		mpz_t xpowi;
		mpz_init(xpowi);

		mpz_pow_ui(xpowi, x, i);

		mpz_addmul(ret, p[i], xpowi);

		mpz_clear(xpowi);
	}

}



void eval5(mpz_t p[5], mpz_t x, mpz_t ret){
	mpz_init(ret);

	for (int i = 0; i < 5; i++){
		mpz_t xpowi;
		mpz_init(xpowi);

		mpz_pow_ui(xpowi, x, i);

		mpz_addmul(ret, p[i], xpowi);

		mpz_clear(xpowi);
	}

}

void points(mpz_t p[3], mpz_t ps[5]){
	mpz_inits(ps[0], ps[1], ps[2], ps[3], ps[4], 0);
	mpz_set(ps[0], p[2]);


	mpz_t two, mone, one, zero;
	mpz_inits(two, mone, one, zero, 0);
	mpz_set_ui(two, 2);
	mpz_set_si(mone, -1);
	mpz_set_ui(one, 1);
	mpz_set_ui(zero, 0);


	eval(p, two, ps[1]);
	eval(p, mone, ps[2]);
	eval(p, one, ps[3]);
	eval(p, zero, ps[4]);

	mpz_clears(two, mone, one, zero, 0);

}

void clearPoints(mpz_t ps[5]){
	for (int i = 0; i < 5; i++)
	{
		mpz_clear(ps[i]);
	}
}

void calcAB(mpz_t ABpts[5], mpz_t AB[5]){
	mpz_inits(AB[0], AB[1], AB[2], AB[3], AB[4], 0);

	mpz_t tmp, tmp2, tmp3;
	mpz_t diff1, diff2, diff3;
	mpz_inits(tmp, tmp2, tmp3, 0);
	mpz_inits(diff1, diff2, diff3, 0);

	mpz_set(AB[0], ABpts[4]);
	mpz_set(tmp, ABpts[4]);

	//gmp_printf("AB0 : \n %Zd \n %Zd \n\n", AB[0], tmp);
	mpz_set_ui(tmp, 0);
	mpz_set_ui(tmp2, 0);
	mpz_set_ui(tmp3, 0);


	mpz_addmul_ui(AB[1], ABpts[0], 12);
	mpz_submul_ui(AB[1], ABpts[1], 1);
	mpz_submul_ui(AB[1], ABpts[2], 2);
	mpz_addmul_ui(AB[1], ABpts[3], 6);
	mpz_submul_ui(AB[1], ABpts[4], 3);
	mpz_tdiv_q_ui(AB[1], AB[1], 6);



	mpz_addmul_ui(tmp2, ABpts[0], 2);
	mpz_add(tmp2, tmp2, ABpts[3]); //add 2P(inf) and P(1)

	mpz_set(tmp3, ABpts[4]);
	mpz_tdiv_q_ui(tmp3, tmp3, 2);

	mpz_sub(tmp2, tmp2, tmp3); //sub 1/2P(0)

	mpz_set_ui(tmp3, 0);

	mpz_sub(tmp3, tmp3, ABpts[1]);
	mpz_tdiv_q_ui(tmp3, tmp3, 2);
	mpz_sub(tmp3, tmp3, ABpts[2]);
	mpz_tdiv_q_ui(tmp3, tmp3, 3);

	mpz_add(tmp2, tmp2, tmp3); //add (-P(2)/2 -B)/P(-1)

	if (!mpz_divisible_ui_p(ABpts[4], 2)){
		mpz_sub_ui(tmp2, tmp2, 1); //sub 1 if P(0) is odd
	}

	mpz_sub(diff1, AB[1], tmp2);


	//gmp_printf("tmp2 : \n %Zd \n AB1 : %Zd \n\n", tmp2, AB[1]);
	mpz_set_ui(tmp, 0);
	mpz_set_ui(tmp2, 0);
	mpz_set_ui(tmp3, 0);







	mpz_submul_ui(AB[2], ABpts[0], 2);
	mpz_addmul_ui(AB[2], ABpts[2], 1);
	mpz_addmul_ui(AB[2], ABpts[3], 1);
	mpz_submul_ui(AB[2], ABpts[4], 2);
	mpz_tdiv_q_ui(AB[2], AB[2], 2);

	

	mpz_sub(tmp2, tmp2, ABpts[0]);
	mpz_sub(tmp2, tmp2, ABpts[4]);
	mpz_add(tmp3, tmp3, ABpts[2]);
	mpz_add(tmp3, tmp3, ABpts[3]);
	mpz_tdiv_q_ui(tmp3, tmp3, 2);
	mpz_add(tmp2, tmp2, tmp3);

	mpz_sub(diff2, AB[2], tmp2);

	//gmp_printf("AB2 : \n %Zd \n %Zd \n %Zd \n\n", AB[2], tmp, tmp2);
	mpz_set_ui(tmp, 0);
	mpz_set_ui(tmp2, 0);
	mpz_set_ui(tmp3, 0);







	mpz_submul_ui(AB[3], ABpts[0], 12);
	mpz_addmul_ui(AB[3], ABpts[1], 1);
	mpz_submul_ui(AB[3], ABpts[2], 1);
	mpz_submul_ui(AB[3], ABpts[3], 3);
	mpz_addmul_ui(AB[3], ABpts[4], 3);
	mpz_tdiv_q_ui(AB[3], AB[3], 6);

	mpz_submul_ui(tmp2, ABpts[0], 2);

	mpz_set(tmp3, ABpts[1]);
	mpz_sub(tmp3, tmp3, ABpts[2]);

	mpz_tdiv_q_ui(tmp3, tmp3, 6);

	mpz_add(tmp2, tmp2, tmp3);

	mpz_set(tmp3, ABpts[4]);
	mpz_sub(tmp3, tmp3, ABpts[3]);

	mpz_tdiv_q_ui(tmp3, tmp3, 2);

	mpz_add(tmp2, tmp2, tmp3);

	mpz_sub(diff3, AB[3], tmp2);

	//gmp_printf("AB3 : \n %Zd \n %Zd \n %Zd \n\n", AB[3], tmp, tmp2);
	mpz_set_ui(tmp, 0);
	mpz_set_ui(tmp2, 0);
	mpz_set_ui(tmp3, 0);



	mpz_set(AB[4], ABpts[0]);
	mpz_set(tmp, ABpts[0]);

	//gmp_printf("AB4 : \n %Zd \n %Zd \n\n", AB[4], tmp);
	mpz_set_ui(tmp, 0);
	mpz_set_ui(tmp2, 0);
	mpz_set_ui(tmp3, 0);

	//gmp_printf("diff : %Zd %Zd %Zd\n", diff1, diff2, diff3);

	mpz_clears(tmp, tmp2, tmp3, 0);
	mpz_clears(diff1, diff2, diff3, 0);

}

void toom3_mpn(mp_limb_t* a, mp_limb_t* b, int nb_limbs, mp_limb_t* ab, mp_limb_t* scratch){
	
	int nb_block_coeff = (nb_limbs + 2) /  3;
	int nb_block_last_coeff = nb_limbs - 2 * nb_block_coeff;

	#define a0 a
	#define a1 a + nb_block_coeff
	#define a2 a + 2*nb_block_coeff

	#define b0 b
	#define b1 b + nb_block_coeff
	#define b2 b + 2*nb_block_coeff

	#define bpts_offset  				nb_block_coeff * 3 + 3
	#define abpts_offset bpts_offset  + nb_block_coeff * 3 + 3
	#define ab_offset 	 abpts_offset + 10*nb_block_coeff + 6
	#define aux_offset   ab_offset 	  + 10*nb_block_coeff + 6




	#define Apts_inf  a2 									 //	nb_block_last_coeff
	#define Apts_zero a0 									 // nb_block_coeff
	mp_limb_t* Apts_one =  scratch; 						 // nb_block_coeff + 1
	mp_limb_t* Apts_mone = scratch + nb_block_coeff + 1; 	 // nb_block_coeff + 1
	mp_limb_t* Apts_two =  scratch + nb_block_coeff * 2 + 2; // nb_block_coeff + 1
	
	
	

	#define Bpts_inf  b2 									 				// nb_block_last_coeff
	#define Bpts_zero b0 													// nb_block_coeff
	mp_limb_t* Bpts_one =  scratch + bpts_offset;							// nb_block_coeff + 1
	mp_limb_t* Bpts_mone = scratch + bpts_offset + nb_block_coeff + 1; 		// nb_block_coeff + 1
	mp_limb_t* Bpts_two =  scratch + bpts_offset + nb_block_coeff * 2 + 2;	// nb_block_coeff + 1

	int Apts_mone_sign = 1;
	int Bpts_mone_sign = 1;

	int Apts_cmp_skip = 0;
	int Bpts_cmp_skip = 0;

	mp_limb_t* ABpts_inf =  scratch + abpts_offset; 							// nb_block_coeff * 2
	mp_limb_t* ABpts_zero = scratch + abpts_offset + nb_block_coeff * 2; 		// nb_block_coeff * 2
	mp_limb_t* ABpts_one =  scratch + abpts_offset + nb_block_coeff * 4; 		// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_mone = scratch + abpts_offset + nb_block_coeff * 6 + 2; 	// nb_block_coeff * 2 + 1?
	mp_limb_t* ABpts_two =  scratch + abpts_offset + nb_block_coeff * 8 + 4; 	// nb_block_coeff * 2 + 2

	int ABpts_mone_sign;

	mp_limb_t* ab0 = scratch + ab_offset;							// nb_block_coeff * 2
	mp_limb_t* ab1 = scratch + ab_offset + nb_block_coeff * 2;		// nb_block_coeff * 2 + 2
	mp_limb_t* ab2 = scratch + ab_offset + nb_block_coeff * 4 + 2; 	// nb_block_coeff * 2 + 2
	mp_limb_t* ab3 = scratch + ab_offset + nb_block_coeff * 6 + 4;	// nb_block_coeff * 2 + 2
	mp_limb_t* ab4 = scratch + ab_offset + nb_block_coeff * 8 + 6; 	// nb_block_coeff * 2

	mp_limb_t* aux_inter_6 = scratch + aux_offset;		// nb_block_coeff * 2
	//mp_limb_t* ab3_aux = scratch + aux_offset;


	printf("mpn version\n");


	printf("%d %d %d\n",nb_limbs, nb_block_coeff, nb_block_last_coeff);

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





	/*
		Bpts calc, B is always positive so B0, B1, B2 are also always positive
		
		Bpts_mone is the only one that can be negative, Bpts_mone_sign will track it's sign

		1 : pos, 0 : neg
	*/

	//Bpts_inf <- b2 :

	mpn_copyd(Bpts_inf, b2, nb_block_last_coeff);



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





	//ABpts calc

	mpn_mul_n(ABpts_inf, Apts_inf, Bpts_inf, nb_block_last_coeff);
	mpn_mul_n(ABpts_two, Apts_two, Bpts_two, nb_block_coeff + 1);

	mpn_mul_n(ABpts_mone, Apts_mone, Bpts_mone, nb_block_coeff + 1); //nb_block_coeff * 2 + 1
	ABpts_mone_sign = Apts_mone_sign == Bpts_mone_sign;

	mpn_mul_n(ABpts_one, Apts_one, Bpts_one, nb_block_coeff + 1); //nb_block_coeff * 2 + 1
	mpn_mul_n(ABpts_zero, Apts_zero, Bpts_zero, nb_block_coeff);





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





	// trying to recreate steps from interpolate (mpn version):

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

	//gmp_printf("(6) v2 : %Nd\n\n", ABpts_two, nb_block_coeff * 2 + 2);




	/*
		(7) v1 <- v1 - vinf : (1 0 1 0 0) - (1 0 0 0 0) = (0 0 1 0 0)
	*/

	mpn_sub(ab2, ABpts_one, nb_block_coeff * 2 + 1, ABpts_inf, nb_block_coeff * 2);

	//gmp_printf("(7) v1 : %Nd\n\n", ABpts_one, nb_block_coeff * 2 + 2);




	/*
		(8) vm1 <- vm1-v2 : (0 1 0 1 0) - (0 1 0 0 0) = (0 0 0 1 0)
	*/

	mpn_sub_n(ab1, ABpts_mone, ab3, nb_block_coeff * 2 + 1);

	//gmp_printf("(8) vm1 : %Nd\n\n", ABpts_mone, nb_block_coeff * 2 + 2);


	//AB calc

	//AB0 <- P(0)
	mpn_copyd(ab0, ABpts_zero, nb_block_coeff * 2);


	



	//AB4 <- P(inf)
	mpn_copyd(ab4, ABpts_inf, nb_block_coeff * 2);



	gmp_printf("AB0 : %Nx \n\n", ab0, 2 * nb_block_coeff);
	gmp_printf("AB1 : %Nx \n\n", ab1, 2 * nb_block_coeff + 1);
	gmp_printf("AB2 : %Nx \n\n", ab2, 2 * nb_block_coeff + 1);
	gmp_printf("AB3 : %Nx \n\n", ab3, 2 * nb_block_coeff + 1);
	gmp_printf("AB4 : %Nx \n\n", ab4, 2 * nb_block_coeff);

	gmp_printf("C    : %Nx \n\n", ab0, nb_block_coeff * 10 + 4);

	mp_limb_t* ab_ = calloc(nb_block_coeff * 10 + 4, sizeof(mp_limb_t));

	mp_limb_t* ret = calloc(nb_block_coeff * 10 + 4, sizeof(mp_limb_t));
	
	mp_limb_t* acc = calloc(nb_block_coeff * 10 + 4, sizeof(mp_limb_t));

	

	int size_table[5] = {nb_block_coeff * 2, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2};
	mp_limb_t* ab_table[5] = {ab0, ab1, ab2, ab3, ab4};

	#define current ab_table[i]
	#define current_size size_table[i]

	for (int i = 0; i < 5; i++){

		mpn_zero(acc, nb_block_coeff * 10 + 4);
		
		mpn_copyd(acc + 43 * i, current, current_size);

		

		/*for (int j = 0; j < i; j++){
			for (int k = 0; k < 43; k++){
				mpn_lshift(acc, acc, nb_block_coeff * 10 + 4, 63);
				mpn_lshift(acc, acc, nb_block_coeff * 10 + 4, 1);

			}
		}*/


		mpn_add_n(ret, ret, acc, nb_block_coeff * 10 + 4);

					

	}

	gmp_printf("C : %Nx \n\n", 	ret, nb_limbs * 2);

	free(acc);
	free(ret);
	free(ab_);

	//cheating version

	/*
		
	*/

	//mpn_zero(ab2, nb_block_coeff * 2 + 1);
	//mpn_zero(ab3, nb_block_coeff * 2 + 1);
	//mpn_zero(ab4, nb_block_coeff * 2);




	/*mp_limb_t carry = 0;

	int global_nb_shift_block = 0;

	int size_table[5] = {nb_block_coeff * 2, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2 + 1, nb_block_coeff * 2};
	mp_limb_t* ab_table[5] = {ab0, ab1, ab2, ab3, ab4};

	#define current ab_table[i]
	#define next ab_table[i+1]
	
	for (int i = 0; i < 4; i++){


		int current_size = size_table[i];
		int next_size = size_table[i+1];

		int nb_shift_block = 0;
		int nb_shift = 0;

		mp_limb_t mask = 1ULL << 63;



		
		

		

		while (current[current_size - nb_shift_block - 1] == 0){
			nb_shift_block++;
		}


		while (((current[current_size - nb_shift_block - 1] & mask) == 0)){
			mask >>= 1;
			nb_shift++;
		}

		
		


		int shift_block_offset = global_nb_shift_block + nb_shift_block;

		printf("%d %d %d %d\n", global_nb_shift_block, nb_shift_block, shift_block_offset, nb_shift);



		if (shift_block_offset){
			next -= shift_block_offset; 
			//next is now on the first avaiable block after current

			mpn_rshift(next, next, next_size + shift_block_offset, shift_block_offset * 64);
			//shifts all of next blocks to align them with current

		}

		if (nb_shift){

			gmp_printf("%p %Nx\n", next, next, 1);

			carry = mpn_rshift(next, next, next_size, nb_shift);



			current[current_size - 1] |= carry;
			//shift the nb_shift bits out of next into current to make them a continuous block

			next -= 1; // next starts one block earlier, merged with current's end block

		}

		global_nb_shift_block += nb_shift_block; //update global shift for next loop
		

		printf("\n");

	}
	printf("%d\n", nb_block_coeff * 10 + 4);
	gmp_printf("C       : %Nx \n\n", ab0, 10 * nb_block_coeff + 4);*/


	

}

void toom3(mpz_t a, mpz_t b, mpz_t ab){
	#define inf 0
	#define two 1
	#define mone 2
	#define one 3
	#define zero 4

	mpz_t A[3];
	mpz_t B[3];

	/*gmp_printf("A : : %Zd\n\n", a);
	gmp_printf("B : : %Zd\n\n", b);*/

	printf("\n \n");
	
	mp_limb_t* scratch = calloc(10000, sizeof(mp_limb_t));

	toom3_mpn(a->_mp_d, b->_mp_d, b->_mp_size, ab->_mp_d, scratch);

	free(scratch);

	printf("\n \n");

	printf("mpz version\n");

	//gmp_printf("A : %Zd\n\n", a);
	//gmp_printf("B : %Zd\n\n", b);

	toPoly(a, A);
	toPoly(b, B);

	

	for (int i = 0; i < 3; i++){
		//gmp_printf("A%d %d : %Zd\n\n", i, A[i]->_mp_size, A[i]);
	}

	for (int i = 0; i < 3; i++){
		//gmp_printf("B%d %d : %Zd\n\n", i, B[i]->_mp_size, B[i]);
	}

	mpz_t evalA, evalB;

	eval(A, powl, evalA);
	eval(B, powl, evalB);

	//gmp_printf("evalA : %Zd\n\n", evalA);
	//gmp_printf("evalB : %Zd\n\n", evalB);

	mpz_clear(evalA);
	mpz_clear(evalB);

	mpz_t Apts[5];
	mpz_t Bpts[5];

	points(A, Apts);
	points(B, Bpts);

	for (int i = 0; i < 5; i++){
		//gmp_printf("Apts%d : %Zx\n\n", i, Apts[i]);
	}
	for (int i = 0; i < 5; i++){
		//gmp_printf("Bpts%d : %Zx\n\n", i, Bpts[i]);
	}

	mpz_t ABpts[5];

	for (int i = 0; i < 5; i++)
	{
		mpz_init(ABpts[i]);
		mpz_mul(ABpts[i], Apts[i], Bpts[i]); //this should be changed to kara once it's done
		//gmp_printf("ABpts%d : %Zx\n\n", i, ABpts[i]);
	}

	mpz_t AB[5];
	calcAB(ABpts, AB);

	

	//trying to recreate steps from interpolate :

	/* 
		(1) v2 <- v2-vm1 < v2+|vm1| : (16 8 4 2 1) - (1 -1 1 -1  1) = (15 9 3  3  0)
		we need can do v2 <- v2 - |vm1|

		then

		v2 <- v2 / 3 : (15 9 3  3  0)/3 =  (5 3 1 1 0)
	*/

	mpz_t aux;
	mpz_init(aux);

	mpz_abs(aux, ABpts[mone]);

	mpz_sub(ABpts[two], ABpts[two], aux);

	mpz_tdiv_q_ui(ABpts[two], ABpts[two], 3);

	//test (1)

	mpz_set(aux, AB[1]);
	mpz_add(aux, aux, AB[2]);
	mpz_addmul_ui(aux, AB[3], 3);
	mpz_addmul_ui(aux, AB[4], 5);

	//gmp_printf("(1)tv2 : %Zd\n\n", 	aux);

	//gmp_printf("(1) v2 : %Zd\n\n", ABpts[two]);
	mpz_set_ui(aux, 0);




	/*
		(2) vm1 <- tm1 := (v1 - vm1) / 2 : [(1 1 1 1 1) - (1 -1 1 -1 1)] / 2 = (0  1 0  1 0)
		we can do v1 - |vm1|
	*/

	mpz_abs(aux, ABpts[mone]);

	mpz_sub(ABpts[mone], ABpts[one], aux);

	mpz_tdiv_q_ui(ABpts[mone], ABpts[mone], 2);

	//test (2)

	mpz_set(aux, AB[1]);
	mpz_add(aux, aux, AB[3]);

	//gmp_printf("(2)tvm1 : %Zd\n\n", 	aux);

	//gmp_printf("(2) vm1 : %Zd\n\n", ABpts[mone]);
	mpz_set_ui(aux, 0);




	/*
		(3) v1 <- t1 := v1 - v0 : (1 1 1 1 1) - (0 0 0 0 1) = (1 1 1 1 0)
	*/

	mpz_sub(ABpts[one], ABpts[one], ABpts[zero]);

	//test (3)

	mpz_set(aux, AB[1]);
	mpz_add(aux, aux, AB[2]);
	mpz_add(aux, aux, AB[3]);
	mpz_add(aux, aux, AB[4]);

	//gmp_printf("(3)tv1 : %Zd\n\n", 	aux);

	//gmp_printf("(3) v1 : %Zd\n\n", ABpts[one]);
	mpz_set_ui(aux, 0);




	/*
		(4) v2 <- t2 := ((v2-vm1)/3-t1)/2 = (v2-vm1-3*t1)/6 : [(5 3 1 1 0) - (1 1 1 1 0)]/2 = (2 1 0 0 0)
		(v2 - vm1)/3 is v2 currently
		t1 is v1 currently
		so we do v2 <- (v2 - v1) / 2
	*/

	mpz_sub(ABpts[two], ABpts[two], ABpts[one]);
	mpz_tdiv_q_ui(ABpts[two], ABpts[two], 2);

	//test (4)

	mpz_set(aux, AB[3]);
	mpz_addmul_ui(aux, AB[4], 2);

	//gmp_printf("(4)tv2 : %Zd\n\n", 	aux);

	//gmp_printf("(4) v2 : %Zd\n\n", ABpts[two]);
	mpz_set_ui(aux, 0);




	/*
		(5) v1 <- t1-tm1 : (1 1 1 1 0) - (0 1 0 1 0) = (1 0 1 0 0)
		t1 is v1 currently
		tm1 is vm1 currently
		so we do v1 <- v1 - vm1
	*/

	mpz_sub(ABpts[one], ABpts[one], ABpts[mone]);

	//test (5)

	mpz_set(aux, AB[2]);
	mpz_add(aux, aux, AB[4]);
	//gmp_printf("(5)tv1 : %Zd\n\n", 	aux);

	//gmp_printf("(5) v1 : %Zd\n\n", ABpts[one]);
	mpz_set_ui(aux, 0);



	/*
		(6) v2 <- v2 - 2*vinf : (2 1 0 0 0) - 2*(1 0 0 0 0) = (0 1 0 0 0)
	*/

	mpz_submul_ui(ABpts[two], ABpts[inf], 2);

	//test (6)

	mpz_set(aux, AB[3]);
	//gmp_printf("(6)tv2 : %Zd\n\n", 	aux);

	//gmp_printf("(6) v2 : %Zd\n\n", ABpts[two]);
	mpz_set_ui(aux, 0);




	/*
		(7) v1 <- v1 - vinf : (1 0 1 0 0) - (1 0 0 0 0) = (0 0 1 0 0)
	*/

	mpz_sub(ABpts[one], ABpts[one], ABpts[inf]);

	//test (7)

	mpz_set(aux, AB[2]);
	//gmp_printf("(7)tv1 : %Zd\n\n", 	aux);

	//gmp_printf("(7) v1 : %Zd\n\n", ABpts[one]);
	mpz_set_ui(aux, 0);




	/*
		(8) vm1 <- vm1-v2 : (0 1 0 1 0) - (0 1 0 0 0) = (0 0 0 1 0)
	*/

	mpz_sub(ABpts[mone], ABpts[mone], ABpts[two]);

	mpz_set(aux, AB[1]);
	//gmp_printf("(8)tvm1 : %Zd\n\n", 	aux);

	//gmp_printf("(8) vm1 : %Zd\n\n", ABpts[mone]);
	mpz_set_ui(aux, 0);


	mpz_clear(aux);

	for (int i = 0; i < 5; ++i)
	{
		//gmp_printf("AB%d : %Zx\n\n", i, AB[i]);
	}

	eval5(AB, powl, ab);

	gmp_printf("C : %Zx\n\n", ab);


	clearPoly3(A);
	clearPoly3(B);
	clearPoly5(AB);

	clearPoints(Apts);
	clearPoints(Bpts);
	clearPoints(ABpts);

}

int main(int argc, char const *argv[])
{

	if (argc != 3){
		printf("use : ./main a b\n");
		return -1;
	}


	mpz_t a, b;

	mpz_init_set_str(a, argv[1], 10);
	mpz_init_set_str(b, argv[2], 10);

	mpz_abs(a, a);
	mpz_abs(b, b);

	//gmp_printf("a : %Zd\n\n", a);
	//gmp_printf("b : %Zd\n\n", b);


	

	mpz_init(powl);

	//toom3_block_size = (8192/3)+1;

	toom3_block_size = 86*32; //43*64

	mpz_ui_pow_ui(powl, 2, toom3_block_size);
	gmp_printf("%Zx\n\n", powl);

	mpz_t ab, ab_;

	mpz_init(ab);
	

	mpz_mul(ab, a, b);


		

	toom3(a, b, ab_);

	//gmp_printf("ab : %Zd\n\n", ab);


	//gmp_printf("ab_ : %Zd\n\n", ab_);

	/*mpz_t x, y, z;

	mpz_t expo;

	mpz_init(expo);

	mpz_ui_pow_ui(expo, 2, size-1);

	unsigned long seed = time(NULL);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);

	for (int i = 0; i < 1000; i++){
		mpz_inits(x, y, z, 0);

		mpz_urandomb(x, state, size-1);
		mpz_urandomb(y, state, size-1);

		mpz_add(x, x, expo);
		mpz_add(y, y, expo);




		toom3(x, y, z);

		mpz_clears(x, y, z, 0);
	}

	mpz_clear(expo);*/



	

	mpz_clear(a);
	mpz_clear(b);

	mpz_clear(powl);

	mpz_clear(ab);
	mpz_clear(ab_);

	//825071401847911350104434016317510780278067283201566754357130842744461195005906395199943718051043673223466730721751949974449248930827451854619414374697619260438961946911749222259780365386291751814038486373465100001426528819046274070279571501847524807503058829803292949581020312542244512310826381340725147604205504421489781205421741398803728752958064593794237603834697854344675329642780388241249095996086652843251798333913223849757378836474373313709743648650234857926624244596951216675142492484734159681684027468402207806205439865347969018612132747592285828382232580799515586663970404103616260646140491130336474241553520642020525516089660137732716330093125235058568988907834750945317769702575635573206226618491699434251943292774344986722857013298412092723076929179307523912203128114090534314456523272023692884978468594858168466441251073816641214891282204410271065824121943315600314597749227889573226269102256807459458734770833396972179130110520264774827679636735786580446623637598159822709594268849412087335528825773591001341459368429154598701743523091608321042347244119166311797831282266269359216738748446338689877553411067459022750494469035170847673664139376175016704738727113860401649828339318679003054736839277206870283192227831617539598383681659913506027472064659239184916246405287069652026544516606976680852350358887886206450202813665604368556089115665756191850387879202714813769550164656889996790902837958128457465941454508790612946593532755196596172284379325906458537834759016596544236993315228152015361031327422453823347007209987459925728347755977454611604649086058581101577344837613218863165735163293232182683983328352703931948934018941328642629738690331539423670935877924451188873671016602730573554563515446477428747122842429919747738097859733596295890237973718047460047408090635470616266779559457882525693610773409953808315167993488749629373478589355791418072020121684946874867922541121111773785091209812728948633266300346422009914729246105964646426829796911062733042949990258625136306447831184696881439764159871493270186815550008827682405281004158733425515557787170576476464934246505105516394418797543153382468729141414796208872453491051467183673152632331727482671336127896268657700582305093675403324150839387303020629188315382088652129003636647961202147261937023281222481237476524670557355735034515267765254502192363955049952594249904232853763728215835537753171258068562061271398805505638547424578150686558691852996512332614361427686993333226212556141099
	//1040724227190522569630755652129160827694143677356231111679162613310384812513826156920001337951731733252439451773942329711576331031075388181073205387002647911964870675318313779624546744891941923582662768200643923498087042155090747334264605591055719099054660179442985304183674856345222456203141013132783158026392261784546092663565284360810095080963414517973759428807308480748749405939264543234807629877959771010008500986026755612237214127065065682213247280231436672239382742108687805769246228790675623889650913518897320187281438414324171794721486035631943422988638315706516913755000989519985732988847636467288332068423945396568728774797004690730240937488833477950515961712397149827148532326268744543555414874562383957385212908921553443988347351243724205595393037428798028298795262033319686849523225316197092611669333529369527041131055806735322820589597691019490296965075074611949214671632253212442437287895749875048294779169292260278408281807600604044787786420843200490828998206663519657220986689098469015608296121555106687463466098507046333038459548527290675084624902018225330540015079208569239410860651729169600026538021093536101759269323191745744034127619151337392042938701581789047534566094111233826083128003120300724490166544545950082498896525316271652063113289130459104100910949515260533921618243600991898837745565419037787907192023975388902936573255341128088053177352079631846450634708443576941765118933866878403682378191342830102706575174538032007822713765608715084301489584481017441779354389943140962536827489130424381719446333456166080391539908818006693735322743003117149244960001644831801398197504218550494006074302257405572658535866738243049409014343138492166755614524621114462740338014386532142248008103530459299929759857087174009531577072856373904234103023908393413898530001516513598816961835269182444015792097128483864528151553600805230496890804015818146570100338461348221776040194599242394232236042591196200639620572280360991597776303554303777030005076450781316150291689845789273612793759884849769552756321274974171373506844634834324714247461365069063218212236649029391095026045486640845355130518837871259383539779725472099076767051496243018212383305637648419449797609209843065498424316352760733617983521631688238693766111166460703047175290778080866278130923400472373840994367371755395951869925957342405088746693753778792727622856685251783816873445207541785986821548839203566447870397667029305331773079406720033511085015551163060248885137965664152471611

	//546621712088986378140988503513700063756618322632475755909680051766442271373996947093447644294171452181949046770187525921922483109102314949207233522484198296782117880902700689707347336132164160263833884373799402438967265339569283294320356713621659685650555397275382693637272959633305537109722612345637993260273163610186416024645677555919407541789240373149301322220076951150265772904404210449418047180318821014543452743750869493593141708896971046376384099304566127862794040586162819897130504021627532549316102201016604389005230223076033170243721573135415155716688734549382721708934179191143388627159830136340413822753955106623966438069362055202674447223598658975680676765124064130235176204956272675265813345582869958053203500920698275262721978546674936469431520899437398551116152858699339354100737184285201748242022900277787202348179820315121101395176131713910274116607345502038118707916632260145954667678458134233991455496151761051482892726427552822776886180454554738764755866768039394184632990612147979805751504339437234196667702869993602542113318849629941246578496414682708408225243710314181326251817202356437996084667923572967399617743403423070965249642074113204713811539080997737289973242384084921718742897907332420681112924651487092228562609669059099715439466569470780438573768217740621879182261626411688479255732542805819234008138165340907783764782187113927809908340847020768273040349429107367520717817156018337520542337726398611339934574368899071146569425180094075730959592740302632322522335705861301754917651530432148485915951357669172886041911744070964935815983924869459736453201846794154084639097327196350887525302191143274537858752967325686224823016666554567470116504330780337437986962583420340571161041764447167460871683496816726293637238971446732639914589534052338272248206973194345867704218892457058677737610811459394339424020335840108221391708556395265786509755994864957583336020163352925570390640511994483368159560012315909090559136733248524233900704042019174358044978723468601693286181071680416690186080850932488256530356601104428255980019575583987994926507149775950965940110196646593786680383521590714662382660815671045552662105088133718418379004513998334265566275037326371521267516726018916309568554615424214268941302966151410781268377828995242130098509699353168987917434875089833763841992722626041396821522768863218147747135625636666260065293244125196663512156608923924853154603482014681412980756472031161796389573200826918558270910114612332794219
	//308871815054961287395965713975868099810479797434657009935421017179443561414976429706996986711551897466104440153581882403748883968420131532265143682878711150458656334889881646248461050593599404911768274976095641239046756171428540270188675990107084574338547259261340886551201378287010447969822788433643036241036650931445272390295283828976639226771635483259495839427193205279401856812091063385906386643601914048817305258819912702242993486980915533760072471179827512109249806636950033217889404288577883635859561474247655680814545473656347791476180782771864935545986185885565836994649995977212849105210844177725626919369459289924982281917992939746399086948424050213504645921210982411244873624632260645810758480994203263933858250001646693918321976090642446425910959382771117270519405292806821076428834039577264798810350321949500305591628751581943992753332001959401832599744616610979618214489685706089865454281894789673143503393295430942628931627212780085244806245165788416413709962561685302711995065159755732836642516268359486950576927395932701487968997157333507213244537902482128073280231300569649531466234353968570673399033017007070641716593918373955594784012907882550845297784139178276374793057398633121272358586812825648795100172283137849399371470054453309919768496717126233814803989983424333103182925952247874858105908110381299091653869797233208180486171563423687510831362906771093227014522730587169112628527834170413215577470394526227070566271239634498527558657062314567655121219104505640094787797790434389677305353850468640679815863023303383300978468780944864934921400153392533600213190241791234484893773888199299205298569205081594282105678121216287153926135564096592077325638335782592510562642832664648528430056409783559109179752967801006050042172150738225538907785886405776096115749903414116903524606544319648431378011800513135892726447972330386588980263004004612278813724290681754403342530677860780579244711202651026476782701322399070193634101429219162828009370349368187351697100092714638764079916849747867690669621891495378046836998628961171768380209284499318721286011810232411381786813082536471550463816747197731908409316376347533142717196348609466303253448802993493949006018300898892909760272629356915924616255824717183143875700369772993277792252070723643074057920793190144775542603103351951174164876298133603552765074753922369643965980022163990704913851029805702575600796639708992867890876281808008804870777085181438127405223045634601487736809239379315684100

	return 0;
}