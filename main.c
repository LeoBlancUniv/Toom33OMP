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

	#define bpts_offset 5*nb_block_coeff + 3
	#define abpts_offset bpts_offset + 5*nb_block_coeff + 3
	#define ab_offset abpts_offset + 10*nb_block_coeff + 6
	#define aux_offset ab_offset + 15*nb_block_coeff

	mp_limb_t* Apts_inf =  scratch; 						 // nb_block_last_coeff 
	mp_limb_t* Apts_two =  scratch + nb_block_coeff * 4 + 2; // nb_block_coeff + 1
	mp_limb_t* Apts_mone = scratch + nb_block_coeff * 2; 	 // nb_block_coeff + 1
	mp_limb_t* Apts_one =  scratch + nb_block_coeff * 3 + 1; // nb_block_coeff + 1
	mp_limb_t* Apts_zero = scratch + nb_block_coeff; 		 // nb_block_coeff

	mp_limb_t* Bpts_inf =  scratch + bpts_offset; 							// nb_block_last_coeff
	mp_limb_t* Bpts_two =  scratch + bpts_offset + nb_block_coeff * 4 + 2; // nb_block_coeff + 1
	mp_limb_t* Bpts_mone = scratch + bpts_offset + nb_block_coeff * 2; 	// nb_block_coeff + 1
	mp_limb_t* Bpts_one =  scratch + bpts_offset + nb_block_coeff * 3 + 1; // nb_block_coeff + 1
	mp_limb_t* Bpts_zero = scratch + bpts_offset + nb_block_coeff; 		// nb_block_coeff

	int Apts_mone_sign = 1;
	int Bpts_mone_sign = 1;

	mp_limb_t* ABpts_inf =  scratch + abpts_offset; 							// nb_block_coeff * 2
	mp_limb_t* ABpts_two =  scratch + abpts_offset + nb_block_coeff * 4; 		// nb_block_coeff * 2 + 2
	mp_limb_t* ABpts_mone = scratch + abpts_offset + nb_block_coeff * 6 + 2; 	// nb_block_coeff * 2 + 2
	mp_limb_t* ABpts_one =  scratch + abpts_offset + nb_block_coeff * 8 + 4; 	// nb_block_coeff * 2 + 2
	mp_limb_t* ABpts_zero = scratch + abpts_offset + nb_block_coeff * 2; 		// nb_block_coeff * 2

	mp_limb_t* ab0 = scratch + ab_offset;
	mp_limb_t* ab1 = scratch + ab_offset + nb_block_coeff *15;
	mp_limb_t* ab2 = scratch + ab_offset + nb_block_coeff * 4; 		// nb_block_coeff * 2 + 2
	mp_limb_t* ab3 = scratch + ab_offset + nb_block_coeff * 6 + 2;	// nb_block_coeff * 2 + 2
	mp_limb_t* ab4 = scratch + ab_offset + nb_block_coeff * 2; 		// nb_block_coeff * 2

	mp_limb_t* ab3_aux = scratch + aux_offset;


	printf("%d %d\n", nb_block_coeff, nb_block_last_coeff);

	
	/*gmp_printf("%Nd \n", a, nb_limbs);
	
	gmp_printf("%Nd \n", a0, nb_block_coeff);
	gmp_printf("%Nd \n", a1, nb_block_coeff);
	gmp_printf("%Nd \n", a2, nb_block_coeff);*/


	//Apts calc

	//Apts_inf <- a2
	mpn_copyd(Apts_inf, a2, nb_block_last_coeff);


	//Apts_two <- 4a2 + 2a1 + a0
	mpn_copyd(Apts_two, a2, nb_block_last_coeff); //nb_block_last_coeff

	if (mpn_lshift(Apts_two, Apts_two, nb_block_last_coeff, 1)){ //nb_block_last_coeff +1
		Apts_two[nb_block_last_coeff] = 1;
	}
	
	if (mpn_add(Apts_two, a1, nb_block_coeff, Apts_two, nb_block_last_coeff + 1)){ //nb_block_coeff+1
		Apts_two[nb_block_coeff] = 1;
	}

	mpn_lshift(Apts_two, Apts_two, nb_block_coeff + 1, 1); // nb_block_coeff + 1
		

	mpn_add(Apts_two, Apts_two, nb_block_coeff + 1, a0, nb_block_coeff); // nb_block_coeff + 1


	//Apts_mone <- a2 + a0 - a1
	if (mpn_add(Apts_mone, a0, nb_block_coeff, a2, nb_block_last_coeff)){ // nb_block_coeff + 1
		Apts_mone[nb_block_coeff] = 1;
	}

	

	if (mpn_sub(Apts_mone, Apts_mone, nb_block_coeff + 1, a1, nb_block_coeff)){
		mpn_neg(Apts_mone, Apts_mone, nb_block_coeff + 1);
	}


	//Apts_one <- a2 + a1 + a0
	if (mpn_add_n(Apts_one, a0, a1, nb_block_coeff)){
		Apts_one[nb_block_coeff] = 1;
	}


	mpn_add(Apts_one, Apts_one, nb_block_coeff + 1, a2, nb_block_last_coeff);

	
	//Apts_zero <- a0
	mpn_copyd(Apts_zero, a0, nb_block_coeff);


	gmp_printf("%Nd \n\n", Apts_inf,  nb_block_last_coeff);
	gmp_printf("%Nd \n\n", Apts_two,  nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Apts_mone, nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Apts_one,  nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Apts_zero, nb_block_coeff);

	//Bpts calc

	//Bpts_inf <- b2
	mpn_copyd(Bpts_inf, b2, nb_block_last_coeff);


	//Bpts_two <- 4b2 + 2b1 + b0
	mpn_copyd(Bpts_two, b2, nb_block_last_coeff);

	if (mpn_lshift(Bpts_two, Bpts_two, nb_block_last_coeff, 1)){
		Bpts_two[nb_block_last_coeff] = 1;
	}
	
	if (mpn_add(Bpts_two, b1, nb_block_coeff, Bpts_two, nb_block_last_coeff + 1)){
		Bpts_two[nb_block_coeff] = 1;
	}

	if(mpn_lshift(Bpts_two, Bpts_two, nb_block_coeff, 1)){
		Bpts_two[nb_block_coeff] = 1;
	}

	mpn_add(Bpts_two, Bpts_two, nb_block_coeff + 1, b0, nb_block_coeff);


	//Bpts_mone <- b2 + b0 - b1
	if (mpn_add(Bpts_mone, b0, nb_block_coeff, b2, nb_block_last_coeff)){
		Bpts_mone[nb_block_coeff] = 1;
	}

	if (mpn_sub(Bpts_mone, Bpts_mone, nb_block_coeff + 1, b1, nb_block_coeff)){
		mpn_neg(Bpts_mone, Bpts_mone, nb_block_coeff + 1);
	}


	//Bpts_one <- b2 + b1 + b0
	if (mpn_add_n(Bpts_one, b0, b1, nb_block_coeff)){
		Bpts_one[nb_block_coeff] = 1;
	}


	mpn_add(Bpts_one, Bpts_one, nb_block_coeff + 1, b2, nb_block_last_coeff);

	
	//Bpts_zero <- b0
	mpn_copyd(Bpts_zero, b0, nb_block_coeff);


	/*gmp_printf("%Nd \n\n", Bpts_inf,  nb_block_last_coeff);
	gmp_printf("%Nd \n\n", Bpts_two,  nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Bpts_mone, nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Bpts_one,  nb_block_coeff + 1);
	gmp_printf("%Nd \n\n", Bpts_zero, nb_block_coeff);*/



	//ABpts calc

	mpn_mul_n(ABpts_inf, Apts_inf, Bpts_inf, nb_block_coeff);
	mpn_mul_n(ABpts_two, Apts_two, Bpts_two, nb_block_coeff + 1);
	mpn_mul_n(ABpts_mone, Apts_mone, Bpts_mone, nb_block_coeff + 1);
	mpn_mul_n(ABpts_one, Apts_one, Bpts_one, nb_block_coeff + 1);
	mpn_mul_n(ABpts_zero, Apts_zero, Bpts_zero, nb_block_coeff);

	/*gmp_printf("%Nd \n\n", ABpts_inf,  2 * nb_block_coeff);
	gmp_printf("%Nd \n\n", ABpts_two,  2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ABpts_mone, 2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ABpts_one,  2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ABpts_zero, 2 * nb_block_coeff);*/



	//AB calc

	//AB0 <- P(0)
	mpn_copyd(ab0, ABpts_zero, nb_block_coeff * 2);


	//AB2 <- -P(0) -P(inf) + (P(1) + P(-1))/2, no aux var req

	mpn_copyd(ab2, ABpts_one, nb_block_coeff * 2 + 2); // P(1)
	if (mpn_add_n(ab2, ab2, ABpts_mone, nb_block_coeff * 2 + 2)){ // P(1) + P(-1)
		ab2[nb_block_coeff * 2 + 2] = 1;
	}
	mpn_rshift(ab2, ab2, nb_block_coeff * 2 + 2, 1); // (P(1) + P(-1))/2

	mpn_sub(ab2, ab2, nb_block_coeff * 2 + 2, ABpts_zero, 2 * nb_block_coeff); // (P(1) + P(-1))/2 - P(0)
	
	mpn_sub(ab2, ab2, nb_block_coeff * 2 + 2, ABpts_inf, 2 * nb_block_coeff); // (P(1) + P(-1))/2 - P(0) - P(inf)


	//AB3 <- -2P(inf) + (P(2) - P(-1))/6 + (P(0) - P(1))/2, need one aux var
	mpn_copyd(ab3, ABpts_inf, nb_block_coeff * 2); // P(inf)
	
	if (mpn_lshift(ab3, ab3, nb_block_coeff * 2, 1)){ // 2P(inf)
		ab3[nb_block_coeff * 2] = 1;
	}

	mpn_neg(ab3, ab3, nb_block_coeff * 2+1); // -2P(inf)

	
	

	mpn_copyd(ab3_aux, ABpts_one, nb_block_coeff * 2 + 2); // aux : P(1)
	mpn_neg(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2); // aux : -P(1)
	if (mpn_add(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2, ABpts_zero, nb_block_coeff * 2)){ // aux : -P(1) + P(0)
		ab3_aux[nb_block_coeff * 2  + 2] = 1;
	}
	mpn_rshift(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2, 1); // aux : (-P(1) + P(0)) / 2
	
	//mpn_neg(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2);

	//gmp_printf("%Nx \n\n", ab3_aux,  2 * nb_block_coeff + 2);


	if (mpn_add(ab3, ab3_aux, nb_block_coeff * 2 + 2, ab3, nb_block_coeff * 2)){ // -2P(inf) +  (-P(1) + P(0)) / 2
		ab3_aux[nb_block_coeff * 2 + 2] = 1;
	}

	mpn_zero(ab3_aux, nb_block_coeff * 2 + 2); // aux : 0

	mpn_copyd(ab3_aux, ABpts_mone, nb_block_coeff * 2 + 2); // aux : P(-1)
	mpn_neg(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2); // aux : -P(-1)

	mpn_add_n(ab3_aux, ab3_aux, ABpts_two, nb_block_coeff * 2 + 2); // aux : -P(-1) + P(2)
	
	mpn_rshift(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2, 1); // aux : (-P(-1) + P(2)) / 2
	mpn_divexact_by3(ab3_aux, ab3_aux, nb_block_coeff * 2 + 2); // aux : (-P(-1) + P(2)) / 6

	mpn_add_n(ab3, ab3, ab3_aux, nb_block_coeff * 2 + 2); // -2P(inf) + (P(2) - P(-1))/6 + (P(0) - P(1))/2




	//AB4 <- P(inf)
	mpn_copyd(ab4, ABpts_inf, nb_block_coeff * 2);



	/*gmp_printf("%Nd \n\n", ab0,  2 * nb_block_coeff);
	//gmp_printf("%Nd \n\n", ab1,  2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ab2, 2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ab3,  2 * nb_block_coeff + 2);
	gmp_printf("%Nd \n\n", ab4, 2 * nb_block_coeff);*/


}

void toom3(mpz_t a, mpz_t b, mpz_t ab){
	mpz_t A[3];
	mpz_t B[3];

	/*gmp_printf("A : : %Zd\n\n", a);
	gmp_printf("B : : %Zd\n\n", b);*/

	printf("\n \n");
	
	mp_limb_t* scratch = calloc(10000, sizeof(mp_limb_t));

	toom3_mpn(a->_mp_d, b->_mp_d, b->_mp_size, ab->_mp_d, scratch);

	free(scratch);

	printf("\n \n");

	toPoly(a, A);
	toPoly(b, B);

	

	for (int i = 0; i < 3; i++){
		//gmp_printf("A%d %d : %Zd\n\n", i, A[i]->_mp_size, A[i]);
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
		gmp_printf("Apts%d : %Zd\n\n", i, Apts[i]);
	}
	for (int i = 0; i < 5; i++){
		//gmp_printf("Bpts%d : %Zd\n\n", i, Bpts[i]);
	}

	mpz_t ABpts[5];

	for (int i = 0; i < 5; i++)
	{
		mpz_init(ABpts[i]);
		mpz_mul(ABpts[i], Apts[i], Bpts[i]); //this should be changed to kara once it's done
		//gmp_printf("ABpts%d : %Zd\n\n", i, ABpts[i]);
	}

	mpz_t AB[5];
	calcAB(ABpts, AB);

	for (int i = 0; i < 5; ++i)
	{
		//gmp_printf("AB%d : %Zd\n\n", i, AB[i]);

	}

	mpz_t tmp, tmp2;
	mpz_init(tmp);
	mpz_init(tmp2);

	mpz_set(tmp, ABpts[0]);
	mpz_mul_ui(tmp, tmp, 2);
	mpz_neg(tmp, tmp);

//	gmp_printf("%Zx\n\n", tmp);

	mpz_set(tmp2, ABpts[3]);
	mpz_neg(tmp2, tmp2);
	mpz_add(tmp2, tmp2, ABpts[4]);
	mpz_divexact_ui(tmp2, tmp2, 2);

	//gmp_printf("%Zx\n\n", tmp2);


	mpz_clear(tmp);
	mpz_clear(tmp2);

	eval5(AB, powl, ab);

	//gmp_printf("ab_ : %Zd\n\n", ab);


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


	mpz_t ab, ab_;

	mpz_init(ab);
	

	mpz_mul(ab, a, b);
	//gmp_printf("ab : %Zd\n\n", ab);


		

	toom3(a, b, ab_);

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


	return 0;
}