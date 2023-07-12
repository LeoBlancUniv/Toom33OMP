#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

//#include "structs.h"

#include <gmp.h>
#include "gmp_stuff.c"

#include "toom33_mul.c"
#include "mul_mpn.h"

#define NTEST 501
#define NSAMPLES 1001

#define UNUSED(X) (void)(X)

int64_t randomint64(void)
{
	return (((int64_t)rand() ^ rand()) << 32) | ((int64_t)rand() ^ rand());
}

mp_ptr scratch, scratch2;
int skip;


// NTEST*NSAMPLES must be odd
// it's easier to compute median value


/**** Measurements procedures according to INTEL white paper

	"How to benchmark code execution times on INTEL IA-32 and IA-64"

*****/

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

static inline void gmp_lomidhi32_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *soak, mp_limb_t *c_limbs, mp_limb_t *soak2, int nb_limbs){

		UNUSED(soak);
		UNUSED(soak2);

		lomidhi32(a_limbs, b_limbs, nb_limbs, c_limbs);

		/*mp_ptr x = calloc(nb_limbs * 2, sizeof(mp_limb_t));

		mpn_mul_n(x, a_limbs, b_limbs, nb_limbs);


		
		if (skip == 0 && mpn_cmp(x, c_limbs, 2 * nb_limbs) != 0){
			skip = 1;
			mp_ptr xor = calloc(nb_limbs * 2, sizeof(mp_limb_t));
			mpn_xor_n(xor, x, c_limbs, nb_limbs);
			
			
			free(xor);

			int hd = mpn_hamdist(x, c_limbs, nb_limbs);

			

			gmp_printf("a : %Nd\n\n", a_limbs, nb_limbs);
			gmp_printf("b : %Nd\n\n", b_limbs, nb_limbs);

			gmp_printf("c : %Nx\n\n", c_limbs, nb_limbs * 2);
			gmp_printf("x : %Nx\n\n", x, nb_limbs * 2);
			gmp_printf("d : %Nx\n\n", xor, nb_limbs * 2);
			printf("%d\n", hd);

		}

		

		free(x);*/
}

static inline void gmp_lohi22_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *soak, mp_limb_t *c_limbs, mp_limb_t *soak2, int nb_limbs){

		UNUSED(soak);
		UNUSED(soak2);

		lohi22(a_limbs, b_limbs, nb_limbs, c_limbs);

		/*mp_ptr x = calloc(nb_limbs * 2, sizeof(mp_limb_t));

		mpn_mul_n(x, a_limbs, b_limbs, nb_limbs);


		
		if (skip == 0 && mpn_cmp(x, c_limbs, 2 * nb_limbs) != 0){
			skip = 1;
			mp_ptr xor = calloc(nb_limbs * 2, sizeof(mp_limb_t));
			mpn_xor_n(xor, x, c_limbs, nb_limbs);
			
			
			free(xor);

			int hd = mpn_hamdist(x, c_limbs, nb_limbs);

			

			gmp_printf("a : %Nd\n\n", a_limbs, nb_limbs);
			gmp_printf("b : %Nd\n\n", b_limbs, nb_limbs);

			gmp_printf("c : %Nx\n\n", c_limbs, nb_limbs * 2);
			gmp_printf("x : %Nx\n\n", x, nb_limbs * 2);
			gmp_printf("d : %Nx\n\n", xor, nb_limbs * 2);
			printf("%d\n", hd);

		}

		

		free(x);*/
}

static inline void gmp_toom3_mpn_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *soak, mp_limb_t *c_limbs, mp_limb_t *soak2, int nb_limbs){

		UNUSED(soak);
		UNUSED(soak2);

		scratch = scratch2;

		toom3_mpn(a_limbs, b_limbs, nb_limbs, c_limbs, scratch, skip);

		/*mp_ptr x = calloc(nb_limbs * 2, sizeof(mp_limb_t));

		mpn_mul_n(x, a_limbs, b_limbs, nb_limbs);


		
		if (skip == 0 && mpn_cmp(x, c_limbs, 2 * nb_limbs) != 0){
			skip = 1;
			mp_ptr xor = calloc(nb_limbs * 2, sizeof(mp_limb_t));
			mpn_xor_n(xor, x, c_limbs, nb_limbs);
			
			
			free(xor);

			int hd = mpn_hamdist(x, c_limbs, nb_limbs);

			

			gmp_printf("a : %Nd\n\n", a_limbs, nb_limbs);
			gmp_printf("b : %Nd\n\n", b_limbs, nb_limbs);

			gmp_printf("c : %Nx\n\n", c_limbs, nb_limbs * 2);
			gmp_printf("x : %Nx\n\n", x, nb_limbs * 2);
			gmp_printf("d : %Nx\n\n", xor, nb_limbs * 2);
			printf("%d\n", hd);

		}

		

		free(x);*/

}

static inline void gmp_toom3_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *soak, mp_limb_t *c_limbs, mp_limb_t *soak2, int nb_limbs)
{
	UNUSED(soak);
	UNUSED(soak2);

	scratch = scratch2;

	mpn_toom33_mul(c_limbs, a_limbs, nb_limbs, b_limbs, nb_limbs, scratch);

	/*mp_ptr x = calloc(nb_limbs * 2, sizeof(mp_limb_t));

		mpn_mul_n(x, a_limbs, b_limbs, nb_limbs);


		
		if (skip == 0 && mpn_cmp(x, c_limbs, 2 * nb_limbs) != 0){
			skip = 1;
			mp_ptr xor = calloc(nb_limbs * 2, sizeof(mp_limb_t));
			mpn_xor_n(xor, x, c_limbs, nb_limbs);
			
			
			free(xor);

			int hd = mpn_hamdist(x, c_limbs, nb_limbs);

			

			gmp_printf("a : %Nx\n\n", a_limbs, nb_limbs);
			gmp_printf("b : %Nx\n\n", b_limbs, nb_limbs);

			gmp_printf("c : %Nx\n\n", c_limbs, nb_limbs * 2);
			gmp_printf("x : %Nx\n\n", x, nb_limbs * 2);
			gmp_printf("d : %Nx\n\n", xor, nb_limbs * 2);
			printf("%d\n", hd);

		}

		

		free(x);*/



	
}



static inline void gmp_lowlevel_wrapper(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *c_limbs, mp_limb_t *q_limbs, int nb_limbs)
{
	mpn_mul_n(c_limbs, a_limbs, b_limbs, nb_limbs); // compute: z = y*x
	//mpn_tdiv_qr(q_limbs, a_limbs, 0, c_limbs, (nb_limbs*2), p_limbs, nb_limbs); // compute: y = z%p
}



static inline uint64_t gmpbench(mpz_t A, mpz_t B, mpz_t modul_p, gmp_randstate_t r, mp_limb_t *c_limbs, mp_limb_t *q_limbs, uint8_t W, void (*gmp_wrapper)(mp_limb_t *a_limbs, mp_limb_t *b_limbs,
		mp_limb_t *p_limbs, mp_limb_t *c_limbs, mp_limb_t *q_limbs, int nb_limbs))
{

	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0,
	meanTimermax = 0, t1, t2, diff_t, *statTimer;
	uint64_t *cycles = (uint64_t *)calloc(NTEST,sizeof(uint64_t));
	mp_limb_t *p_limbs, *a_limbs, *b_limbs;
	int nb_limbs = mpz_size(modul_p);
	for(int i=0;i<NTEST;i++)
	{
	// Here we "heat" the cache memory.
		mpz_urandomm(A, r, modul_p);
		mpz_urandomm(B, r, modul_p);
		
		p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
		a_limbs = mpz_limbs_modify (A, nb_limbs);
		b_limbs = mpz_limbs_modify (B, nb_limbs);
		gmp_wrapper(a_limbs, b_limbs, p_limbs, c_limbs, q_limbs, nb_limbs);
		mpn_zero(scratch2, 10000);
		mpn_zero(c_limbs, 2 * nb_limbs);

	}

		

	
	for(int i=0;i<NSAMPLES;i++)
	{
		// Here we generate a random dataset to use for our test each iteration.
		mpz_urandomm(A, r, modul_p);
		mpz_urandomm(B, r, modul_p);
	
		p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
		a_limbs = mpz_limbs_modify (A, nb_limbs);
		b_limbs = mpz_limbs_modify (B, nb_limbs);
		timermin = (uint64_t)0x1<<63;
		timermax = 0;
		memset(cycles,0,NTEST*sizeof(uint64_t));
		for(int j=0;j<NTEST;j++)
		{
			//printf("j : %d\n", j);
			t1 = cpucyclesStart();
			// We call the function W times to get an accurate measurement.
			for(int soak=0; soak < W; soak++){
				gmp_wrapper(a_limbs, b_limbs, p_limbs, c_limbs, q_limbs, nb_limbs);
			}
			t2 = cpucyclesStop();

			mpn_zero(scratch2, 10000);
			mpn_zero(c_limbs, 2 * nb_limbs);

			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			if(timermin > diff_t) timermin = diff_t;
			else if(timermax < diff_t) timermax = diff_t;
			cycles[j]=diff_t;
		}
		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles,NTEST);
		medianTimer += statTimer[1];
		free(statTimer);

	}
	free(cycles);
	return medianTimer/NSAMPLES/W; // We divide by W since we measured W calls.
}

void do_benchgmp(uint64_t retcycles[3], const char* pstr, const uint8_t W)
{
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	mpz_t modul_p;
	mpz_init(modul_p);
	
	mpz_set_str(modul_p, pstr, 0);
	
	int nb_limbs = mpz_size(modul_p);
	
	mpz_t A, B;
	mpz_inits(A, B, NULL);
	
	mp_limb_t *p_limbs, *c_limbs, *q_limbs;
	
	q_limbs = (mp_limb_t*) calloc ((nb_limbs+1), sizeof(mp_limb_t));
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	
	retcycles[0] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_lowlevel_wrapper);
	
	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));

	retcycles[1] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_toom3_wrapper);

	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));

	retcycles[2] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_toom3_mpn_wrapper);
	
	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));

	retcycles[3] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_lohi22_wrapper);

	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));

	retcycles[4] = gmpbench(A, B, modul_p, r, c_limbs, q_limbs, W, gmp_lomidhi32_wrapper);


	/*p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	free(c_limbs);
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));
	mpn_binvert(mip_limbs, p_limbs, nb_limbs, c_limbs);
	
	retcycles[1] = gmpbench(A, B, modul_p, r, mip_limbs, q_limbs, W, gmp_montgomery_wrapper);
	
	mp_limb_t mip0;
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	retcycles[2] = gmpbench(A, B, modul_p, r, &mip0, q_limbs, W, gmp_montgomeryCIOS_wrapper);
	*/
	mpz_clears(modul_p, A, B, NULL);
	free(c_limbs);
	free(q_limbs);
	gmp_randclear(r);
}

int main(void)
{
	scratch2 = calloc(10000, sizeof(mp_limb_t));
	skip = 0;
	
	const char *p8192 = "825071401847911350104434016317510780278067283201566754357130842744461195005906395199943718051043673223466730721751949974449248930827451854619414374697619260438961946911749222259780365386291751814038486373465100001426528819046274070279571501847524807503058829803292949581020312542244512310826381340725147604205504421489781205421741398803728752958064593794237603834697854344675329642780388241249095996086652843251798333913223849757378836474373313709743648650234857926624244596951216675142492484734159681684027468402207806205439865347969018612132747592285828382232580799515586663970404103616260646140491130336474241553520642020525516089660137732716330093125235058568988907834750945317769702575635573206226618491699434251943292774344986722857013298412092723076929179307523912203128114090534314456523272023692884978468594858168466441251073816641214891282204410271065824121943315600314597749227889573226269102256807459458734770833396972179130110520264774827679636735786580446623637598159822709594268849412087335528825773591001341459368429154598701743523091608321042347244119166311797831282266269359216738748446338689877553411067459022750494469035170847673664139376175016704738727113860401649828339318679003054736839277206870283192227831617539598383681659913506027472064659239184916246405287069652026544516606976680852350358887886206450202813665604368556089115665756191850387879202714813769550164656889996790902837958128457465941454508790612946593532755196596172284379325906458537834759016596544236993315228152015361031327422453823347007209987459925728347755977454611604649086058581101577344837613218863165735163293232182683983328352703931948934018941328642629738690331539423670935877924451188873671016602730573554563515446477428747122842429919747738097859733596295890237973718047460047408090635470616266779559457882525693610773409953808315167993488749629373478589355791418072020121684946874867922541121111773785091209812728948633266300346422009914729246105964646426829796911062733042949990258625136306447831184696881439764159871493270186815550008827682405281004158733425515557787170576476464934246505105516394418797543153382468729141414796208872453491051467183673152632331727482671336127896268657700582305093675403324150839387303020629188315382088652129003636647961202147261937023281222481237476524670557355735034515267765254502192363955049952594249904232853763728215835537753171258068562061271398805505638547424578150686558691852996512332614361427686993333226212556141099";

	uint64_t cyclesGMP8192[5];

	do_benchgmp(cyclesGMP8192, p8192, 1);


	printf("low : %lu \n", cyclesGMP8192[0]);
	printf("toom33_original : %lu \n", cyclesGMP8192[1]);
	printf("toom33_mpn : %lu \n", cyclesGMP8192[2]);
	printf("lohi22 : %lu \n", cyclesGMP8192[3]);
	printf("lomidhi32 : %lu \n", cyclesGMP8192[4]);

	free(scratch2);


	return 0;
}
