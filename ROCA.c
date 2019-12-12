// gcc -I/usr/local/include ROCA.c -o ROCA -lgmp -L/usr/local/lib -std=c99
#define _XOPEN_SOURCE 500
#include <stdint.h>
#include <stdio.h>			/* printf() */
#include <stdlib.h>                     /* abort() */
#include <unistd.h>                     /* getopt() */
#include <gmp.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "mpi.h"
#include "mpprime.h"

/* Constants */

//----------------------- create primorial ------------------------
const int prime_list[200] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223};

/* Command Line Parsing  */

#define DEFAULT_TARGET_BITSIZE 512
#define DEFAULT_FUDGE_FACTOR 1
#define DEFAULT_PRIME_COUNT 69
#define DEFAULT_SEED 0x1337
#define DEFAULT_COUNT 32
#define DEFAULT_VERBOSE 1
#define DEFAULT_OFFSET 0

struct _cmdline_params_struct {
  uint64_t target_bitsize; //< target bitsize of x
   int64_t fudge_factor;   //< we want 2x+1 * 4x+1 to be 1024 bits, so use this to fudge.
  uint32_t prime_count;    //< number of primes
  uint32_t count;
  unsigned long int seed;  //< seed
  int32_t verbose; //< print passes
  uint32_t offset; //< how much we want to offset the starting value of k by
};

typedef struct _cmdline_params_struct cmdline_params_t[1];

static inline void print_help_and_exit() {
  printf("-b    0 < target bit-size <= 8192 (default: %d)\n", DEFAULT_TARGET_BITSIZE);
  printf("-f    fudge factor (default: %d)\n", DEFAULT_FUDGE_FACTOR);
  printf("-p    0 < prime count <= 200 (default: %d)\n", DEFAULT_PRIME_COUNT);
  printf("-c    log2 of number of trials (default: %d)\n", DEFAULT_COUNT);
  printf("-s    random seed in hex (default: %x)\n", DEFAULT_SEED);
  printf("-v    record this number of passes (default: False)");
  printf("-o    offset on starting k value, where offset*c (default: %d)\n", DEFAULT_OFFSET);
  abort();
}

static inline void parse_cmdline(cmdline_params_t params, int argc, char *argv[]) {
  params->target_bitsize = DEFAULT_TARGET_BITSIZE;
  params->fudge_factor = DEFAULT_FUDGE_FACTOR;
  params->prime_count = DEFAULT_PRIME_COUNT;
  params->seed = DEFAULT_SEED;
  params->count = DEFAULT_COUNT;
  params->verbose = DEFAULT_VERBOSE;
  params->offset = DEFAULT_OFFSET;

  int c;
  while ((c = getopt(argc, argv, "b:f:c:p:v:s:o:")) != -1) {
    switch(c) {
    case 'b':
      params->target_bitsize = (uint64_t)strtoul(optarg,NULL,10);
      break;
    case 'f':
      params->fudge_factor = (uint64_t)strtoul(optarg,NULL,10);
      break;
    case 'p':
      params->prime_count = (uint32_t)strtoul(optarg,NULL, 10);
      break;
    case 'c':
      params->count = (uint32_t)strtoul(optarg, NULL, 10);
      break;
    case 's':
      params->seed = (unsigned long int)strtoul(optarg, NULL, 16);
      break;
    case 'v':
      params->verbose = (int32_t)strtoul(optarg,NULL, 10);
      break;
    case 'o':
      params->offset = (int32_t)strtoul(optarg,NULL, 10);
      break;
    case ':':  /* without operand */
      print_help_and_exit();
    case '?':
      print_help_and_exit();
    }
  }

  printf("-b %zu -f %zu -p %d -s %zx -c %d -v %d -o %d\n", params->target_bitsize,
         params->fudge_factor, params->prime_count, params->seed, params->count, params->verbose, params->offset);

  if (params->target_bitsize > 8192)
    print_help_and_exit();
  if (params->prime_count > 200)
    print_help_and_exit();
}

/* Logging */

void logit(int res, mpz_t x, mpz_t n) {

  char tmp[2000];

  snprintf(tmp, 2000, "%02d:0x%s:0x%s", res, mpz_get_str(NULL, 16, x), mpz_get_str (NULL, 16, n));
  FILE *fh = fopen("ROCA.log", "a");
  fprintf(fh, "%s\n", tmp);
  fclose(fh);
}

/* Primality Testing */

int millerrabin(mpz_t n, mpz_t nm1, mpz_t x, mpz_t y, mpz_t q, unsigned long int k)
//
// Performs Miller rabin on n, to base x. y is a redundant parameter, but i have
// mimicked the call and arguments GMP uses, for better integration.
//
{
  mpz_powm(y,x,q,n);  //mpz_powm(y,x,q,n) = y = x^q mod n
  if (mpz_cmp_ui (y, 1L) == 0 || mpz_cmp (y, nm1) == 0)
    return 1;

  unsigned long int i;
  for (i = 1; i < k; i++)
    {
      mpz_powm_ui (y, y, 2L, n);
      if (mpz_cmp (y, nm1) == 0)
	return 1;
      if (mpz_cmp_ui (y, 1L) <= 0)
	return 0;
    }
  return 0;
}

int singlemrtest(mpz_t n, gmp_randstate_t rstate)
//
// Takes the mpz_t n the number we test and rstate the randomness for base selection
// this function used to take the second parameter as "a" the base to perform the test on
// i just moved the base selection inside of here so it isn't defined in the massive for loop
//
{
  mpz_t nm1,q,y,x,nm3;
  mpz_init(q);
  mpz_init(nm1);
  mpz_init(nm3);
  mpz_init(x);
  mpz_init(y);

  mpz_sub_ui (nm3, n, 3L);
  mpz_urandomm (x, rstate, nm3);
  mpz_add_ui (x, x, 2L);

  mpz_sub_ui(nm1,n,1L); // n-1
  unsigned long int k;
  k = mpz_scan1 (nm1, 0L);
  mpz_tdiv_q_2exp (q, nm1, k);

  int is_prime = millerrabin(n,nm1,x,y,q,k);
  mpz_clear(nm1);
  mpz_clear(q);
  mpz_clear(y);
  mpz_clear(nm3);
  mpz_clear(x);
  return is_prime;
}

void mpz_to_mpint(mp_int *out, mpz_t in)
{
  char* in_str = mpz_get_str(NULL, 16, in);
  mp_init(out);
  mp_read_radix(out, in_str, 16);
  printf("n: %s\n", in_str);
  free(in_str);
}

int test(mpz_t n){
  //
  // This function is used to test the maximum rounds of mpz_probab_prime_p n passes
  // with the returning valyes -1,0 understood as "no mr passed, was proven prime by trial
  // division or fermat test done first".
  //
  int count = 0; // we can start this at say 5 or 10 (no point)
  int countorig = count;
  // Convert n to NSS mp_int
  mp_int n_mp_int;
  mpz_to_mpint(&n_mp_int, n);
  srand(0);
  int res = mpp_pprime(&n_mp_int, count) == MP_YES ? 1 : 0;
  // int res= mpz_probab_prime_p (n, count);

  // fix NSS rand output
  while (res!=0){ // if the inital test passed, we want to see how far we can go
    if (count >= 10) {
      printf(" >>> passed MR test (%d) with n := ", count);
      exit(0);
    }
    count += 1;
    // res= mpz_probab_prime_p (n, count);
    srand(0);
    res = mpp_pprime(&n_mp_int, count) == MP_YES ? 1 : 0;
  }
  mp_clear(&n_mp_int);
  if (count-1 <countorig){
    return -1;
  }
  else{
    return count-1;
  }
}

int main(int argc, char *argv[]) {
 mpz_t M,k,x,n,p1,p2,one,two,four,off;
 mpz_init(x);
 mpz_init(p1);
 mpz_init(p2);
 mpz_init(n);
 mpz_init(off);
 mpz_init_set_ui(M, 1);
 mpz_init_set_ui(k, 1);
 mpz_init_set_ui(one, 1);
 mpz_init_set_ui(two, 2);
 mpz_init_set_ui(four, 4);

 // fix NSS rand output
 srand(0);

 cmdline_params_t params;
 parse_cmdline(params, argc, argv);
 //mpz_init_set_ui(off, (1ULL)<<params->offset); // set the offset from user input
 mpz_init_set_ui(off,params->offset);
 mpz_mul_2exp(off,off,params->count);
 // ---------------------- seed base selection for MR -------------
 //The randomness is seeded on the time, but this doesn't really matter as we
 //are only using this source of randomness to choose bases for Miller-Rabin.
 //We could exclude the use of randomness entirely, and call Miller-Rabin to
 // base 2 each time. This shouldn't make a difference to the accuracy.

 gmp_randstate_t rstate;
 gmp_randinit_default (rstate);
 gmp_randseed_ui(rstate, params->seed);

 for (uint32_t i = 0; i < params->prime_count; i++) {
   mpz_mul_si(M, M, prime_list[i]);
 }
 size_t Mbits = mpz_sizeinbase (M, 2); // get bitsize of M

 //-------------------- generate k ---------------------------
 int kbits = params->target_bitsize - Mbits - params->fudge_factor;
 mpz_mul_2exp (k, k, kbits);

 // ---------- add the offset value (always even) -----------
 if (params->offset != 0) {
   mpz_add(k,k,off); //
 }
// ------------ begin BIG loop on k ---------------------
 for (uint64_t i = 0; i <= (1ULL)<<params->count; i++) {
   mpz_add_ui(k,k,1);
   mpz_mul(x,k,M);
   mpz_add_ui(x,x,189); // offset to beat Fermat test.
   mpz_mul_2exp(p1,x,1);
   mpz_add_ui(p1,p1,1);
   if (singlemrtest(p1,rstate)==0) { // this first test takes up the majority of computation time
       continue;
   }
   mpz_mul_2exp(p2,x,2);
   mpz_add_ui(p2,p2,1);
   if (singlemrtest(p2,rstate)==0) {
       continue;
   }
   mpz_mul(n,p1,p2); // n = p1*p2
   int res = test(n);
   if (res >= params->verbose) {
     logit(res, x, n);
   }
 }
 mpz_clear(M);
 mpz_clear(k);
 mpz_clear(p1);
 mpz_clear(p2);
 mpz_clear(n);
 mpz_clear(x);
 mpz_clear(off);
 mpz_clear(one);
 mpz_clear(two);
 mpz_clear(four);
 return 0;
}
