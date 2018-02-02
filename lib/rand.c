
#include <assert.h>
#include <stdio.h>

unsigned long int rand_seed() {
  /*
    returns a random unsigned long integer read from the standard unix
    pseudorandom device /dev/urandom
  */

  FILE *f = fopen("/dev/urandom", "r");
  assert(f != NULL);

  unsigned long int seed;
  fread(&seed, sizeof(unsigned long int), 1, f);

  fclose(f);

  return seed;
}
