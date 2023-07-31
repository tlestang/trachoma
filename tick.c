#include <stdint.h>

void tick(uint8_t *inf,
	  uint8_t *dis,
	  int *clock,
	  double *bactld,
	  unsigned int *ages,
	  unsigned int *count,
	  int n) {
  int i, nblocks;
  uint8_t new_s, new_d, new_id, trans;

  nblocks = b / 8 + 1;
  for (i = 0; i < nblocks; ++i) {
    trans = 0;
    for (j=0; j < 8; ++j)
      trans |= (((uint8_t)(!clock[j + i * 8])) << (7 - j));

    new_s = *(dis+i) & ~*(inf+i) & trans;
    new_d = *(inf+i) & *(dis+i) & trans;
    new_id = *(inf+i) & ~*(dis+i) & trans;

    new_i = get_new_infections()
  }
    
