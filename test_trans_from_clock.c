#include <stdio.h>
#include <stdint.h>

int main() {
  int i, j;
  int nblocks = 1;
  uint8_t trans, mask;
  int clock[8] = {14, 0, 0, 1, 2, 0, 7, 18};

  for (i = 0; i < nblocks; ++i) {
    trans = 0;
    for (j=0; j < 8; ++j) {
      trans |= (((uint8_t)(!clock[j + i * 8])) << (7 - j));
    }
    printf("trans = %u", trans);
  }
}
