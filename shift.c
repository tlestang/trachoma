#include <stdio.h>
#include <stdint.h>

void printbytearray(uint8_t *b, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < 8; ++j)
      printf("%d", !!((char)(b[i] << j) & '\x80'));
    printf("%c", (i < n - 1) ? ' ' : '\n');
  }
}

void shift(uint8_t * A, int shift, int size) {
  int i;
  uint8_t buf[2], mask;
  for (i = 0, mask = 1; i < shift; ++i)
    mask *= 2;
  mask =- 1;
  buf[0] = 0;
  for (i = 0; i < size; ++i) {
    buf[1] = buf[0] << (8 - shift); buf[0] = 0; // make space in buffer
    buf[0] |= (A[i] & mask); // record
    A[i] >>= shift; // shift
    A[i] |= buf[1]; // apply
  }
}

void bgd_death(uint8_t *A, int idx, int size) {
  uint8_t buf, mask;
  int n, m, i;
  m = idx % 8; // local index in byte block
  n = (idx / 8) + 1; // size of subarray
  for (mask = 1, i = 0; i < 8 - m - 1; ++i)
    mask *= 2;
  mask -= 1;
  buf = A[n-1] & mask;
  shift(A, 1, n);
  A[n-1] = (A[n-1] & ~mask) | buf;
}

void rotate(uint8_t *A, int n, int size) {
  int i, k;
  k = n / 8;
  shift(A, n%8, size - k);
  printbytearray(A, size);
  for (i = size-1; i >= k; --i)
    A[i] = A[i - k];
  for (i = 0; i < k; ++i)
    A[i] = 0;
}

int main() {
  uint8_t b[] = {72, 98, 13, 56, 112};
  int n = 5, i;
  uint8_t buf[2];
  printbytearray(b, n);
  //rotate(b, 17, n);
  bgd_death(b, 23, n);
  // shift(b, 3, 5);
  printbytearray(b, n);
}
