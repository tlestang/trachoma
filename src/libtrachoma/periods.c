#include <math.h>

#define MIN_ID 11
#define MIN_D 1
#define INF_RED 0.45
#define AG 0.00179
#define AQ 0.0368

int setidtime(int base, int count, int age) {
  double b = base - MIN_ID;
  double a = INF_RED * (1 - count);
  return (int)floor(b * exp(a) + MIN_ID);
}

int setdtime(int base, int count, int age) {
  double b = base - MIN_D;
  double a = AQ * (1 - count) - AG * age;
  return (int)floor(b * exp(a) + MIN_D);
}

int setlatenttime(int base, int count) {
  return (int)base;
}
