#include <math.h>

#define MIN_ID 11
#define MIN_D 1
#define IND_RED 0.45
#define AQ 0.00179
#define AG 0.0368

int setidtime(double base, int count, int age) {
  double b = base - MIN_ID;
  double a = INF_RED * (1 - count);
  return (int)floor(b * exp(a) + MIN_ID);
}

int setdtime(double base, int count, int age) {
  double b = base - MIN_D;
  double a = AQ * (1 - count) - AG * age;
  return (int)floor(b * exp(a) + MIN_D);
}

int setlatenttime(double base, int count, int age) {
  return base;
}
