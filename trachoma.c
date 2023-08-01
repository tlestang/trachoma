#include <stdint.h>
#include <math.h>
#include <stdio.h>

void apply_rules(uint8_t *inf,
		 uint8_t *dis,
		 uint8_t *lat,
		 int *clock, int n) {
  int i, nblocks, j;
  uint8_t new_s, new_d, clearinf, trans;

  nblocks = n / 8;

  for (i = 0; i < nblocks; ++i) {
    trans = 0;
    for (j=0; j < 8; ++j)
      trans |= (((uint8_t)(!clock[j + i * 8])) << (7 - j));

    new_s = *(dis+i) & ~*(inf+i) & trans;
    new_d = *(inf+i) & ~*(lat+i) & trans;
    clearinf = *(lat+i) & trans;

    dis[i] = dis[i] & ~new_s | clearinf;
    inf[i] = inf[i] & ~new_d; // | new_i;
    lat[i] = lat[i] & ~clearinf; // | new_i;
  }
}

#define BETA 0.21
#define V_1 1.
#define V_2 2.6
#define PHI 1.4
#define EPSILON 0.5

void get_infection_prob(int *ages, double *ld, double *prob,
		int *groups, int ngroups, int n) {
  int n_prev, n_ingroup, igroup, k;
  double lam[3], pop_ratio[3], one_over_popsize, meanld, sum_ld;
  // unsigned int groups[] = {468, 780, 3121};
  double epsm;

  k = 0;
  one_over_popsize = 1. / n;
  epsm = 1. - EPSILON;
  
  for (k = 0, igroup = 0; igroup < ngroups; ++igroup) {
    n_prev = k;
    for(sum_ld = 0; (ages[k] < groups[igroup]) &&  k < n; ++k)
      sum_ld += ld[k];

    n_ingroup = (k - n_prev);
    meanld = sum_ld / n_ingroup;
    lam[igroup] = BETA * V_1 * meanld + BETA * V_2 * pow(meanld, PHI + 1);
    pop_ratio[igroup] = (double)n_ingroup * one_over_popsize;
  }

  A = (double *) malloc(sizeof(double) * ngroups);
  A[0] = -lam[0]*pop_ratio[0] - lam[1]*epsm*pop_ratio[1] - lam[2]*epsm*pop_ratio[2];
  A[1] = -lam[0]*pop_ratio[0]*epsm - lam[1]*pop_ratio[1] - lam[2]*epsm*pop_ratio[2];
  A[2] = -lam[0]*pop_ratio[0]*epsm - lam[1]*epsm*pop_ratio[1] - lam[2]*pop_ratio[2];

  for (k = 0, igroup = 0; igroup < ngroups; ++igroup) {
    while ((ages[k] < groups[igroup]) &&  k < n)
      prob[k++] = 1. - exp(A[igroup]);
  }
  free(A);
}
