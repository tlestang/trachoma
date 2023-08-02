#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "periods.h"

void get_infection_prob(int*, double*, double*,
			int*, int, int);

void apply_rules(uint8_t *inf,
		 uint8_t *dis,
		 uint8_t *lat,
		 int *clock, int *ages, double *bactload, int n) {
  int i, nblocks, j;
  uint8_t new_s, new_d, clearinf, trans, isinf, infect;
  uint8_t *new_i;
  int groups[] = {468, 780, 3121}; int ngroups;
  double *prob;

  nblocks = n / 8;

  ngroups = 3;

  // TODO: Sort ages and bactload here
  prob = (double *) malloc(sizeof(double) * n);
  get_infection_prob(ages, bactload, prob, groups, ngroups, n);

  new_i = (uint8_t *) malloc(sizeof(uint8_t) * nblocks);
  for (i = 0; i < nblocks; ++i) {
    trans = 0;
    for (j=0; j < 8; ++j) {
      trans |= (((uint8_t)(!clock[j + i * 8])) << (7 - j));
      isinf = (inf[i] << j) & '\x80';
      infect = !isinf && ((rand() / RAND_MAX) < prob[j + i * 8]);
      new_i[i] |= infect << (7 - j);
    }

    new_s = *(dis+i) & ~*(inf+i) & trans;
    new_d = *(inf+i) & ~*(lat+i) & trans;
    clearinf = *(lat+i) & trans;

    dis[i] = dis[i] & ~new_s | clearinf;
    inf[i] = inf[i] & ~new_d | new_i[i];
    lat[i] = lat[i] & ~clearinf | new_i[i];

    for (j=0; j < 8; ++j) {

      k = j + i * 8;
      // TODO: set MAXAGE and define bgd_death()
      if (ages[k] == MAX_AGE || bgd_death()) {
	clock[k] = -1;
	count[k] = 0;
	ages[k] = -1;
	mask = ~(1 << (7 - j));
	dis[i] &= mask;
	inf[i] &= mask;
	lat[i] &= mask;
	continue;
      }

      isnewd = (new_d << j) & '\x80';
      isclearinf = (clearinf << j) & '\x80';
      isnewi = (new_i << j) & '\x80';

      if (!(isnewd || isclearinf || isnewi))
	continue;

      // TODO: Implement time functions
      if (isnewd) {
	clock[k] = setdtime(D_base[k], count[k], ages[k]);
	bactload[k] = 0.;
	continue;
      } else if (isclearinf) {
	clock[k] = setidtime(ID_base[k], count[k], ages[k]);
	bactload[k] = get_load(count[k]);
	continue;
      } else if (isnewi) {
	clock[k] = setlatenttime(latent_base[k], count[k], ages[k]);
	count[k]++;
      }
    }
    for (j=0; j < 8; ++j)
      ages[k] += 1;
  }
  free(prob);
  free(new_i);
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
  double epsm, *A;

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
