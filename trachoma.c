#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "periods.h"
#include "shift.h"

#define MAX_AGE 3120
#define TAU 1. / (40. * 52.)

int *D_base, *ID_base, *latent_base;

uint8_t *inf, *dis, *lat, *new_i;
int *clockm, *ages, *count;
double *bactload, *prob;

int groups[] = {468, 780, 3121}; int ngroups = 3;

const double BGD_DEATH_RATE = 1. - exp(-TAU);

void get_infection_prob(int);
double get_load(int);

void apply_rules(int n, int nblocks, int times) {
  int i, j, k, t;

  for (t = 0; t < times; ++t) {

    get_infection_prob(n);

    // if first indiv in byte is at MAX_AGE, then all indivs after are
    // as well. No need to process these blocks.
    for (i = 0; i < nblocks && (ages[i * 8] < MAX_AGE); ++i) {
      uint8_t trans = 0;
      uint8_t isinf, infect;
      for (j=0; j < 8; ++j) {
	k = j + i * 8;
	trans |= (((uint8_t)(!clockm[k])) << (7 - j));
	isinf = (inf[i] << j) & '\x80';
	double draw = ((double)rand() / RAND_MAX);
        infect = !isinf && ( draw < prob[k]);
	new_i[i] |= infect << (7 - j);
      }

      uint8_t new_s = *(dis+i) & ~*(inf+i) & trans;
      uint8_t new_d = *(inf+i) & ~*(lat+i) & trans;
      uint8_t clearinf = *(lat+i) & trans;

      dis[i] = dis[i] & ~new_s | clearinf;
      inf[i] = inf[i] & ~new_d | new_i[i];
      lat[i] = lat[i] & ~clearinf | new_i[i];

      for (j=0; j < 8; ++j) {
	k = j + i * 8;
	uint8_t isnewd = (new_d << j) & '\x80';
	uint8_t isclearinf = (clearinf << j) & '\x80';
	uint8_t isnewi = (new_i[i] << j) & '\x80';
	if (isnewd) {
	  clockm[k] = setdtime(D_base[k], count[k], ages[k]);
	  bactload[k] = 0.;
	  continue;
	} else if (isclearinf) {
	  clockm[k] = setidtime(ID_base[k], count[k], ages[k]);
	  bactload[k] = get_load(count[k]);
	  continue;
	} else if (isnewi) {
	  clockm[k] = setlatenttime(latent_base[k], count[k], ages[k]);
	  count[k]++;
	}
      }
    } // nblocks

    // Background mortality
    for (i = 0; ages[i] < MAX_AGE && i < n; ++i) {
      if (rand() / (double)RAND_MAX < BGD_DEATH_RATE) {
	bgd_death_bitarray(inf, i, nblocks);
	bgd_death_bitarray(dis, i, nblocks);
	bgd_death_bitarray(lat, i, nblocks);
	rotate(clockm, 1, i + 1, -1);
	rotate(ages, 1, i + 1, 0);
	rotate(count, 1, i + 1, 0);
	rotate_double(bactload, 1, i + 1, 0.);
      }
      ages[i] += 1;
    }

    // Natural death for people aged > MAX_AGE
    rotate(clockm, n - i, n, -1);
    rotate(count, n - i, n, 0);
    rotate(ages, n - i, n, 0);
    rotate_double(bactload, n - i, n, 0.);
    rotate_bitarray(inf, n - i, nblocks);
    rotate_bitarray(dis, n - i, nblocks);
    rotate_bitarray(lat, n - i, nblocks);
  }
}



#define BETA 0.21
#define V_1 1.
#define V_2 2.6
#define PHI 1.4
#define EPSILON 0.5

void get_infection_prob(int n) {
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
      sum_ld += bactload[k];

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

double get_load(int ninf) {
  double b1 = 1., ep2 =  - 0.114;

  return b1 * exp((ninf - 1) * ep2);
}

int main() {
  int n = 16;
  int nblocks = n / 8;
  prob = (double *) malloc(sizeof(double) * n);
  new_i = (uint8_t *) malloc(sizeof(uint8_t) * nblocks);

  inf = (uint8_t *) malloc(sizeof(uint8_t) * nblocks);
  dis = (uint8_t *) malloc(sizeof(uint8_t) * nblocks);
  lat = (uint8_t *) malloc(sizeof(uint8_t) * nblocks);
  clockm = (int *) malloc(sizeof(int) * n);
  ages = (int *) malloc(sizeof(int) * n);
  count = (int *) malloc(sizeof(int) * n);
  bactload = (double *) malloc(sizeof(double) * n);

  ID_base = (int *) malloc(sizeof(int) * n);
  D_base = (int *) malloc(sizeof(int) * n);
  latent_base = (int *) malloc(sizeof(int) * n);

  int ages_l[] = {203,  388,  586,  846, 1323,
		  1327, 1400, 1446, 1824, 1917,
		  2032, 2222, 2650, 3120, 3117,
		  3119};
  double bactload_l[] = {
    0.54589385, 0.64449697, 0.94586293, 0.89743933, 0.74142636,
    0.42041218, 0.30355981, 0.76096054, 0.62761568, 0.88280128,
    0.16847746, 0.41943052, 0.56054159, 0.09329043, 0.67442931,
    0.16027857
  };

  inf[0] =139; inf[1] = 239;
  dis[0] = 62; dis[1] = 26;
  lat[0] = 129; lat[1] = 229;
  int i;
  for (i = 0; i < n; ++i) {
    count[i] = 0;
    clockm[i] = 2;
    ages[i] = ages_l[i];
    bactload[i] = bactload_l[i];
    ID_base[i] = 2;
    latent_base[i] = 2;
    D_base[i] = 2;
  }
  clockm[0] = 0; clockm[1] = -1;
  clockm[5] = 0; clockm[10] = 0;
  clockm[12] = 0;

  // for (i = 0; i < 8; ++i)
    // printf("%c", clock[i] == 0 ? 'X' : 'O');
  // printf(" ");
  // for (i = 0; i < 8; ++i)
    // printf("%c", clock[8 + i] == 0 ? 'X' : 'O');
  // printf("\n");
  /* printbytearray(dis, nblocks); */
  /* printbytearray(inf, nblocks); */
  /* printbytearray(lat, nblocks); */

  apply_rules(n, nblocks, 1);
  // printf("\n");

  // printf("%c", '\n');

  /* printbytearray(dis, nblocks); */
  /* printbytearray(inf, nblocks); */
  /* printbytearray(lat, nblocks); */
}

void set_arrays(uint8_t *inf_m, uint8_t *dis_m, uint8_t *lat_m,
	 uint8_t *new_i_m, int *clock_m, int *ages_m, int *count_m,
	 double *bactload_m, double *prob_m) {
  ages = ages_m;
  inf = inf_m;
  dis = dis_m;
  lat = lat_m;
  new_i = new_i_m;
  clockm = clock_m;
  count = count_m;
  bactload = bactload_m;
  prob = prob_m;
}

void set_times(int *D_base_m, int *ID_base_m, int *latent_base_m) {
  D_base = D_base_m;
  ID_base = ID_base_m;
  latent_base = latent_base_m;
}
