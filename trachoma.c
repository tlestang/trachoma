#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "periods.h"
#include "shift.h"

int *D_base, *ID_base, *latent_base;
double *prob;

int *groups, ngroups;

double BGD_DEATH_RATE;

void get_infection_prob(int);
double get_load(int);

struct state {
  int n;
  uint8_t *inf;
  uint8_t *dis;
  uint8_t *lat;
  int *clockm;
  int *ages;
  int *count;
  double *bactload;
};

void apply_rules(struct state st, int times) {
  int i, j, k, t;

  int nblocks = st.n / 8;

  for (t = 0; t < times; ++t) {

    get_infection_prob(st.n);

    // if first indiv in byte is at max_age, then all indivs after are
    // as well. No need to process these blocks.
    int max_age = groups[ngroups - 1];
    for (i = 0; i < nblocks && (st.ages[i * 8] < max_age); ++i) {
      uint8_t trans = 0, new_i = 0;
      uint8_t isinf, infect;
      for (j=0; j < 8; ++j) {
	k = j + i * 8;
	trans |= (((uint8_t)(!st.clockm[k])) << (7 - j));
	isinf = (st.inf[i] << j) & '\x80';
	double draw = ((double)rand() / RAND_MAX);
        infect = !isinf && ( draw < prob[k]);
	new_i |= infect << (7 - j);
      }

      uint8_t new_s = *(st.dis+i) & ~*(st.inf+i) & trans;
      uint8_t new_d = *(st.inf+i) & ~*(st.lat+i) & trans;
      uint8_t clearinf = *(st.lat+i) & trans;

      st.dis[i] = st.dis[i] & ~new_s | clearinf;
      st.inf[i] = st.inf[i] & ~new_d | new_i;
      st.lat[i] = st.lat[i] & ~clearinf | new_i;

      update_indivs(new_i, clearinf, new_d,
		    st.bactload, st.clockm, st.count);
    } // nblocks

    int n_old = background_mortality(st, BGD_DEATH_RATE, max_age);
    old_age_mortality(st, n_old);
  }
}

int background_mortality(struct state st, double rate, int max_age) {
  int i;
  for (i = 0; st.ages[i] < max_age && i < st.n; ++i) {
    if (rand() / (double)RAND_MAX < rate) {
      bgd_death_bitarray(st.inf, i, nblocks);
      bgd_death_bitarray(st.dis, i, nblocks);
      bgd_death_bitarray(st.lat, i, nblocks);
      rotate(st.clockm, 1, i + 1, -1);
      rotate(st.ages, 1, i + 1, 0);
      rotate(st.count, 1, i + 1, 0);
      rotate_double(st.bactload, 1, i + 1, 0.);
    }
    st.ages[i] += 1;
  }
  return st.n - i;
}

void old_age_mortality(struct state st, int nold) {
  int nblocks = st.n / 8;
  rotate(st.clockm, nold, st.n, -1);
  rotate(st.count, nold, st.n, 0);
  rotate(st.ages, nold, st.n, 0);
  rotate_double(st.bactload, nold, st.n, 0.);
  rotate_bitarray(st.inf, nold, nblocks);
  rotate_bitarray(st.dis, nold, nblocks);
  rotate_bitarray(st.lat, nold, nblocks);
  }
}

void update_indivs(uint8_t new_i, uint8_t new_d,
		   uint8_t clearinf, double bactload,
		   int clockm, int count, int block_id) {
  int j, k;
  for (j=0; j < 8; ++j) {
    k = j + block_id * 8;
    uint8_t isnewd = (new_d << j) & '\x80';
    uint8_t isclearinf = (clearinf << j) & '\x80';
    uint8_t isnewi = (new_i << j) & '\x80';
    if ((new_d << j) & '\x80') { // is new D
      st.clockm[k] = setdtime(D_base[k], st.count[k], st.ages[k]);
      st.bactload[k] = 0.;
      continue;
    } else if ((clearinf << j) & '\x80') { // is new clearinf
      st.clockm[k] = setidtime(ID_base[k], st.count[k], st.ages[k]);
      st.bactload[k] = get_load(count[k]);
      continue;
    } else if ((new_i << j) & '\x80') { // is new I
      st.clockm[k] = setlatenttime(latent_base[k], st.count[k], st.ages[k]);
      st.count[k]++;
    }
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
	 int *clock_m, int *ages_m, int *count_m,
	 double *bactload_m, double *prob_m) {
  ages = ages_m;
  inf = inf_m;
  dis = dis_m;
  lat = lat_m;
  clockm = clock_m;
  count = count_m;
  bactload = bactload_m;
  prob = prob_m;
}

void set_base_periods(int *latent_base_m, *ID_base_m, int *D_base_m) {
  D_base = D_base_m;
  ID_base = ID_base_m;
  latent_base = latent_base_m;
}

void set_infection_parameters(double v1_m, double v2_m,
			      double phi_m, double eps_m) {
  V_1 = v1_m;
  V_2 = v2_m;
  PHI = phi_m;
  EPSILON = eps_m;
}

void set_background_mortality(double rate) {
  BGD_DEATH_RATE = rate
}

void set_groups(int groups_m[], int ngroups_m) {
  groups = groups_m;
  ngroups = ngroups_m;
}
