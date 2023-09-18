#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "periods.h"
#include "shift.h"
#include "trachoma.h"

int *D_base, *ID_base, *latent_base;
int *groups, ngroups;
double BGD_DEATH_RATE;

void get_infection_prob(int*, double*, int, double, double*);
double get_bact_load(int);
void spread(struct state, uint8_t, int);
void remove_indiv(struct state, int);
void old_age_mortality(struct state, int);

/**
 *  Step the trachoma model over multiple iterations.
 *
 *  The ``step`` function is the main entry point in the
 *  ``libtrachoma`` library.  One iteration consits of the following
 *  stages:
 *
 *  - Compute infection probability for each individual
 *  - Select individuals for infection
 *  - Apply transition rules
 *  - Apply background mortality
 *  - Apply old age mortality
 *
 *  :param st: The population to evolve in time
 *  :param: out: The output container
 *  :param niter: The number of iterations to perform
 *  :param beta: The value of the beta parameter to simulate for
 */
void step(struct state st, struct output *out, int niter, double beta) {
  int i, j, t;

  int nblocks = st.n / 8;
  double *prob = (double *) malloc(st.n * sizeof(double));

  for (t = 0; t < niter; ++t) {

    if (out != NULL) {
      int write_offset = t * nblocks;
      memcpy(out->lat+write_offset, st.lat, nblocks);
      memcpy(out->inf+write_offset, st.inf, nblocks);
      memcpy(out->dis+write_offset, st.dis, nblocks);
      memcpy(out->ages+(t * st.n), st.ages, st.n * sizeof(int));
      *(out->nrecords) += 1;
    }

    get_infection_prob(st.ages, st.bactload, st.n, beta, prob);

    // if first indiv in byte is at max_age, then all indivs after are
    // as well. No need to process these blocks.
    int max_age = groups[ngroups - 1];
    for (i = 0; i < nblocks && (st.ages[i * 8] < max_age); ++i) {
      uint8_t new_i = 0, isinf, infect;
      for (j = 0; j < 8; ++j) {
	isinf = (st.inf[i] << j) & '\x80';
	int k = j + i * 8;
	infect = !isinf && ((double)rand() / RAND_MAX < prob[k]);
	new_i |= infect << (7 - j);
      }
      spread(st, new_i, i);
    } // nblocks

    /* Apply background mortality */
    for (i = 0; st.ages[i] < max_age && i < st.n; ++i) {
      if (rand() / (double)RAND_MAX < BGD_DEATH_RATE)
	remove_indiv(st, i);
      st.ages[i]++;
    }

    old_age_mortality(st, st.n - i);
  }
  free(prob);
}

void update_indivs(struct state, uint8_t, uint8_t, uint8_t, int);

void spread(struct state st, uint8_t new_i, int blk) {
  uint8_t trans = 0;
  int j, k;
  for (j=0; j < 8; ++j) {
    k = j + blk * 8;
    trans |= (((uint8_t)(!st.clockm[k])) << (7 - j));
  }
  uint8_t new_s = st.dis[blk] & ~st.inf[blk] & trans;
  uint8_t new_d = st.inf[blk] & ~st.lat[blk] & trans;
  uint8_t clearinf = st.lat[blk] & trans;

  update_indivs(st, new_i, new_d, clearinf, blk);

  st.dis[blk] = st.dis[blk] & ~new_s | clearinf;
  st.inf[blk] = st.inf[blk] & ~new_d | new_i;
  st.lat[blk] = st.lat[blk] & ~clearinf | new_i;
}

void remove_indiv(struct state st, int idx) {
  int nblocks = st.n / 8;
  bgd_death_bitarray(st.inf, idx, nblocks);
  bgd_death_bitarray(st.dis, idx, nblocks);
  bgd_death_bitarray(st.lat, idx, nblocks);
  rotate(st.clockm, 1, idx + 1, -1);
  rotate(st.ages, 1, idx + 1, 0);
  rotate(st.count, 1, idx + 1, 0);
  rotate_double(st.bactload, 1, idx + 1, 0.);
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

void update_indivs(struct state st, uint8_t new_i, uint8_t new_d,
		   uint8_t clearinf, int block_id) {
  int j, k;
  for (j=0; j < 8; ++j) {
    k = j + block_id * 8;
    if ((new_d << j) & 0x80) { // is new D
      st.clockm[k] = setdtime(D_base[k], st.count[k], st.ages[k]);
      st.bactload[k] = 0.;
      continue;
    } else if ((clearinf << j) & 0x80) { // is new clearinf
      st.clockm[k] = setidtime(ID_base[k], st.count[k], st.ages[k]);
      continue;
    } else if ((new_i << j) & 0x80) { // is new I
      st.clockm[k] = setlatenttime(latent_base[k], st.count[k]);
      st.bactload[k] = get_bact_load(st.count[k]);
      st.count[k]++;
    }
    else
      st.clockm[k]--;
  }
}



double V_1;
double V_2;
double phi;
double epsilon;

void get_infection_prob(int *ages, double* bactload, int n,
			double beta, double *prob) {
  int n_prev, n_ingroup, igroup, k;
  double lam[3], pop_ratio[3], one_over_popsize, meanld, sum_ld;
  // unsigned int groups[] = {468, 780, 3121};
  double epsm, *A;

  k = 0;
  one_over_popsize = 1. / n;
  epsm = 1. - epsilon;

  for (k = 0, igroup = 0; igroup < ngroups; ++igroup) {
    n_prev = k;
    for(sum_ld = 0; (ages[k] < groups[igroup]) &&  k < n; ++k)
      sum_ld += bactload[k];

    n_ingroup = (k - n_prev);
    meanld = sum_ld / n_ingroup;
    lam[igroup] = beta * V_1 * meanld + beta * V_2 * pow(meanld, phi + 1);
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

double get_bact_load(int ninf) {
  double b1 = 1., ep2 =  - 0.114;

  return b1 * exp((ninf - 1) * ep2);
}

/**
 *  Set global pointers to memory describing the base infection periods.
 *
 *  Each individual in the population is assigned a base value for the
 *  time spent in each of the infection stages.  This value will
 *  typically be used for setting the period of time the individual
 *  will stay in that stage following the first infection.  Period
 *  values for subsequent infection will be function of the base
 *  values.
 *
 *  Note that, as oppsesed to infected state arrays, base periods
 *  arrays are not members of the State structure.  Instead they are
 *  set globally.
 *
 *  This function is meant to be used by Python code to set pointers
 *  to data handle by NumPy arrays:
 *
 *  .. code-block:: python
 *
 *      from numpy.ctypeslib import ndpointer
 *
 *      clib = ctypes.CDLL("/path/to/trachoma/c/lib.so")
 *      clib.set_base_periods.argtypes = [
 *              ndpointer(dtype=np.int32, ndim=1),
 *              ndpointer(dtype=np.int32, ndim=1),
 *              ndpointer(dtype=np.int32, ndim=1),
 *          ]
 *      clib.set_base_periods.restype = None
 *      latent_base = numpy.array([2] * 8)  # let NumPy handle memory
 *      ID_base = numpy.array([2] * 8)
 *      D_base = numpy.array([2] * 8)
 *      clib.set_base_periods(latent_base, ID_base, D_base)
 *
 *  :param latent_base_m: A pointer to a memory location containing
 *    the base value for the latent period of each individual in the
 *    population
 *  :param ID_base_m: A pointer to a memory location containing
 *    the base value for the *ID* period of each individual in the
 *    population
 *  :param D_base_m: A pointer to a memory location containing
 *    the base value for the *D* period of each individual in the
 *    population
 */
void set_base_periods(int *latent_base_m, int *ID_base_m, int *D_base_m) {
  D_base = D_base_m;
  ID_base = ID_base_m;
  latent_base = latent_base_m;
}

/**
 *  Set global variales to values related to infection probability.
 *
 *  :param v1_m: Value for the :math:`v_1` infection parameter.
 *  :param v1_m: Value for the :math:`v_2` infection parameter.
 *  :param v1_m: Value for the :math:`\phi` infection parameter.
 *  :param v1_m: Value for the :math:`\epsilon` infection parameter.
 */
void set_infection_parameters(double v1_m, double v2_m,
			      double phi_m, double eps_m) {
  // TODO: Call parameters by their actual name in docstring
  V_1 = v1_m;
  V_2 = v2_m;
  phi = phi_m;
  epsilon = eps_m;
}

/**
 *  Set value of background mortality probability.
 *
 *  The background mortality probability is expressed as
 *
 *  .. math::
 *
 *     r = 1 - e^{1 / \tau}
 *
 *  where :math:`\tau` is the background mortality rate.
 *
 *  :param prob: Mortality probability value
 */
void set_background_mortality(double prob) {
  // TODO: Take rate as a parameter instead of probablity.
  BGD_DEATH_RATE = prob;
}

/**
 *  Set global pointer to memory describing the age groups upper
 *  boundaries.
 *
 *  .. note::
 *
 *     The minimum age is ``0``.
 *
 *  :param groups_m: pointer to an array describing the age grouops
 *    upper boundaries.  The last element of these array should be the
 *    maximum age any individual in the population can get to.
 *  :param ngroups_m: the number of age groups.
 */
void set_groups(int groups_m[], int ngroups_m) {
  groups = groups_m;
  ngroups = ngroups_m;
}
