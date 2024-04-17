#ifndef _TRACHOMA
#define _TRACHOMA

/**
 *  Represent the state of a population of individuals.
 */
struct state {
  /** Population size */
  int n;
  /** Infected state array */
  uint8_t *inf;
  /** Diseased state array */
  uint8_t *dis;
  /** latent state array */
  uint8_t *lat;
  /** Clock array named 'clockm' to avoid clash with built-in clock()
   function */
  int *clockm;
  /** ages array */
  int *ages;
  /** Infection count array */
  int *count;
  /** Bacterial load array */
  double *bactload;
};

/**
 *  Hold pointers to records of various population properties.
 *
 *  This structure holds information about the number of records made
 *  so far.  New records are appended to previous one. For instance,
 *  after K records have been made, the ``ages`` array will be filled
 *  with K * N integers, where N is the population size.
 *
 *  .. seealso::
 *
 *     :py:class:`ntdmc_trachoma.output.Output`.
 */
struct output {
  /** Infected states */
  uint8_t *inf;
  /** Diseased states */
  uint8_t *dis;
  /** Latent states */
  uint8_t *lat;
  /** Individuals' age */
  int *ages;
  /** Current number of records */
  int *nrecords;
};
#endif
