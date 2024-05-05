"""Provide access to Python classes representing different types of
events, currently:

- MDA intervention
- Simulation start
- Simulation end

Events are represented as Callable classes, that is regular Python
classes for which a __call__ method is implemented. This method *must*
take a `Population` object as only argument. It is expected to modify
this object in-place as part of modelling the impact of the event.

Instance of event classes are supposed to represent all occurences of
the event for a given intervention. For example, the events list for a
annual MDA campaign from 2030 to 2040 will contain 10 references to
the *same* MDA instance.

.. note::

    By construct, all MDA in a given intervention (campaign) have the
    same parameters.

Event-specific ``Population`` attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Event objects often need to keep memory of their previous impact on a
population.  For example a MDA object, to work on a population, must
be aware of the number of times an individual received the treatment
in the past.  A vaccination event would likely must be able to tell
whether an individual in the population has received one or more doses
already.  Events objects do not keep any state related any Population
they work on. Instead, they 'attach' information to the `Population`
instance as a event object specific attribute. The attribute is named
according to the format:

    <event type>_<event object id>_<generic attribute name>

where `<event object id>` is the underlying identifier of the Python
object, obtained with `id()`. Therefore the attribute is specific to
the event *instance*, not the class itself. It means two different
instances of the same event class can attach independant information
to the same population object. for instance two different MDA
campaigns. Event attributes are created the first time the event
object processes the `Population` object, and only update subsequent
times.

Example:
~~~~~~~~

>>> rng = numpy.random.default_rng()
>>> pop = Population(ages=ranges(1,17), latent_base=[16] * 2, rng)
>>> mdaevent = MDA(coverage=0.8, efficacy=0.85, rho=0.78)
>>> pop.mda_treatment_count
AttributeError: 'MDA' object has no attribute 'mda_treatment_count'
>>> mdaevent(pop)
>>> pop.mda_treatment_count
np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 ,0, 0, 0, 1, 0], dtype=np.int32)

"""
# TODO: Rename 'mda' module as 'events'
from dataclasses import dataclass

import numpy as np

from state import Population


# When an event doesn't need to hold any state, we can just define as
# a function. After all, a function is a callable class.
def start_sim(pop):
    pass


def end_sim(pop):
    pass


@dataclass
class MDA_legacy_data:
    ages: np.array
    treatment_count: np.array


class MDA_legacy:
    r"""Represent a MDA event.

    MDA events attach a `treatment_count` attribute to Population
    objets they are applied to.  This attribute is a NumPy integer
    array keeping track of how many times each individual in the
    population was administered the treatment.

    The probability of treatment is function of how many times this
    individual received the treatment in the past. Denoting coverage
    and compliance levels by :math:`C` and :math:`\rho` respectively,
    it is expressed as

    .. math::

       p = \frac{C(1 - \rho) + m\rho}{1 + \rho (k - 1)}

    where :math:`m` and `k` represent the individual's current
    treatment count and how many this the MDA event was performed in
    the past, respectively.

    To minimise the number of computations, we recast the above expressions into

    .. math::

       p = \frac{b + m}{a + k}

    where :math:`a = \rho^{-1} - 1` and :math:`b = a * C`.

    :param coverage: Coverage level
    :type coverage: number between 0 and 1
    :param efficacy: Treatment efficacy value
    :type efficacy: number between 0 and 1
    :param rho: Compliance level
    :type rho: number between 0 and 1
    """
    def __init__(self, coverage, efficacy, rho):
        self.a = (1. / rho) - 1
        self.coverage = coverage
        self.b = self.a * coverage
        self.rng = np.random.default_rng()
        self.count = 0
        self.efficacy = efficacy

    def distribute_treatment(self, popsize, treatment_count):
        p = (
            self.b + treatment_count
        ) / (self.a + self.count)
        return (self.rng.uniform(size=popsize) < p)

    def __call__(self, pop):
        # Start by checking whether or not the MDA instance registered
        # itself to the Population object passed to __call__. If it,
        # distribute treatment and update attribute. If not, set the
        # attribute.
        if hasattr(pop, "mda_legacy_data"):
            treatment_count = pop.mda_legacy_data.treatment_count
            # We need to determine which individuals were reset since
            # the last __call__, and fix the current value of their
            # treatment count. I.e. it should be set to 0.

            # A way to check whether or not an individual was reset is
            # to check if they are younger than they were last time!
            treatment_count[pop.ages - pop.mda_ages < 0] = 0
            t = self.distribute_treatment(pop.size, treatment_count)
            treatment_count[t] += 1
        else:
            # First MDA application is actually a special case where
            # treatment probablity is he coverage value itself.
            t = self.rng.uniform(size=pop.size) < self.coverage
            setattr(pop, attrname, t.astype(np.int32))
        pop.mda_legacy_data.ages = pop.ages.copy()

        ntreated = np.count_nonzero(t)
        cured = np.zeros(pop.size, dtype=np.bool_)
        cured[t] = self.rng.uniform(size=ntreated) < self.efficacy
        # At this point it's tempting to write
        #
        # pop.inf[cured] = True
        #
        # but remember that pop.inf is in fact a call to the `inf`
        # property getter so it won't actually modify the internal
        # _inf attribute.
        inf = pop.inf
        inf[cured] = False
        pop.inf = inf
        pop.bact_load[cured] = 0

        self.count += 1

        return ntreated, np.count_nonzero(cured)


@dataclass
class MDA_data:
    ages: np.array
    treatment_probability: np.array
    beta_dist_params: tuple[float]


class MDA:
    def __init__(self, coverage: float,
                 correlation: float, efficacy: float,
        ):
        self.rng = np.random.default_rng()
        self.a = coverage * (1. - correlation) / correlation
        self.b = (1. - coverage) * (1. - correlation) / correlation
        self.efficacy = efficacy

    def __call__(self, pop: Population):
        if not hasattr(pop, "MDA_data"):
            pop.MDA_data = MDA_data(
                ages=pop.ages,
                treatment_probability = self.tment_prob,
                beta_dist_params = (self.a, self.b)
            )

        was_reset = pop.MDA_data.ages < pop.ages
        tment_prob_for_resets = np.beta(
            *MDA_data.beta_dist_params,
            size=np.count_nonzero(was_reset),
        )
        if pop.MDA_data.beta_dist_params != (self.a, self.b):
            cur_prob = pop.MDA_data.tment_prob.copy()
            cur_prob[was_reset] = tment_prob_for_resets
            self.tment_prob = self.tment_prob[np.argsort(cur_prob)]
        else:
            self.tment_prob[was_reset] = tment_prob_for_resets
        pop.MDA_data = MDA_data(
            ages=pop.ages,
            treatment_probability = self.tment_prob,
            beta_dist_params = (self.a, self.b)
        )

        t = self.rng.uniform(size=pop.size) < self.tment_prob
        ntreated = np.count_nonzero(t)
        cured = np.zeros(pop.size, dtype=np.bool_)
        cured[t] = self.rng.uniform(size=ntreated) < self.efficacy
        # At this point it's tempting to write
        #
        # pop.inf[cured] = True
        #
        # but remember that pop.inf is in fact a call to the `inf`
        # property getter so it won't actually modify the internal
        # _inf attribute.
        inf = pop.inf
        inf[cured] = False
        pop.inf = inf
        pop.bact_load[cured] = 0

        return ntreated, np.count_nonzero(cured)
