import numpy as np
from parameters import MIN_ID, MIN_D, INF_RED, AQ, AG


def id_time(base_period, infection_count):
    return np.round(
        (base_period - MIN_ID)
        * np.exp(-INF_RED * (infection_count - 1))
        + MIN_ID
    )


def d_time(base_period, infection_count, ages):
    return np.round(
        (base_period - MIN_D)
        * np.exp(AQ * (1 - infection_count) - AG * ages)
        + MIN_D
        )

def latent_time(base_period):
    return base_period
        
