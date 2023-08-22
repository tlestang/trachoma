from dataclasses import dataclass
from collections import namedtuple
import json


@dataclass
class SimulationParameters:
    population_size: int
    max_age: int
    mean_age: int
    v1: float
    v2: float
    phi: float
    epsilon: float
    av_I_duration: float
    av_ID_duration: float
    av_D_duration: float
    inf_red: float
    min_d: int
    min_id: int
    ag: float
    aq: float
    tau: int
    groups: list


def get_params(path):
    with open(path, 'r') as f:
        d = json.load(f)
    return SimulationParameters(**d)

AverageDurations = namedtuple(
    typename="AverageDurations", field_names=["I", "ID", "D"],
)
