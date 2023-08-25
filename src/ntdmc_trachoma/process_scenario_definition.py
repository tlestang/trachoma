import json
from pathlib import Path
import sys

import ntdmc_trachoma.mda as mda


def process_intervention(inter: dict, gen_dates):
    t = inter.pop("type")
    try:
        cls = getattr(mda, t)
    except AttributeError:
        sys.stderr.write(
            f"Error: could not find event type {t}.\n"
        )
        return []
    interval = inter.pop("interval")
    event = cls(**inter)
    return [
        (date, event) for date in gen_dates(interval)
    ]


def process_scenario_definition(path: Path):
    with path.open() as f:
        d = json.load(f)
    all_events = [
        (d["start"], mda.start_sim),
        (d["end"], mda.end_sim),
    ]
    for p in d["programs"]:
        start = p["start"]
        end = p["end"]

        def gen_dates(interval):
            date = start
            while (date < end):
                yield date
                date += interval

        for events in map(
                lambda x: process_intervention(x, gen_dates),
                p["interventions"]
        ):
            all_events.extend(events)
    # TODO: What if first event is at the same time as
    # start_sim event?
    return sorted(all_events, key=lambda x: x[0])
