from datetime import datetime, timedelta
import sys

import mda


def process_intervention(inter: dict, gen_dates):
    t = inter.pop("type")
    try:
        cls = getattr(mda, t)
    except AttributeError:
        sys.stderr.write(
            f"Error: could not find event type {t}.\n"
        )
        return []
    interval = timedelta(weeks=inter.pop("interval"))
    event = cls(**inter)
    return [
        (event, date) for date in gen_dates(interval)
    ]


def process_scenario_definition(d: dict):
    all_events = []
    for p in d["programs"]:
        start = datetime.fromisoformat(p["start"])
        end = datetime.fromisoformat(p["end"])

        def gen_dates(interval):
            date = start
            while(date <= end):
                yield date
                date += interval

        for events in map(
                lambda x: process_intervention(x, gen_dates),
                p["interventions"]
        ):
            all_events.extend(events)
    return sorted(all_events, key=lambda x: x[1])
