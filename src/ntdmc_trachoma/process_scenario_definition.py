import json
from pathlib import Path
import sys

import ntdmc_trachoma.mda as mda


def process_scenario_definition(path: Path):
    """Read a scenario definition file and build a list of events.

    A scenario defines a list of programs, themselves defining
    one or several lists of interventions.  For instance, a program
    could define a bi-annual MDA intervention between 2024 and 2034,
    as well as a annual vaccination campaign between 2030 and 2034.

    The list of events is returned sorted according to increasing
    dates.  The final list always start with a "simulation start"
    event and ends with a "simulation end" event. See
    :py:func:`events.start_sim <ntdmc_trachoma.mda.start_sim>`.
    
    :param path: Location of a valid scenario definition file.
    :type path: pathlib.Path

    """
    with path.open() as f:
        d = json.load(f)
    all_events = [
        (d["start"], mda.start_sim),
        (d["end"], mda.end_sim),
    ]
    for p in d["programs"]:
        start = p["start"]
        end = p["end"]

        # Generate list of dates for events before end of program is
        # reached. Not so useful when 'dates' are only integers, but
        # ultimately we might want to work with actual datetime.date
        # objects.
        def gen_dates(interval):
            date = start
            while (date < end):
                yield date
                date += interval

        # For each intervention in the program dict, generate a list
        # of event and append them in the list of all events.
        for events in map(
                lambda x: create_events_for_intervention(x, gen_dates),
                p["interventions"]
        ):
            all_events.extend(events)
    # TODO: What if first event is at the same time as
    # start_sim event?
    return sorted(all_events, key=lambda x: x[0])


def create_events_for_intervention(inter: dict, gen_dates):
    # Start by checking that the 'type' entry dictionary corresponds
    # to a class impelemetend in the 'mda' module.
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
