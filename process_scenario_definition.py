from datetime import datetime, timedelta
import json
import sys

import mda

with open("scenario.json", "r") as f:
    s = json.load(f)

events = []
for p in s["programs"]:
    start = datetime.fromisoformat(s["start"])
    end = datetime.fromisoformat(s["end"])
    for inter in s["interventions"]:
        t = inter.pop("type")
        try:
            cls = getattr(t, mda)
        except AttributeError:
            sys.stderr.write(
                f"Error: could not find event type {t}.\n"
            )
            sys.exit(1)
        d = inter.pop("interval")
        interval = timedelta(weeks=d["weeks"])
        event = cls(**inter)
        events.append([
            (event, start + i * interval)
            for i in range(0, (start-end).days, interval.days)
        ])
events_flattened = []
for e in events:
    events_flattened.extend(e)
events_flattened.sort(key=lambda x: x[1])
            
        
