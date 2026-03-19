"""
generate_synthetic_run_file.py

Generates a synthetic .run file for testing bEPIC, by simulating P-wave
trigger times at real-ish seismic stations for a known earthquake location.

Synthetic event: M4.0 near Berkeley, CA
"""

import numpy as np
import pandas as pd
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Synthetic event parameters
# ---------------------------------------------------------------------------
EQ_LON      = -122.27   # degrees
EQ_LAT      =   37.87   # degrees
EQ_DEPTH_KM =    8.0    # km (matches assumed depth in calculate_likelihood)
EQ_ORIGIN   = datetime(2024, 1, 15, 10, 30, 0, tzinfo=timezone.utc).timestamp()
VELOCITY    =    6.0    # km/s (constant velocity model)

# ---------------------------------------------------------------------------
# Synthetic stations (real-ish NCEDC stations in Northern California)
# ---------------------------------------------------------------------------
stations = [
    {'station': 'BKS',  'channel': 'BHZ', 'network': 'BK', 'location': '--', 'longitude': -122.2352, 'latitude': 37.8762},
    {'station': 'OAKD', 'channel': 'HHZ', 'network': 'BK', 'location': '--', 'longitude': -122.2153, 'latitude': 37.7719},
    {'station': 'MCCM', 'channel': 'HHZ', 'network': 'NC', 'location': '--', 'longitude': -122.4130, 'latitude': 37.9938},
    {'station': 'SUTB', 'channel': 'HHZ', 'network': 'NC', 'location': '--', 'longitude': -122.4786, 'latitude': 37.7607},
    {'station': 'FARB', 'channel': 'HHZ', 'network': 'BK', 'location': '--', 'longitude': -123.0018, 'latitude': 37.6978},
    {'station': 'WENL', 'channel': 'HHZ', 'network': 'NC', 'location': '--', 'longitude': -121.9617, 'latitude': 37.6142},
]

# ---------------------------------------------------------------------------
# Compute synthetic trigger times
# ---------------------------------------------------------------------------
def haversine_km(lon1, lat1, lon2, lat2):
    """Horizontal distance in km between two lon/lat points."""
    R = 6371.0
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(a))

rows = []
for version, sta in enumerate(stations):
    horiz_dist = haversine_km(EQ_LON, EQ_LAT, sta['longitude'], sta['latitude'])
    total_dist = np.sqrt(horiz_dist**2 + EQ_DEPTH_KM**2)  # slant distance
    travel_time = total_dist / VELOCITY

    # Add small random noise to simulate real pick uncertainty (~0.1s)
    noise = np.random.normal(0, 0.1)
    trigger_time = EQ_ORIGIN + travel_time + noise
    trigger_time_str = datetime.utcfromtimestamp(trigger_time).strftime('%Y-%m-%dT%H:%M:%S.%f')

    rows.append({
        'version':      version,        # each station arrives as a new version
        'order':        version,
        'station':      sta['station'],
        'channel':      sta['channel'],
        'network':      sta['network'],
        'location':     sta['location'],
        'longitude':    sta['longitude'],
        'latitude':     sta['latitude'],
        'trigger time': trigger_time,
        'tterr':        noise,           # residual — noise we added
        'logPd':        np.random.uniform(-1.5, -0.5),  # synthetic log peak displacement
    })

df = pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
postgres_id        = '999999'
output_dir         = f'/tmp/bEPIC_test/{postgres_id}/'
output_filename    = output_dir + postgres_id + '.run'

import os
os.makedirs(output_dir, exist_ok=True)
df.to_csv(output_filename, index=False)

print(f'Synthetic .run file written to {output_filename}')
print(f'\nTrue epicenter: {EQ_LON}, {EQ_LAT}')
print(f'Origin time:    {datetime.utcfromtimestamp(EQ_ORIGIN)}')
print(f'\n{df.to_string(index=False)}')