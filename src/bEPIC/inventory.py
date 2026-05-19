import pandas as pd


def get_fdsn_station_inventory(origin_time, networks='*', client='IRIS'):
    """
    Query an FDSN web service for stations active at origin_time.

    Intended for case-study use: call once at the start of a mainshock/
    aftershock sequence and reuse the returned DataFrame for all subsequent
    events.  For large benchmark catalogs spanning years use
    seismic_benchmark.benchmark.runner.get_unique_stations() instead.

    Parameters
    ----------
    origin_time : obspy.UTCDateTime, datetime, or str
        The epoch at which station availability is evaluated.
    networks : str
        FDSN network wildcard (default '*' = all networks).
    client : str
        FDSN client identifier or URL (default 'IRIS').

    Returns
    -------
    pandas.DataFrame with columns: station, network, longitude, latitude.
        Suitable for passing directly to EPIC_PARAMS.station_inventory.
    """
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime

    t = UTCDateTime(origin_time)
    c = Client(client)
    inv = c.get_stations(
        network=networks,
        station='*',
        starttime=t,
        endtime=t,
        level='station',
    )
    rows = []
    for net in inv:
        for sta in net:
            rows.append({
                'station':   sta.code,
                'network':   net.code,
                'longitude': sta.longitude,
                'latitude':  sta.latitude,
            })
    return pd.DataFrame(rows, columns=['station', 'network', 'longitude', 'latitude'])
