#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 10:15:50 2022

@author: amy

"""

from bEPIC import data_util
from bEPIC import locate,likelihood,magnitude,prior,posterior
import os
import shutil
import numpy as np
import pandas as pd
bepic = "/home/daned/2024_NEHRP/bEPIC"

# initialize by seeing if everything is in the correct directory structure and if a .run file
# exists
 
def initialize_bEPIC_event(project_parent_directory, postgres_id):
    """
    Prepares the directory structure and input files for a bEPIC event.

    - Creates required output subdirectories (EPIC, plots, bEPIC, USGS) if missing
    - Moves any existing EPIC log files into the EPIC/ subdirectory
    - Generates a .run file if one doesn't exist
    - Attempts to fetch a USGS event summary if one isn't present

    Args:
        project_parent_directory (str): Root path containing per-event project folders.
        postgres_id (str|int):          Unique event identifier (zero-padded to 6 digits if int).
    """

    if type(postgres_id) != str:
        postgres_id = str(postgres_id).zfill(6)

    event_dir = project_parent_directory + postgres_id

    if not os.path.exists(event_dir):
        print('directory does not exist')
        return

    # --- Create required output subdirectories if missing ---
    for subdir in ['EPIC', 'plots', 'bEPIC', 'USGS']:
        subdir_path = f"{event_dir}/{subdir}/"
        if not os.path.exists(subdir_path):
            print(f' ... making a {subdir} output directory')
            os.mkdir(subdir_path)

    # --- Move any existing EPIC log files into the EPIC/ subdirectory ---
    epic_log_files = [
        f"{postgres_id}_epic_location_log.txt",
        f"{postgres_id}_event_summary_log.txt",
        f"{postgres_id}_station_summary_log.txt",
        f"{postgres_id}_station_counts_log.txt",
        f"{postgres_id}_event_triggers_log.txt",
        f"{postgres_id}_location_triggers_log.txt",
        f"{postgres_id}_misc_log.txt",
    ]

    for filename in epic_log_files:
        src  = f"{event_dir}/{filename}"
        dest = f"{event_dir}/EPIC/{filename}"
        if os.path.exists(src):
            shutil.move(src, dest)

    # --- Generate a .run file if one doesn't exist ---
    run_file = f"{event_dir}/{postgres_id}.run"
    if not os.path.exists(run_file):
        print(f' ... event needs a .run file')
        print(f' ... creating a .run file for event id {postgres_id}')
        data_util.generate_run_file(project_parent_directory, postgres_id)

    # --- Fetch USGS event summary if not already present ---
    usgs_summary = f"{event_dir}/USGS/usgs_event_summary.txt"
    if not os.path.exists(usgs_summary):
        print(' ... no USGS event found ... attempting to find one')
        try:
            data_util.search_for_USGS_event(project_parent_directory, postgres_id)
        except Exception as e:
            print(f' ... USGS search failed: {e}')


def run_bEPIC(project_parent_directory, postgres_id, velocity_model, GridSize, GridSpacing):
    """
    Main entry point for the bEPIC earthquake location pipeline.

    Iterates over successive versions (station configurations) of an event,
    computing a likelihood, prior, and posterior location estimate at each step.
    Results are saved to a tab-separated log file.

    Args:
        project_parent_directory (str): Root path containing per-event project folders.
        postgres_id (str|int):          Unique event identifier (zero-padded to 6 digits if int).
        velocity_model:                 Velocity model passed to the likelihood calculator.
        GridSize (float):               Spatial extent of the search grid.
        GridSpacing (float):            Resolution of the search grid.
    """

    # Ensure postgres_id is a zero-padded 6-character string (e.g. 42 -> '000042')
    if type(postgres_id) != str:
        postgres_id = str(postgres_id).zfill(6)

    # Load the .run file — contains per-station arrival time data across all versions
    run_df = pd.read_csv(project_parent_directory + postgres_id + '/' + postgres_id + '.run')
    run_df['sigma'] = np.ones(len(run_df))  # initialize uncertainty weights to 1

    # ---------------------------------------------------------------------------
    # Main loop: iterate over each version (each represents a new station coming online)
    # ---------------------------------------------------------------------------
    version = 0 # Version is the current iteration - current station
    while version <= np.max(np.unique(run_df['version'])):
        print('postgres id: ', postgres_id, '| version: ', version, '|')

        relocate = True  # default to recomputing location at each version

        # Get the subset of stations that were active at this version, excluding bad picks (tterr > -999)
        idx    = np.where((run_df['version'] == version) &
                          (run_df['tterr'] > -999))[0]
        
        # Get the df of just these subset of stations
        sta_df = run_df.iloc[idx].reset_index(drop=True)

        if version == 0:
            # For the first version, estimate initial location from the first two stations
            CenterPoint = locate.get_two_station_location(sta_df)

            # Initialize the results dataframe
            columns = ['version', 'num stations', 'likelihood lon', 'likelihood lat',
                       'likelihood mag', 'prior lon', 'prior lat',
                       'posterior lon', 'posterior lat', 'posterior mag']
            bEPIC_df = pd.DataFrame(columns=columns)
        else:
            # For subsequent versions, center the grid on the previous posterior estimate
            # Centerpoint - [float,float]
            CenterPoint = [bEPIC_df['posterior lon'].iloc[-1], bEPIC_df['posterior lat'].iloc[-1]]

            # Skip relocation if the station count hasn't changed since the last version
            if len(sta_df) == bEPIC_df['num stations'].iloc[-1]:
                relocate = False

        if relocate:
            # --- Likelihood: compute the data-driven location probability surface ---
            (likelihood_function, misfit,
             likelihood_lon, likelihood_lat) = likelihood.calculate_likelihood(
                CenterPoint, sta_df, velocity_model, GridSize, GridSpacing)

            # --- Prior: compute the seismicity-based prior probability surface ---
            prior_function, prior_lon, prior_lat = prior.compute_prior(
                CenterPoint, GridSize, GridSpacing, ANSS_timestamp=None)

            # --- Posterior: combine likelihood and prior to get the best location estimate ---
            post, posterior_lon, posterior_lat = posterior.compute_posterior(
                CenterPoint, GridSize, GridSpacing, prior_function, likelihood_function)

        else:
            # Reuse previous version's locations if station configuration hasn't changed
            likelihood_lon = bEPIC_df['likelihood lon'].iloc[-1]
            likelihood_lat = bEPIC_df['likelihood lat'].iloc[-1]
            posterior_lon  = bEPIC_df['posterior lon'].iloc[-1]
            posterior_lat  = bEPIC_df['posterior lat'].iloc[-1]

        # Compute magnitude estimates at both the likelihood and posterior locations
        likelihood_mag = magnitude.compute_magnitude(run_df, version, [likelihood_lon, likelihood_lat])
        posterior_mag  = magnitude.compute_magnitude(run_df, version, [posterior_lon, posterior_lat])

        # Append this version's results to the log dataframe
        bEPIC_df.loc[len(bEPIC_df.index)] = [
            version, int(len(sta_df)),
            likelihood_lon, likelihood_lat, likelihood_mag,
            prior_lon, prior_lat,
            posterior_lon, posterior_lat, posterior_mag
        ]

        version += 1

    # Save the full results log to disk as a tab-separated file
    bEPIC_df.to_csv(
        project_parent_directory + postgres_id + '/bEPIC/' + postgres_id + '_bEPIC_log.txt',
        sep='\t', index=False
    )