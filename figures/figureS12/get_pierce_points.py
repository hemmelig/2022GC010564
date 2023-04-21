# -*- coding: utf-8 -*-
"""
Short script to evaluate the 100-km depth pierce points for the measurements
presented in:

    Bacon, C. A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D. G., & Tongkul, F. (2022). The signature of lithospheric
    anisotropy at post-subduction continental margins: New insight from
    XKS splitting analysis in northern Borneo. Geochemistry, Geophysics,
    Geosystems, 23, e2022GC010564. https://doi.org/10.1029/2022GC010564

"""

import subprocess
from typing import Tuple

import pandas as pd


def parse_pierce_point_from_arrival(arrival: str) -> Tuple[float, float]:
    """
    Extracts geographic information about the 100 km depth pierce point from
    the output of the TauP pierce program.

    Parameters
    ----------
    arrival: Output from TauP pierce dumped as a string.

    Returns
    -------
    pierce_lon, pierce_lat: Longitude and latitude of the pierce point.

    """

    # Parse out 100 km pierce points for each phase path
    arrivals = {}
    for line in arrival.split("\n"):
        try:
            if line[0] == ">":
                arrival = f"{line.split()[1]}_{line.split()[3]}"
                arrivals.setdefault(arrival, {})
                continue
            line = line.split()
            if line[1] == "100.0":
                arrivals[arrival].update({"traveltime": line[0]})
                arrivals[arrival].update({"lon": line[4]})
                arrivals[arrival].update({"lat": line[3]})
        except IndexError:
            break

    # Select the fastest phase path and the receiver-side pierce point
    keys = list(arrivals.keys())
    phase = keys[0].split("_")[0]
    key = f"{phase}_{min(key.split('_')[1] for key in keys)}"

    return arrivals[key]["lon"], arrivals[key]["lat"]


def main():
    measurements = pd.read_csv(
        "../data/splitting_results/supp_file_6_all_splitting_measurements.csv"
    )
    stations = pd.concat(
        [
            pd.read_csv("../data/station_files/supp_file_1_YC.sta", header=None),
            pd.read_csv("../data/station_files/supp_file_2_MY.sta", header=None),
        ]
    )

    pierce_lons, pierce_lats = [], []
    for _, measurement in measurements.iterrows():
        # Get source information
        source_loc = (measurement.lon, measurement.lat, measurement.depth)
        source_lon, source_lat, source_depth = source_loc

        # Get station information
        station_name = measurement.station
        station_name = station_name if station_name != "MALB" else "SBB5a"
        station_info = stations[stations[1] == station_name].values[0]
        receiver_lon, receiver_lat = station_info[3], station_info[2]

        # Get phase information
        phase = measurement.phase

        cmd = (
            "taup pierce "
            f"-evt {source_lat} {source_lon} "
            f"-sta {receiver_lat} {receiver_lon} "
            f"-h {source_depth} -ph {phase}, -pierce 100.0"
        )
        arrival = subprocess.check_output(cmd, shell=True, universal_newlines=True)

        pierce_lon, pierce_lat = parse_pierce_point_from_arrival(arrival)
        pierce_lons.append(pierce_lon)
        pierce_lats.append(pierce_lat)

    measurements["pierce_lon"] = pierce_lons
    measurements["pierce_lat"] = pierce_lats

    measurements.to_csv("measurements_with_pierce_points.csv", index=False)


if __name__ == "__main__":
    main()
