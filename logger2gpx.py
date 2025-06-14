# logger2gpx.py

"""Extract GPX track data from a sqlite file containing position
data. The sqlite file should contain a table named `location_data` with
columns `timestamp`, `latitude`, `longitude`.

The program offers the opportunity to filter the data based on a minimum distance
between points.
"""

import argparse
import math
import os
import sqlite3
from datetime import datetime
from typing import Iterator
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, tostring

import yaml


def fetch_data(db_path: str,
               start_time: float | None = None,
               end_time: float | None = None) -> Iterator[tuple[float|int, float, float]]:
    """
    Generator function that fetches track data from the database within an optional time range,
    ordered by timestamp.

    This function queries the database table `location_data` to retrieve data including
    timestamp, latitude, and longitude. Optional parameters `start_time` and
    `end_time` can be supplied to filter the results within a specific time range. If no
    time range is provided, all available records in the table are retrieved.

    Parameters
        db_path: Path to the SQLite database file
        start_time: Optional Unix epoch time for start of range. All returned records will be
            after this time
        end_time: Optional Unix epoch time for end of range. All returned records will before or
            at this time

    Yields
        Tuples (timestamp, latitude, longitude) ordered by timestamp
    """

    params = []
    conditions = []

    # Build query with optional time constraints
    query = "SELECT timestamp, latitude, longitude FROM location_data"

    if start_time:
        conditions.append("timestamp > ?")
        params.append(start_time)
    if end_time:
        conditions.append("timestamp <= ?")
        params.append(end_time)

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    query += " ORDER BY timestamp"

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(query, params)
        for row in cursor:
            yield row



def filter_data(data: Iterator[tuple[float|int, float, float]],
                min_distance: float) -> Iterator[tuple[float | int, float, float]]:
    """
    Generator function that filters a dataset of geographic coordinates such
    that consecutive points meet a minimum distance criterion.

    Parameters
        data: An iterator that returns tuples, where each tuple consists of
            (timestamp, latitude, longitude). None
            for invalid or missing coordinates.
        min_distance: The minimum distance in meters.  Zero (0) to include all points.
`
    Yields
        tuples (timestamp, latitude, longitude) ordered by timestamp, containing only points
        that are either the first valid point or exceed the specified minimum
        distance from the previous valid point.
    """

    last_lat, last_lon = None, None

    for row in data:
        # Unpack the row
        timestamp, lat, lon = row

        # Skip points with NULL coordinates
        if lat is None or lon is None:
            continue

        # Apply minimum distance filter
        if min_distance > 0 and last_lat is not None and last_lon is not None:
            distance = calculate_distance(last_lat, last_lon, lat, lon)
            if distance < min_distance:
                continue

        yield row
        last_lat, last_lon = lat, lon


def create_gpx(config_dict: dict,
               track_data: Iterator[tuple[float | int, float, float]],
               start_time_str: str | None = None,
               end_time_str: str | None = None) -> Element:
    """
    Creates a GPX (GPS Exchange Format) XML structure containing metadata and track segments.

    This function generates a GPX document based on the specified configuration and
    track data. It includes metadata such as the GPX name, description, and track name.
    The GPX document also includes track segments with points (latitude, longitude,
    and timestamp) derived from the provided track data.

    Parameters
        config_dict: A dictionary containing configuration information, including
            details for generating author, GPX name, description, and track name. These
            are used to populate the metadata and track details of the GPX structure.
        track_data: An iterator that returns tuples (timestamp, latitude, longitude).
        start_time_str: An optional string representing the start time. If
            provided, it is incorporated into the metadata description. Defaults to None.
        end_time_str: An optional string representing the end time. If
            provided, it is incorporated into the metadata description. Defaults to None.

    Returns
        The root XML element of the GPX structure, representing the generated GPX document.
    """

    # Get interpolation dictionary:
    interp_dict = config_dict['gpx']

    # Create root GPX element
    gpx = Element('gpx')
    gpx.set('version', '1.1')
    gpx.set('creator', config_dict['gpx']['author'])
    gpx.set('xmlns', 'http://www.topografix.com/GPX/1/1')

    # Add metadata
    metadata = SubElement(gpx, 'metadata')
    name = SubElement(metadata, 'name')
    name.text = config_dict['gpx']['gpx_name'].format_map(interp_dict)

    desc = SubElement(metadata, 'desc')
    desc.text = config_dict['gpx']['description'].format_map(interp_dict)
    if start_time_str and end_time_str:
        desc.text += f" (From {start_time_str} to {end_time_str})"

    # Create track
    trk = SubElement(gpx, 'trk')
    trk_name = SubElement(trk, 'name')
    trk_name.text = config_dict['gpx']['track_name'].format_map(interp_dict)

    # Create track segment
    trkseg = SubElement(trk, 'trkseg')

    # Add track points
    for row in track_data:
        timestamp, lat, lon = row

        if lat is not None and lon is not None:
            trkpt = SubElement(trkseg, 'trkpt')
            trkpt.set('lat', str(f"{lat:.5f}"))
            trkpt.set('lon', str(f"{lon:.5f}"))

            # Add timestamp
            time_elem = SubElement(trkpt, 'time')
            time_elem.text = timestamp_to_iso(timestamp)

    return gpx


def save_gpx(gpx_elements: Element, output_dir: str, filename: str) -> str:
    """
    Saves a GPX (GPS Exchange Format) XML element to a file with pretty-printed formatting.
    If the specified output directory does not exist, it creates the directory.

    Parameters
        gpx_elements: The GPX XML element to be saved.
        output_dir: The directory where the file will be saved.
        filename: The name of the file to save the GPX content.

    Returns
        The full file path of the saved GPX file.
    """
    # Pretty print XML
    rough_string = tostring(gpx_elements, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ", encoding='utf-8')

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)

    # Save to file
    with open(filepath, 'wb') as f:
        f.write(pretty_xml)

    return filepath


# ===============================================================================
#                           UTILITY FUNCTIONS
# ===============================================================================


def timestamp_to_iso(timestamp: float | int) -> str:
    """
    Convert a Unix timestamp to an ISO 8601 formatted string.

    Parameters
        timestamp: The Unix timestamp to be converted.

    Returns
        The ISO 8601 formatted date-time string, with a 'Z' suffix indicating UTC time.
    """
    return datetime.fromtimestamp(timestamp).isoformat() + 'Z'


def parse_time_input(time_str: str | None, time_format: str) -> int | None:
    """
    Parses a time input string into a Unix timestamp in seconds. Returns None if the input
    time string is None.

    Parameters
        time_str: The time string to parse, or None if no input time is given.
            Should be in local time.
        time_format: A format string as defined by the datetime module to interpret
            the given time string.

    Returns
        The time as a Unix timestamp in seconds, or None if the input string is None.
    """
    if not time_str:
        return None
    dt = datetime.strptime(time_str, time_format)
    return int(dt.timestamp())


def calculate_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calculate the great-circle distance in meters between two latitude/longitude points.

    Parameters:
        lat1: Latitude of the first point in degrees. Must not be None.
        lon1: Longitude of the first point in degrees. Must not be None.
        lat2: Latitude of the second point in degrees. Must not be None.
        lon2: Longitude of the second point in degrees. Must not be None.

    Returns:
        The great-circle distance between the two points in meters.
    """
    # Convert to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371000  # Earth radius in meters

    return c * r


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(description='Generate GPX tracks from location database')
    parser.add_argument('--config', default='config.yaml', help='Configuration file path')
    parser.add_argument('--start', help='Start time (YYYY-MM-DD HH:MM:SS)')
    parser.add_argument('--end', help='End time (YYYY-MM-DD HH:MM:SS)')

    args = parser.parse_args()

    # Get the configuration dictionary from a YAML file
    with open(args.config, 'r') as f:
        config_dict = yaml.safe_load(f)

    # Form the filename based on the configuration:
    filename = config_dict['gpx']['filename'].format_map(config_dict['gpx'])
    # Make sure there are no spaces in the filename
    filename = filename.replace(' ', '_')
    # Make sure it ends with suffix .gpx
    if not filename.endswith('.gpx'):
        filename += '.gpx'

    # Calculate the starting and ending times
    start_time_str = args.start or config_dict['default_time_range']['start_time']
    end_time_str = args.end or config_dict['default_time_range']['end_time']
    start_timestamp = parse_time_input(start_time_str, config_dict['time_format'])
    end_timestamp = parse_time_input(end_time_str, config_dict['time_format'])

    # It's time to process the data! Start by fetching the position data...
    print(f"Fetching data from {start_time_str or 'beginning'} to "
          f"{end_time_str or 'end'} local time...")
    raw_data = fetch_data(config_dict['database']['path'],
                          start_timestamp, end_timestamp)

    # ...then filter it...
    filtered_data = filter_data(raw_data,
                                config_dict['gpx']['filters']['min_distance_between_points'])

    # ... then create the GPX XML elements...
    gpx_elements = create_gpx(config_dict, filtered_data, start_time_str, end_time_str)

    # ... then finally save it to a file.
    filepath = save_gpx(gpx_elements, config_dict['gpx']['output_directory'], filename)
    print(f"GPX track saved to: {filepath}")


if __name__ == "__main__":
    main()
