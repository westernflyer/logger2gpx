# logger2gpx.py

"""Extract GPX track data from a sqlite file containing position
data. The sqlite file should contain a table named `location_data` with
columns `timestamp`, `latitude`, `longitude`.

The program offers the opportunity to filter the data based on a minimum distance
between points.
"""

import sqlite3
import yaml
import argparse
import os
from datetime import datetime
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom
import math


class GPXGenerator:
    def __init__(self, config_path="config.yaml"):
        """Initialize the GPX generator with configuration."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        self.db_path = self.config['database']['path']

    def connect_db(self):
        """Connect to the SQLite database."""
        return sqlite3.connect(self.db_path)

    @staticmethod
    def timestamp_to_iso(timestamp):
        """Convert Unix timestamp to ISO format for GPX."""
        return datetime.fromtimestamp(timestamp).isoformat() + 'Z'

    def parse_time_input(self, time_str):
        """Parse time string. The string should be in local time. Returns unix epoch time"""
        if time_str:
            dt = datetime.strptime(time_str, self.config['time_format'])
            return int(dt.timestamp())
        return None

    @staticmethod
    def calculate_distance(lat1, lon1, lat2, lon2):
        """Calculate distance between two points using Haversine formula."""
        if None in [lat1, lon1, lat2, lon2]:
            return float('inf')

        # Convert to radians
        lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

        # Haversine formula
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
        c = 2 * math.asin(math.sqrt(a))
        r = 6371000  # Earth radius in meters

        return c * r

    def fetch_track_data(self, start_time=None, end_time=None):
        """
        Fetches track data from the database within an optional time range, ordered by timestamp.

        This method queries the database table `location_data` to retrieve data including
        timestamp, latitude, and longitude. Optional parameters `start_time` and
        `end_time` can be supplied to filter the results within a specific time range. If no
        time range is provided, all available records in the table are retrieved. The query
        is executed using a database connection established via `connect_db`.

        :param start_time: Optional; the start of the time range for filtering records.
            Only records with a timestamp greater than this value will be included.
        :param end_time: Optional; the end of the time range for filtering records.
            Only records with a timestamp equal to or less than this value will be included.
        :return: A list of tuples representing the fetched data, where each tuple consists
            of a timestamp, latitude, longitude, and heading.
        :rtype: list[tuple]
        """
        with self.connect_db() as conn:
            cursor = conn.cursor()

            # Build query with optional time constraints
            query = "SELECT timestamp, latitude, longitude FROM location_data"
            params = []

            conditions = []
            if start_time:
                conditions.append("timestamp > ?")
                params.append(start_time)
            if end_time:
                conditions.append("timestamp <= ?")
                params.append(end_time)

            if conditions:
                query += " WHERE " + " AND ".join(conditions)

            query += " ORDER BY timestamp"

            cursor.execute(query, params)
            data = cursor.fetchall()

        return data


    def filter_track_points(self, data):
        """Apply filters to track points."""
        filtered_data = []
        last_lat, last_lon = None, None
        min_distance = self.config['gpx']['filters']['min_distance_between_points']
        skip_null = self.config['gpx']['filters']['skip_null_coordinates']

        for row in data:
            timestamp, lat, lon = row

            # Skip points with NULL coordinates if configured
            if skip_null and (lat is None or lon is None):
                continue

            # Apply minimum distance filter
            if min_distance > 0 and last_lat is not None and last_lon is not None:
                distance = self.calculate_distance(last_lat, last_lon, lat, lon)
                if distance < min_distance:
                    continue

            filtered_data.append(row)
            last_lat, last_lon = lat, lon

        return filtered_data

    def create_gpx(self, track_data, start_time_str=None, end_time_str=None):
        """Create GPX XML from track data."""

        # Get interpolation dictionary:
        interp_dict = self.config['gpx']

        # Create root GPX element
        gpx = Element('gpx')
        gpx.set('version', '1.1')
        gpx.set('creator', self.config['gpx']['author'])
        gpx.set('xmlns', 'http://www.topografix.com/GPX/1/1')

        # Add metadata
        metadata = SubElement(gpx, 'metadata')
        name = SubElement(metadata, 'name')
        name.text = self.config['gpx']['gpx_name'].format_map(interp_dict)

        desc = SubElement(metadata, 'desc')
        desc.text = self.config['gpx']['description'].format_map(interp_dict)
        if start_time_str and end_time_str:
            desc.text += f" (From {start_time_str} to {end_time_str})"

        # Create track
        trk = SubElement(gpx, 'trk')
        trk_name = SubElement(trk, 'name')
        trk_name.text = self.config['gpx']['track_name'].format_map(interp_dict)

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
                time_elem.text = GPXGenerator.timestamp_to_iso(timestamp)

        return gpx

    def save_gpx(self, gpx_element, filename):
        """Save GPX to file with pretty formatting."""
        # Create output directory if it doesn't exist
        output_dir = self.config['gpx']['output_directory']
        os.makedirs(output_dir, exist_ok=True)

        # Pretty print XML
        rough_string = tostring(gpx_element, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ", encoding='utf-8')

        # Save to file
        filepath = os.path.join(output_dir, filename)
        with open(filepath, 'wb') as f:
            f.write(pretty_xml)

        return filepath


    def generate_gpx_track(self, start_time_str=None, end_time_str=None):
        """Main method to generate GPX track."""
        # Parse time inputs
        start_timestamp = self.parse_time_input(start_time_str)
        end_timestamp = self.parse_time_input(end_time_str)

        # Use default times from config if not provided
        if not start_time_str and 'start_time' in self.config.get('default_time_range', {}):
            start_time_str = self.config['default_time_range']['start_time']
            start_timestamp = self.parse_time_input(start_time_str)

        if not end_time_str and 'end_time' in self.config.get('default_time_range', {}):
            end_time_str = self.config['default_time_range']['end_time']
            end_timestamp = self.parse_time_input(end_time_str)

        # Fetch data
        print(f"Fetching data from {start_time_str or 'beginning'} to {end_time_str or 'end'} local time...")
        raw_data = self.fetch_track_data(start_timestamp, end_timestamp)

        if not raw_data:
            print("No data found for the specified time range.")
            return None

        print(f"Found {len(raw_data)} raw data points.")

        # Filter data
        filtered_data = self.filter_track_points(raw_data)
        print(f"Filtered to {len(filtered_data)} track points.")

        if not filtered_data:
            print("No valid track points after filtering.")
            return None

        # Create GPX
        gpx_element = self.create_gpx(filtered_data, start_time_str, end_time_str)

        # Generate filename and save
        filename = self.config['gpx']['filename'].format_map(self.config['gpx'])
        # Make sure there are no spaces in the filename
        filename = filename.replace(' ', '_')
        if not filename.endswith('.gpx'):
            filename += '.gpx'
        filepath = self.save_gpx(gpx_element, filename)

        print(f"GPX track saved to: {filepath}")
        return filepath


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(description='Generate GPX tracks from location database')
    parser.add_argument('--config', default='config.yaml', help='Configuration file path')
    parser.add_argument('--start', help='Start time (YYYY-MM-DD HH:MM:SS)')
    parser.add_argument('--end', help='End time (YYYY-MM-DD HH:MM:SS)')

    args = parser.parse_args()

    generator = GPXGenerator(args.config)

    generator.generate_gpx_track(args.start, args.end)


if __name__ == "__main__":
    main()
