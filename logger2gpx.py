# logger2gpx.py

"""Extract GPX track data from a sqlite file containing position
data. The sqlite file should contain a table named `location_data` with
columns `timestamp`, `latitude`, `longitude`."""

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

    def timestamp_to_iso(self, timestamp):
        """Convert Unix timestamp to ISO format for GPX."""
        return datetime.fromtimestamp(timestamp).isoformat() + 'Z'

    def parse_time_input(self, time_str):
        """Parse time string. The string should be in local time. Returns unix epoch time"""
        if time_str:
            dt = datetime.strptime(time_str, self.config['time_format'])
            return int(dt.timestamp())
        return None

    def calculate_distance(self, lat1, lon1, lat2, lon2):
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

    def subsample_data(self, data):
        """Subsample data according to a specified interval."""
        if not self.config['gpx']['subsampling']['enabled']:
            return data

        if not data:
            return data

        interval = self.config['gpx']['subsampling']['interval_seconds']
        method = self.config['gpx']['subsampling']['method']

        print(f"Subsampling data with {interval}s interval using '{method}' method...")

        subsampled = []

        if method == 'first':
            subsampled = self._subsample_first(data, interval)
        elif method == 'last':
            subsampled = self._subsample_last(data, interval)
        elif method == 'closest_to_interval':
            subsampled = self._subsample_closest(data, interval)
        elif method == 'average':
            subsampled = self._subsample_average(data, interval)
        else:
            raise ValueError(f"Unknown subsampling method: {method}")

        print(f"Subsampled from {len(data)} to {len(subsampled)} points")
        return subsampled

    def _subsample_first(self, data, interval):
        """Take the first point in each interval."""
        if not data:
            return data

        subsampled = []
        start_time = data[0][0]  # First timestamp
        current_interval_start = start_time

        for row in data:
            timestamp = row[0]

            # If we've moved to a new interval, take this point
            if timestamp >= current_interval_start:
                subsampled.append(row)
                # Move to next interval
                current_interval_start = timestamp + interval

        return subsampled

    def _subsample_last(self, data, interval):
        """Take the last point in each interval."""
        if not data:
            return data

        subsampled = []
        start_time = data[0][0]

        # Group data by intervals
        intervals = {}
        for row in data:
            timestamp = row[0]
            interval_index = (timestamp - start_time) // interval
            intervals[interval_index] = row  # This will keep the last one for each interval

        # Sort by interval index and return
        for interval_index in sorted(intervals.keys()):
            subsampled.append(intervals[interval_index])

        return subsampled

    def _subsample_closest(self, data, interval):
        """Take the point closest to each interval boundary."""
        if not data:
            return data

        subsampled = []
        start_time = data[0][0]

        # Calculate target times for each interval
        current_target = start_time
        end_time = data[-1][0]

        i = 0
        while current_target <= end_time and i < len(data):
            # Find the point closest to current_target
            closest_point = None
            min_diff = float('inf')

            # Look at points around the target time
            while i < len(data):
                timestamp = data[i][0]
                diff = abs(timestamp - current_target)

                if diff < min_diff:
                    min_diff = diff
                    closest_point = data[i]
                    i += 1
                elif diff > min_diff:
                    # We've passed the closest point, break
                    break
                else:
                    i += 1

            if closest_point:
                subsampled.append(closest_point)

            current_target += interval

        return subsampled

    def _subsample_average(self, data, interval):
        """Average all points within each interval."""
        if not data:
            return data

        subsampled = []
        start_time = data[0][0]

        # Group data by intervals
        intervals = {}
        for row in data:
            timestamp, lat, lon = row
            interval_index = (timestamp - start_time) // interval

            if interval_index not in intervals:
                intervals[interval_index] = []
            intervals[interval_index].append(row)

        # Average each interval
        for interval_index in sorted(intervals.keys()):
            points = intervals[interval_index]

            # Calculate averages (skip None values)
            timestamps = [p[0] for p in points]
            latitudes = [p[1] for p in points if p[1] is not None]
            longitudes = [p[2] for p in points if p[2] is not None]

            avg_timestamp = sum(timestamps) / len(timestamps)
            avg_lat = sum(latitudes) / len(latitudes) if latitudes else None
            avg_lon = sum(longitudes) / len(longitudes) if longitudes else None

            subsampled.append((int(avg_timestamp), avg_lat, avg_lon))

        return subsampled

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
        # Create root GPX element
        gpx = Element('gpx')
        gpx.set('version', '1.1')
        gpx.set('creator', self.config['gpx']['author'])
        gpx.set('xmlns', 'http://www.topografix.com/GPX/1/1')

        # Add metadata
        metadata = SubElement(gpx, 'metadata')
        name = SubElement(metadata, 'name')
        name.text = self.config['gpx']['name']

        desc = SubElement(metadata, 'desc')
        desc.text = self.config['gpx']['description']
        if start_time_str and end_time_str:
            desc.text += f" (From {start_time_str} to {end_time_str})"

        # Add subsampling info to description
        if self.config['gpx']['subsampling']['enabled']:
            interval = self.config['gpx']['subsampling']['interval_seconds']
            method = self.config['gpx']['subsampling']['method']
            desc.text += f" - Subsampled every {interval}s using {method} method"

        # Create track
        trk = SubElement(gpx, 'trk')
        trk_name = SubElement(trk, 'name')
        trk_name.text = f"Western Flyer track"

        # Create track segment
        trkseg = SubElement(trk, 'trkseg')

        # Add track points
        for row in track_data:
            timestamp, lat, lon = row

            if lat is not None and lon is not None:
                trkpt = SubElement(trkseg, 'trkpt')
                trkpt.set('lat', str(lat))
                trkpt.set('lon', str(lon))

                # Add timestamp
                time_elem = SubElement(trkpt, 'time')
                time_elem.text = self.timestamp_to_iso(timestamp)

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

        # Apply subsampling
        subsampled_data = self.subsample_data(raw_data)

        # Filter data
        filtered_data = self.filter_track_points(subsampled_data)
        print(f"Filtered to {len(filtered_data)} track points.")

        if not filtered_data:
            print("No valid track points after filtering.")
            return None

        # Create GPX
        gpx_element = self.create_gpx(filtered_data, start_time_str, end_time_str)

        # Generate filename and save
        filename = self.config['gpx']['filename']
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
    parser.add_argument('--no-subsample', action='store_true', help='Disable subsampling for this run')
    parser.add_argument('--interval', type=int, help='Override subsampling interval (seconds)')

    args = parser.parse_args()

    try:
        generator = GPXGenerator(args.config)

        # Override subsampling settings from command line
        if args.no_subsample:
            generator.config['gpx']['subsampling']['enabled'] = False
        if args.interval:
            generator.config['gpx']['subsampling']['interval_seconds'] = args.interval

        generator.generate_gpx_track(args.start, args.end)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
