"""Constants used in the whale project."""

import collections
import os
import sys

# tablename in whaletracks.db
# filename in data
Schema = collections.namedtuple("Schema",
    "tablename columns csv_path")

PROJECT_NAME = "whaletracks"

# PATHS
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
TEST_DIR = os.path.join(PROJECT_DIR, PROJECT_NAME)
TEST_DIR = os.path.join(TEST_DIR, "tests")
PROJECT_CODE = os.path.join(PROJECT_DIR, PROJECT_NAME)
COMMON_CODE = os.path.join(PROJECT_DIR, "common_python")
STATION_FILE = os.path.join(DATA_DIR, "ALL-StationList.csv")
BLUE_DETECTION_FILE = os.path.join(DATA_DIR,
    "ExampleBlueDetections.csv")
DB_PTH = os.path.join(DATA_DIR, "whaletracks.db")

# Add search paths
sys.path.insert(0, PROJECT_CODE)
sys.path.insert(0, COMMON_CODE)

# Columns
AZIMUTH = "azimuth"
CHANNEL = "channel"
CHANNEL_CODE = "channel_code"
CHANNEL_ID = "channel_id"  # NETWORK_CODE.STATION_CODE.CHANNEL_CODE
CHANNEL_TYPES = "channel_types"
CODE = "code"
CREATION_DATE = "creation_date"
DIP = "dip"
DURATION = "duration"
ELEVATION  = 'elevation'
END_DATE  = 'end_date'
END_TIME = 'end_time'
EVENT = "event"
GAIN = "gain"
LATITUDE  = 'latitude'
LONGITUDE  = 'longitude'
MIN_SIGNAL = "min_signal"
NETWORK = "network"
NETWORK_CODE = "network_code"
PEAK_SIGNAL = "peak_signal"
PEAK_TIME = "peak_time"
POLES = "poles"  # semicolon separated values of immaginary numbers
SENSITIVITY_FREQUENCY = "sensitivity_frequency"
SENSITIVITY_VALUE = "sensitivity_value"
SENSOR = "sensor"
START_TIME = "start_time"
START_DATE  = 'start_date'
STATION_CODE = "station_code"
STATION_ID = "station_id"  # NETWORK_CODE.STATION_CODE
THRESHOLD = "threshold"
TERMINATION_DATE  = 'termination_date'
TOTAL_NUMBER_OF_CHANNELS  = 'total_number_of_channels'
VALUE = "value"
ZEROES = "zeroes" # semicolon separated values of immaginary numbers

# Table schemas
SCM_STATION = Schema(
    tablename="stations",
    columns= [CREATION_DATE, ELEVATION, END_DATE, LATITUDE,
    LONGITUDE, NETWORK_CODE, STATION_CODE, START_DATE, TERMINATION_DATE,
    TOTAL_NUMBER_OF_CHANNELS],
    csv_path=None,
    )
SCM_CHANNEL = Schema(tablename="channels",
    columns = [
    AZIMUTH, CHANNEL_ID, CHANNEL_TYPES, DIP, END_DATE, 
    POLES, SENSITIVITY_FREQUENCY, 
    SENSITIVITY_VALUE, SENSOR, START_DATE, STATION_ID, ZEROES],
    csv_path=None,
    )
SCM_DETECTION = Schema(tablename="detections",
    columns=[
    NETWORK_CODE, STATION_CODE, DURATION, END_TIME,
    MIN_SIGNAL, PEAK_SIGNAL, PEAK_TIME,
    START_TIME, STATION_CODE, THRESHOLD],
    csv_path=os.path.join(DATA_DIR, "detections.csv"),
    )
SCM_PEAK = Schema(tablename="peaks",
    columns=[
    NETWORK_CODE, STATION_CODE, START_TIME, END_TIME, VALUE, EVENT],
    csv_path=None,
    )
SCM_STATION_QUALITY = Schema(tablename="station_quality",
    columns=[
    NETWORK_CODE, STATION_CODE, START_TIME, END_TIME],
    csv_path=None,
    )
SCMS = [SCM_STATION, SCM_CHANNEL, SCM_DETECTION]
TABLES = [s.tablename for s in SCMS]
