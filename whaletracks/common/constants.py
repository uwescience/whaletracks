"""Constants used in the whale project."""

import collections
import os
import sys

Schema = collections.namedtuple("Schema", "tablename columns")

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
CHANNEL = "channel"
CODE = "code"
CREATION_DATE = "creation_date"
DURATION = "duration"
ELEVATION  = 'elevation'
END_DATE  = 'end_date'
END_TIME = 'end_time'
LATITUDE  = 'latitude'
LONGITUDE  = 'longitude'
NETWORK = "network"
PEAKTIME = "peaktime"
START_TIME = "StartTime"
START_DATE  = 'start_date'
STATION_CODE = "station_code"
TERMINATION_DATE  = 'termination_date'
TOTAL_NUMBER_OF_CHANNELS  = 'total_number_of_channels'
# Table columns
SCM_STATION = Schema(
    tablename="stations",
    columns= [CREATION_DATE, ELEVATION, END_DATE, LATITUDE,
    LONGITUDE, START_DATE, CODE, TERMINATION_DATE,
    TOTAL_NUMBER_OF_CHANNELS]
    )
CHANNEL_SCM = Schema(
    tablename="channels",
    columns= [CODE]
    )
DETECTION_SCM = Schema(
    tablename="detections",
    columns= [CODE]
    )
SCMS = [SCM_STATION]
TABLES = [s.tablename for s in SCMS]
