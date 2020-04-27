"""Constants used in the whale project."""

import collections
import os
import sys

# tablename in whaletracks.db
# filename in data
Schema = collections.namedtuple("Schema",
    "tablename columns csv_path")

class SchemaContainer(object):
    def __init__(self):
        self.schemas = []
        
    def append(self, tablename, columns, csv_path=None):
        schema = Schema(tablename=tablename, columns=columns, csv_path=csv_path)
        self.schemas.append(schema)
        
SCHEMA = SchemaContainer()
        
        

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
# Columns ended in _EPOCH populated automatically in db_builder
AZIMUTH = "azimuth"
CHANNEL = "channel"
CHANNEL_CODE = "channel_code"
CHANNEL_ID = "channel_id"  # NETWORK_CODE.STATION_CODE.CHANNEL_CODE
CHANNEL_TYPES = "channel_types"
CODE = "code"
CREATION_EPOCH = "creation_epoch" # Float for days since 1-1-1970
CREATION_TIME = "creation_time"  # UTC datetime
DIP = "dip"
DURATION = "duration"
ELEVATION  = 'elevation'
END_EPOCH  = 'end_epoch'  # Float for days since 1-1-1970
END_TIME  = 'end_time'  # UTC datetime
EVENT = "event"
GAIN = "gain"
LATITUDE  = 'latitude'
LONGITUDE  = 'longitude'
MIN_SIGNAL = "min_signal"
NETWORK = "network"
NETWORK_CODE = "network_code"  # str
PEAK_SIGNAL = "peak_signal"
PEAK_EPOCH = "peak_epoch"  # Float for days since 1-1-1970
PEAK_TIME = "peak_time"  # UTC datetime
POLES = "poles"  # semicolon separated values of immaginary numbers
SENSITIVITY_FREQUENCY = "sensitivity_frequency"
SENSITIVITY_VALUE = "sensitivity_value"
SENSOR = "sensor"
SNR = "snr"
START_EPOCH = "start_epoch" # Float for days since 1-1-1970
START_TIME = "start_time"  #UTC datetime
STATION_CODE = "station_code"  # str
STATION_ID = "station_id"  # NETWORK_CODE.STATION_CODE
THRESHOLD = "threshold"
TERMINATION_EPOCH = 'termination_epoch' # Float for days since 1-1-1970
TERMINATION_TIME  = 'termination_time'  # UTC datetime
TOTAL_NUMBER_OF_CHANNELS  = 'total_number_of_channels'
VALUE = "value"
ZEROES = "zeroes" # semicolon separated values of immaginary numbers

# Table schemas
SCM_STATION = Schema(tablename="stations",
    columns=[CREATION_EPOCH, CREATION_TIME, ELEVATION, END_EPOCH,
    END_TIME, LATITUDE, LONGITUDE, NETWORK_CODE, STATION_CODE, 
    START_EPOCH, START_TIME, TERMINATION_EPOCH, TERMINATION_TIME,
    TOTAL_NUMBER_OF_CHANNELS],csv_path=None)
SCHEMA.append(SCM_STATION.tablename, SCM_STATION.columns)

SCM_CHANNEL = Schema(tablename="channels",
    columns=[AZIMUTH, CHANNEL_ID, CHANNEL_TYPES, DIP,
    END_EPOCH, END_TIME, POLES, SENSITIVITY_FREQUENCY, 
    SENSITIVITY_VALUE, SENSOR, START_EPOCH, START_TIME, STATION_ID, ZEROES],
    csv_path=None)
SCHEMA.append(SCM_CHANNEL.tablename, SCM_CHANNEL.columns)

SCM_DETECTION = Schema(tablename="detections",
    columns=[NETWORK_CODE, STATION_CODE, DURATION,
    END_EPOCH, END_TIME, MIN_SIGNAL, PEAK_SIGNAL,
    PEAK_EPOCH, PEAK_TIME,
    START_EPOCH, START_TIME, STATION_CODE, THRESHOLD, SNR],
    csv_path=os.path.join(DATA_DIR, 'detections.csv'))
SCHEMA.append(SCM_DETECTION.tablename, SCM_DETECTION.columns, csv_path=SCM_DETECTION.csv_path)

SCHEMA.append("peaks",
    [NETWORK_CODE, STATION_CODE,
    START_EPOCH, START_TIME, END_EPOCH, END_TIME, VALUE, EVENT])
SCHEMA.append("station_quality",
    [NETWORK_CODE, STATION_CODE, START_EPOCH, START_TIME, END_EPOCH, END_TIME])

SCMS = SCHEMA.schemas
TABLES = [s.tablename for s in SCMS]
