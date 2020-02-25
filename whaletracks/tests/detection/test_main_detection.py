#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:53:17 2020

@author: wader
"""

from whaletracks.detection import main_detection
import whaletracks.common.constants as cn

import numpy as np
import unittest
from obspy import UTCDateTime
COLUMNS =[cn.DURATION, cn.END_TIME, cn.MIN_SIGNAL, cn.PEAK_SIGNAL,
          cn.PEAK_TIME, cn.START_TIME, cn.THRESHOLD, cn.STATION_CODE,
          cn.NETWORK_CODE]

IS_PLOT = False

class TestFunctions(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def testMain1(self):
        return
        startime = ("2011-12-14T12:00:00.000")
        endtime = ("2011-12-14T12:20:00.000")
        station_codes = "J28A,J29A"
        df = main_detection.main(startime, endtime, station_ids=station_codes)
        diff = set(COLUMNS).symmetric_difference(df.columns)
        self.assertEqual(len(diff), 0)
        self.assertGreater(len(df), 0)
        
        
    def testMain2(self):
        startime = ("2011-12-14T12:00:00.000")
        endtime = ("2011-12-14T12:20:00.000")
        station_codes = "J28Z"
        df = main_detection.main(startime, endtime, station_ids=station_codes)
        self.assertEqual(len(df), 0)
        
        
        
if __name__ == "__main__":
    unittest.main()