#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:53:17 2020

@author: wader
"""

from whaletracks.detection import main_detection
import whaletracks.common.constants as cn
import os
import numpy as np
import unittest
from obspy import UTCDateTime

COLUMNS =[cn.DURATION, cn.END_TIME, cn.MIN_SIGNAL, cn.PEAK_SIGNAL,
          cn.PEAK_TIME, cn.START_TIME, cn.THRESHOLD, cn.STATION_CODE,
          cn.NETWORK_CODE]

IS_PLOT = False
IGNORE_TEST = True
TEST_DETECT_PTH=os.path.join(cn.TEST_DIR, "test_detections.csv")
TEST_CHUNK_PTH=os.path.join(cn.TEST_DIR, "test_chunks.csv")

class TestFunctions(unittest.TestCase):
    
    def _cleanUp(self):
        if os.path.isfile(TEST_DETECT_PTH):
            os.remove(TEST_DETECT_PTH)
        if os.path.isfile(TEST_CHUNK_PTH):
            os.remove(TEST_CHUNK_PTH)
    
    def setUp(self):
        self._cleanUp()
        
    def tearDown(self):
        self._cleanUp()
    
    def testMain1(self):
        self.assertFalse(os.path.isfile(TEST_DETECT_PTH))
        self.assertFalse(os.path.isfile(TEST_CHUNK_PTH))
        
        startime = ("2011-12-14T12:00:00.000")
        endtime = ("2011-12-14T12:20:00.000")
        station_ids = "J28A,J29A"
        #import pdb; pdb.set_trace()
        df = main_detection.main(startime, endtime,
            station_ids=station_ids,detection_pth=TEST_DETECT_PTH,
            chunk_pth=TEST_CHUNK_PTH)
        diff = set(COLUMNS).symmetric_difference(df.columns)
        self.assertEqual(len(diff), 0)
        self.assertGreater(len(df), 0)
        self.assertTrue(os.path.isfile(TEST_DETECT_PTH))
        self.assertTrue(os.path.isfile(TEST_CHUNK_PTH))
        
    def testMain2(self):
        self.assertFalse(os.path.isfile(TEST_DETECT_PTH))
        self.assertFalse(os.path.isfile(TEST_CHUNK_PTH))
        if IGNORE_TEST:
          return
        startime = ("2011-12-14T12:00:00.000")
        endtime = ("2011-12-14T12:20:00.000")
        station_ids = "J28Z"
        df = main_detection.main(startime, endtime,
            station_ids=station_ids,detection_pth=TEST_DETECT_PTH,
            chunk_pth=TEST_CHUNK_PTH)
        self.assertEqual(len(df), 0)
        self.assertFalse(os.path.isfile(TEST_DETECT_PTH))
        self.assertFalse(os.path.isfile(TEST_CHUNK_PTH))
        
        
if __name__ == "__main__":
    unittest.main()
