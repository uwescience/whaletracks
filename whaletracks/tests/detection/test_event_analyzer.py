#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:53:17 2020

@author: wader
"""

from whaletracks.detection.event_analyzer import EventAnalyzer
from whaletracks.detection import event_analyzer 
import whaletracks.common.constants as cn
import numpy as np
import unittest
from random import seed
from random import random
import matplotlib.pyplot as plt
from obspy import UTCDateTime

IS_PLOT = False

class TestEventAnalyzer(unittest.TestCase):
    
    def setUp(self):
        seed(1)
        times = range(100)
        values = np.random.randint(-100,100,len(times))*random()
        inds=values<0
        values[inds] = 0
        #import pdb; pdb.set_trace()
        self.analyzer = EventAnalyzer(times, values, UTCDateTime("2011-12-14T12:00:00.000"),
                                 dur=np.random.randint(2,4)*(times[1]-times[0]), 
                                 prominence=np.mean(values))
     
    def testConstructorSimple(self):
        times = range(10)
        values = [0, 0, 1, 2, 3, 3, 2, 1, 0, 0]
        analyzer = EventAnalyzer(times, values, UTCDateTime("2011-12-14T12:00:00.000"),
                                 dur=4, 
                                 prominence=1.5)
        diff = set(analyzer.df.columns).difference(cn.SCM_DETECTION.columns)
        self.assertEqual(len(diff), 0)
    
    def testConstructorSophisticated(self):
        self.assertGreater(len(self.analyzer.df), 1)
        #import pdb; pdb.set_trace()     
    
        
    def testPlot(self):
        self.analyzer.plot(is_plot=IS_PLOT)
        
        
        
        
if __name__ == "__main__":
    unittest.main()