import whaletracks.common.constants as cn
from whaletracks.common import util
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import sys
import unittest
from obspy import UTCDateTime

IGNORE_TEST = False


class TestFunctions (unittest.TestCase):

  def testComplexifyString(self):
    stg = "(-0.03691, 0.03712);(-0.03691, -0.03712)"
    result = util.complexifyString(stg)
    self.assertEqual(len(result), 2)
    self.assertTrue(all([isinstance(c, complex) for c in result]))
    self.assertEqual(result[0], result[1].conjugate())
    
    result = util.complexifyString("")
    self.assertEqual(len(result), 0)

  def testDatestrToEpoch(self):
     
      detlist = ['2011-12-24T00:15:28.452000Z',
                 '2011-12-24T00:15:56.612000Z',
                 '2011-12-24T00:22:44.164000Z',
                 '2011-12-24T00:24:57.028000Z',
                 '2011-12-24T00:27:11.428000Z']
      
      result = util.datestrToEpoch(detlist)
      self.assertEqual(len(detlist),len(result))
      sort_result=np.sort(result)
       
      for j in range(1,len(result)+1):
          res=j-1
          
          self.assertEqual(sort_result[res],result[res])
          
          
  def testDatetimeToEpoch(self):
     
      detlist = ['2011-12-24T00:15:28.452000Z',
                 '2011-12-24T00:15:56.612000Z',
                 '2011-12-24T00:22:44.164000Z',
                 '2011-12-24T00:24:57.028000Z',
                 '2011-12-24T00:27:11.428000Z']
      
      datetime_list=UTCDateTime.
      result = util.datestrToEpoch(detlist)
      self.assertEqual(len(detlist),len(result))
      sort_result=np.sort(result)
       
      for j in range(1,len(result)+1):
          res=j-1
          
          self.assertEqual(sort_result[res],result[res])        
          

  def testAddEpochColumns(self):
       
      test_dct = {'START_TIME','END_TIME','THRESH','STATION'}
      
      expected_result = ['START_TIME','END_TIME','THRESH',
                         'STATION','START_EPOCH','END_EPOCH']
      test_df = pd.DataFrame([],columns=test_dct)
      
      epoch_df = util.addEpochColumns(test_df)
      #import pdb; pdb.set_trace()
      self.assertEqual(len(epoch_df.columns),len(expected_result))
      
      for j in range(1,len(expected_result)+1):
          res=j-1
          self.assertTrue(epoch_df.columns[res],expected_result[res])
   

if __name__ == '__main__':
  unittest.main()
