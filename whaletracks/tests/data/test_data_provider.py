import whaletracks.common.constants as cn
from whaletracks.data.data_provider import DataProvider
from common_python.testing import helpers

import numpy as np
import pandas as pd
import sys
import unittest

IGNORE_TEST = False

class TestDataProvider (unittest.TestCase):

  def setUp(self):
    self.provider = DataProvider()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(helpers.isValidDataFrame(
        self.provider.df_station,
        expected_columns=cn.STATION_SCM.columns))

  def testUpdateDB(self):
    if IGNORE_TEST:
      return
    #self.provider.updateDB()



if __name__ == '__main__':
  unittest.main()
