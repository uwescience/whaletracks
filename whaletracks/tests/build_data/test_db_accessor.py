import whaletracks.common.constants as cn
from whaletracks.build_data.db_accessor import DBAccessor
from whaletracks.build_data.db_builder import DBBuilder
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import sys
import unittest

IGNORE_TEST = False
TEST_DB_PTH = os.path.join(cn.TEST_DIR, "test_db_accessor.db")


class TestDBAccessor (unittest.TestCase):

  def _cleanUp(self):
    if os.path.isfile(TEST_DB_PTH):
      os.remove(TEST_DB_PTH)

  def _init(self):
    builder = DBBuilder(db_pth=TEST_DB_PTH)
    builder.build()
    self.accessor = DBAccessor(db_pth=TEST_DB_PTH)

  def setUp(self):
    self._init()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(os.path.isfile(TEST_DB_PTH))
    self.assertEqual(len(self.accessor._tables), 0)

  def testDfStation(self):
    if IGNORE_TEST:
      return
    df = self.accessor.df_station
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_STATION.columns))

  def testDfChannel(self):
    if IGNORE_TEST:
      return
    df = self.accessor.df_channel
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_CHANNEL.columns))

  def testDfDetection(self):
    if IGNORE_TEST:
      return
    df = self.accessor.df_detection
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_DETECTION.columns))
   

if __name__ == '__main__':
  unittest.main()
