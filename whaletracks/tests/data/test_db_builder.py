import whaletracks.common.constants as cn
from whaletracks.data.db_builder import DBBuilder
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import sys
import unittest

IGNORE_TEST = False
TEST_DB_PTH = os.path.join(cn.TEST_DIR, "test_db_builder.db")


class TestDBBuilder (unittest.TestCase):

  def _cleanUp(self):
    if os.path.isfile(TEST_DB_PTH):
      os.remove(TEST_DB_PTH)

  def setUp(self):
    self._cleanUp()
    self.builder = DBBuilder(db_pth=TEST_DB_PTH)

  def tearDown(self):
    self._cleanUp()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.builder._db_pth, TEST_DB_PTH)

  def testGetStationDF(self):
    if IGNORE_TEST:
      return
    df = self.builder._getStationDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_STATION.columns))

  def testBuild(self):
    if IGNORE_TEST:
      return
    self.assertFalse(os.path.isfile(TEST_DB_PTH))
    self.builder.build()
    self.assertTrue(os.path.isfile(TEST_DB_PTH))



if __name__ == '__main__':
  unittest.main()
