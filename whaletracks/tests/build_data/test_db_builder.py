import whaletracks.common.constants as cn
from whaletracks.build_data.db_builder import DBBuilder
from common_python.testing import helpers

import copy
import os
from obspy.core.utcdatetime import UTCDateTime
import pickle
import unittest

IGNORE_TEST = False
TEST_DB_PTH = os.path.join(cn.TEST_DIR, "test_db_builder.db")
TEST_DB_PCL = os.path.join(cn.TEST_DIR, "test_db_builder.pcl")
# Used save values of the network if they exist
if os.path.isfile(TEST_DB_PCL):
  BUILDER = pickle.load(open(TEST_DB_PCL, "rb"))
else:
  BUILDER = DBBuilder(db_pth=TEST_DB_PTH,
      start_time=UTCDateTime("2015-01-01"),
      end_time=UTCDateTime("2015-01-10"))
  pickle.dump(BUILDER, open(TEST_DB_PCL, "wb"))


class TestDBBuilder (unittest.TestCase):

  def _cleanUp(self):
    if os.path.isfile(TEST_DB_PTH):
      os.remove(TEST_DB_PTH)

  def setUp(self):
    self.builder = copy.deepcopy(BUILDER)
    self._cleanUp()

  def tearDown(self):
    self._cleanUp()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.builder._db_pth, TEST_DB_PTH)

  def testMakeStationDF(self):
    if IGNORE_TEST:
      return
    df = self.builder._makeStationDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_STATION.columns))

  def testBuild(self):
    if IGNORE_TEST:
      return
    self.assertFalse(os.path.isfile(TEST_DB_PTH))
    self.builder.build()
    self.assertTrue(os.path.isfile(TEST_DB_PTH))

  def testMakeStationDF(self):
    if IGNORE_TEST:
      return
    df = self.builder._makeChannelDF()
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=cn.SCM_CHANNEL.columns))



if __name__ == '__main__':
  unittest.main()
