import whaletracks.common.constants as cn
from whaletracks.data.db_access import DBAccessor
from common_python.testing import helpers

import numpy as np
import pandas as pd
import sys
import unittest

IGNORE_TEST = False

class TestDBAccessor (unittest.TestCase):

  def setUp(self):
    self.provider = DBAccessor()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(len(self._db_tables), 0)



if __name__ == '__main__':
  unittest.main()
