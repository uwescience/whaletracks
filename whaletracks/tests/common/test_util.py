import whaletracks.common.constants as cn
from whaletracks.common import util
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import sys
import unittest

IGNORE_TEST = False


class TestFunctions (unittest.TestCase):

  def testComplexifyString(self):
    stg = "(-0.03691, 0.03712);(-0.03691, -0.03712)"
    result = util.complexifyString(stg)
    self.assertEqual(len(result), 2)
    self.assertTrue(all([isinstance(c, complex) for c in result]))
    self.assertEqual(result[0], result[1].conjugate())
    #
    result = util.complexifyString("")
    self.assertEqual(len(result), 0)

   

if __name__ == '__main__':
  unittest.main()
