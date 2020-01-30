#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:24:15 2020

@author: wader
"""

import whaletracks.detect_calls as detect
import numpy as np
import unittest


class TestFunctions(unittest.TestCase):
    
    def testDefaultScaleFunction(self):
        Sxx = np.array( [ [1, 1], [2,2]] )
        [vmin, vmax] = detect.defaultScaleFunction(Sxx)
        #import pdb; pdb.set_trace()
        self.assertLess(vmin, vmax)
        
    def testPlotwav(self):
        samp=50
        data=np.repeat(.5,10000)
        
        
        [f1, t1, Sxx1] = detect.plotwav(samp, data)
        
        #import pdb; pdb.set_trace()
        self.assertEqual(len(f1), Sxx1.shape[0])
        self.assertEqual(len(t1), Sxx1.shape[1])
        
        [f, t, Sxx] = detect.plotwav(samp, data, plotflag=False)
        
        self.assertEqual(len(f), Sxx.shape[0])
        self.assertEqual(len(t), Sxx.shape[1])
        
    def testBuildkernel(self):
        f0 = 16 
        f1 = 14.6 
        bdwdth = 1 
        dur = 10
        samp=50
        f=np.array(range(1, 25))
        t=np.array(range(1,6000))
        
        
        [tvec, fvec, BlueKernel] = detect.buildkernel(f0, f1, bdwdth, dur,
        f, t, samp, plotflag=True)
        
        #import pdb; pdb.set_trace()
        self.assertEqual(len(tvec),BlueKernel.shape[1])
        self.assertEqual(len(fvec),BlueKernel.shape[0])
        
        
    def testXcorr(self):
        samp=50
        data=np.repeat(.5,10001)
        f0 = 16 
        f1 = 14.6 
        bdwdth = 1 
        dur = 10
        [f, t, Sxx] = detect.plotwav(samp, data, plotflag=False)
        [tvec, fvec, BlueKernel] = detect.buildkernel(f0, f1, bdwdth, dur,
        f, t, samp, plotflag=False)
        
        
        [t_scale, CorrVal_scale] = detect.xcorr(t,f,Sxx,tvec,fvec,BlueKernel,
        plotflag=True)
        
        #import pdb; pdb.set_trace()
        self.assertEqual(len(t_scale),len(CorrVal_scale))

        
        
        
if __name__ == "__main__":
    unittest.main()