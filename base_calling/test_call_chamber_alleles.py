#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for the compare_mappings.py module

Created on Tue Jun 17 16:44:51 2014

@author: Peter Edge
"""

import unittest
import call_chamber_alleles as cca
from math import log10

class test_safenext(unittest.TestCase):    

    def test1(self):
        lst = [1,2,3,4,5]
        it  = iter(lst)
        
        cur = cca.safenext(it)
        i = 0
        while(cur):
            assert(cur == lst[i])
            cur = cca.safenext(it)
            i += 1
            

class test_format_chrom(unittest.TestCase):    

    def test1(self):
        chrnum = '1'
        assert(cca.format_chrom(chrnum) == 'chr1')   

    def test1rev(self):
        chrnum = 'chr1'
        assert(cca.format_chrom(chrnum) == 'chr1')           

    def testX(self):
        chrnum = 'X'
        assert(cca.format_chrom(chrnum) == 'chrX')   

    def testXrev(self):
        chrnum = 'chrX'
        assert(cca.format_chrom(chrnum) == 'chrX') 


class test_prod(unittest.TestCase):    

    def test1(self):
        lst = [4,5]
        assert(cca.prod(lst) == 20)   

    def test2(self):
        lst = [1,2,3,4,5]
        assert(cca.prod(lst) == 120)            

    def test3(self):
        lst = [0.5,0.5,0.2,0.2,0.1]
        print(cca.prod(lst))
        self.assertAlmostEqual(cca.prod(lst),0.001)    

    def test4(self):
        chrnum = 'chrX'
        assert(cca.format_chrom(chrnum) == 'chrX') 
        
class test_addlogs(unittest.TestCase):    

    def f(self, p1, p2):
        p3 = log10(p1 + p2)
        res = cca.addlogs(log10(p1),log10(p2))
        self.assertAlmostEqual(p3, res)
    
    def test1(self):
        self.f(0.001,0.000001)
        
    def test2(self):
        self.f(5e-20,0.001)

    def test3(self):
        self.f(0.001,5e-20)
        
    def test4(self):
        self.f(10,20)
        
    def test5(self):
        self.f(0.01,0.232823832)
    
        
class test_subtractlogs(unittest.TestCase):    

    def f(self, p1, p2):
        p3 = log10(p1 - p2)
        res = cca.subtractlogs(log10(p1),log10(p2))
        self.assertAlmostEqual(p3, res)
    
    def test1(self):
        self.f(0.001,0.000001)
        
    def test2(self):
        self.f(0.001,1e-10)

    def test3(self):
        self.f(0.001,1e-20)
        
    def test4(self):
        self.f(20,10)
        
    def test5(self):
        self.f(1,0.5)
                
if __name__ == '__main__':

    # this try/catch block allows test suite to be run in spyder IDE interpreter
    # without hanging forever afterwards
    #via: http://stackoverflow.com/questions/9202772/tests-succeed-still-get-traceback
    try:
        unittest.main()
    except SystemExit as inst:
        pass
    