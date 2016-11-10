#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for the compare_mappings.py module

Created on Tue Jun 17 16:44:51 2014

@author: Peter Edge
"""

import unittest
import sissorhands 
from math import log10

class test_addlogs(unittest.TestCase):

    def f(self, p1, p2):
        p3 = log10(p1 + p2)
        res = sissorhands.addlogs(log10(p1),log10(p2))
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

class test_hom_het_config_probs(unittest.TestCase):

    # there is the potential for more sophisticated tests but this should suffice
    # for the time being
    def test_hom_probs_sum_to_1(self):
        total = -1e100
        for v in sissorhands.hom_config_probs.values():
            total = sissorhands.addlogs(total, v)
        self.assertAlmostEqual(10**total,1.0)

    def test_het_probs_sum_to_1(self):
        total = -1e100
        for v in sissorhands.het_config_probs.values():
            total = sissorhands.addlogs(total, v)
        self.assertAlmostEqual(10**total,1.0)


class test_singlecell_config(unittest.TestCase):

    def test0(self):
        cfg = sissorhands.singlecell_config(0,[24])
        self.assertTrue(len(cfg) == 1)
        cfg = sissorhands.singlecell_config(1,[24])
        self.assertTrue(len(cfg) == 1)

F = 0.00019952623149688788 # value corresponing to phred F value

class test_parse_mpileup_base_qual(unittest.TestCase):

    def test0(self):
        raw_bd = '..,+2AA,,'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)

    def test1(self):
        raw_bd = '..,,,+2AA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test2(self):
        raw_bd = '..,,,+25AAAAAAAAAAAAAAAAAAAAAAAAA'
        raw_qd = 'FFFFF'
        ref_base = 'G'
        bd_exp = list('GGGGG')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test3(self):
        raw_bd = '+25AAAAAAAAAAAAAAAAAAAAAAAAA..,,,'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        

    def test4(self):
        raw_bd = '..,-2AA,,'
        raw_qd = 'FFFFF'
        ref_base = 'A'
        bd_exp = list('AAAAA')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)

    def test5(self):
        raw_bd = '..,,,-2AA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test6(self):
        raw_bd = '..,,,-25AAAAAAAAAAAAAAAAAAAAAAAAA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test7(self):
        raw_bd = '-25AAAAAAAAAAAAAAAAAAAAAAAAA..,,,'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test8(self):
        raw_bd = 'AAaaA'
        raw_qd = '!!!!!'
        ref_base = 'C'
        bd_exp = list('AAAAA')
        qd_exp = [1.0]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test9(self):
        raw_bd = 'A.,aA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('ACCAA')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test10(self):
        raw_bd = 'A.,aA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('ACCAA')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
        
    def test11(self):
        raw_bd = 'A.,$$^.AA'
        raw_qd = 'FFFFF'
        ref_base = 'C'
        bd_exp = list('ACCAA')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
     
    def test12(self):
        raw_bd = '..$$$$^C,,,-25AAAAAAAAAAAAAAAAAAAAAAAAA'
        raw_qd = '!!!!!'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [1.0]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)
 
    def test13(self):
        raw_bd = '.*..<>..'
        raw_qd = 'F!FF!!FF'
        ref_base = 'C'
        bd_exp = list('CCCCC')
        qd_exp = [F]*5
        bd, qd = sissorhands.parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)
        self.assertEqual(bd, bd_exp)
        self.assertEqual(qd, qd_exp)                       
    
###########################################################################
# TODO
###########################################################################

class test_multicell_config(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_compute_mixed_allele_priors(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_pr_one_chamber_data(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_precompute_pr_one_chamber(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_pr_genotype(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)
        
if __name__ == '__main__':

    # this try/catch block allows test suite to be run in spyder IDE interpreter
    # without hanging forever afterwards
    #via: http://stackoverflow.com/questions/9202772/tests-succeed-still-get-traceback
    try:
        unittest.main()
    except SystemExit as inst:
        pass
