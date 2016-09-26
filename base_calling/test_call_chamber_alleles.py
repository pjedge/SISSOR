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
            self.assertEqual(cur,lst[i])
            cur = cca.safenext(it)
            i += 1


class test_format_chrom(unittest.TestCase):

    def test1(self):
        chrnum = '1'
        self.assertEqual(cca.format_chrom(chrnum),'chr1')

    def test1rev(self):
        chrnum = 'chr1'
        self.assertEqual(cca.format_chrom(chrnum),'chr1')

    def testX(self):
        chrnum = 'X'
        self.assertEqual(cca.format_chrom(chrnum),'chrX')

    def testXrev(self):
        chrnum = 'chrX'
        self.assertEqual(cca.format_chrom(chrnum),'chrX')


class test_prod(unittest.TestCase):

    def test1(self):
        lst = [4,5]
        self.assertEqual(cca.prod(lst),20)

    def test2(self):
        lst = [1,2,3,4,5]
        self.assertEqual(cca.prod(lst),120)

    def test3(self):
        lst = [0.5,0.5,0.2,0.2,0.1]
        self.assertAlmostEqual(cca.prod(lst),0.001)

    def test4(self):
        chrnum = 'chrX'
        self.assertEqual(cca.format_chrom(chrnum),'chrX')

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

class test_remove_substrings(unittest.TestCase):

    def test1(self):
        cur_string = 'she sells sea shells by the sea shore'
        replace_list = ['se','ore','by']
        expected_result = 'she lls a shells  the a sh'
        res = cca.remove_multiple_strings(cur_string, replace_list)
        self.assertEqual(expected_result, res)

    def test2(self):
        cur_string = 'she sells sea shells by the sea shore'
        replace_list = []
        expected_result = 'she sells sea shells by the sea shore'
        res = cca.remove_multiple_strings(cur_string, replace_list)
        self.assertEqual(expected_result, res)

    def test3(self):
        cur_string = 'she sells sea shells by the sea shore'
        replace_list = ['se']
        expected_result = 'she lls a shells by the a shore'
        res = cca.remove_multiple_strings(cur_string, replace_list)
        self.assertEqual(expected_result, res)

    def test4(self):
        cur_string = 'blahblahblahblahblather'
        replace_list = ['blather','blah']
        expected_result = ''
        res = cca.remove_multiple_strings(cur_string, replace_list)
        self.assertEqual(expected_result, res)

class test_hom_het_config_probs(unittest.TestCase):

    # there is the potential for more sophisticated tests but this should suffice
    # for the time being
    def test_hom_probs_sum_to_1(self):
        total = -1e100
        for v in cca.hom_config_probs.values():
            total = cca.addlogs(total, v)
        self.assertAlmostEqual(10**total,1.0)

    def test_het_probs_sum_to_1(self):
        total = -1e100
        for v in cca.het_config_probs.values():
            total = cca.addlogs(total, v)
        self.assertAlmostEqual(10**total,1.0)


class test_singlecell_config(unittest.TestCase):

    def test0(self):
        cfg = cca.singlecell_config(0,[24])
        self.assertTrue(len(cfg) == 1)
        cfg = cca.singlecell_config(1,[24])
        self.assertTrue(len(cfg) == 1)


class test_multicell_config(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)


class test_compute_mixed_allele_priors(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_pr_one_chamber_data(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)


class test_precompute_pr_one_chamber_data(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)


class test_pr_all_chamber_data(unittest.TestCase):

    def test1(self):
        self.assertTrue(True)

class test_pr_allele(unittest.TestCase):

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
