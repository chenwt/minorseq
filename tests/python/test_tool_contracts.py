
import tempfile
import unittest
import re
import os.path

import pbcommand.testkit

TEST_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TEST_DIR, "data")

SIV_DATA = "/pbi/dept/secondary/siv/testdata/minorseq-test"
skip_if_no_testdata = unittest.skipUnless(os.path.exists(SIV_DATA),
                                          "Test data not available")


@skip_if_no_testdata
class TestFuse(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "fuse"
    REQUIRES_PBCORE = False
    INPUT_FILES = [
        "/pbi/dept/secondary/siv/testdata/minorseq-test/mix_hxb2.consensusalignmentset.xml"
    ]


@skip_if_no_testdata
class TestCleric(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "cleric"
    REQUIRES_PBCORE = False
    INPUT_FILES = [
        "/pbi/dept/secondary/siv/testdata/minorseq-test/mix_hxb2.consensusalignmentset.xml",
        "/pbi/dept/secondary/siv/testdata/minorseq-test/consensus.referenceset.xml",
        "/pbi/dept/secondary/siv/testdata/minorseq-test/hxb2.referenceset.xml"
    ]


@skip_if_no_testdata
class TestJuliet(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "juliet"
    REQUIRES_PBCORE = False
    INPUT_FILES = [
        "/pbi/dept/secondary/siv/testdata/minorseq-test/cleric.consensusalignmentset.xml"
    ]
