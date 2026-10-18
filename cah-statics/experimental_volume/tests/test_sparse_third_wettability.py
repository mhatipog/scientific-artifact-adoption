"""Offline, zero-network checks for sparse_third_wettability.py."""
import struct
import unittest
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sparse_third_wettability import (
    BYTES_PER_FRAME, EXPECTED_FRAMES, analyze, crossing,
    find_pore, parse_first_ifd, parse_nd2_info,
)


def synthetic_pore():
    first=np.full((512,512),3500,dtype=np.uint16)
    first[69,116:439]=500
    first[391,116:439]=500
    first[69:392,116]=500
    first[69:392,438]=500
    return first


class SparseHoldoutOffline(unittest.TestCase):
    def test_known_nikon_path_and_julian_conversion(self):
        info = parse_nd2_info(
            "Location = /media/Theta_95/Q1/drainage/test1035.nd2\n"
            "dTimeAbsolute = 2459529.1673781597\n"
            "dCalibration = 2.485922277608956")
        self.assertTrue(info["julian_date_interpretation_utc"].startswith(
            "2021-11-10T16:01:01"))
        with self.assertRaises(ValueError):
            parse_nd2_info("Location = /Theta_90/Q1/a.nd2\n"
                           "dTimeAbsolute = 2459529.1673781597")

    def test_relative_crossings_no_time_shift(self):
        t=[0.,1.,2.,3.]
        y=[0.,.6,.5,1.]
        self.assertAlmostEqual(crossing(t,y,.3),.5)
        self.assertAlmostEqual(crossing(t,y,.8),2.5)
        self.assertIsNone(crossing(t,y,1.1))

    def test_geometry_paired_walls(self):
        box=find_pore(synthetic_pore())
        for got,want in zip(box,(116,69,438,391)):
            self.assertLessEqual(abs(got-want),1)

    def test_header_proves_contiguity(self):
        desc=b"ImageJ=1.53f\nimages=4848\nfinterval=0.033328503370285034\n\x00"
        off=2_090_691
        header=bytearray(4096)
        header[:4]=b"MM\x00*"
        struct.pack_into(">I",header,4,8)
        struct.pack_into(">H",header,8,15)
        data=[
            (254,4,1,0),(256,4,1,512),(257,4,1,512),
            (258,3,1,16<<16),(262,3,1,1<<16),
            (270,2,len(desc),194),(273,4,1,off),
            (277,3,1,1<<16),(278,3,1,512<<16),
            (279,4,1,BYTES_PER_FRAME),
            (282,5,1,305),(283,5,1,313),
            (296,3,1,1<<16),
            (50838,4,4850,321),(50839,1,2070970,19721)]
        for k,row in enumerate(data):
            struct.pack_into(">HHII",header,10+12*k,*row)
        struct.pack_into(">I",header,190,off+EXPECTED_FRAMES*BYTES_PER_FRAME)
        header[194:194+len(desc)]=desc
        tags,start=parse_first_ifd(bytes(header))
        self.assertEqual(start,off)
        self.assertEqual(tags[50839][2],19721)
        struct.pack_into(">I",header,190,off+BYTES_PER_FRAME*3)
        with self.assertRaises(ValueError):
            parse_first_ifd(bytes(header))

    def test_synthetic_monotone_physical_front(self):
        first=synthetic_pore()
        images=np.repeat(first[None,:,:],48,axis=0)
        for j in range(48):
            cut=123+min(310,max(0,(j-10)*11))
            images[j,73:387,121:cut]=2800
        report,rows=analyze(first,images)
        summary=report["threshold_summary"]["180_0.5"]
        self.assertGreater(summary["duration_s"],0)
        self.assertLess(summary["duration_s"],1.2)
        self.assertGreater(summary["final_mask"],summary["start_mask"]+.65)
        self.assertEqual(len(rows),48)
        self.assertTrue(report["front_speed_mm_s"]>0)


if __name__ == "__main__":
    unittest.main()
