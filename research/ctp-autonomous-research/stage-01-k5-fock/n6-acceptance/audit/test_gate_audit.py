"""Failure-path unit tests using locked native traces, no physics computation."""
import copy
import tempfile
import unittest
from pathlib import Path
from gate_audit import (ALG_ATOL, EPS_GATE, EXPECTED, F2, FD, HS, Invalid, LABELS, M, MOMENTS, P, S,
    audit_cf_structure, change, epsilon, extract_centers, extract_direct, extract_dop, get_rows, json_read,
    parse_cf, strict_identity, symmetries)

HERE=Path(__file__).resolve().parent
DATA=HERE/"input"

class NativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cf=json_read(DATA/"k5-n6-cf4.json")
        cls.direct=json_read(DATA/"k5-n6-direct-moments.json")
        cls.dop=json_read(DATA/"k5-n6.json")
        cls.case=cls.cf["runs"][-1]

    def test_no_nan_duplicate(self):
        for invalid in ['{"a":NaN}','{"a":Infinity}','{"a":0,"a":1}']:
            with self.subTest(value=invalid), tempfile.TemporaryDirectory() as d:
                p=Path(d)/"broken.json";p.write_text(invalid)
                with self.assertRaises(Invalid):json_read(p)

    def test_complete_direct_and_dop(self):
        p=DATA
        a,_=extract_direct(self.direct,p/"k5-n6-direct-moments.json",p/"direct_moments.py",p/"sd_residual.py",
                           p/"k5-n6.json",p/"k5-n6-segmented.json")
        ca,h,_=extract_dop(self.dop,p/"k5-n6.json",p/"scan.py")
        self.assertEqual(len(a),16);self.assertEqual(tuple(h),HS)
        self.assertTrue(all(e["pass_strict_less"] for e in epsilon(change(ca,a,(P,S)),"native").values()))

    def test_config_and_float_failure(self):
        d=copy.deepcopy(self.cf)
        d["K"]=6
        with self.assertRaisesRegex(Invalid,"K lock mismatch"):strict_identity(d,"CF4")
        d=copy.deepcopy(self.cf)
        d["events"][0][0]=float("inf")
        with self.assertRaisesRegex(Invalid,"finite number"):strict_identity(d,"CF4")
        d=copy.deepcopy(self.cf)
        d["t_init"]=-11.999999999
        with self.assertRaisesRegex(Invalid,"t_init lock mismatch"):strict_identity(d,"CF4")

    def test_provenance_mismatch(self):
        p=DATA
        with self.assertRaisesRegex(Invalid,"SHA-256 mismatch"):
            extract_dop(self.dop,p/"k5-n6.json",p/"direct_moments.py")

    def test_stencil_h_lock(self):
        c=extract_centers(self.case["central_observables"],"CF4",moments=True,need_an=True)
        rr=copy.deepcopy(self.case["rows"]);rr[-1]["h"]*=2
        with self.assertRaisesRegex(Invalid,"stencil h lock mismatch"):get_rows(rr,c,"CF4")

    def test_nonfinite_or_bool_component(self):
        for k in (float('nan'),float('inf'),False):
            with self.subTest(field=k):
                s=copy.deepcopy(self.direct["components"])
                s[LABEl][P][0]=k
                with self.assertRaisesRegex(Invalid,"finite number"):extract_centers(s,"direct",moments=True)

    def test_raw_moment_missing_or_inconsistent(self):
        s=copy.deepcopy(self.direct["components"])
        s[LABEl]["M31"].pop("01")
        with self.assertRaisesRegex(Invalid,"pair set"):extract_centers(s,"direct",moments=True)
        s=copy.deepcopy(self.direct["components"])
        s[LABEl]["M31_F2_terms"][0][1]+=1e-7
        with self.assertRaisesRegex(Invalid,"product mismatch"):extract_centers(s,"direct",moments=True)

    def test_symmetry_signed_flip_corruption(self):
        c=extract_centers(self.direct["components"],"direct",moments=True)
        symmetries(c,(P,S)+MOMENTS,"direct")
        c[LABEl]["F2.23"]+=1e-9j
        with self.assertRaisesRegex(Invalid,"contour conjugation fails"):
            symmetries(c,("F2.23",),"direct")

    def test_stencil_raw_perturbation_detected(self):
        g=copy.deepcopy(self.case)
        c=extract_centers(g["central_observables"],"CF4",moments=True,need_an=True)
        r=get_rows(g["rows"],c,"CF4")
        g["rows"][0]["components"][LABEl]["F4c_t_plus_h"][0]+=1e-11
        with self.assertRaisesRegex(Invalid,"reconstructed FD dtt"):
            audit_cf_structure(g,c,r,"CF4")

    def test_threshold_is_strict_and_undefined_fails(self):
        info={q:{"n_nonzero_reference":16,"max_relative":EPS_GATE} for q in (P,S)}
        self.assertFalse(epsilon(info,"threshold")[P]["pass_strict_less"])
        info[P]["n_nonzero_reference"]=15
        with self.assertRaisesRegex(Invalid,"undefined"):epsilon(info,"threshold")

    def test_two_cf_dt_both_required(self):
        p=DATA
        d=copy.deepcopy(self.cf);d["runs"].pop()
        with self.assertRaisesRegex(Invalid,"incomplete integration dt series"):
            parse_cf(d,p/"k5-n6-cf4.json",p/"gate_cf4.py",p/"sd_residual.py",p/"phase7_fast_ops.py",
                    p/"phase7-cutoff-scan-001.tar.xz",p/"k5-n6.json",(.1,.05))

LABEl=LABELS[0]
if __name__=="__main__":
    unittest.main()
