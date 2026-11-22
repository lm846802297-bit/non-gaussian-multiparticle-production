"""Independent signed audit of baseline and high K data. Optional fallback input."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

LABELS=tuple("".join(i) for i in itertools.product("+-",repeat=4))
TR=str.maketrans("+-","-+")
COLS=("F4c","M33","connected_composite_source","SD_residual","analytic_residual")
NATIVE_CASES={
 "K5N6L20":(5,6,20., "k5-n6-segmented.json"),
 "K5N7L20":(5,7,20., "k5-n7.json"),
 "K6N6L20":(6,6,20., "k6-n6-A.json"),
 "K6N7L20":(6,7,20., "k6-n7-A.json"),
 "K6N6L24":(6,6,24., "k6-n6-A-L24.json"),
 "K6N7L24":(6,7,24., "k6-n7-A-L24.json"),
}
def z(x,k):
    a,b=x[k]
    z=complex(a,b)
    assert math.isfinite(abs(z)),(k,x)
    return z
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def audit(path,meta):
    d=json.loads(path.read_text())
    K,N,L,_=meta
    assert (d["K"],d["Nmax"],d["L"])==(K,N,L)
    assert d["source_commit"]=="eab1f5b3198d03eafc3f162299c9038956beaa0b"
    assert d["events"]==[[1.,0.],[.5,1.],[-.5,2.],[-1.,3.]]
    assert d["basis_dimension"]==math.comb(2*K+N+1,N)
    assert d["m"]==d["T"]==1. and d["lambda0"]==.1
    assert d["hermiticity_error"]<1.e-14
    e=[]
    first=d["rows"][0]["components"]
    assert [r["h"] for r in d["rows"]]==[.00625,.003125,.0015625]
    for row in d["rows"]:
        h=row["h"]; c=row["components"]
        assert set(c)==set(LABELS)
        meas={}
        for col in COLS:
            abs_max=max(abs(z(c[a],col)) for a in LABELS)
            ce=max(abs(z(c[a],col)-z(c[a.translate(TR)],col).conjugate()) for a in LABELS)
            lt=max(abs(z(c[a],col)-z(c[a[0].translate(TR)+a[1:]],col)) for a in LABELS)
            assert ce<1.e-15
            assert lt<(2.e-15/h**2 if col=="SD_residual" else 3.e-14)
            meas[col]={"max_abs":abs_max,"sign_flip_conjugation_abs_error":ce,"latest_time_abs_error":lt}
            if col!="SD_residual":
                assert max(abs(z(c[a],col)-z(first[a],col)) for a in LABELS)<1.e-15
        for a in LABELS:
            v,s,an=z(c[a],"SD_residual"),z(c[a],"connected_composite_source"),z(c[a],"analytic_residual")
            src=.1/(math.cosh(d["events"][0][0])**2)*s/6.
            exp=abs(v)/(abs(v-src)+abs(src)+1e-30)
            expa=abs(an)/(abs(an-src)+abs(src)+1e-30)
            assert abs(exp-c[a]["SD_relative_residual"])<2e-16
            assert abs(expa-c[a]["analytic_relative_residual"])<2e-16
            assert abs(abs(v-an)-c[a]["fd_error_vs_analytic"])<1.e-20
        x=max(c[a]["SD_relative_residual"] for a in LABELS)
        assert abs(x-row["max_relative"])<1.e-16
        e.append({"h":h,"signed_checks":meas,"max_relative":x,
                  "limiting_branch":max(LABELS,key=lambda a:c[a]["SD_relative_residual"])})
    series=[]
    for label in LABELS:
        v=[z(r["components"][label],"SD_residual") for r in d["rows"]]
        an=z(first[label],"analytic_residual")
        series.append({"branch":label,"FD_error":[abs(w-an) for w in v],
               "complex_Richardson_abs_error":[abs((4*v[i+1]-v[i])/3-an) for i in (0,1)]})
    summary={"file":path.name,"sha256":sha(path),"runner_sha256":d["runner_sha256"],
                "V_nnz":d["V_nnz"],"basis_dim":d["basis_dimension"],"rows":e,
                "analytic_max_normalized":max(v["analytic_relative_residual"] for v in first.values()),
                "signed_Richardson":series}
    return summary,first

def compare(x,y):
    out={}
    for c in COLS[:3]:
        val=[(a,z(x[a],c),z(y[a],c)) for a in LABELS]
        a,old,new=max(val,key=lambda t:abs((t[2]-t[1])/t[1]))
        b,oldu,newu=max(val,key=lambda t:abs(t[2]-t[1]))
        ft,gt=z(x["++++"],c),z(y["++++"],c)
        out[c]={"max_relative":abs(new-old)/abs(old),"limiting_branch_relative":a,
                "max_absolute":abs(newu-oldu),"limiting_branch_absolute":b,
                "min_old_magnitude":min(abs(t[1]) for t in val),
                "allforward_complex_delta":[(gt-ft).real,(gt-ft).imag],
                "allforward_relative":abs(gt-ft)/abs(ft)}
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",type=Path,default=Path("/tmp/ctpproject/run"))
    p.add_argument("--main",type=Path,default=Path("/tmp/ctp_next/run"))
    p.add_argument("--fallback",type=Path)
    p.add_argument("--out",type=Path,required=True)
    q=p.parse_args()
    record={}
    for k,v in NATIVE_CASES.items():
        f=(q.base if k.startswith("K5") else q.main)/v[-1]
        if not f.exists() and q.fallback is not None:f=q.fallback/v[-1]
        if f.exists():record[k]=audit(f,v)
    comparisons={}
    for name,pair in {
        "fixed_L20_K5_to_K6_N7":("K5N7L20","K6N7L20"),
        "matched_Lambda_volume_N7":("K5N7L20","K6N7L24"),
        "fixed_K6_L20_N6_to_N7":("K6N6L20","K6N7L20"),
        "fixed_K6_L24_N6_to_N7":("K6N6L24","K6N7L24"),
        "fixed_K6_N7_L20_to_L24":("K6N7L20","K6N7L24"),
        "fixed_L20_K5_to_K6_N6":("K5N6L20","K6N6L20"),
        "matched_Lambda_volume_N6":("K5N6L20","K6N6L24"),
    }.items():
        if set(pair)<=record.keys():comparisons[name]=compare(*(record[x][1] for x in pair))
    result={"completed":sorted(record),"not_found":sorted(NATIVE_CASES.keys()-record.keys()),
            "signed_audit":{k:v[0] for k,v in record.items()},"comparisons":comparisons}
    q.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("completed",*result["completed"],"; missing",*result["not_found"])
