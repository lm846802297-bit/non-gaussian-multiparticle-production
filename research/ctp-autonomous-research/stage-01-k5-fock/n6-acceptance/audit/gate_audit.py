"""K5/N6 complex-valued CTP gate: raw inputs, algebra, and numerical refinements."""
import argparse
import hashlib
from itertools import product
import json
import math
from pathlib import Path

PIN="eab1f5b3198d03eafc3f162299c9038956beaa0b"
EXPECTED={"K":5,"Nmax":6,"L":20.0,"m":1.0,"T":1.0,"lambda0":0.1,
          "t_init":-12.0,"events":[[1.0,0.0],[0.5,1.0],[-0.5,2.0],[-1.0,3.0]]}
HS=(0.00625,0.003125,0.0015625)
P="F4c"
S="connected_composite_source"
M="M33"
FD="SD_residual"
AN="analytic_residual"
M31=("01","02","03")
F2=("23","13","12")
LABELS=tuple("".join(s) for s in product("+-",repeat=4))
CONJ=str.maketrans("+-","-+")
MOMENTS=(M,)+tuple("M31."+k for k in M31)+tuple("F2."+k for k in F2)
MOMENT_AXES={M:(0,1,2,3),**{g+q:tuple(map(int,q)) for g,ks in (("M31.",M31),("F2.",F2)) for q in ks}}
ALG_ATOL=2e-13
EPS_GATE=1e-3
SOURCE_COEFF=EXPECTED["lambda0"]/math.cosh(EXPECTED["events"][0][0]/EXPECTED["T"])**2/6.

class Invalid(ValueError):pass

def require(ok,msg):
    if not ok:raise Invalid(msg)

def finite(x,key):
    require(type(x) in (float,int) and math.isfinite(x),f"{key}: expected finite number")
    return float(x)

def as_complex(v,where):
    require(type(v) in (tuple,list) and len(v)==2,f"{where}: expected [Re,Im]")
    return complex(finite(v[0],where+".real"),finite(v[1],where+".imag"))

def zpack(c):return [c.real,c.imag]

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()

def json_no_dup(kvs):
    result={}
    for k,v in kvs:
        require(k not in result,f"JSON has duplicate key: {k}")
        result[k]=v
    return result

def json_no_nonfinite(s):raise Invalid("JSON non-finite constant: "+s)

def json_read(path):
    with path.open() as handle:
        return json.load(handle,object_pairs_hook=json_no_dup,parse_constant=json_no_nonfinite)

def expected_binary_labels(obj,where):
    require(type(obj) is dict and set(obj)==set(LABELS),f"{where}: not exactly 16 branches")

def strict_identity(d,where,t_init_inferred=False):
    cfg=d.get("inputs",d.get("parameters",d.get("configuration",d.get("config",d))))
    require(type(cfg) is dict,where+": wrong configuration type")
    for k,v in EXPECTED.items():
        aliases={"t_init":("t_init","t0","initial_time","t_initial","tinit"),
                "K":("K","Kcut","kcut"),"Nmax":("Nmax","N","nmax"),
                "lambda0":("lambda0","lam","lam0"),"m":("m","mass"),
                "L":("L","length"),"T":("T","duration")}.get(k,(k,))
        values=[src[x] for src in (cfg,d) for x in aliases if x in src]
        if not values and k=="t_init" and t_init_inferred:continue
        require(bool(values),where+": missing "+k)
        if k=="events":
            require(all(type(a) in (tuple,list) and len(a)==len(v) for a in values),where+": events shape mismatch")
            for a in values:
                require(all(type(p) in (tuple,list) and len(p)==2 for p in a),where+": event pairs malformed")
                for p,e in zip(a,v):
                    require([finite(y,where+": event coordinate") for y in p]==e,where+": event lock mismatch")
        else:
            require(all(finite(y,where+": "+k)==v for y in values),where+": "+k+" lock mismatch")
    dim=d.get("basis_dimension",cfg.get("basis_dimension",None))
    require(finite(dim,where+": basis_dimension")==12376,where+": basis dimension mismatch")
    commit=d.get("source_commit",cfg.get("source_commit",d.get("phase7_git_base")))
    require(commit==PIN,where+": source commit mismatch")
    if not t_init_inferred:
        convention=d.get("branch_convention",cfg.get("branch_convention",
                 d.get("branch_label_convention",d.get("label_convention",d.get("branch_order")))))
        if isinstance(convention,str):
            require("+" in convention and "ket" in convention.lower(),where+": branch sign/ket unconfirmed")
    # The equation concerns latest insertion x0, so θ derivatives never hit other sites in all three stencils.
    require(EXPECTED["events"][0][0]-max(HS)>max(t for t,x in EXPECTED["events"][1:]),"invalid locked derivative geometry")

def binary_field(obj,k,where):
    return as_complex(obj.get(k),f"{where}.{k}")

def extract_centers(obj,where,moments=False,need_an=False):
    expected_binary_labels(obj,where)
    book={}
    for label in LABELS:
        r=obj[label]; tag=where+"."+label
        require(isinstance(r,dict),tag+": not a record")
        q={k:binary_field(r,k,tag) for k in (P,S,M)}
        if need_an:q[AN]=binary_field(r,AN,tag)
        if moments:
            for group,keys in (("M31",M31),("F2",F2)):
                require(isinstance(r.get(group),dict) and set(r[group])==set(keys),tag+" wrong "+group+" pair set")
                for key in keys:q[group+"."+key]=binary_field(r[group],key,tag+"."+group)
            products=[q["M31."+a]*q["F2."+b] for a,b in zip(M31,F2)]
            raw=r.get("M31_F2_terms")
            require(type(raw) is list and len(raw)==3,tag+": expected 3 M31_F2_terms")
            for i,v in enumerate(raw):
                require(abs(as_complex(v,tag+": subtraction term")-products[i])<=ALG_ATOL,tag+": product mismatch")
            require(abs(q[M]-sum(products)-q[S])<=ALG_ATOL,tag+": connected source recomposition mismatch")
            for i,a in enumerate(M31):
                alias="M31_j"+str(i+1)
                if alias in r:require(abs(binary_field(r,alias,tag)-q["M31."+a])<=ALG_ATOL,tag+": M31 alias mismatch")
            for b in F2:
                alias="F2_kl_"+b
                if alias in r:require(abs(binary_field(r,alias,tag)-q["F2."+b])<=ALG_ATOL,tag+": F2 alias mismatch")
        book[label]=q
    return book

def flip_at(label,i):return label[:i]+label[i].translate(CONJ)+label[i+1:]

def symmetries(table,fields,where,h=None):
    summary={}
    for q in fields:
        axes=MOMENT_AXES.get(q,(0,1,2,3))
        conj=lambda a:abs(table[a][q]-table[a.translate(CONJ)][q].conjugate())
        latest=lambda a:abs(table[a][q]-table[flip_at(a,min(axes))][q])
        e_conj=max(conj(a) for a in LABELS)
        e_latest=max(latest(a) for a in LABELS)
        e_ignored=max((abs(table[a][q]-table[flip_at(a,i)][q])
                         for i in range(4) if i not in axes for a in LABELS),default=0.)
        latest_tol=2e-15/h**2 if q==FD and h is not None else ALG_ATOL
        require(e_conj<=ALG_ATOL,f"{where}: {q} contour conjugation fails: {e_conj:.4g}")
        require(e_latest<=latest_tol,f"{where}: {q} latest-time invariance fails: {e_latest:.4g}")
        require(e_ignored<=ALG_ATOL,f"{where}: {q} spectator sign invariance fails: {e_ignored:.4g}")
        summary[q]={"conjugation_abs_max":e_conj,"latest_time_abs_max":e_latest,
                    "spectator_flip_abs_max":e_ignored,"latest_test_atol":latest_tol}
    return summary

def named_close(x,y,where):
    require(abs(x-y)<=ALG_ATOL,where+": invariant mismatch")

def optional_signed_normalization(raw,k,v,src,where):
    if k not in raw:return
    expected=abs(v)/(abs(v-src)+abs(src)+1e-30)
    require(abs(finite(raw[k],where+"."+k)-expected)<2e-11,where+" normalized residual inconsistent with signed value")

def get_rows(rows,central,where,fd_alias=FD):
    require(type(rows) is list and len(rows)==3,where+": missing three FD h rows")
    require(tuple(finite(r.get("h"),where+": h") for r in rows)==HS,where+": stencil h lock mismatch")
    run={}
    for row in rows:
        h=row["h"]; a=row.get("components"); tag=where+f".h={h}"
        expected_binary_labels(a,tag)
        run[h]={}
        for k in LABELS:
            w=tag+"."+k;b=a[k]
            rf=binary_field(b,fd_alias,w)
            qa=central[k]
            r={FD:rf}
            for q in (P,S,M,AN):
                if q in b:named_close(binary_field(b,q,w),qa[q],w+"."+q)
            residual_norm=SOURCE_COEFF*qa[S]
            optional_signed_normalization(b,"SD_relative_residual",rf,residual_norm,w)
            optional_signed_normalization(b,"analytic_relative_residual",qa[AN],residual_norm,w)
            run[h][k]=r
    return run

def change(ref,new,fields):
    out={}
    for q in fields:
        rows={}
        for a in LABELS:
            r,n=ref[a][q],new[a][q]
            da=abs(n-r)
            rows[a]={"reference":zpack(r),"candidate":zpack(n),"signed_delta":zpack(n-r),
                     "absolute_delta":da,"relative_to_reference":da/abs(r) if r else None}
        max_abs=max(rows,key=lambda k:rows[k]["absolute_delta"])
        finite_relative=[k for k,x in rows.items() if x["relative_to_reference"] is not None]
        max_rel=max(finite_relative,key=lambda k:rows[k]["relative_to_reference"]) if finite_relative else None
        out[q]={"min_reference_branch_abs":min(abs(ref[a][q]) for a in LABELS),
               "min_candidate_branch_abs":min(abs(new[a][q]) for a in LABELS),
               "max_absolute":rows[max_abs]["absolute_delta"],"limiting_branch_absolute":max_abs,
               "max_relative":rows[max_rel]["relative_to_reference"] if max_rel is not None else None,
               "limiting_branch_relative":max_rel,"n_nonzero_reference":len(finite_relative),
               "all_16":rows}
    return out

def stencil(central,rows):
    table={}
    for a in LABELS:
        an=central[a][AN]
        zz=[rows[h][a][FD] for h in HS]
        ee=[r-an for r in zz]
        c=[e/h**2 for e,h in zip(ee,HS)]
        power=[math.log2(abs(v)/abs(w)) if v and w else None for v,w in zip(ee,ee[1:])]
        ri=[(4*zz[i+1]-zz[i])/3 for i in (0,1)]
        phase=lambda e,f:(f*e.conjugate()).real/(abs(e)*abs(f)) if e and f else None
        table[a]={"signed_R_SD":[zpack(r) for r in zz],"signed_R_exact_derivative":zpack(an),
            "signed_E_stencil_only":[zpack(e) for e in ee],"abs_E_stencil":[abs(e) for e in ee],
            "observed_order_h_of_E":power,"phase_alignment_of_E":[phase(x,y) for x,y in zip(ee,ee[1:])],
            "signed_E_over_h2":[zpack(x) for x in c],"signed_Richardson":[zpack(r) for r in ri],
            "abs_Richardson_minus_exact_derivative":[abs(r-an) for r in ri]}
    return table

def epsilon(record,where):
    result={}
    for q in (P,S):
        row=record[q]
        require(row["n_nonzero_reference"]==16,where+": undefined relative error for "+q)
        value=row["max_relative"]
        result[q]={"observed":value,"threshold":EPS_GATE,"pass_strict_less":value<EPS_GATE}
    return result

def require_hash(recorded,f,where):
    require(type(recorded) is str and recorded==sha(f),where+": SHA-256 mismatch")

def abs_scalar_matches(r,key,v,w):
    require(abs(finite(r.get(key),w+"."+key)-v)<=ALG_ATOL,w+"."+key+": mismatch")

def audit_cf_structure(run,central,rows,where):
    # Independently reconstruct the equation and both signed differential insertions.
    result={}
    for tag in LABELS:
        s=run["central_observables"][tag]
        F=central[tag][P];S6=SOURCE_COEFF*central[tag][S]
        named_close(binary_field(s,"source_term",where+"."+tag),S6,where+"."+tag+".source_term")
        Dxx=binary_field(s,"dxx_F4c",where+"."+tag)
        Dtt=binary_field(s,"analytic_dtt_F4c",where+"."+tag)
        K=Dtt-Dxx+EXPECTED["m"]**2*F
        named_close(binary_field(s,"analytic_box_plus_mass_F4c",where+"."+tag),K,where+"."+tag+".analytic_box_plus_mass_F4c")
        named_close(central[tag][AN],K+S6,where+"."+tag+".analytic EOM")
        abs_scalar_matches(s,"analytic_scale",abs(K)+abs(S6)+1e-30,where+"."+tag)
        optional_signed_normalization(s,"analytic_relative_residual",central[tag][AN],S6,where+"."+tag)
    for row in run["rows"]:
        h=row["h"]
        step={}
        for tag in LABELS:
            s=row["components"][tag]; w=where+f".h={h}."+tag
            F=central[tag][P]
            mp,mm=(binary_field(s,key,w) for key in ("F4c_t_plus_h","F4c_t_minus_h"))
            dtt=((mp-F)-(F-mm))/h**2
            Dxx=binary_field(run["central_observables"][tag],"dxx_F4c",w)
            B=dtt-Dxx+EXPECTED["m"]**2*F
            S6=SOURCE_COEFF*central[tag][S]
            named_close(binary_field(s,"dtt_F4c",w),dtt,w+" reconstructed FD dtt")
            named_close(binary_field(s,"box_plus_mass_F4c",w),B,w+" reconstructed FD kinetic+mass")
            named_close(rows[h][tag][FD],B+S6,w+" reconstructed signed EOM")
            abs_scalar_matches(s,"SD_scale",abs(B)+abs(S6)+1e-30,w)
            abs_scalar_matches(s,"SD_residual_abs",abs(rows[h][tag][FD]),w)
            abs_scalar_matches(s,"fd_error_vs_analytic",abs(rows[h][tag][FD]-central[tag][AN]),w)
            step[tag]={"signed_rebuilt_FD_residual":zpack(B+S6),"delta_from_native":abs(B+S6-rows[h][tag][FD])}
        result[str(h)]=step
    return result

def extract_dop(d,p,runner):
    require_hash(d.get("runner_sha256"),runner,"DOP pinned segmentation runner")
    strict_identity(d,"DOP",t_init_inferred=True)
    for k,v in (("rtol",2e-10),("atol",2e-12),("max_step",.1)):
        require(finite(d.get(k),"DOP."+k)==v,"DOP integration settings mismatch")
    require(isinstance(d.get("rows"),list) and len(d["rows"])==3,"DOP no 3 stencil rows")
    c=extract_centers(d["rows"][0].get("components"),"DOP.central",need_an=True)
    return c,get_rows(d["rows"],c,"DOP"),{"sha256":sha(p),"runner_sha256":sha(runner),
          "t_init_and_contour_convention":"locked in provenance-bound archived/pinned DOP runner and operators; not independently JSON-asserted"}

def extract_direct(d,p,runner,pin,dop_path,dop_other=None):
    require_hash(d.get("runner_sha256"),runner,"independent signed moment runner")
    require_hash(d.get("operator_sha256"),pin,"independent signed moment operator")
    strict_identity(d,"direct")
    require(d.get("source_order")==[M,"M31:01","M31:02","M31:03","F2:23","F2:13","F2:12"],"direct source seven-factor order lock mismatch")
    require(d.get("independently_evaluated_contour_labels")==list(LABELS),"direct signed label order mismatch")
    require("independent" in d.get("solver","").lower(),"direct calculation route is not identified")
    require(d.get("operator_git_blob")=="29139f6f8a35d3d874c401d240318ed2c8cadc40","direct SHA-1 source pin not matched")
    for k,v in (("rtol",2e-10),("atol",2e-12),("max_step",.1)):
        require(finite(d.get(k),"direct."+k)==v,"direct integration settings mismatch")
    prov={"sha256":sha(p),"runner_sha256":sha(runner),"operator_sha256":sha(pin),
          "contour_convention":"verified via the provenance-bound direct insertion-path implementation: + ket"}
    if "archived_comparison" in d:
        alt=dop_other if dop_other is not None else dop_path
        require_hash(d["archived_comparison"].get("archived_sha256"),alt,"direct archived reference")
        prov["native_archived_comparison_sha256"]=sha(alt)
    c=extract_centers(d.get("components"),"direct.central",moments=True)
    for tag,r in d["components"].items():
        w="direct.central."+tag
        extra=r.get("F2_0j_for_F4")
        require(isinstance(extra,dict) and set(extra)==set(M31),w+" missing extra 2-point functions for F4c")
        F4raw=binary_field(r,"F4raw",w)
        terms=[binary_field(extra,a,w)*c[tag]["F2."+b] for a,b in zip(M31,F2)]
        named_close(c[tag][P],F4raw-sum(terms),w+" direct F4c recomposition")
    return c,prov

def parse_cf(d,p,runner,pin,fast,archive,dop_path,steps):
    strict_identity(d,"CF4")
    convention=d.get("branch_convention")
    require(isinstance(convention,str) and "+ chronological ket" in convention and "event order 0,1,2,3" in convention,"CF4 ordered branch lock mismatch")
    require_hash(d.get("phase7_original_archive_sha256"),archive,"CF4 restored Phase VII archive")
    require_hash(d.get("baseline_dop853_sha256"),dop_path,"CF4 DOP baseline native trace")
    hashes=d.get("native_source_sha256",{})
    for key,f in (("gate_cf4.py",runner),("sd_residual.py",pin),("phase7_fast_ops.py",fast)):
        require_hash(hashes.get(key),f,"CF4 "+key)
    src=d.get("runs")
    require(isinstance(src,list) and len(src)==len(steps),"CF4 incomplete integration dt series")
    require(tuple(finite(s.get("integrator_max_step"),"CF4 max step") for s in src)==steps,"CF4 integration dt lock mismatch")
    require(len(steps)>=2 and all(u==2*v for u,v in zip(steps,steps[1:])),"expected geometric CF4 refinement")
    data={}
    for dt,r in zip(steps,src):
        w=f"CF4.dt={dt}"
        require(r.get("engine",{}).get("order")==4,w+": integrator order metadata missing/invalid")
        c=extract_centers(r.get("central_observables"),w,moments=True,need_an=True)
        rows=get_rows(r.get("rows"),c,w)
        data[dt]={"central":c,"rows":rows,"reconstructed":audit_cf_structure(r,c,rows,w),
                  "vacuum_norm_t1":finite(r.get("vacuum_norm_t1"),w+".vacuum_norm"),"engine":r["engine"]}
    return data,{"sha256":sha(p),"native_source_sha256":hashes,"phase7_restored_archive":sha(archive)}

def signed_report(central,rows,w):
    sem={str(h):symmetries(rows[h],(FD,),w,h) for h in HS}
    return {"h_grid":HS,"signed_stencil_convergence":stencil(central,rows),
            "central_symmetries":symmetries(central,tuple(central[LABELS[0]]),w),"stencil_symmetries":sem,
            "max_relative_analytic":max(abs(x[AN])/(abs(x[AN]-SOURCE_COEFF*x[S])+abs(SOURCE_COEFF*x[S])+1e-30)
                                         for x in central.values())}

def separation(a,b):
    c=tuple(a["central"][LABELS[0]])
    report={"central_at_fixed_configuration":change(a["central"],b["central"],c),
            "fixed_h_integration_refinement":{}}
    for h in HS:
        u={z:{FD:a["rows"][h][z][FD],
              "signed_E_stencil_only":a["rows"][h][z][FD]-a["central"][z][AN]} for z in LABELS}
        v={z:{FD:b["rows"][h][z][FD],
              "signed_E_stencil_only":b["rows"][h][z][FD]-b["central"][z][AN]} for z in LABELS}
        report["fixed_h_integration_refinement"][str(h)]=change(u,v,(FD,"signed_E_stencil_only"))
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cf","direct","dop","cf-runner","direct-runner","dop-runner","pinned-module","restored-fast-ops","restored-archive","out"):
        p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--direct-dop-reference",type=Path)
    p.add_argument("--expected-dt",nargs="+",type=float,default=[.1,.05])
    p.add_argument("--fock-context",type=Path)
    a=p.parse_args()
    require(a.pinned_module.is_file(),"pinned SD module unavailable")
    pd=a.pinned_module.read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(pd)).encode()+b"\0"+pd).hexdigest()
    require(actual=="29139f6f8a35d3d874c401d240318ed2c8cadc40","archived SD module Git blob mismatch")
    dop,dr,dpro=extract_dop(json_read(a.dop),a.dop,a.dop_runner)
    dc,mpro=extract_direct(json_read(a.direct),a.direct,a.direct_runner,a.pinned_module,a.dop,a.direct_dop_reference)
    c,cfpro=parse_cf(json_read(a.cf),a.cf,a.cf_runner,a.pinned_module,a.restored_fast_ops,a.restored_archive,a.dop,tuple(a.expected_dt))
    steps=tuple(c);fine=c[steps[-1]]["central"]
    cm=change(dc,fine,(P,S)+MOMENTS)
    cd=change(dop,fine,(P,S,M,AN))
    md=change(dop,dc,(P,S,M))
    refinement={f"{u:g}_to_{v:g}":separation(c[u],c[v]) for u,v in zip(steps,steps[1:])}
    gates={"CF4_fine_vs_DOP":epsilon(cd,"CF4 fine vs DOP"),
           "CF4_fine_vs_independent_direct":epsilon(cm,"CF4 fine vs independent direct"),
           "independent_direct_vs_DOP":epsilon(md,"independent direct vs DOP"),
           "CF4_last_refinement":epsilon(list(refinement.values())[-1]["central_at_fixed_configuration"],"CF4 final refinement")}
    allg=all(v["pass_strict_less"] for r in gates.values() for v in r.values())
    report={"configuration_locked":EXPECTED,"basis_dimension":12376,"core_source_commit":PIN,
            "auditor_sha256":sha(Path(__file__)),"input_provenance":{"DOP":dpro,"CF4":cfpro,"direct":mpro},
            "signed_branch_and_stencil":{"DOP":signed_report(dop,dr,"DOP"),
                "CF4":{str(d):{
                    **signed_report(x["central"],x["rows"],f"CF4.dt={d}"),
                    "vacuum_norm_t1":x["vacuum_norm_t1"],"engine":x["engine"],
                    "native_EOM_reconstruction":x["reconstructed"]} for d,x in c.items()},
                "direct_central_symmetries":symmetries(dc,(P,S)+MOMENTS,"direct")},
            "integration_refinement_separate_from_FD":refinement,
            "comparison_CF4_fine_vs_DOP":cd,
            "comparison_CF4_fine_vs_DOP_same_FD_h":{
                str(h):change(dr[h],c[steps[-1]]["rows"][h],(FD,)) for h in HS},
            "comparison_CF4_fine_vs_direct_before_and_after_source_subtraction":cm,
            "comparison_direct_vs_DOP":md,
            "fixed_cutoff_numeric_gates":gates,"fixed_cutoff_numeric_gates_all_pass":allg,
            "no_physics_pass_threshold":"Conjugation and algebra tolerances guard arithmetic consistency; no analytic or finite-stencil SD residual physics threshold, observed order threshold, or Fock threshold is imposed."}
    if a.fock_context:
        f=json_read(a.fock_context)
        require("physical Fock cutoff sensitivity" in f.get("scope",""),"Fock sensitivity object mislabeled")
        require(finite(f.get("observables",{}).get(P,{}).get("relative_to_N6"),"archived Fock F4c")>0,
                 "no N5→N6 Fock displacement supplied")
        report["separate_physical_N5_to_N6_displacement"]={"sha256":sha(a.fock_context),"context":f,
                "excluded_from_fixed_Nmax_acceptance":True}
    a.out.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"Signed identity, provenance, and reconciliation complete; numeric gates pass={allg}.")
    for r,v in gates.items():
        print(r,*[(m,z["observed"]) for m,z in v.items()])
    return 0 if allg else 2

if __name__=="__main__":
    raise SystemExit(main())
