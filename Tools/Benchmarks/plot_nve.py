#!/usr/bin/env python3
"""Plot measured NVE refinement without converting missing cases to results.

Optional plotting dependency: matplotlib==3.10.7. This is not needed by the
reference generator or native benchmark runner.
"""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--qualification",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();data=json.loads(a.qualification.read_text())
    assert data["schema"]=="numivivo.org/md-nve-refinement-result/v1"
    assert not a.out.exists(),"plot output exists"
    fig,ax=plt.subplots(figsize=(9,6.3),layout="constrained")
    labels={"alanine-smooth":"Alanine", "protein-smooth":"DHFR", "complex-vacuum":"Vacuum complex",
            "dna-smooth":"DNA", "membrane-smooth":"Membrane", "ions-smooth":"Salt water",
            "water-orthogonal-smooth":"Orthogonal water"}
    failures=[]
    for row in data["cases"]:
        if not row.get("measurements"):
            name="Original water input" if row["identifier"]=="water" else row["identifier"]
            failures.append(name+": "+row["outcome"]);continue
        values=sorted(row["measurements"],key=lambda m:m["timeStepPS"])
        ax.plot([m["timeStepPS"]*1000 for m in values],[m["rmsDeviationPerParticle"] for m in values],
            "o-",linewidth=1.8,markersize=5,label=labels.get(row["identifier"],row["identifier"])+(" (failed)" if row["outcome"]!="passed" else ""))
    ax.set_xscale("log",base=2);ax.set_yscale("log")
    ax.set_xticks([0.25,0.5,1],labels=["0.25","0.5","1.0"])
    ax.set_xlabel("Time step (femtoseconds)")
    ax.set_ylabel("RMS energy deviation (kJ/mol per particle)")
    ax.set_title("Smaller time steps reduce energy error",loc="left",fontsize=17,pad=20)
    ax.grid(True,which="major",alpha=0.2);ax.legend(frameon=False,ncol=2,fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    fig.supxlabel("0.1 ps constrained NVE · smooth periodic LJ; unchanged vacuum model\n"
                  "Short conservation study; equilibrium and speed not assessed.\n"+"; ".join(failures),fontsize=9)
    fig.savefig(a.out,metadata={"Creator":"NumiVivo benchmark plot; matplotlib "+matplotlib.__version__,"Date":None} if a.out.suffix==".svg" else None,dpi=160)
    if a.out.suffix==".svg":
        a.out.write_text("\n".join(line.rstrip() for line in a.out.read_text().splitlines())+"\n")
    print(json.dumps(dict(qualificationSHA256=hashlib.sha256(a.qualification.read_bytes()).hexdigest(),
        plotSHA256=hashlib.sha256(a.out.read_bytes()).hexdigest(),matplotlib=matplotlib.__version__,retainedFailures=failures)))


if __name__=="__main__":main()
