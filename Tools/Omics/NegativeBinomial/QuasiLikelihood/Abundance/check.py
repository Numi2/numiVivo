#!/usr/bin/env python3
"""Controlled abundance checks against pinned edgeR and 80-digit score roots."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import mpmath as mp

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--binary', required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
sha = lambda data: hashlib.sha256(data).hexdigest()
rows, labels = [], []
patterns = [
    ('equal', [0,2,7,11], [math.log(1000)]*4),
    ('unequal', [0,3,7,0], [0,math.log(10),math.log(100),math.log(1000)]),
    ('zero', [0,0,0], [0,math.log(10),math.log(100)]),
    ('singleton', [31], [math.log(1000)]),
    ('large', [0,10**6,1,10**5], [math.log(1e3),math.log(1e6),math.log(1e4),math.log(1e7)]),
]
for name, counts, offsets in patterns:
    for phi in [0,1e-8,0.05,4,100]:
        for prior in [0,0.5,2,10]:
            if name == 'zero' and prior == 0:
                continue
            rows.append(dict(counts=counts, offsets=offsets, dispersion=phi, priorCount=prior))
            labels.append(dict(id=f'{name}-phi{phi}-prior{prior}', expected='finite', challenging=False))
for name, counts, offsets, phi in [
        ('large-exact-poisson', [2**53,2**52,0], [0,math.log(10),math.log(100)], 0),
        ('large-exact-nb', [2**53,2**52,0], [0,math.log(10),math.log(100)], .05),
        ('wide-libraries', [0,1000,0], [-300,0,300], .2),
        ('small-libraries', [0,3,7], [-600,-599,-598], .2),
        ('large-libraries', [0,3,7], [598,599,600], .2)]:
    rows.append(dict(counts=counts,offsets=offsets,dispersion=phi,priorCount=2))
    labels.append(dict(id=name,expected='finite',challenging=True))
for name, row in [
        ('zero-without-prior', dict(counts=[0,0],offsets=[0,0],dispersion=.2,priorCount=0)),
        ('inexact-count', dict(counts=[2**64-1],offsets=[0],dispersion=.2)),
        ('offset-outside-domain', dict(counts=[1],offsets=[701],dispersion=.2)),
        ('scaled-prior-underflow', dict(counts=[1,1],offsets=[-700,700],dispersion=.2)),
        ('work-bound', dict(counts=[0,3,7,0],offsets=[0,1,2,3],dispersion=.2,maximumScoreEvaluations=1))]:
    rows.append(row)
    labels.append(dict(id=name,expected='error',challenging=False))
payload = json.dumps(dict(rows=rows),separators=(',',':')).encode()
assert not (args.out/'native.json').exists(), 'Retain prior attempts before rerunning'
(args.out/'input.json').write_bytes(payload)
(args.out/'labels.json').write_text(json.dumps(labels,indent=2)+'\n')
native_run = subprocess.run(['ssh','macmini',args.binary],input=payload,capture_output=True)
(args.out/'native.log').write_bytes(native_run.stderr)
(args.out/'native.json').write_bytes(native_run.stdout)
assert native_run.returncode == 0
source = Path(__file__).parent
with (args.out/'reference.log').open('wb') as log:
    reference_run = subprocess.run(['/opt/homebrew/bin/Rscript',str(source/'reference.R'),
                                   str(args.out/'input.json'),str(args.out/'reference.json')],
        stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
assert reference_run.returncode == 0
native = json.loads(native_run.stdout)['rows']
reference = json.loads((args.out/'reference.json').read_text())['results']
mp.mp.dps = 80
checks = []
for label, row, result, ref in zip(labels,rows,native,reference,strict=True):
    check = dict(**label,native=result,reference=ref)
    if label['expected'] == 'error':
        check['status'] = 'passed' if 'error' in result else 'failed-expected-error'
    elif 'fit' not in result:
        check['status'] = 'failed-native'
    else:
        offsets = [mp.mpf(v) for v in row['offsets']]
        libraries = [mp.exp(o) for o in offsets]
        average = mp.fsum(libraries)/len(libraries)
        prior = mp.mpf(row.get('priorCount',2))
        z = [mp.mpf(y)+prior*l/average for y,l in zip(row['counts'],libraries)]
        adjusted = [l+2*prior*l/average for l in libraries]
        phi = mp.mpf(row['dispersion'])
        def score(beta):
            mu = [l*mp.exp(beta) for l in adjusted]
            return mp.fsum((y-m)/(1+phi*m) for y,m in zip(z,mu))
        start = mp.log(mp.fsum(z)/mp.fsum(adjusted))
        low, high, step = start, start, mp.mpf(1)
        while score(low) < 0:
            low = start-step;step *= 2
        step = mp.mpf(1)
        while score(high) > 0:
            high = start+step;step *= 2
        for _ in range(240):
            middle = (low+high)/2
            if score(middle) > 0: low = middle
            else: high = middle
        answer = ((low+high)/2+mp.log(10**6))/mp.log(2)
        error = abs(result['fit']['log2CountsPerMillion']-float(answer))
        reference_errors = {k:abs(ref[k]-float(answer)) if isinstance(ref.get(k),(int,float)) and math.isfinite(ref[k]) else None for k in ['value','tighterValue']}
        reference_pass = reference_errors['value'] is not None and reference_errors['value'] <= 2e-8
        check.update(status='passed' if error<=2e-8 and (reference_pass or label['challenging']) else 'failed',
                     independentLog2CPM=mp.nstr(answer,75),independentAbsoluteError=error,
                     referenceAbsoluteErrors=reference_errors,referenceDefaultPassed=reference_pass)
    checks.append(check)
summary = dict(status='passed' if all(c['status']=='passed' for c in checks) else 'failed',
               cases=len(checks),checks=checks,maximumIndependentAbsoluteError=max(c.get('independentAbsoluteError',0) for c in checks),
               referenceDisagreements=[c['id'] for c in checks if c.get('referenceDefaultPassed') is False],
               binarySHA256=subprocess.check_output(['ssh','macmini','shasum','-a','256',args.binary]).decode().split()[0],
               protocolSHA256=sha((source/'PROTOCOL.md').read_bytes()),inputSHA256=sha(payload),
               nativeSHA256=sha(native_run.stdout),referenceSHA256=sha((args.out/'reference.json').read_bytes()),
               referenceScriptSHA256=sha((source/'reference.R').read_bytes()))
(args.out/'checks.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='checks'},indent=2))
raise SystemExit(0 if summary['status']=='passed' else 1)
