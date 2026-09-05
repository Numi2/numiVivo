"""Record the sources and actual completed reports; never manufacture pass logs."""
import datetime
import hashlib
import json
import platform
import sys
from pathlib import Path

out, rng_source, rng_fragment = map(Path, sys.argv[1:4])
mode = sys.argv[4]
files = sys.argv[5:] + [
    'Tools/Posterior/PortableSupport.swift', 'Tools/Posterior/PortableChecks.swift',
    'Tools/Posterior/DomainChecks.swift', 'Tools/Posterior/run_portable_checks.sh',
    'Tools/Posterior/write_manifest.py', 'Tools/Posterior/reference-data.json',
    'Tools/Posterior/domain-reference.json',
    'Examples/target-engagement/synthetic-assay-inference.json',
    'Examples/btk-aggregate-benchmark/observations.json',
    'Examples/target-engagement/finite-drug.json', 'Examples/target-engagement/design-candidates.json',
]
def identity(path):
    data=Path(path).read_bytes()
    return dict(path=str(path),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
                gitBlobSHA1=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest())
reports={name:json.loads((out/f'{name}.json').read_text()) for name in ['PortableChecks','DomainChecks']}
manifest=dict(schemaVersion=1,recordedAtUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    scope='Selected production Swift numeric sources with harness artifact-identity adapters; not a full package build',
    platform=platform.platform(),toolchain=(out/'toolchain.log').read_text().strip(),
    productionFlags=['-O','-whole-module-optimization','-swift-version','6','-warnings-as-errors','-enable-testing'],
    harnessFlags=['-Onone','-parse-as-library','-swift-version','6','-warnings-as-errors'],
    assertionsPassed=sum(x['checksPassed'] for x in reports.values()),
    fullApplePackageBuilt=False,metalExecuted=False,clinicalPredictionValidated=False,
    randomSourceMode=mode,randomSource=identity(rng_source),compiledRandomExcerpt=identity(rng_fragment),
    sources=[identity(p) for p in files],reports={name:identity(out/f'{name}.json') for name in reports},
    limitations=['Artifact identity is a harness adapter, not production CryptoKit/store validation.',
                 'SBC-style ranks are descriptive with dependent SMC particles; 24 trials do not establish universal calibration.',
                 'Finite-pool conformance is for declared analytic and one independent eight-state ODE fixture.',
                 'Published clinical cohort facts are inspected as aggregates; test cohort predictions are explicitly synthetic.'])
# Only source and output basenames are published for temporary paths.
for key in ['randomSource','compiledRandomExcerpt']:
    manifest[key]['path']=Path(manifest[key]['path']).name
for entry in manifest['reports'].values(): entry['path']=Path(entry['path']).name
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
