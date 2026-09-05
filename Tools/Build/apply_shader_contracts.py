#!/usr/bin/env python3
"""One-time source migration; the published Metal files remain authority."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
resources = root / 'Sources/NumiVivoShaders/Resources'
changes, entry_count = {}, 0
for path in sorted(resources.glob('*.metal')):
    text = path.read_text()
    if not re.search(r'^namespace \w+\s*\{', text, re.M):
        continue
    names = re.findall(r'^kernel void (\w+)\(', text, re.M)
    if not names:
        continue
    if len(names) != len(set(names)):
        raise RuntimeError('duplicate kernel host name within ' + path.name)
    updated, count = re.subn(r'^kernel void (\w+)\(',
        lambda m: '[[host_name("' + m[1] + '")]] kernel void ' + m[1] + '(', text, flags=re.M)
    if count != len(names): raise RuntimeError('incomplete kernel migration')
    entry_count += count
    changes[path] = updated

path = resources / 'NumiVivoProgramPackRuntime.metal'
text = changes.get(path, path.read_text())
# All eight exact tokens are a kinetic local, not half-precision scalar types.
uses = len(re.findall(r'\bhalf\b', text))
if uses not in (0, 8): raise RuntimeError('ambiguous half identifier migration')
changes[path] = re.sub(r'\bhalf\b', 'halfSaturation', text)
path = root / 'Sources/NumiVivoKit/Kinetics/VivoTargetEngagementProgramSource.swift'
text = path.read_text()
for before, after in [('[targetSpecies[0], drugInput]', '[targetSpecies[0]]'),
                      ('[targetSpecies[0], competitorInput]', '[targetSpecies[0]]')]:
    if text.count(before) > 1: raise RuntimeError('ambiguous reservoir replacement')
    text = text.replace(before, after)
changes[path] = text
# All checks finish before writing. This changes the externally maintained
# reservoir model only; finite-drug closed-system kinetics are not altered.
for path, text in changes.items(): path.write_text(text)
print('Explicit namespace-independent kernel host names:', entry_count)
