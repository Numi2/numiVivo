#!/usr/bin/env python3
"""One-time checked source migration; no runtime rewriting or solver fallback."""
from pathlib import Path
import json
import re
root=Path(__file__).resolve().parents[2]
changes={}
def replace(path,before,after):
    p=root/path;text=changes.get(p,p.read_text())
    if text.count(before)==1:text=text.replace(before,after)
    elif text.count(after)!=1:raise RuntimeError('ambiguous source context: '+path)
    changes[p]=text
for name,new in [('NumiVivoMDPME','halfSpan'),('NumiVivoHybridExecution','halfSaturation'),('NumiVivoKernels','halfSaturation')]:
    p=root/('Sources/NumiVivoShaders/Resources/'+name+'.metal')
    changes[p]=re.sub(r'\bhalf\b',new,p.read_text())
include='#include "NumiVivoErrorFunctions.metalh"'
for name in ['NumiVivoMDPME','NumiVivoMDPMECorrections']:
    p=root/('Sources/NumiVivoShaders/Resources/'+name+'.metal');text=changes.get(p,p.read_text())
    if include not in text:text=text.replace('#include <metal_stdlib>','#include <metal_stdlib>\n'+include,1)
    text=re.sub(r'(?<!::)\berfc\(', 'nvivo_math::erfc(',text)
    text=re.sub(r'(?<!::)\berf\(', 'nvivo_math::erf(',text)
    text=text.replace('erfv*ir*ir2-(2.0f*beta*0.5641895835477563f)*gaussian*ir2','nvivo_math::erfMinusGaussian(br)*ir*ir2')
    changes[p]=text
replace('Package.swift','.copy("Resources/NumiVivoMetalABI.h"),','.copy("Resources/NumiVivoMetalABI.h"),\n            .copy("Resources/NumiVivoErrorFunctions.metalh"),')
loader='var source = try String(contentsOf: url, encoding: .utf8)\n            let mathInclude = '+json.dumps(include)+'''
            if source.contains(mathInclude) {
                guard let header = Bundle.module.url(forResource: "NumiVivoErrorFunctions", withExtension: "metalh")
                    ?? Bundle.module.url(forResource: "NumiVivoErrorFunctions", withExtension: "metalh", subdirectory: "Resources") else {
                    throw NumiVivoShaderError.sourceResourceMissing
                }
                source = source.replacingOccurrences(of: mathInclude, with: try String(contentsOf: header, encoding: .utf8))
            }'''
replace('Sources/NumiVivoShaders/MetalShaderLibrary.swift','let source = try String(contentsOf: url, encoding: .utf8)',loader)
replace('Sources/NumiVivoKit/Kinetics/VivoTargetEngagementProgramSource.swift','"expression": ["any": [','"expression": ["all": [')
replace('Sources/NumiVivoKit/Kinetics/VivoTargetEngagementProgramSource.swift','["gt": [total,','["lte": [total,')
replace('Sources/NumiVivoKit/Kinetics/VivoTargetEngagementProgramSource.swift','["lt": [total,','["gte": [total,')
replace('Tests/NumiVivoIntegrationTests/AppleExecutionTests.swift','#expect(pending.canCommit)','try #require(pending.canCommit)')
for p,text in changes.items():p.write_text(text)
print('Updated',len(changes),'source files with explicit math and invariant contracts')
