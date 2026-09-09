import Foundation

/// Self-contained, offline research review client. Untrusted input is carried as
/// base64 and rendered only with textContent. No analytics, remote fonts or uploads.
public enum VivoNeoantigenHTML {
    public static func render(_ report: VivoNeoantigenReport) throws -> Data {
        let bytes = try VivoCanonicalJSON.encode(report)
        let fingerprint = try VivoNeoantigenWorkbench.digest(bytes)
        let payload = bytes.base64EncodedString()
        let html = #"""
        <!doctype html><html lang="en"><head><meta charset="utf-8">
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; connect-src 'none'; base-uri 'none'; form-action 'none'">
        <title>NumiVivo · Neoantigen research review</title>
        <style>
        :root{color-scheme:light;--ink:#142c3b;--muted:#4f626d;--line:#d5dfe4;--paper:#f4f7f8;--accent:#145b5b}
        *{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 system-ui,-apple-system,sans-serif}
        main{max-width:1160px;margin:auto;padding:40px 24px 64px}header{border-bottom:1px solid var(--line);padding-bottom:24px}
        .eyebrow{font-size:12px;text-transform:uppercase;letter-spacing:.13em;font-weight:750;color:var(--accent)}h1{font-size:clamp(28px,5vw,42px);line-height:1.15;letter-spacing:-.025em;margin:12px 0}h2{font-size:21px;margin:0 0 12px}h3{margin:0;font-size:18px}p{margin:8px 0}.muted{color:var(--muted)}
        .notice{background:#fff5da;border:1px solid #d8c08c;border-radius:10px;padding:16px 20px;margin:24px 0}.grid{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}.metric,.panel,.candidate{background:white;border:1px solid var(--line);border-radius:12px;padding:20px}.metric strong{display:block;font-size:30px}.panel{margin-top:22px}
        .toolbar{display:flex;gap:14px;flex-wrap:wrap;align-items:end;margin:18px 0}.toolbar label{flex:1;min-width:180px}label{display:block;font-size:14px;font-weight:650}input,select,textarea,button{font:inherit}input,select,textarea{width:100%;padding:10px 12px;border:1px solid #aebdc5;border-radius:7px;background:white;color:var(--ink)}textarea{min-height:80px;resize:vertical}button{border:0;border-radius:8px;padding:11px 18px;background:var(--accent);color:white;font-weight:650;cursor:pointer}button:disabled{opacity:.45;cursor:not-allowed}button:focus-visible,input:focus-visible,select:focus-visible,textarea:focus-visible{outline:3px solid #78a6bf;outline-offset:2px}
        .candidate{margin:14px 0}.candidate-top{display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap}.pill{background:#edf2f5;border-radius:6px;padding:3px 8px;font-size:12px;font-weight:700}.facts{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}.facts dt{font-size:12px;color:var(--muted)}.facts dd{margin:4px 0;font-weight:600;overflow-wrap:anywhere}.gaps{padding:12px 16px;background:#f6f8f9;border-radius:7px}.gap{margin:6px 0;font-size:14px}.review-controls{display:grid;grid-template-columns:220px 1fr;gap:14px;margin-top:16px}details{margin-top:14px}summary{cursor:pointer;font-weight:650}table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13px}td{padding:7px;border-top:1px solid var(--line);text-align:left;vertical-align:top;overflow-wrap:anywhere}td:first-child{width:40%;font-weight:600}code{font-size:12px;overflow-wrap:anywhere}.finding{border-left:3px solid #d5b05c;padding:6px 14px;margin:16px 0}footer{margin-top:32px;color:var(--muted);font-size:13px}#status{min-height:24px}a{color:var(--accent)}
        @media(max-width:680px){main{padding:24px 16px}.grid{gap:8px}.metric{padding:12px}.metric strong{font-size:24px}.metric span{font-size:12px}.facts{grid-template-columns:1fr 1fr}.review-controls{grid-template-columns:1fr}.candidate,.panel{padding:16px}}
        </style></head><body><main>
        <header><div class="eyebrow">NumiVivo / research workbench</div><h1>Review the evidence.<br>Keep uncertainty visible.</h1><p class="muted">A local review of imported neoantigen predictions. No data leave this page.</p></header>
        <div class="notice"><strong id="classification"></strong><p>Not a vaccine prescription or a clinical result. Retaining a candidate means retaining it for research only.</p></div>
        <div class="grid"><div class="metric"><strong id="count"></strong><span>Imported candidate rows</span></div><div class="metric"><strong id="missing"></strong><span>Rows with additional evidence gaps</span></div><div class="metric"><strong id="decisionCount">0</strong><span>Research decisions drafted</span></div></div>
        <section class="panel"><h2>Case and next actions</h2><p id="case"></p><p id="source" class="muted"></p><div id="findings"></div><p id="integrity" role="status">Checking report byte integrity…</p><details><summary>Provenance and scope</summary><div id="provenance"></div><div id="limitations"></div></details></section>
        <section class="panel"><h2>Candidate evidence</h2><p class="muted">Rows remain in source order. No clinical ranking or automated acceptance is applied.</p><div class="toolbar"><label>Search gene, HLA or peptide<input id="search" type="search" placeholder="Search candidates" autocomplete="off"></label><label>Show<select id="filter"><option value="all">All rows</option><option value="gaps">Additional evidence gaps</option><option value="drafted">Drafted decisions</option></select></label></div><p id="empty" class="muted"></p><div id="candidates"></div></section>
        <section class="panel"><h2>Export research review</h2><p class="muted">The export binds your decisions to this exact report. Your reviewer ID is self-declared, not an authenticated signature.</p><div class="toolbar"><label>Reviewer identifier<input id="reviewer" placeholder="researcher-01" maxlength="128" autocomplete="off"></label><button id="export" disabled>Export review JSON</button></div><p id="status" role="status" aria-live="polite"></p><p class="muted">Exports stay on your device. Keep source data, reports and review files in institution-approved storage.</p></section>
        <footer><p>Research-use reference-data release · no patient intake, treatment authorization or manufacturing functions.</p><p>Report SHA-256: <code id="digest"></code></p></footer>
        </main><script>
        'use strict';
        const encoded = '\#(payload)', expected = '\#(fingerprint)';
        const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
        const report=JSON.parse(new TextDecoder().decode(bytes));
        const $=id=>document.getElementById(id), drafts=new Map(); let integrity=false;
        const text=(tag,value,cls)=>{const n=document.createElement(tag);n.textContent=value;if(cls)n.className=cls;return n};
        const value=x=>x===null||x===undefined?'Not supplied':String(x);
        $('classification').textContent=report.manifest.dataClass==='synthetic'?'SYNTHETIC DEMONSTRATION — no biological measurements':'PUBLIC REFERENCE DATA — research use only';
        $('case').textContent=report.manifest.caseID+' · '+report.manifest.assembly+' · '+report.state;
        $('source').textContent=report.manifest.sourceCitation;
        $('count').textContent=report.candidates.length;
        $('missing').textContent=report.candidates.filter(c=>c.evidenceGaps.length>1).length;
        $('digest').textContent=expected;
        for(const f of report.findings){const n=text('div','','finding');n.append(text('strong',f.message),text('p','Next action: '+f.nextAction),text('small','Responsible: '+f.owner));$('findings').append(n)}
        for(const l of report.limitations)$('limitations').append(text('p',l,'gap'));
        const provenance=[['Adapter',report.manifest.externalRun.format],['Declared pVACseq version',report.manifest.externalRun.version],['Annotation',report.manifest.annotationRelease],['HLA nomenclature',report.manifest.hlaNomenclatureRelease],['Implementation SHA-256',report.implementationSHA256],['Case fingerprint',report.caseFingerprint],['Source TSV SHA-256',report.sourceSHA256],['Declared predictor versions',JSON.stringify(report.manifest.externalRun.modelVersions)]];
        function table(entries){const t=document.createElement('table');for(const [k,v] of entries){const tr=document.createElement('tr');tr.append(text('td',k),text('td',String(v)));t.append(tr)}return t}
        $('provenance').append(table(provenance));
        function updateCount(){$('decisionCount').textContent=drafts.size}
        function draw(){const q=$('search').value.trim().toLowerCase(),f=$('filter').value;$('candidates').replaceChildren();let count=0;
          for(const c of report.candidates){if(![c.gene,c.transcript,c.hlaAllele,c.mutantPeptide].join(' ').toLowerCase().includes(q)||(f==='gaps'&&c.evidenceGaps.length<=1)||(f==='drafted'&&!drafts.has(c.id)))continue;count++;
            const card=text('article','','candidate'),top=text('div','','candidate-top');top.append(text('h3',c.gene),text('span','Source line '+c.sourceLine+' · not clinically qualified','pill'));card.append(top);
            const facts=text('dl','','facts');for(const [k,v] of [['Mutant peptide',c.mutantPeptide],['HLA allele',c.hlaAllele],['Predicted median binding (nM)',value(c.predictedMedianBindingNM)],['Gene expression (source units)',value(c.geneExpression)],['Tumor DNA VAF',value(c.tumorDNAFraction)],['Tumor RNA VAF',value(c.tumorRNAFraction)],['Normal DNA VAF',value(c.normalDNAFraction)],['Variant interval (0-based)',c.chromosome+':'+c.start+'–'+c.stop]]){const d=document.createElement('div');d.append(text('dt',k),text('dd',v));facts.append(d)}card.append(facts);
            const gaps=text('div','','gaps');gaps.append(text('strong','Evidence requiring review'));for(const gap of c.evidenceGaps)gaps.append(text('p',gap,'gap'));card.append(gaps);
            const details=document.createElement('details');details.append(text('summary','Inspect original source fields'),table(Object.entries(c.sourceFields)));card.append(details);
            const controls=text('div','','review-controls'),label=text('label','Research disposition'),select=document.createElement('select');select.setAttribute('aria-label','Disposition for '+c.gene);
            for(const [v,t] of [['','Not reviewed'],['retainForResearch','Retain for research'],['exclude','Exclude'],['deferReview','Defer pending evidence']]){const o=text('option',t);o.value=v;select.append(o)}select.value=drafts.get(c.id)?.disposition||'';label.append(select);
            const reasonLabel=text('label','Rationale (required for export)'),reason=document.createElement('textarea');reason.maxLength=4096;reason.setAttribute('aria-label','Rationale for '+c.gene);reason.value=drafts.get(c.id)?.rationale||'';reasonLabel.append(reason);controls.append(label,reasonLabel);card.append(controls);
            const change=()=>{if(select.value)drafts.set(c.id,{candidateID:c.id,disposition:select.value,rationale:reason.value});else drafts.delete(c.id);updateCount()};select.addEventListener('change',change);reason.addEventListener('input',change);$('candidates').append(card)
          }$('empty').textContent=count?'Showing '+count+' row(s).':report.state==='blocked'?'Import blocked. Resolve the case findings before candidate review.':report.state==='noCandidates'?'The supplied report contains no candidate rows. This does not establish absence of neoantigens.':'No rows match this view.'
        }
        $('search').addEventListener('input',draw);$('filter').addEventListener('change',draw);draw();
        (async()=>{try{if(!globalThis.crypto?.subtle)throw Error('unavailable');const actual=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),b=>b.toString(16).padStart(2,'0')).join('');integrity=actual===expected;if(!integrity)throw Error('mismatch');$('integrity').textContent='Report byte integrity verified locally. Scientific correctness and authorship are not verified.';$('export').disabled=report.state==='blocked'||!report.candidates.length}catch{$('integrity').textContent='Report integrity could not be verified in this browser. Export is disabled; use native verification.'}})();
        $('export').addEventListener('click',()=>{const id=$('reviewer').value.trim();if(!integrity||!/^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$/.test(id)){$('status').textContent='Enter a bounded reviewer identifier before exporting.';return}const decisions=report.candidates.filter(c=>drafts.has(c.id)).map(c=>drafts.get(c.id));if(!decisions.length||decisions.some(d=>!d.rationale.trim()||new TextEncoder().encode(d.rationale).length>4096)){$('status').textContent='Choose at least one disposition and provide a rationale of at most 4,096 UTF-8 bytes for every drafted decision.';return}const review={schema:'numivivo.org/neoantigen-review/v1',reportSHA256:expected,reviewerID:id,reviewedAt:new Date().toISOString().replace(/\.\d{3}Z$/,'Z'),decisions};const blob=new Blob([JSON.stringify(review,null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='neoantigen-research-review.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);$('status').textContent='Research review exported. No treatment or manufacturing action was authorized.'});
        </script></body></html>
        """#
        return Data(html.utf8)
    }
}
