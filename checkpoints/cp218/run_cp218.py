#!/usr/bin/env python3
import hashlib, json, math, random, re, urllib.request
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

SOURCE_URL='https://www.voynich.nu/data/IT2a-n.txt'
EXPECTED_SHA='7f27a8b0feed8f6de0a99900df6bf912dd1d295c38e5f830bac8b41c3f536fb5'
OUT=Path(__file__).resolve().parent/'results'; OUT.mkdir(exist_ok=True)
NREP=8192; SEED=21820260824
GROUPS=('Bio','HerbalA','HerbalB','PharmaA','Stars')
REPS=('TOKEN','PREFIX16','COARSE6')
STAGES=('FIRST','LATER')
METRICS=('NEXT1','INV_GAP')
BIO=set(range(75,85)); STARS=set(range(103,117))
FINE=("SH","CH","QX","OX","D","OL","A","L","O","Y","K","T","P","S","R","OTHER")
COARSE_MAP={
 'QX':'Q','OX':'Q','SH':'SHCH','CH':'SHCH','D':'DOL','OL':'DOL',
 'A':'VOWELLIKE','L':'VOWELLIKE','O':'VOWELLIKE','Y':'VOWELLIKE',
 'K':'GALLOWS_OTHER_PREFIX','T':'GALLOWS_OTHER_PREFIX','P':'GALLOWS_OTHER_PREFIX',
 'S':'GALLOWS_OTHER_PREFIX','R':'GALLOWS_OTHER_PREFIX','OTHER':'OTHER'
}
FIRST_TARGET={
 'NEXT1':{'Bio':0.08497311143193834,'HerbalA':0.098819476050479,'HerbalB':0.08667321600357315,'PharmaA':0.08963058800756074,'Stars':0.10458723949298963},
 'INV_GAP':{'Bio':0.06362337453848818,'HerbalA':0.06383892797159722,'HerbalB':0.06444343232996992,'PharmaA':0.055598134036812574,'Stars':0.06156018471963578}
}
LATER_TARGET={
 'NEXT1':{'Bio':0.12028686122038759,'HerbalA':0.09521606839867496,'HerbalB':0.039329216535098875,'PharmaA':0.0857583774250441,'Stars':0.08633375235981074},
 'INV_GAP':{'Bio':0.07306261637411987,'HerbalA':0.06132236633716184,'HerbalB':0.0344660904134381,'PharmaA':0.06272907743344253,'Stars':0.05477165896469922}
}
PAGE_RE=re.compile(r'^<f(\d+)([rv]\d*)>\s+<!\s*(.*?)>')
LOCUS_RE=re.compile(r'^<f(\d+)([rv]\d*)\.([^,>]+),([^>]+)>\s*(.*)$')
META_RE=re.compile(r'\$([A-Z])=([^\s>]+)')
INLINE_RE=re.compile(r'<[^>]*>'); CURLY_RE=re.compile(r'\{[^}]*\}')

def fetch_source():
    req=urllib.request.Request(SOURCE_URL,headers={'User-Agent':'voynich-research-cp218/1.0'})
    data=urllib.request.urlopen(req,timeout=30).read(); sha=hashlib.sha256(data).hexdigest()
    if sha!=EXPECTED_SHA: raise ValueError(f'checksum mismatch {sha}')
    return data.decode('utf-8','strict')

def tokenize(s):
    s=CURLY_RE.sub('.',s); s=INLINE_RE.sub('.',s); s=re.sub(r'[^A-Za-z?]+','.',s)
    return [x.lower() for x in s.split('.') if x]

def classify(tok):
    t=tok.lower()
    if '?' in t or not t.isalpha(): return 'OTHER'
    if t.startswith('sh'): return 'SH'
    if t.startswith('ch'): return 'CH'
    if t.startswith('qok') or t.startswith('qot'): return 'QX'
    if t.startswith('ok') or t.startswith('ot'): return 'OX'
    if t.startswith('d'): return 'D'
    if t.startswith('ol'): return 'OL'
    if t.startswith('a'): return 'A'
    if t.startswith('l'): return 'L'
    if t.startswith('o'): return 'O'
    if t.startswith('y'): return 'Y'
    if t.startswith('k'): return 'K'
    if t.startswith('t'): return 'T'
    if t.startswith('p'): return 'P'
    if t.startswith('s'): return 'S'
    if t.startswith('r'): return 'R'
    return 'OTHER'

def represent(tok,rep):
    if rep=='TOKEN': return tok.lower()
    f=classify(tok)
    if rep=='PREFIX16': return f
    if rep=='COARSE6': return COARSE_MAP[f]
    raise ValueError(rep)

def page_group(f,md):
    if f in BIO: return 'Bio'
    if f in STARS: return 'Stars'
    if md.get('I')=='H' and md.get('L')=='A': return 'HerbalA'
    if md.get('I')=='H' and md.get('L')=='B': return 'HerbalB'
    if md.get('I')=='P' and md.get('L')=='A': return 'PharmaA'
    return None

def parse(text):
    pars=[]; current=None; cur=None; pid=0
    for raw in text.splitlines():
        pm=PAGE_RE.match(raw)
        if pm:
            f=int(pm.group(1)); side=pm.group(2); md=dict(META_RE.findall(pm.group(3)))
            current={'folio':f,'side':side,'group':page_group(f,md)}; cur=None
            continue
        m=LOCUS_RE.match(raw)
        if not m or current is None: continue
        f=int(m.group(1)); side=m.group(2); tag=m.group(4); body=m.group(5)
        if f!=current['folio'] or side!=current['side']: cur=None; continue
        g=current['group']; p0=g is not None and tag.endswith('P0') and tag[:1] in {'@','+','=','*'}
        if not p0: cur=None; continue
        if '<%>' in body:
            cur={'id':pid,'group':g,'folio':f,'page':f'f{f}{side}','lines':[]}; pid+=1
        if cur is None: continue
        toks=tokenize(body)
        if len(toks)>=2: cur['lines'].append(toks)
        if '<$>' in body:
            if len(cur['lines'])>=3: pars.append(cur)
            cur=None
    return pars

@lru_cache(maxsize=None)
def expected_inv_gap(n,k):
    den=math.comb(n,k)
    return sum((1.0/d)*(math.comb(n-d,k-1)/den) for d in range(1,n-k+2))

def rep_lines(lines,rep):
    return [set(represent(t,rep) for t in line) for line in lines]

def stage_scores(lines,stage,order=None):
    ordered=lines if order is None else [lines[i] for i in order]
    where=defaultdict(list)
    for i,line in enumerate(ordered):
        for x in line: where[x].append(i)
    n=len(ordered); next1=[]; inv=[]
    for inds in where.values():
        k=len(inds)
        if stage=='FIRST':
            if k<2: continue
            gaps=[inds[1]-inds[0]]
        else:
            if k<3: continue
            gaps=[inds[j+1]-inds[j] for j in range(1,k-1)]
        en=k/n; ei=expected_inv_gap(n,k)
        for gap in gaps:
            next1.append((1.0 if gap==1 else 0.0)-en)
            inv.append((1.0/gap)-ei)
    if not next1: return None
    return {'NEXT1':sum(next1)/len(next1),'INV_GAP':sum(inv)/len(inv),'events':len(next1)}

def build(pars,rep,stage):
    rows=[]
    for p in pars:
        lines=rep_lines(p['lines'],rep); s=stage_scores(lines,stage)
        if s is not None:
            rows.append({'id':p['id'],'group':p['group'],'folio':p['folio'],'lines':lines,'obs':s})
    return rows

def section_means(rows,vals):
    sums=defaultdict(float); nums=Counter()
    for r,v in zip(rows,vals): sums[r['group']]+=v; nums[r['group']]+=1
    return {g:(sums[g]/nums[g] if nums[g] else None) for g in GROUPS}

def aggregate(sec):
    return sum(sec[g] for g in GROUPS)/len(GROUPS)

def observed(rows,metric):
    sec=section_means(rows,[r['obs'][metric] for r in rows]); return {'sections':sec,'global':aggregate(sec)}

def validity(rows):
    pc=Counter(r['group'] for r in rows)
    fol={g:len({r['folio'] for r in rows if r['group']==g}) for g in GROUPS}
    ok=all(pc[g]>=15 and fol[g]>=5 for g in GROUPS)
    return {'ok':ok,'paragraphs':dict(pc),'folios':fol}

def leave_one_section_out(sec):
    return {drop:sum(sec[g] for g in GROUPS if g!=drop)/(len(GROUPS)-1) for drop in GROUPS}

def randorder(n,rng):
    o=list(range(n)); rng.shuffle(o); return o

def main():
    pars=parse(fetch_source())
    data={(rep,stage):build(pars,rep,stage) for rep in REPS for stage in STAGES}
    val={(rep,stage):validity(data[(rep,stage)]) for rep in REPS for stage in STAGES}
    obs={(rep,stage,m):observed(data[(rep,stage)],m) for rep in REPS for stage in STAGES for m in METRICS}

    fp_first=max(abs(obs[('TOKEN','FIRST',m)]['sections'][g]-FIRST_TARGET[m][g]) for m in METRICS for g in GROUPS)
    fp_later=max(abs(obs[('TOKEN','LATER',m)]['sections'][g]-LATER_TARGET[m][g]) for m in METRICS for g in GROUPS)
    repro=(fp_first<=1e-9 and fp_later<=1e-9)

    ge={(rep,stage,m):{'global':0,**{g:0 for g in GROUPS}} for rep in REPS for stage in STAGES for m in METRICS}
    if repro:
        rng=random.Random(SEED)
        # Precompute representation lines by paragraph id for synchronized orders.
        rep_by_id={(rep,p['id']):rep_lines(p['lines'],rep) for rep in REPS for p in pars}
        for _ in range(NREP):
            orders={p['id']:randorder(len(p['lines']),rng) for p in pars}
            for rep in REPS:
                for stage in STAGES:
                    rows=data[(rep,stage)]; values={m:[] for m in METRICS}
                    for r in rows:
                        s=stage_scores(rep_by_id[(rep,r['id'])],stage,orders[r['id']])
                        if s is None: raise ValueError('eligibility changed under line permutation')
                        for m in METRICS: values[m].append(s[m])
                    for m in METRICS:
                        sec=section_means(rows,values[m]); glob=aggregate(sec); target=obs[(rep,stage,m)]
                        if glob>=target['global']: ge[(rep,stage,m)]['global']+=1
                        for g in GROUPS:
                            if sec[g]>=target['sections'][g]: ge[(rep,stage,m)][g]+=1

    den=NREP+1
    pvals={(rep,stage,m):{q:((ge[(rep,stage,m)][q]+1)/den if repro else None) for q in ('global',)+GROUPS} for rep in REPS for stage in STAGES for m in METRICS}

    combo={}
    rep_support={}
    for rep in REPS:
        valid_rep=all(val[(rep,stage)]['ok'] for stage in STAGES)
        combo_pass=True
        for stage in STAGES:
            for m in METRICS:
                o=obs[(rep,stage,m)]; p=pvals[(rep,stage,m)]; loo=leave_one_section_out(o['sections'])
                pos_sections=sum(o['sections'][g]>0 for g in GROUPS)
                sig_sections=sum(p[g]<.05 for g in GROUPS) if repro else 0
                passed=(valid_rep and o['global']>0 and p['global']<.01 and pos_sections>=4 and all(v>0 for v in loo.values()) and sig_sections>=2)
                combo[(rep,stage,m)]={'pass':passed,'positive_sections':pos_sections,'significant_sections':sig_sections,'loo':loo}
                combo_pass=combo_pass and passed
        rep_support[rep]=valid_rep and combo_pass

    if not repro:
        formal='CP218_REPRODUCIBILITY_GATE_FAILED'
    elif rep_support['COARSE6']:
        formal='CP218_COMMON_RECURRENCE_PROCESS_SURVIVES_COARSE6'
    elif rep_support['PREFIX16']:
        formal='CP218_COMMON_RECURRENCE_PROCESS_SURVIVES_PREFIX16_NOT_COARSE6'
    elif rep_support['TOKEN']:
        formal='CP218_COMMON_RECURRENCE_PROCESS_REQUIRES_TOKEN_LEVEL_UNDER_FROZEN_LADDER'
    else:
        formal='CP218_COMMON_RECURRENCE_PROCESS_NOT_RECONFIRMED'

    result={
      'formal':formal,'source_sha256':EXPECTED_SHA,'n_replicates':NREP,'seed':SEED,
      'fingerprint_first_max_abs_error':fp_first,'fingerprint_later_max_abs_error':fp_later,'reproducibility_pass':repro,
      'validity':{f'{rep}:{stage}':val[(rep,stage)] for rep in REPS for stage in STAGES},
      'observed':{f'{rep}:{stage}:{m}':obs[(rep,stage,m)] for rep in REPS for stage in STAGES for m in METRICS},
      'pvalues':{f'{rep}:{stage}:{m}':pvals[(rep,stage,m)] for rep in REPS for stage in STAGES for m in METRICS},
      'combo_gates':{f'{rep}:{stage}:{m}':combo[(rep,stage,m)] for rep in REPS for stage in STAGES for m in METRICS},
      'representation_support':rep_support
    }
    (OUT/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True),encoding='utf-8')
    L=['# CP218 Results','',f'Formal: `{formal}`','',f'FIRST fingerprint max abs error={fp_first}',f'LATER fingerprint max abs error={fp_later}',f'Reproducibility={"PASS" if repro else "FAIL"}','']
    for rep in REPS:
        L += [f'## {rep}',f'Representation support: **{"PASS" if rep_support[rep] else "FAIL"}**']
        for stage in STAGES:
            L.append(f'- {stage} validity: {val[(rep,stage)]}')
            for m in METRICS:
                o=obs[(rep,stage,m)]; p=pvals[(rep,stage,m)]; c=combo[(rep,stage,m)]
                L.append(f'  - {m}: global={o["global"]:.9f}, p={p["global"]:.9f}, sections={o["sections"]}, positive_sections={c["positive_sections"]}, significant_sections={c["significant_sections"]}, gate={"PASS" if c["pass"] else "FAIL"}')
        L.append('')
    (OUT/'SUMMARY.md').write_text('\n'.join(L),encoding='utf-8')
    print('CP218_RESULT_BEGIN'); print((OUT/'SUMMARY.md').read_text()); print('CP218_RESULT_END')

if __name__=='__main__': main()
