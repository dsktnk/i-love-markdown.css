#!/usr/bin/env python3
import hashlib, importlib.util, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'checkpoints'/'cp441'/'stage_b_result.json'
PREREG_BLOB_SHA='6951050bfe923f68e1f60a520ad30f4a302a1cc5'
SEED=int(PREREG_BLOB_SHA[:16],16)
N_PERM=9999
GROUPS=('Bio','HerbalA','HerbalB','PharmaA','Stars')
CATS=('SH','CH','QX','OX','D','OL','A','L','O','Y','K','T','P','S','R','OTHER')
K=len(CATS)
SMOOTH=0.5
LO=-6.0; HI=12.0; TOL=1e-8; MAXITER=512

CP218_PATH=ROOT/'checkpoints'/'cp218'/'run_cp218.py'
spec=importlib.util.spec_from_file_location('cp218_frozen',CP218_PATH)
cp=importlib.util.module_from_spec(spec); spec.loader.exec_module(cp)

def canonical_paragraph(p):
    n=len(p['lines']); where={}
    for t,line in enumerate(p['lines']):
        for tok in set(line): where.setdefault(tok,set()).add(t)
    ids={tok:sorted(s) for tok,s in where.items() if len(s)>=3}
    if not ids:return None
    return {'p':p,'risk':sum(n-1-v[0] for v in ids.values()),'pos':sum(len(v)-1 for v in ids.values())}

def fold_of(folio):
    return int(hashlib.sha256(str(folio).encode()).hexdigest(),16)%5

def build_lines():
    pars=cp.parse(cp.fetch_source())
    rows=[]; risk=pos=npars=0
    for p in pars:
        z=canonical_paragraph(p)
        if not z:continue
        npars+=1; risk+=z['risk']; pos+=z['pos']
        for line in p['lines']:
            seq=[cp.represent(tok,'PREFIX16') for tok in line]
            cnt=Counter(seq)
            rows.append({'group':p['group'],'folio':str(p['folio']),'fold':fold_of(p['folio']),
                         'counts':tuple(cnt[c] for c in CATS),'length':len(seq)})
    if (npars,risk,pos)!=(247,6154,1827):raise RuntimeError('canonical paragraph mismatch')
    if len(rows)!=1954 or sum(r['length'] for r in rows)!=16148 or len({r['folio'] for r in rows})!=82:
        raise RuntimeError('canonical line/token/folio mismatch')
    return rows

def log_mult_coeff(n):
    L=sum(n)
    return math.lgamma(L+1)-sum(math.lgamma(x+1) for x in n)

def train_p(train):
    C=[0]*K
    for r in train:
        for j,x in enumerate(r['counts']):C[j]+=x
    den=sum(C)+K*SMOOTH
    return tuple((x+SMOOTH)/den for x in C)

def ll_m0(n,p):
    return log_mult_coeff(n)+sum(x*math.log(q) for x,q in zip(n,p) if x)

def ll_dm(n,p,kappa):
    L=sum(n); out=log_mult_coeff(n)+math.lgamma(kappa)-math.lgamma(kappa+L)
    for x,q in zip(n,p):
        a=kappa*q
        out+=math.lgamma(a+x)-math.lgamma(a)
    return out

def train_dm_ll(train,p,logk):
    kappa=math.exp(logk)
    return sum(ll_dm(r['counts'],p,kappa) for r in train)

def golden_max(fn,lo=LO,hi=HI,tol=TOL,maxiter=MAXITER):
    phi=(1+5**0.5)/2
    inv=1/phi
    c=hi-(hi-lo)*inv; d=lo+(hi-lo)*inv
    fc=fn(c); fd=fn(d)
    it=0
    while (hi-lo)>tol and it<maxiter:
        if fc>fd:
            hi=d; d=c; fd=fc; c=hi-(hi-lo)*inv; fc=fn(c)
        else:
            lo=c; c=d; fc=fd; d=lo+(hi-lo)*inv; fd=fn(d)
        it+=1
    x=(lo+hi)/2
    return x,fn(x),it

def stage_a_geometry(rows,allowed):
    sec={}
    for g in allowed:
        rr=[r for r in rows if r['group']==g]; folds={}
        for f in range(5):
            test=[r for r in rr if r['fold']==f]; train=[r for r in rr if r['fold']!=f]
            toks=sum(r['length'] for r in train)
            nz=sum(sum(r['counts'][j] for r in train)>0 for j in range(K))
            folds[str(f)]={'test_lines':len(test),'train_lines':len(train),'train_tokens':toks,'nonzero_categories':nz}
        sec[g]={'lines':len(rr),'tokens':sum(r['length'] for r in rr),'folios':len({r['folio'] for r in rr}),
                'nonempty_folds':sum(folds[str(f)]['test_lines']>0 for f in range(5)),'folds':folds}
    return sec

def check_floors(sec):
    return (all(sec[g]['nonempty_folds']>=3 for g in sec) and
            all(d['train_lines']>=20 and d['train_tokens']>=100
                for g in sec for d in sec[g]['folds'].values() if d['test_lines']>0))

def run_universe(rows,allowed,archive=False):
    sub=[r for r in rows if r['group'] in allowed]
    geom=stage_a_geometry(sub,allowed)
    if not check_floors(geom):raise RuntimeError('Stage-A structural floor failed')
    sec_delta={g:0.0 for g in allowed}; folio=defaultdict(float); fits={}
    for g in allowed:
        rr=[r for r in sub if r['group']==g]
        fits[g]={}
        for f in range(5):
            test=[r for r in rr if r['fold']==f]
            if not test:continue
            train=[r for r in rr if r['fold']!=f]
            p=train_p(train)
            logk,obj,it=golden_max(lambda z:train_dm_ll(train,p,z))
            kappa=math.exp(logk)
            boundary=(abs(logk-LO)<1e-5 or abs(logk-HI)<1e-5)
            fits[g][str(f)]={'kappa':kappa,'log_kappa':logk,'iterations':it,'boundary':boundary,'train_lines':len(train)}
            for r in test:
                d=ll_dm(r['counts'],p,kappa)-ll_m0(r['counts'],p)
                sec_delta[g]+=d; folio[r['folio']]+=d
    pooled=sum(sec_delta.values())
    return {'pooled':pooled,'sections':sec_delta,'folio':dict(folio),'fits':fits,'geometry':geom}

def main():
    rows=build_lines()
    full=run_universe(rows,GROUPS,archive=True)
    vals=[full['folio'][k] for k in sorted(full['folio'],key=lambda x:int(x))]
    if abs(sum(vals)-full['pooled'])>1e-8:raise RuntimeError('folio/pooled mismatch')
    rng=random.Random(SEED); ge=0
    for _ in range(N_PERM):
        sim=sum(v if rng.getrandbits(1) else -v for v in vals)
        if sim>=full['pooled']:ge+=1
    p=(1+ge)/(N_PERM+1)
    loo={}
    for omit in GROUPS:
        allowed=tuple(g for g in GROUPS if g!=omit)
        loo[omit]=run_universe(rows,allowed)['pooled']
    pos=sum(full['sections'][g]>0 for g in GROUPS)
    gates={'pooled_delta_ll_positive':full['pooled']>0,
           'positive_sections_at_least_4of5':pos>=4,
           'folio_signflip_p_lt_0_01':p<.01,
           'all_full_loo_positive':all(loo[g]>0 for g in GROUPS)}
    passed=all(gates.values())
    result={'formal':('PREFIX16_COMPOSITION_OVERDISPERSION_PREDICTIVELY_CONFIRMED' if passed else
                      'NO_CONFIRMATORY_PREFIX16_COMPOSITION_OVERDISPERSION_GAIN'),
            'pass':passed,'prereg_blob_sha':PREREG_BLOB_SHA,'seed':SEED,'n_perm':N_PERM,
            'source_sha256':cp.EXPECTED_SHA,'pooled_delta_ll':full['pooled'],
            'section_delta_ll':full['sections'],'positive_sections':pos,
            'folio_blocks':len(vals),'folio_block_sums':full['folio'],
            'permutations_ge_observed':ge,'p_one_sided_plus_one':p,
            'fold_training_kappa':full['fits'],'full_leave_one_section_out_delta_ll':loo,
            'gates':gates,'stage_a_geometry':full['geometry'],
            'interpretation_ceiling':'PREFIX16 composition overdispersion conditional on observed line length only; no exact-token generator, grammar, syntax, language, cipher, plaintext, semantics, or unique causal generator.'}
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True),encoding='utf-8')
    print('CP441_STAGE_B_BEGIN');print(json.dumps(result,indent=2,sort_keys=True));print('CP441_STAGE_B_END')

if __name__=='__main__':main()
