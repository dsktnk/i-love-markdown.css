#!/usr/bin/env python3
import hashlib, importlib.util, json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'checkpoints'/'cp441'/'stage_a_counts.json'
PREREG_BLOB_SHA='6951050bfe923f68e1f60a520ad30f4a302a1cc5'
GROUPS=('Bio','HerbalA','HerbalB','PharmaA','Stars')
CATS=('SH','CH','QX','OX','D','OL','A','L','O','Y','K','T','P','S','R','OTHER')

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
    rows=[]; risk=pos=0; npars=0
    for p in pars:
        z=canonical_paragraph(p)
        if not z:continue
        npars+=1; risk+=z['risk']; pos+=z['pos']
        for line in p['lines']:
            seq=[cp.represent(tok,'PREFIX16') for tok in line]
            rows.append({'group':p['group'],'folio':str(p['folio']),'fold':fold_of(p['folio']),'seq':seq})
    if (npars,risk,pos)!=(247,6154,1827):
        raise RuntimeError(f'canonical mismatch {(npars,risk,pos)}')
    if len(rows)!=1954:raise RuntimeError(f'line mismatch {len(rows)}')
    toks=sum(len(r['seq']) for r in rows)
    if toks!=16148:raise RuntimeError(f'token mismatch {toks}')
    if len({r['folio'] for r in rows})!=82:raise RuntimeError('folio mismatch')
    return rows

def summarize(rows):
    out={}
    for g in GROUPS:
        sec=[r for r in rows if r['group']==g]
        folds={}
        for f in range(5):
            test=[r for r in sec if r['fold']==f]
            train=[r for r in sec if r['fold']!=f]
            cnt=Counter(x for r in train for x in r['seq'])
            lengths=[len(r['seq']) for r in train]
            folds[str(f)]={
                'test_lines':len(test),
                'test_folios':len({r['folio'] for r in test}),
                'train_lines':len(train),
                'train_folios':len({r['folio'] for r in train}),
                'train_tokens':sum(cnt.values()),
                'nonzero_categories':sum(cnt[c]>0 for c in CATS),
                'train_min_length':min(lengths) if lengths else None,
                'train_max_length':max(lengths) if lengths else None,
            }
        out[g]={
            'lines':len(sec),
            'tokens':sum(len(r['seq']) for r in sec),
            'folios':len({r['folio'] for r in sec}),
            'nonempty_folds':sum(folds[str(f)]['test_lines']>0 for f in range(5)),
            'folds':folds,
        }
    return out

def main():
    rows=build_lines(); sec=summarize(rows)
    floors={
      'all_sections_at_least_3_nonempty_folds':all(sec[g]['nonempty_folds']>=3 for g in GROUPS),
      'every_nonempty_test_fold_train_lines_at_least_20':all(
          d['train_lines']>=20 for g in GROUPS for d in sec[g]['folds'].values() if d['test_lines']>0),
      'every_nonempty_test_fold_train_tokens_at_least_100':all(
          d['train_tokens']>=100 for g in GROUPS for d in sec[g]['folds'].values() if d['test_lines']>0),
    }
    licensed=all(floors.values())
    result={
      'formal':'CP441_STAGE_A_LICENSED' if licensed else 'CP441_STAGE_A_NOT_FEASIBLE',
      'scientific_effect_inspected':False,
      'prereg_blob_sha':PREREG_BLOB_SHA,
      'source_sha256':cp.EXPECTED_SHA,
      'categories':CATS,
      'canonical':{'paragraphs':247,'lines':1954,'tokens':16148,'folios':82,'risk_rows':6154,'recurrence_events':1827},
      'sections':sec,'floors':floors
    }
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True),encoding='utf-8')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
