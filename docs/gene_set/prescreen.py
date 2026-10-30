import sys,re,csv,collections,json
want=set(open(sys.argv[1]).read().split())
out=sys.argv[2]
PLP={'Pathogenic','Likely pathogenic','Pathogenic/Likely pathogenic'}
BLB={'Benign','Likely benign','Benign/Likely benign'}
AA3={'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val'}
PRE=re.compile(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)')
plp=collections.defaultdict(set); blb=collections.defaultdict(set); mx=collections.defaultdict(int)
r=csv.reader(sys.stdin,delimiter='\t'); h={c:i for i,c in enumerate(next(r))}
gi,ni,si,ai,ti=h['GeneSymbol'],h['Name'],h['ClinicalSignificance'],h['Assembly'],h['Type']
for row in r:
    try:
        g=row[gi]
        if g not in want: continue
        if row[ai]!='GRCh38' or row[ti]!='single nucleotide variant': continue
        m=PRE.search(row[ni])
        if not m: continue
        ref,pos,alt=m.group(1),int(m.group(2)),m.group(3)
        if ref not in AA3 or alt not in AA3 or ref==alt: continue
        if pos>mx[g]: mx[g]=pos
        s=row[si].strip()
        if s in PLP: plp[g].add(pos)
        elif s in BLB: blb[g].add(pos)
    except (IndexError,ValueError): continue
buyuk,ok,tek,yok=[],[],[],[]
for g in sorted(want):
    p,b,L=len(plp[g]),len(blb[g]),mx[g]
    if L==0: yok.append(g)
    elif L>2700: buyuk.append((g,p,b,L))
    elif p>=1 and b>=1: ok.append((g,p,b,L,round(p/L,4)))
    else: tek.append((g,p,b,L))
print(f"  ACMG SF v3.2 — {len(want)} gen")
print(f"  ├── AlphaFold reddeder (>2700 aa) : {len(buyuk):>2}  {', '.join(x[0] for x in buyuk)}")
print(f"  ├── bir sınıf boş → BLOCKED       : {len(tek):>2}  {', '.join(x[0] for x in tek)}")
print(f"  ├── ClinVar'da missense yok       : {len(yok):>2}  {', '.join(yok) if yok else '—'}")
print(f"  └── ÇALIŞTIRILABİLİR              : {len(ok):>2}")
json.dump(ok, open(out+'/acmg_runnable.json','w'))
