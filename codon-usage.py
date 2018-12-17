from Bio import SeqIO,Entrez
from matplotlib import pyplot as plt
# %matplotlib inline (uncomment if you want inline graphics in a Jupyter notebook)

colors=['r','b']
Entrez.email="carlton.petermark.3v@kyoto-u.ac.jp"

block=1000
codontable={}

handle=Entrez.esearch(retmax=100000, 
                      db="nucleotide", 
                      term="158[BioProject]", 
                      idtype="acc")
gs=Entrez.read(handle)
genelist=gs['IdList']
handle.close()

runs=int(len(genelist)/block)+1
recorded_proteins=0
recorded_codons=0

for run in range(runs):
    print("run",run+1," of ",runs)
    h2=Entrez.efetch(db="nucleotide", 
                     id=genelist[(run*block):((run+1)*block)], 
                     rettype="gb")
    for t1 in SeqIO.parse(h2,format="genbank"):
        if len(t1.features)>2:
            if 'translation' in t1.features[2].qualifiers:
                recorded_proteins+=1
                dna=str(t1.seq)
                pro=str(t1.seq.translate())
                for i in range(len(pro)):
                    recorded_codons+=1
                    codon=dna[(i*3):(i*3+3)]
                    amino=pro[i]
                    if not(amino in codontable):
                        codontable[amino]={}
                    if not(codon in codontable[amino]):
                        codontable[amino][codon]=0
                    codontable[amino][codon]+=1

print(codontable)
print("proteins examined: %i" % recorded_proteins)
print("codons examined: %i" % recorded_codons)

## we are making a new dictionary now with the relative abundance for each codon
cdict={};
for i in sorted(codontable):
    isum=sum(list(codontable[i].values()))
    for j in codontable[i]:
        cdict[i+"_"+j]=codontable[i][j]/isum

        
## now we can plot the result:
plt.figure(figsize=(18,6))
e=0;ee=0
for aa in sorted(codontable):
    isum=sum(list(codontable[aa].values()))
    ee=ee+1
    for cd in sorted(codontable[aa]):
        plt.bar(e-0.4,codontable[aa][cd]/isum,label=aa,color=colors[ee%2])
        e=e+1

plt.xticks(range(e), sorted(cdict.keys()) , rotation='vertical',family='monospace');
plt.margins(0.01)
plt.grid(axis='x',linewidth='0.5',color='0.1')
plt.savefig("C_elegans_codons_full.png")
plt.close()
