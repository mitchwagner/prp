from utilsPoirel import *
import sys
import numpy
import glob

ppifile = sys.argv[1]
print ppifile
pathwaydir = sys.argv[2]
print pathwaydir

print '\nPPI Table'
lines = readColumns(ppifile,1,2,3,4)
evidences = {'NetPath':set(),'KEGG':set(),'CSBDB':set(),'SPIKE':set()}
weights = {'NetPath':[],'KEGG':[],'CSBDB':[],'SPIKE':[]}
edgeweights = []
multev = set()
edges = set()
for u,v,w,t in lines:
  edges.add((u,v))
  edgeweights.append(float(w))
  for evtype in t.split('|'):
    if 'NetPath' in evtype:
	evidences['NetPath'].add((u,v))
        weights['NetPath'].append(float(w))
    elif 'KEGG' in evtype:
	evidences['KEGG'].add((u,v))
        weights['KEGG'].append(float(w))
    elif 'SPIKE' in evtype:
	evidences['SPIKE'].add((u,v))
        weights['SPIKE'].append(float(w))
    else:
	evidences['CSBDB'].add((u,v))
        weights['CSBDB'].append(float(w))
  if len(t.split('|'))>=2:
    multev.add((u,v))

print 'Type\t#Dir\t#Undir\tTot\tMeanWeight\tStdWeight'
for ev in evidences:
  directed = set([(u,v) for u,v in evidences[ev] if (v,u) not in evidences[ev]])
  undirected = evidences[ev].difference(directed)
  print '%s\t%d\t%d\t%d\t%.2f\t%.2f' % (ev,len(directed),len(undirected)/2.0,len(evidences[ev]),numpy.mean(weights[ev]),numpy.std(weights[ev]))

directed = set([(u,v) for u,v in edges if (v,u) not in edges])
undirected = edges.difference(directed)
print 'Total\t%d\t%d\t%d\t%.2f\t%.2f' % (len(directed),len(undirected)/2.0,len(edges),numpy.mean(edgeweights),numpy.std(edgeweights))

print '\nPos:Neg ratio'
lines = readColumns(ppifile,1,2)
edges = set([tuple(sorted([u,v])) for u,v in lines])
nodes = set([u for u,v in edges]).union([v for u,v in edges])

posedges = set()
for f in glob.glob('%s/*-edges.txt' % (pathwaydir)):
  lines = readColumns(f,1,2)
  posedges.update(set([tuple(sorted([u,v])) for u,v in lines]))
posnodes = set([u for u,v in posedges]).union([v for u,v in posedges])

posedges = posedges.intersection(edges)
posnodes = posnodes.intersection(nodes)

negedges = edges.difference(posedges)
negnodes = nodes.difference(posnodes)

print '%d total nodes and %d total edges'  % (len(nodes),len(edges))
print '%d (%.4f) positive nodes and %d (%.4f) positive edges' % (len(posnodes),len(posnodes)/float(len(nodes)),len(posedges),len(posedges)/float(len(edges)))
print '%d (%.4f) negative nodes and %d (%.4f) negative edges' % (len(negnodes),len(negnodes)/float(len(nodes)),len(negedges),len(negedges)/float(len(edges)))
print '%.4f pos:neg node ratio and %.4f pos:neg edge ratio' % (float(len(posnodes))/len(negnodes),float(len(posedges))/len(negedges))
