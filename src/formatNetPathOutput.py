from utilsPoirel import *
import sys

interactome = 'results/formatted/background-interactome-pathlinker-2015.txt'
npdir = 'results/pathlinker-signaling-children-reg/weighted/netpath/pathlinker'

pathways = readItemSet('data/netpath-analyzed-pathways.txt')

lines = readColumns(interactome,1,2,3,4)
ppiedges = {}
for u,v,w,t in lines:
    ppiedges[tuple([u,v])] = [w,t]

for pathway in pathways:
    pledges = readColumns('%s/%s-k_20000-paths.txt' % (npdir,pathway),1,3)
    edges = {}
    edgeorder = []
    for k,path in pledges:
        p = path.split('|')
        for i in range(len(p)-1):
            e = tuple([p[i],p[i+1]])
            if e not in edges:
                edges[e] = []
                edgeorder.append(e)
            edges[e].append(int(k))
    
    out = open('results/formatted/%s-pathlinker-20000_paths.txt'% (pathway),'w')
    out.write('#tail\thead\tedge_weight\tedge_type\tpaths\n')
    for e in edgeorder:
        out.write('%s\t%s\t%s\t%s\t%s\n' % (e[0],e[1],ppiedges[e][0],ppiedges[e][1],'|'.join([str(a) for a in edges[e]])))
    out.close()
    print 'Wrote %s' % (pathway)
