## prints # edges, precision, recall for different algorithms
## for wntforexperiments runs.
from utilsPoirel import *

for experiment in ['wnt-all-receptors','netpath']:
    if experiment == 'wnt-all-receptors':
        negtypes = ['adjacent']
    else:
        negtypes = ['adjacent','file']

    for negtype in negtypes:

        print '\n\nEXPERIMENT: ',experiment,'NegType for last two columns:',negtype,'\n'

        if negtype == 'adjacent':
            print '|Method| #Edges| Thres| NonePrec| NoneRec| AdjPrec |AdjRec|'
        else:
            print '|Method |#Edges |Thres| NonePrec |NoneRec| FilePrec| FileRec|'
        ## PathLinker
        ktotest = [200,300,800,1000]
        numedges = []
        nonepr = []
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pathlinker/Wnt-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        lines = readColumns(infile,3,5,6)
        for k in ktotest:
            index = min([i for i in range(len(lines)) if float(lines[i][0]) > k])
            numedges.append(index+1)
            nonepr.append((float(lines[index][1]),float(lines[index][2]),float(lines[index][0])))

        adjpr = []
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pathlinker/Wnt-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        lines = readColumns(infile,3,5,6)
        for k in ktotest:
            index = min([i for i in range(len(lines)) if float(lines[i][0]) > k])
            adjpr.append((float(lines[index][1]),float(lines[index][2]),float(lines[index][0])))
        for i in range(len(ktotest)):
            print '|PathLinker (top %d paths) | %d | %d | %.3f | %.3f | %.3f | %.3f|' % \
                (ktotest[i],numedges[i],nonepr[i][2],nonepr[i][0],nonepr[i][1],adjpr[i][0],adjpr[i][1])


        ### PageRank
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pagerank/Wnt-q_0.50-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        lines = readColumns(infile,3,5,6)
        nonepr = [(float(lines[i][1]),float(lines[i][2]),float(lines[i][0])) for i in numedges]
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pagerank/Wnt-q_0.50-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        lines = readColumns(infile,5,6)
        adjpr = [(float(lines[i][0]),float(lines[i][1])) for i in numedges]
        for i in range(len(numedges)):
            print '|PageRank | %d | %e | %.3f | %.3f | %.3f | %.3f|' % \
                (numedges[i],nonepr[i][2]*-1,nonepr[i][0],nonepr[i][1],adjpr[i][0],adjpr[i][1])

        ### IPA
        nmaxes = [5,10,15]
        for nmax in nmaxes:
            infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/ipa/Wnt-nmax%d-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment,nmax)
            nonelines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
            infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/ipa/Wnt-nmax%d-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,nmax,negtype)
            adjlines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
            print '|IPA (nmax=%d) | %d |N/A | %.3f | %.3f | %.3f | %.3f|' % \
                (nmax,len(nonelines),float(nonelines[0][2]),\
                 float(nonelines[0][3]),float(adjlines[0][2]),float(adjlines[0][3]))

        ### Shortest Paths
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/shortestpaths/Wnt-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        nonelines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/shortestpaths/Wnt-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        adjlines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        print '|ShortestPaths | %d |N/A | %.3f | %.3f | %.3f | %.3f|' % \
                (len(nonelines),float(nonelines[0][2]),\
                 float(nonelines[0][3]),float(adjlines[0][2]),float(adjlines[0][3]))

        ### PCSF
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pcsf/Wnt-prize5-omega0.01-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        nonelines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/pcsf/Wnt-prize5-omega0.01-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        adjlines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        print '|PCSF | %d | N/A| %.3f | %.3f | %.3f | %.3f|' % \
                (len(nonelines),float(nonelines[0][2]),\
                 float(nonelines[0][3]),float(adjlines[0][2]),float(adjlines[0][3]))


        ### RESPONSENET
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/responsenet/Wnt-gamma_20-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        nonelines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/responsenet/Wnt-gamma_20-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        adjlines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        print '|ResponseNet | %d | N/A | %.3f | %.3f | %.3f | %.3f|' % \
                (len(nonelines),float(nonelines[0][2]),\
                 float(nonelines[0][3]),float(adjlines[0][2]),float(adjlines[0][3]))

        ### anat
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/anat/Wnt-alpha0.00-exclude_none-sample_50X-edge-precision-recall.txt' % (experiment)
        nonelines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        infile = 'results/pathlinker-signaling-children-reg/weighted/%s/precision-recall/anat/Wnt-alpha0.00-exclude_%s-sample_50X-edge-precision-recall.txt' % (experiment,negtype)
        adjlines = [(u,v,p,r) for u,v,val,p,r in readColumns(infile,1,2,3,5,6) if val != 'Inf']
        print '|ANAT | %d | N/A | %.3f | %.3f | %.3f | %.3f|' % \
                (len(nonelines),float(nonelines[0][2]),\
                 float(nonelines[0][3]),float(adjlines[0][2]),float(adjlines[0][3]))
