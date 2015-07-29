#!/usr/bin/python

## Code to post linker pathway results to GraphSpace
## (http://holmes.cs.vt.edu/graphs)
##
## Anna Ritz, Aug 2013
## Majority of code taken from post-predicted-pathway.py, post-example-pathway.py
## Required Modules: 
##  utilsPoirel: svn checkout file:///home/murali/svnrepo/src/python/scripts/trunk/utilsPoirel.py
##  Net: contained in trunk/net-graphspace/lib/ (add to PYTHONPATH)

## Import Statements
import sys
import os
import re
import urllib
import json
import glob
from  collections import defaultdict
from utilsPoirel import *
from optparse import OptionParser
import graphspace_interface as graphspace
import networkx as nx

# VERSIONS: see src/python/Interactomes/Human/generate-human-interactome.py.
## dictionary of {id : filename}
## for every version.  Old (unused) files:
## 'PANTHER':'/data/annaritz/datasets/DBParsers/trunk/PANTHER/edgefiles/',
VERSIONS = {
    '2013linker': { # w/ vinayagam
        'CSBDB':'%s/csbdb/2013-04-19-human_ppi-r.txt',
        'vinayagam': '%s/interactions/human-direction-vinayagam/2011-vinayagam-wanker-sci-signaling-human-directed-pin.txt',
        'NetPath':'%s/interactions/netpath/pathways/',
        'KEGG':'%s/interactions/kegg/2015-03-23/hsa/edge-files/',
        'SPIKE':'%s/interactions/spike/edge-files/',       
    },
    '2015pathlinker': { # no vinayagam
        'CSBDB':'%s/csbdb/2013-04-19-human_ppi-r.txt',
        'NetPath':'%s/interactions/netpath/pathways/',
        'KEGG':'%s/interactions/kegg/2015-03-23/hsa/edge-files/',
        'SPIKE':'%s/interactions/spike/edge-files/',       
    },   
}
#### Add versions that differ by KEGG only
# KEGG points to original files from 2013.
VERSIONS['pathlinker-old-kegg-buggy'] = VERSIONS['2015pathlinker'].copy()
VERSIONS['pathlinker-old-kegg-buggy']['KEGG'] = '/data/annaritz/datasets/SignalingDatabases/KEGG/hsa_pathways/'
# KEGG points to old, correct files from 2013.
VERSIONS['pathlinker-old-kegg-fixed'] = VERSIONS['2015pathlinker'].copy()
VERSIONS['pathlinker-old-kegg-fixed']['KEGG'] = '/data/annaritz/datasets/SignalingDatabases/KEGG/hsa_pathways_bidirected_edges/'
# KEGG points newly-parsed files from 2013.
VERSIONS['pathlinker-old-kegg-new-parse'] = VERSIONS['2015pathlinker'].copy()
VERSIONS['pathlinker-old-kegg-new-parse']['KEGG'] = '%s/interactions/kegg/2013-10-22/hsa/edge-files/'

#### Add versions with identical copies of the same set of evidence sources.
VERSIONS['pathlinker-old-kegg-buggy-oldannotations'] = VERSIONS['pathlinker-old-kegg-buggy'].copy()
VERSIONS['pathlinker-signaling'] = VERSIONS['2015pathlinker'].copy()
VERSIONS['pathlinker-signaling-children'] = VERSIONS['2015pathlinker'].copy()
VERSIONS['pathlinker-signaling-children-reg'] = VERSIONS['2015pathlinker'].copy()

# Filtered detection method (from parse-human-ppi.py)
METHODS_TO_IGNORE = ('MI:0036', # domain fusion
                       'MI:0046', # experimental knowledge based
                       'MI:0063', # interaction prediction
                       'MI:0064', # interologs mapping
                       'MI:0085', # phylogenetic profile
                       'MI:0087', # predictive text mining
                       'MI:0105', # structure based prediction
                       'MI:0363', # inferred by author
                       'MI:0364', # inferred by curator
                       'MI:0686', # unspecified method coexpression
                       'MI:0045', # experimental interaction detection
)

# Dictionaries of node and edge properties
COLORS = { 'inner':'gray',
           #'inner':'#5CD65C', #green
           'source':'#8CE1DC', #blue
           'target':'#FFFF60', #yellow
           'norm_edge':'black',
           'no_evid_edge':'red', #NOTE: might want to pick a diff color
           'red':'#FE2E2E', 
           'gray':'#C8C6C6',
           'pink':'#F5A9A9',
           'kegg':'#AC58FA',#'#FAAC58',
           'netpath':'#01DF01',#'#31B404',
           'both':'#CC2EFA',
           'neither':'#D8D8D8',#'#848484',
           'white':'#FFFFFF',
           'darkgray':'#6E6E6E',
           'crosstalk':'#FACC2E',#'#D0A9F5',
}
NODESHAPES = { 'target':'rectangle',
               'source':'triangle',
               'inner':'ellipse',
               'ligand':'hexagon',
}

# global variables will be populated by main()
DATADIR=''
MAPFILE=''
NAME2UNIPROT={}
UNIPROT2NAME={}
NAME2ENTREZ={}
ENTREZ2NAME={}
PPI = ''
PPIEDGES = []
PPIWEIGHTS = {}
EVIDENCEFILES={}

USERNAME = 'annaritz@vt.edu'
PASSWORD = 'platypus'
GROUP='2015PathLinkerReconstructions'

##########################################################
def getPredictedEdges(predfile,increase,decrease,thres):
    edges = {}
    if increase:
        lines = readColumns(predfile,1,2,3)
        edges = {(UNIPROT2NAME.get(u,u),UNIPROT2NAME.get(v,v)):float(t) for u,v,t in lines if float(t) <= thres}
    elif decrease:
        lines = readColumns(predfile,1,2,3)
        edges = {(UNIPROT2NAME.get(u,u),UNIPROT2NAME.get(v,v)):float(t) for u,v,t in lines if float(t) >= thres}
    else:
        lines = readColumns(predfile,1,2)
        edges = {(UNIPROT2NAME.get(u,u),UNIPROT2NAME.get(v,v)):0 for u,v in lines}

    print '%d edges in reconstruction' % (len(edges))

    undiredges = set([tuple(sorted([u,v])) for u,v in edges])
    print '%d edges in reconstruction excluding direction'  %(len(undiredges))

    # make sure all nodes are in PPI
    toremove = set()
    for u,v in edges.keys():
        if (u,v) not in PPIEDGES and (v,u) not in PPIEDGES:
            print 'WARNING: removing (%s,%s) from predicted edges.' % (u,v)
            toremove.add((u,v))
    
    edges = {e:edges[e] for e in edges if e not in toremove}
    print '%d edges remain after removing ones not in PPI.' % (len(edges))
    return edges

##########################################################
def getLigandInformation(ligandfile,prededges,receptors):
    uniprotligands = readItemSet(ligandfile,1)
    ligandnodes = [UNIPROT2NAME[u] for u in uniprotligands if u in UNIPROT2NAME]
    print '%d ligands (%s)' % (len(ligandnodes),','.join(ligandnodes))

    # only get edges that connect ligands to receptors in the predicted set.
    prednodes = set([u for u,v in prededges]).union(set([v for u,v in prededges]))
    predreceptors = receptors.intersection(prednodes)
    
    ligandedges = set([tuple(sorted([u,v])) for u,v in PPIEDGES if \
                       (u in ligandnodes and v in predreceptors) or \
                       (u in predreceptors and v in ligandnodes)])
    print '%d ligand edges' % (len(ligandedges))
    return ligandedges    

##########################################################
def constructGraph(receptors,tfs,prededges,increase,decrease,thres,netpath,kegg,prededgefile,undirected,nolabels,ligandedges):
    evidence = getEvidence(prededges)

    # NetworkX object
    G = nx.DiGraph(directed=True)

    # get metadata
    desc = getGraphDescription(increase,decrease,thres,prededgefile,netpath,kegg)
    metadata = {'description':desc,'tags':['pathlinker-paper']}

    prednodes = set([t for t,h in prededges]).union(set([h for t,h in prededges]))
    if netpath == None:
        netpathnodes = set()
        netpathedges = set()
    else:
         tpedges = readColumns('%s/%s-edges.txt' % (EVIDENCEFILES['NetPath'],netpath),1,2)
         netpathedges = [tuple(sorted((UNIPROT2NAME[u],UNIPROT2NAME[v]))) for u,v in tpedges if u in UNIPROT2NAME and v in UNIPROT2NAME]
         netpathnodes = set([u for u,v in netpathedges]).union(set([v for u,v in netpathedges]))

    if kegg == None:
        keggnodes = set()
        keggedges = set()
    else:
        keggids = {'BCR':'hsa04662',
                   'EGFR1':'hsa04012',
                   'Hedgehog':'hsa04340',
                   'TCR':'hsa04660',
                   'TGF_beta_Receptor':'hsa04350',
                   'Wnt':'hsa04310'
               }
        tpedges = readColumns('%s/%s-edges.txt' % (EVIDENCEFILES['KEGG'],keggids[kegg]),1,2)
        keggedges = [tuple(sorted((UNIPROT2NAME[u],UNIPROT2NAME[v]))) for u,v in tpedges if u in UNIPROT2NAME and v in UNIPROT2NAME]
        keggnodes = set([u for u,v in keggedges]).union(set([v for u,v in keggedges]))

    
    ## re-evaluate ligand edges
    # make sure that this goes from a ligand to a receptor in the pathway 
    # (not a random receptor):
    ligandedges = set([(u,v) for u,v in ligandedges if (u in netpathnodes and v in netpathnodes.intersection(prednodes)) or\
                       (v in netpathnodes and u in netpathnodes.intersection(prednodes))])
    ligandnodes = set([u for u,v in ligandedges]).union(set([v for u,v in ligandedges]))
     
    addednodes = set()
    for n in prednodes:
        if nolabels:
            name = ''
        elif increase or decrease:
            val = min([prededges[(u,v)] for (u,v) in prededges if n==u or n==v])
            if int(val)==val:
                name = '%s\n%d' % (n,val)
            else:
                name = '%s\n%.2e' % (n,val)
        else:
            name = n

        #determine node shape by receptors/tfs
        if n in receptors:
            nodeshape = NODESHAPES['source']
        elif n in tfs:
            nodeshape = NODESHAPES['target']
        elif n in ligandnodes:
            nodeshape = NODESHAPES['ligand']
        else: #internal node
            nodeshape = NODESHAPES['inner']

        #crosstalk nodes
        crosstalknodes = ['SMAD','NOTCH','MAPK','PIK3','EGFR','SRC']

        # determine node color:
        if nolabels:
            if n in receptors and n in netpathnodes:
                htmlcolor = COLORS['source']
            elif n in tfs and n in netpathnodes:
                htmlcolor = COLORS['target']
            else:
                htmlcolor = COLORS['white'] # white nodes for nolabels
        else:
            if n in netpathnodes:
                htmlcolor = COLORS['netpath']
            elif n in keggnodes:
                htmlcolor = COLORS['kegg']
            elif any([cn in n for cn in crosstalknodes]):
                htmlcolor = COLORS['crosstalk']
            else:
                htmlcolor = COLORS['neither']
        
        edgeswithnode = set([(t,h) for t,h in prededges if t==n or h==n])
        pathswithnode = set([prededges[e] for e in edgeswithnode])
        annotation = getNodeAnnotation(n,pathswithnode)

        ## add node
        if nolabels:
            graphspace.add_node(G,n,label=name,color=htmlcolor,shape=nodeshape,popup=annotation,k=int(min(pathswithnode)),width=20,height=20)
        else:
            graphspace.add_node(G,n,label=name,color=htmlcolor,shape=nodeshape,popup=annotation,k=int(min(pathswithnode)),width=45,height=45,textoutlinecolor=htmlcolor)
        addednodes.add(n)

    # add ligand nodes if necessary
    # these are white hexagons.
    for n in ligandnodes:
        if n not in addednodes:
            graphspace.add_node(G,n,label=n,color=COLORS['white'],shape=NODESHAPES['ligand'],popup='Ligand (not in paths)',k=0,textoutlinecolor=htmlcolor)
            print G.get_node_data(n)
            addednodes.add(n)
            print 'Adding ligand node %s..' % (n)

    # Add edges to graph
    seen = set()
    for (tail,head) in prededges:       
        if (head,tail) in seen:
            continue
        seen.add((tail,head))

        # get edge direction
        if (head,tail) in PPIEDGES:
            edgedir = False
        else:
            edgedir = True

        # determine edge color:
        edgewidth = 3
        htmlcolor = COLORS['neither']            
        if tuple(sorted((tail,head))) in netpathedges:
            htmlcolor = COLORS['netpath']
            #edgewidth=3
        elif tuple(sorted((tail,head))) in keggedges:
            htmlcolor = COLORS['kegg']
            #edgewidth=3

        annotation = getEdgeAnnotation(tail,head, prededges[(tail,head)],evidence,undirected)
        graphspace.add_edge(G,tail,head,color=htmlcolor,directed=edgedir,width=edgewidth,popup=annotation,k=int(prededges[(tail,head)]))

    # add ligand edges if necessary
    # these are thin dashed black lines
    for (tail,head) in ligandedges:
        if (tail,head) in seen or (head,tail) in seen: # already written; skip
            continue
        seen.add((tail,head))
        if (head,tail) in PPIEDGES:
            edgedir = False
        else:
            edgedir = True
        graphspace.add_edge(G,tail,head,color=COLORS['darkgray'],directed=edgedir,width=0.25,popup='Ligand-receptor edge (not in paths)',k=0)
    return G,metadata

##########################################################
def getGraphDescription(increase,decrease,thres,predfile,netpath,kegg):
    if increase:
        desc = 'Edges ranked less than or equal to %.4f' % (thres)
    elif decrease:
        desc = 'Edges ranked greater than or equal to %.4f' % (thres)
    else:
        desc = 'Edge set (no ranking)'
    desc += '<hr />Edge file: %s<br>' % (predfile)
    if netpath != None and kegg != None:
        desc += 'Coloring netpath %s nodes/edges %s; coloring other nodes/edges %s if they are in KEGG %s' % (netpath,COLORS['netpath'],COLORS['kegg'],kegg)
    elif netpath != None:
        desc += 'Coloring netpath %s nodes/edges %s' % (netpath,COLORS['netpath'])
    elif kegg != None:
        desc += 'Coloring kegg %s nodes/edges %s' % (kegg,COLORS['kegg'])

    desc += '<hr /><b>Sources of Evidence</b><br> Sources are black if they were used to construct the interactome and <FONT COLOR="gray">gray</FONT> if they were not used to construct the interactome.' 
    desc+='<ul>'

    desc+='<li><FONT COLOR="black"><b>CSBDB</b> are protein-protein interactions from the CSBDB database</FONT></li>'
    desc+='<li><FONT COLOR="black"><b>KEGG</b> are protein interactions from the KEGG database</FONT></li>'
    desc+='<li><FONT COLOR="black"><b>NetPath</b> are protein interactions from the NetPath database</FONT></li>'
    desc+='<li><FONT COLOR="black"><b>SPIKE</b> are protein interactions from the SPIKE database</FONT></li>'
    #desc+='<li><FONT COLOR="black"><b>Vinayagam et. al., 2011</b> is a publication that assigns confidence scores to the edge directions</FONT></li>'
    desc+='</ul>'
    return desc

##########################################################
def getNodeAnnotation(name, pathset):
    pathlist = [p for p in pathset]

    #fullname = csbdb.get_description(NAME2UNIPROT[name])
    #annotation = '<b>Recommended name:</b>%s' % (fullname)
    annotation = '<b>%s</b><hr />' % (name)
    if len(pathlist)>0:
        annotation += '<b>Edge Rankings</b>: %s' % (','.join(str(i) for i in sorted(pathlist)))

    #List Uniprot accession number
    uid = NAME2UNIPROT[name]
    uniproturl = 'http://www.uniprot.org/uniprot/%s' % (uid)
    annotation += '<hr /><b>Uniprot ID:</b> <a style="color:blue" href="%s" target="UniProtKB">%s</a>' % (uniproturl,uid)

    #List EntrezGene IDs:
    for e in [i for i in ENTREZ2NAME if ENTREZ2NAME[i] == name]:
        entrezurl = 'http://www.ncbi.nlm.nih.gov/gene/%s' % (e)
        annotation += '<br><b>Entrez ID:</b> <a style="color:blue" href="%s" target="EntrezGene">%s</a>' % (entrezurl,e)

    #List IHOP link
    ihopurl = 'http://www.ihop-net.org/UniPub/iHOP/?search=%s&field=UNIPROT__AC&ncbi_tax_id=9606' % (uid)
    annotation += '<br><b>Search in </b> <a style="color:blue" href="%s" target="IHOP">IHOP</a>' % (ihopurl)

    return annotation

##########################################################
def getEdgeAnnotation(tail, head, pathset,evidence,undirected):
    if undirected and (tail,head) not in PPIWEIGHTS:
        annotation='Weight: %.4f' % (PPIWEIGHTS[(head,tail)])
    else:
        annotation='Weight: %.4f' % (PPIWEIGHTS[(tail,head)])
    annotation+='<hr /><b>Edge Ranking</b>: %.4f' % (pathset)

    annotation+='<hr /><h><b>Sources of Evidence</b></h>'
    evidencelist = evidence[tuple(sorted((tail,head)))]
    annotation+='<dl>'

    empty = True
    for e in sorted(evidencelist):
        if len(evidencelist[e]) == 0:
            continue
        if e == 'vinayagam':
            vinayagamurl = 'http://stke.sciencemag.org/cgi/content/abstract/4/189/rs8'
            vinayagamname = '<a style="color:blue" href="%s" target="Vinayagam">Vinayagam et. al., 2011</a>' % (vinayagamurl)
            annotation+='<dt>%s</dt>' % (vinayagamname)
        else:
            annotation+='<dt>%s</dt>' % (e)

        annotation+='<ul>'
        for desc in evidencelist[e]:
            if desc in annotation:
                continue
            filtered = False
            for m in METHODS_TO_IGNORE: # CSBDB
                if m in desc:
                    filtered = True
                    break
            # new June 18: keep indirect-effect.
            #if 'indirect-effect' in desc: # KEGG
            #    filtered = True
            if filtered:
                annotation+='<li><FONT COLOR="gray">%s</FONT></li>' % (desc)
            else:
                empty = False
                annotation+='<li><FONT COLOR="black">%s</FONT></li>' % (desc)
        annotation+='</ul>'
    annotation+='<dl>'
    if empty:
            sys.exit('ERROR: Unsorted edge (%s,%s) only has ignored sources of evidence.' % (e[0],e[1]))
    return annotation

##########################################################
def getEvidence(edges):

    evidence = {} 
    sortededges = set([tuple(sorted((t,h))) for t,h in edges])
    for e in sortededges:
        evidence[e] = {'CSBDB':[],'vinayagam':[],'SPIKE':[],'NetPath':[],'KEGG':[]}
    
    keggnamefile = '%s/interactions/kegg/2015-03-23/hsa/HSA_PATHWAY_LIST_FORMATTED.txt' % (DATADIR)
    keggnames = {l[0]:l[1].replace('_',' ') for l in readColumns(keggnamefile,1,2)}

    # CSBDB
    if 'CSBDB' in EVIDENCEFILES:
        dbfile = open(EVIDENCEFILES['CSBDB'],'r')
        for line in dbfile.readlines():
            row = line.strip().split('\t')
            if row[2] not in UNIPROT2NAME or row[3] not in UNIPROT2NAME:
                continue
            sortede = tuple(sorted((UNIPROT2NAME[row[2]],UNIPROT2NAME[row[3]])))
            if sortede not in sortededges:
                continue

            # interactiontype
            interactiontype = ''
            matchObj = re.search('(MI:\d+).*(\(.*\))',row[5])
            if matchObj:
                interactiontype = '%s %s' % (matchObj.group(1),matchObj.group(2))
            else:
                matchObj = re.search('(MI:\d+)',row[5])
                if matchObj:
                    interactiontype = matchObj.group(1)

            # detection method
            detectionmethod = ''
            matchObj = re.search('(MI:\d+).*(\(.*\))',row[6])
            if matchObj:
                detectionmethod = '%s %s' % (matchObj.group(1),matchObj.group(2))
            else:
                matchObj = re.search('(MI:\d+)',row[5])
                if matchObj:
                    detectionmethod = matchObj.group(1)

            # publication ids
            pubids = parseCSBDBpubs(row[7])

            description = 'from %s: %s %s %s' % (row[8],interactiontype,detectionmethod,pubids)
            evidence[sortede]['CSBDB'].append(description)
        dbfile.close()

    # Vinayagam
    if 'vinayagam' in EVIDENCEFILES:
        lines = readColumns(EVIDENCEFILES['vinayagam'],2,4,5)
        for t,h,s in lines:
           if t not in entrez2uniprot or h not in entrez2uniprot:
               continue
           if entrez2uniprot[t] not in UNIPROT2NAME or entrez2uniprot[h] not in UNIPROT2NAME:
               continue
           t = UNIPROT2NAME[entrez2uniprot[t]]
           h = UNIPROT2NAME[entrez2uniprot[h]]
           sortede = tuple(sorted((t,h)))
           if sortede not in sortededges:
               continue
           description = '%s -> %s Score = %s' % (t,h,s)
           evidence[sortede]['vinayagam'].append(description)

    # NetPath
    if 'NetPath' in EVIDENCEFILES:
        pathwayfiles = glob.glob('%s*-edges.txt' % (EVIDENCEFILES['NetPath']))
        for f in pathwayfiles:
            matchObj = re.search('%s(.*)-edges.txt'  % (EVIDENCEFILES['NetPath']),f)
            if matchObj:
                pathwayname = matchObj.group(1)
            else:
                sys.exit('ERROR! Cannot get pathway name from %s' % f)
            lines = readColumns(f,1,2)
            for t,h in lines:

                if t not in UNIPROT2NAME or h not in UNIPROT2NAME:
                    continue

                t = UNIPROT2NAME[t]
                h = UNIPROT2NAME[h]
                sortede = tuple(sorted((t,h)))
                if sortede not in sortededges:
                    continue
                description = '%s' % (pathwayname)
                evidence[sortede]['NetPath'].append(description)
            
    # SPIKE
    if 'SPIKE' in EVIDENCEFILES:
        pathwayfiles = glob.glob('%s*-edges.txt' % (EVIDENCEFILES['SPIKE']))
        for f in pathwayfiles:
            matchObj = re.search('%s(.*)_\d+_\d+.spike-edges.txt' % (EVIDENCEFILES['SPIKE']),f)
            if matchObj:
                pathwayname = matchObj.group(1)
                pathwayname = pathwayname.replace('_',' ')
            else:
                matchObj = re.search('%s(.*).spike-edges.txt' % (EVIDENCEFILES['SPIKE']),f)
                if matchObj:
                    pathwayname = matchObj.group(1)
                    pathwayname = pathwayname.replace('_',' ')
                else:
                    matchObj = re.search('%s(.*)-edges.txt' % (EVIDENCEFILES['SPIKE']),f)
                    if matchObj:
                        pathwayname = matchObj.group(1)
                        pathwayname = pathwayname.replace('_',' ')
                    else:
                        sys.exit('ERROR! Cannot get pathway name from %s' % f)
            lines = readColumns(f,1,2,5,6)
            for t,h,integrity,pubmedids in lines:
                if t not in UNIPROT2NAME or h not in UNIPROT2NAME:
                    continue
                t = UNIPROT2NAME[t]
                h = UNIPROT2NAME[h]
                sortede = tuple(sorted((t,h)))
                if sortede not in sortededges:
                    continue
                if integrity == 'Complex':
                    description = '%s (Complex)' % (pathwayname)
                else:
                    description = '%s (%s) ' % (pathwayname,integrity)
                    for p in pubmedids.split(','):
                        description += parseCSBDBpubs('pubmed:%s' % p)+' '
                evidence[sortede]['SPIKE'].append(description)

    # KEGG
    if 'KEGG' in EVIDENCEFILES:
        pathwayfiles = glob.glob('%s*-edges.txt' % (EVIDENCEFILES['KEGG']))
        for f in pathwayfiles:
            matchObj = re.search('%s(.*)-edges.txt' % (EVIDENCEFILES['KEGG']),f)
            if matchObj:
                pathwayname = matchObj.group(1)
                pathwayname = keggnames.get(pathwayname,pathwayname)
            else:
                sys.exit('ERROR! Cannot get pathway name from %s' % f)
            lines = readColumns(f,1,2,5,6)
            for t,h,edgetype,edgesubtype in lines:
                if t not in UNIPROT2NAME or h not in UNIPROT2NAME:
                    continue
                t = UNIPROT2NAME[t]
                h = UNIPROT2NAME[h]
                sortede = tuple(sorted((t,h)))
                if sortede not in sortededges:
                    continue

                description = '%s (%s,%s)' % (pathwayname,edgetype,edgesubtype)
                evidence[sortede]['KEGG'].append(description)

    ## check edges: all should have some evidence!
    for e in sortededges:
        empty = True
        for key in evidence[e].keys():
            if len(evidence[e][key]) != 0:
                empty = False
                break
        if empty:
            for key in evidence[e].keys():
                print key,evidence[e][key]
            sys.exit('ERROR! sorted edge (%s,%s) (%s,%s) does not have any evidence sources.' % (e[0],e[1],NAME2UNIPROT[e[0]],NAME2UNIPROT[e[1]]))

    return evidence
            

##########################################################
def parseCSBDBpubs(line):
    row = line.split(':')
    if row[0] == 'pubmed':
        pubmedurl = 'http://www.ncbi.nlm.nih.gov/pubmed/%s' % (row[1])
        desc = '<a style="color:blue" href="%s" target="PubMed">pmid:%s</a>' % (pubmedurl,row[1])
    elif row[0] == 'doi':
        doiurl = 'http://dx.doi.org/%s' % (row[1])
        desc = '<a style="color:blue" href="%s" target="DOI">doi:%s</a>' % (doiurl,row[1])
    elif row[0] == 'omim':
        omimurl = 'http://omim.org/entry/%d' % (int(row[1]))
        desc = '<a style="color:blue" href="%s" target="OMIM">doi:%s</a>' % (omimurl,row[1])
    elif row[0] == 'imex':
        imexurl = 'http://www.ebi.ac.uk/intact/imex/main.xhtml?query=%s' % (row[1])
        desc = '<a style="color:blue" href="%s" target="IMEX">imex:%s</a>' % (imexurl,row[1])
    else:
        desc = line

    return desc  

##########################################################
def main(args):
    global MAPFILE, NAME2UNIPROT, UNIPROT2NAME, NAME2ENTREZ, ENTREZ2NAME 
    global PPI, PPIEDGES, PPIWEIGHTS
    global EVIDENCEFILES, DATADIR

    # Parse Arguments
    parser = OptionParser(usage='python post_predictions_to_graphspace.py [options]\n')
    parser.add_option('','--ppi',type='str',metavar='STR',\
                      help='PPI network.')
    parser.add_option('','--version',type='str',metavar='STR',\
                      help='Version of interactome, which pulls from different directories.  Options are "%s"' % ('","'.join(VERSIONS.keys())))  
    parser.add_option('','--datadir',type='str',metavar='STR',\
                      help='Data directory from SVN. Required.')
    parser.add_option('','--infile',type='str',metavar='STR',\
                      help='File of edges.  If --increase (resp. --decrease) is specified, looks at a 3rd column and takes all that are <= (resp. >=) thres.') 
    parser.add_option('','--outdir',type='str',metavar='STR',\
                      help='Output directory to write JSON objects to.')
    parser.add_option('','--undirected',action='store_true',default=False,\
                      help='If specified, edges are assumed to be already sorted. Used when reading edge weights.')
    parser.add_option('','--addfzd',action='store_true',default=False,\
                      help='Add FZD4 and FZD6 as receptors (for wntforexperiments)')
    parser.add_option('','--increase',action='store_true',default=False,\
                      help='If specified, input file has a 3rd column that contains the ranking. Edges <= thres are taken.')
    parser.add_option('','--decrease',action='store_true',default=False,\
                      help='If specified, input file has a 3rd column that contains the ranking. Edges <= thres are taken.')
    parser.add_option('','--thres',type='float',default='200',metavar='float',\
                      help='Threshold for ranked edges. Default is 200.')
    parser.add_option('','--netpath',type='str',metavar='STR',\
                      help='Color specified netpath pathway, e.g., --netpath Wnt')
    parser.add_option('','--kegg',type='str',metavar='STR',\
                      help='Color specified kegg pathway, e.g., --kegg Wnt')
    parser.add_option('','--gsid',type='str',metavar='STR',\
                      help='GraphSpace Graph ID')
    parser.add_option('','--nolabels',action='store_true',default=False,\
                      help='Do not show labels.')
    parser.add_option('','--ligandfile',type='str',metavar='STR',\
                      help='pass file of nodes to connect to receptors.  First column contains IDs, second column contains names.')
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()
    if opts.version not in VERSIONS:
        sys.exit('ERROR: Version %s is not one of the possible versions. Possible versions are %s. Exiting.\n' % 
                 (opts.version,','.join(VERSIONS.keys())))
    if opts.infile == None:
        sys.exit('ERROR: Need a file of edges to post. Exiting.')    
    if opts.outdir == None:
        sys.exit('ERROR: Need a directory to output json files to. Exiting')
    if opts.datadir == None:
        sys.exit('ERROR: data directory must be specified. Exiting.')
    if opts.ppi == None:
        sys.exit('ERROR: ppi file must be specified. Exiting.')
    if opts.gsid == None:
        sys.exit('ERROR: graph space ID (gsid) is required. Exiting.')
    if opts.ligandfile and opts.nolabels:
        sys.exit('ERROR: can only have ligand file on an annotated graphspace graph. Exiting.')
    
    # (0) Determine map file, PPI, and versions
    DATADIR = opts.datadir
    # specify sources and pre-pend datadirectory.
    EVIDENCEFILES = VERSIONS[opts.version]
    for s in EVIDENCEFILES:
        if '%' in EVIDENCEFILES[s]:
            EVIDENCEFILES[s] = EVIDENCEFILES[s] % (opts.datadir)

    MAPFILE = '%snamespace-mappers/human-gene-map.txt' % (DATADIR)
    NAME2UNIPROT = readDict(MAPFILE,1,6)
    UNIPROT2NAME = readDict(MAPFILE,6,1)
    NAME2ENTREZ = readDict(MAPFILE,1,5)
    ENTREZ2NAME = readDict(MAPFILE,5,1)
    
    PPI = opts.ppi
    lines = readColumns(PPI,1,2,3)
    PPIEDGES = [(UNIPROT2NAME[u],UNIPROT2NAME[v]) for u,v,w in lines if u in UNIPROT2NAME and v in UNIPROT2NAME]
    PPIWEIGHTS = {(UNIPROT2NAME[u],UNIPROT2NAME[v]):float(w) for u,v,w in lines if u in UNIPROT2NAME and v in UNIPROT2NAME}
    
    # Determine Receptors and TFs
    receptorfile = '%s/receptors/uniprot-target-list.txt' % (DATADIR)
    receptors = set([UNIPROT2NAME[n] for n in readItemSet(receptorfile,1) if n in UNIPROT2NAME])
    if opts.addfzd:    ## ADD FZD4 and FZD6!!
        #print 'ADDING FZD4 and FZD6 AS RECEPTORS'
        receptors.add('FZD4')
        receptors.add('FZD6')
    tffile = '%s/transcription-factors/vaquerizas-ravasi/human-tfs.txt' % (DATADIR)
    tfs = readItemSet(tffile,1)

    # Get Predicted Edges
    prededges = getPredictedEdges(opts.infile,opts.increase,opts.decrease,opts.thres)

    if opts.ligandfile:
        print 'Adding ligand information...'
        ligandedges = getLigandInformation(opts.ligandfile,prededges,receptors)
    else:
        ligandedges = set()

    # Construct NetworkX Graph
    G,metadata = constructGraph(receptors,tfs,prededges,opts.increase,opts.decrease,opts.thres,opts.netpath,opts.kegg,opts.infile,opts.undirected,opts.nolabels,ligandedges)

    # Post to GraphSpace\
    print 'graphID is %s' % (opts.gsid)
    outfile = '%s/%s.json' % (opts.outdir,opts.gsid)
    graphspace.postGraph(G,opts.gsid,outfile=outfile,user=USERNAME,password=PASSWORD,metadata=metadata)
    if GROUP != None:
        graphspace.shareGraph(opts.gsid,user=USERNAME,password=PASSWORD,group=GROUP)
    print 'DONE\n'

if __name__=='__main__':
    main(sys.argv)
