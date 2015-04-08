from utilsPoirel import *
import sys

MAPPINGFILE = '/data/annaritz/datasets/svn-data/namespace-mappers/human-gene-map.txt'

####
def main(args):

    if len(args) != 2:
        sys.exit('error: filename required')
    infile = args[1]
    outfile = args[1]+'.mapped'

    lines = readColumns(MAPPINGFILE,1,6)
    u2n = {}
    n2u = {}
    for n,u in lines:
        if u not in u2n:
            u2n[u] = set()
        if n not in n2u:
            n2u[n] = set()
        u2n[u].add(n)
        n2u[n].add(u)

    out = open(outfile,'w')
    with open(infile) as fin:
        for line in fin:
            if line[0] == '#':
                out.write(line)
                continue
            row = line.strip().split('\t')
            for i in range(len(row)):
                if row[i] in u2n:
                    row[i] = ','.join(u2n[row[i]])
                elif '|' in row[i]:
                    splitrow = row[i].split('|')
                    for j in range(len(splitrow)):
                        if splitrow[j] in u2n:
                            splitrow[j] = ','.join(u2n[splitrow[j]])
                    row[i] = '|'.join(splitrow)
            out.write('\t'.join(row)+'\n')
    out.close()




if __name__=='__main__':
    main(sys.argv)
