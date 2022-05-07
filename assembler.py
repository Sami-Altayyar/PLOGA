# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 11:42:49 2022

@author: sa
"""


# from pyvis.network import Network
# import multiprocessing
# import datetime
# import os
# import threading
# import concurrent.futures

def readf(filename):
    # it could be set to a very -1 to accept all reads
    # mink is the minimum k-mer
    tempset = set()
    s = ""
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":
                    tempset.add(s)
                    s = ""
    if not s == "":
        tempset.add(s)
    f.close()
    reads = list(tempset)
    return reads


def creatdic(reads, mink):
    prefixhash = {}
    r = 0
    for read in reads:
        if (read[0:mink] in prefixhash):
            prefixhash[read[0:mink]].append(r)
        else:
            prefixhash[read[0:mink]] = [r]
        r += 1
    return prefixhash


def bulidgraph(reads, h, mink):
    graph = dict()
    r = 0
    starts = set()
    notstarts = set()
    for read in reads:
        graph[r] = []
        for x in range(1, len(read)):
            if read[x:x + mink] in h:
                l = h[read[x:x + mink]]
                for rd in l:
                    if read[x:] == reads[rd][:len(read[x:])]:
                        graph[r].append((rd, x))
                        notstarts.add(rd)
                        if rd in starts:
                            starts.remove(rd)
                        if not r in notstarts:
                            starts.add(r)
        #        if len (graph[r]) == 0 :
        #            #The read is not connected with any other read
        #            notconnected.append (r)
        #            del graph[r]

        r += 1
    return graph, starts


def findpaths(reads, graph, starts):
    paths = []
    nods = []
    used = set()
    for a in starts:
        temppaths = []
        nods.append([a, "", 0])
        while nods:
            s = nods.pop()
            flag = False
            if not (s[0] in used):
                # net.add_node(s, label=str(s))
                used.add(s[0])
                for n, m in graph[s[0]]:
                    if not n in used:
                        score = (s[2] + (m / len(reads[s[0]]))) / 2
                        nods.append([n, s[1] + reads[s[0]][:m], score])
                        flag = True
            if not flag:
                temppaths.append((s[1] + reads[s[0]][:], int(round(s[2], 2) * 100)))

        temppaths.sort(key=lambda x: len(x[0]))
        for i in range(len(temppaths)):
            flag = False
            for j in range(i, len(temppaths)):
                if temppaths[i][0] in temppaths[j][0] and i != j:
                    flag = True
                    break
            if not flag:
                paths.append(temppaths[i])
    return paths


def PLOGAmain(filename, mink, outputfilename):
    reads = readf(filename)
    h = creatdic(reads, mink)
    g, starts = bulidgraph(reads, h, mink)
    x = findpaths(reads, g, starts)

    of = open(outputfilename, "a")
    for l, s in x:
        of.write(">" + l[:mink] + "_" + str(s) + "%\n")
        of.write(l + "\n")
    of.close()
    del reads
    del h
    del g
    del starts
    del x
    del of

