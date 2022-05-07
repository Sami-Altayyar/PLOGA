# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 16:51:43 2021

@author: sa
"""
import os


def readf(filename):
    # it could be set to a very -1 to accept all reads
    # mink is the minimum k-mer
    reads = []
    s = ""
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":
                    reads.append(s)
                    s = ""
    if not s == "":
        reads.append(s)
    f.close()
    # reads=list(tempset)
    # del (tempset)
    return reads


def creatdic(filename, mink):
    # it could be set to a very -1 to accept all reads
    # mink is the minimum k-mer
    readsdic = {}
    reads = []
    s = ""
    count = 0
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":
                    reads.append(s[:mink] + s[len(s) - mink:])
                    prehash = hash(s[:mink])
                    suffhash = hash(s[len(s) - mink:])
                    if (prehash in readsdic):
                        readsdic[prehash][0].add(count)
                    else:
                        readsdic[prehash] = [set([count]), set(), set()]

                    if (suffhash in readsdic):
                        readsdic[suffhash][1].add(count)
                    else:
                        readsdic[suffhash] = [set(), set([count]), set()]
                    s = ""
                    count += 1
    if not s == "":
        reads.append(s)
        prehash = hash(s[:mink])
        suffhash = hash(s[len(s) - mink:])
        if (prehash in readsdic):
            readsdic[prehash][0].add(count)
        else:
            readsdic[prehash] = [set([count]), set(), set()]

        if (suffhash in readsdic):
            readsdic[suffhash][1].add(count)
        else:
            readsdic[suffhash] = [set(), set([count]), set()]
    f.close()
    # reads=list(tempset)
    # del (tempset)
    return readsdic, reads


def filldic(filename, mink, readsdic):
    s = ""
    count = 0
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":

                    for kmer in range(1, len(s) - mink):
                        thehash = hash(s[kmer:kmer + mink])
                        if (thehash in readsdic):
                            readsdic[thehash][2].add(count)
                    s = ""
                    count += 1
    if not s == "":
        for kmer in range(1, len(s) - mink):
            thehash = hash(s[kmer:kmer + mink])
            if (thehash in readsdic):
                readsdic[thehash][2].add(count)


def findreleations(thereadid, reads, mink):
    global usedreads
    global readsdic
    usedreads.add(thereadid)
    relatedreads = set()
    filnalrelation = set()
    relatedreads.add(thereadid)
    filnalrelation.add(thereadid)
    while (len(relatedreads) > 0):
        readid = relatedreads.pop()
        pre = reads[readid][:mink]
        suf = reads[readid][len(reads[readid]) - mink:]
        pre_reads = readsdic[hash(pre)][2] - usedreads
        suf_reads = readsdic[hash(suf)][2] - usedreads

        for i in pre_reads:
            if readid in readsdic[hash(reads[i][len(reads[i]) - mink:])][2]:
                relatedreads.add(i)
                filnalrelation.add(i)
                usedreads.add(i)

        for i in suf_reads:
            if readid in readsdic[hash(reads[i][:mink])][2]:
                relatedreads.add(i)
                filnalrelation.add(i)
                usedreads.add(i)

    return (filnalrelation)


# def findreleations2(thereadid,reads,mink):
#     #No checking 
#     global usedreads
#     global readsdic
#     usedreads.add(thereadid)
#     relatedreads=set()
#     filnalrelation=set()
#     relatedreads.add(thereadid)
#     filnalrelation.add(thereadid)
#     while (len(relatedreads)>0):
#         readid=relatedreads.pop()
#         pre=reads[readid][:mink]
#         suf=reads[readid][len(reads[readid])-mink:]
#         pre_reads=readsdic[hash(pre)][2].union( readsdic[hash(suf)][2]) - usedreads
#         #suf_reads=readsdic[hash(suf)][2] - usedreads

#         for i in pre_reads:
#             if readid in readsdic[hash(reads[i][len(reads[i])-mink:])][2] or readid in readsdic[hash(reads[i][:mink])][2] :
#                 relatedreads.add(i)
#                 filnalrelation.add(i)
#                 usedreads.add(i)


#     return (filnalrelation)
def cluster(mink, filename):

    usedreads = set()
    readsdic, reads = creatdic(filename, mink)
    filldic(filename, mink, readsdic)
    relations = []
    for i in range(0, len(reads) - 1):
        if i not in usedreads:
            s = findreleations(i, reads, mink)
            relations.append(s)

    print(len(relations))
    print(len(reads))

    # deleting all dictionaries
    del usedreads, readsdic

    reads = readf(filename)
    c = 0
    d = 0
    for i in relations:
        if len(i) == 1:
            c += 1
        else:
            tempfile = filename + str(d) + ".fnn"
            f = open(tempfile, "w")
            for r in i:
                f.write(">" + reads[r][:mink] + "\n")
                f.write(reads[r] + "\n")

            d += 1

    print(c)
    print(d)
