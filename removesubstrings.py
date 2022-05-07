


import multiprocessing
import datetime
import os
import threading
import concurrent.futures

def readf (filename,minOKlength):
    #it could be set to a very -1 to accept all reads
    #mink is the minimum k-mer
    reads=[]
    s=""
    with open(filename,"r") as f:
        for line in f:
            if not line[0]==">":
                s += line.rstrip()
            else:
                if not s=="":
                    reads.append(s)
                    s=""
    if not s=="":
        reads.append(s)
    f.close()
    #reads=list(tempset)
    #del (tempset)
    return reads
def createhash(reads,mink):
    ht={}
    index=0
    for r in reads:
        try:
            ht[r[0:mink]].add(index)
        except:
            ht[r[0:mink]]= set() #[0]
            ht[r[0:mink]].add(index)


        index+=1
    return ht


def findsubstrings(reads,ht,minkmer):
    #global allsubstrings
    #g
    allsubstrings=set()
    nonsupstrings=set()
    readid=0
    readslength = len(reads)
    while  readslength > 0:
        if readid not in allsubstrings:
            #tempreads=set()
            maybesubstringlist=set()
            for kmer in range (0, len(reads[readid])-minkmer+1):
                if (reads[readid][kmer:kmer+minkmer] in ht ):
                    templist=ht[reads[readid][kmer:kmer+minkmer]] #- set(range (readid))
                    for s in templist:
                        if s not in allsubstrings and s > readid:

                            maybesubstringlist.add (s)
    
            for r in maybesubstringlist:
                if (reads[r] in reads[readid]):
                    allsubstrings.add(r)
                    readslength-=1

            nonsupstrings.add(readid)
            readslength-=1

        readid +=1
                    
    return allsubstrings,nonsupstrings


def removered(filename,mink,outfile):

	a = datetime.datetime.now()
	reads=readf(filename,-1)

	print ("total number of reads")
	print (len(reads))


	b = datetime.datetime.now()
	print ("Reading file time")

	print(b-a)


	a = datetime.datetime.now()
	reads.sort(key=len, reverse=True)
	b = datetime.datetime.now()
	print ("Sorting time")
	print(b-a)
	a = datetime.datetime.now()
	ht= createhash(reads,mink)
	b = datetime.datetime.now()
	print ( "finished hash table")
	print(b-a)
	a = datetime.datetime.now()

	allsubstrings,nonsupstrings=findsubstrings(reads,ht,mink)
	#writing the reads back to the file
	f=open(outfile , "w")
	for i in range (0,len(reads)):
	    if i not in allsubstrings:
	        f.write( ">" + reads[i][:10] + "\n")
	        f.write( reads[i] + "\n")
		      
	f.close()
		



	print ("tha substrings length is")
	print (len (allsubstrings))
	b = datetime.datetime.now()
	print ("finding substrings time ")

	print(b-a)


	print ("number of reads withour the substrings")
	print (len (reads) - len (allsubstrings))
	print (len (nonsupstrings)+ len(allsubstrings))


