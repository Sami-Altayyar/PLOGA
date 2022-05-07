
import sys

def removered(filename,mink,outfile):
    # it could be set to a very -1 to accept all reads
    # mink is the minimum k-mer
    # reads = []
    hashs = set()
    s = ""
    count = 0
    subs = 0
    fo = open(outfile , "w")
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":
                    temphash = hash(str(hash(s)) + str(len(s)) + str(s[:mink]) )
                    if not temphash in hashs:
                        hashs.add(temphash)
                        fo.write(">" + s[:mink] + "\n")
                        fo.write(s + "\n")
                        count = count + 1
                    else:
                        subs += 1
                    # reads.append(s)
                    s = ""
    if not s == "":
        temphash = hash(str(hash(s)) + str(len(s)))
        if not temphash in hashs:
            hashs.add(temphash)
            fo.write(">" + s[:mink] + "\n")
            fo.write(s + "\n")
            count = count + 1
        else:
            subs += 1
        # reads.append(s)
    f.close()
    fo.close()
    # reads=list(tempset)
    # del (tempset)
    return


def readf(filename, mink , outfile):
    reads = []
    s = ""
    ht = {}
    index = 0
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                s += line.rstrip()
            else:
                if not s == "":
                    flag = False
                    if s[0:mink] in ht:
                        for i in ht[s[0:mink]]:
                            if s == reads[i]:
                                flag = True
                                break
                    if not flag:
                        reads.append(s)
                        try:
                            ht[s[0:mink]].add(index)
                        except:
                            ht[s[0:mink]] = set()
                            ht[s[0:mink]].add(index)
                        index += 1
                    s = ""
    if not s == "":

        flag = False
        if s[0:mink] in ht:
            for i in ht[s[0:mink]]:
                if s == reads[i]:
                    flag = True
                    break
            if not flag:
                reads.append(s)

        reads.append(s)
    f.close()
    print ("Reading complete ")
    print("Number of reads = " + str (len (reads)))
    print ("Start writing none redundant reads output")
    f = open(outfile , "w")
    for i in range(0, len(reads)):
        f.write(">" + reads[i][:mink] + "\n")
        f.write(reads[i] + "\n")

    f.close()

    print("Writing  complete ")
    return


