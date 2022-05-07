# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
#import subprocess
import sys, getopt
import os
import datetime
import removesubstrings,assembler,cluster_reads

def main(argv):
   inputfile = ''
   outputfile = ''


   try:

      opts, args = getopt.getopt(argv,"hi:o:",[["ifile="],["ofile="]])
   except getopt.GetoptError:
      print ('main.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('main.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile =  arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print ('Input file is "', inputfile)
   print ('Output file is "', outputfile)

   print (opts)
   arguments = len(sys.argv) - 1
   print ("The script is called with %i arguments" % (arguments))


def readarg(argv):
   inputfile = []
   outputfile = ''
   numthreads=4
   mink=6
   train= "complete"
   # arguments = len(sys.argv) - 1
   flag = -1
   for r in sys.argv:

      if not r in ["-i","-o","-k","-p", "-t"]:
         if flag==1:
            inputfile.append(r)

         if flag==0:
            outputfile = r

         if flag == 3:
            numthreads  = int (r)

         if flag == 4:
            mink =  int (r)

         if flag == 5:
            train =  r

      if r == "-o": # Output Folder
         flag = 0
      if r == "-i": # input files
         flag=1
      if r == "-p": # Threads
         flag=3
      if r == "-k": # K-mer
         flag=4
      if r == "-t": # Train type
         flag=5


   return inputfile,outputfile,numthreads,mink,train


def print_hi(name):

    print(f'Hi, {name}')

def remove_red(filename1,outfile1):
# Remove from fastq file paire end
    #hashs = set()
    count = 0
    subs = 0
    f = open(filename1, "r")
    fo = open(outfile1 , "w")
    hashsarr=[set(),set(),set(),set(),set(),set(),set(),set(),set(),set()]
    #fremoved = open(removedfile , "w")
    #t={84: 65, 65: 84, 71: 67, 67: 71}
    while True:
        h1= f.readline()
        s1 = f.readline()
        if h1:
            temphash = hash (s1 )
            ind=temphash % 10
            if temphash in hashsarr[ind]: # hashs:
                subs += 1
            else:
                hashsarr[ind].add (temphash) #hashs.add(temphash)
                fo.write(h1 + s1  )
                count = count + 1
        else:
            break


    f.close()
    fo.close()
    #fremoved.close()
    return count,subs

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi(' Welcom To Protin Level Overlap Graph Assembler (PLOGA)')
    inputfiles , outputfolder,threadnum, mink ,train = readarg(sys.argv[1:])

    print ("Input Files " )
    print (inputfiles)
    print("output Folder ")
    print(outputfolder)
    print("Number of threads = " + str(threadnum))
    print ("minimum Kmer = " + str (mink))
    cmd = "mkdir " + outputfolder
    os.system(cmd)
    print ("Start remove redundant reads")
    a = datetime.datetime.now()
    remove_red(inputfiles[0],outputfolder + "/nored.fasta")
    b = datetime.datetime.now()
    print("Removing Redundant reads time ")
    print(b - a)
    a = datetime.datetime.now()

    print ("Splitting File")
    cmd = "cat "
    cmd += outputfolder + "/nored.fasta | awk "
    cmd += '\'BEGIN{b=1}{if(/^>/&&b){printf"%s\\n",$1;b=0}else if(/^>/){printf"\\n%s\\n", $1} else {printf"%s",$1}}END{printf"\\n"}\''
    cmd += " | split -l 20000 -d -a 10 - "
    cmd += outputfolder + "/x"



    print(cmd)
    #cmd = "./FragGeneScan1.31/run_FragGeneScan.pl -genome=" + outputfolder + "/nored.fasta -out=" + outputfolder + "/fgs -complete=0  -train=complete -thread=" + str(threadnum)
    os.system(cmd)

    print("Splitting Files Completed")
    b = datetime.datetime.now()
    print("Splitting Files time ")
    print(b - a)

    a = datetime.datetime.now()

    print ("Starting FragGeneScan")
    cmd = "ls " + outputfolder + "/x* | xargs -P 32 -i bash -c './FragGeneScan1.31/run_FragGeneScan.pl -genome={} -out={}.fgs -complete=0 -train=" + train + " -thread= " + str (threadnum) + " &> {}.log'"
    print(cmd)
    os.system(cmd)
    print("FragGeneScan Completed")

    #cmd = "cat " + outputfolder +"/x*.faa > " + outputfolder +"/fgs.faa"
    cmd = "ls " + outputfolder +"/x*.faa | xargs -I % sh -c 'cat %; echo "";' > " + outputfolder +"/fgs.faa"
    print (cmd)
    os.system(cmd)
    
    cmd = "rm " + outputfolder +"/x*"
    os.system(cmd)
    
    b = datetime.datetime.now()
    print("FragGeneScan time ")
    print(b - a)

    a = datetime.datetime.now()
    #removreredundent.readf(outputfolder + "/fgs.faa", mink, outputfolder + "/nored.faa")
    print("Start remove redundant ORFs")
    removesubstrings.removered(outputfolder +"/fgs.faa", mink, outputfolder + "/fgsnored.faa")
    b = datetime.datetime.now()
    print("Removing Redundant ORFs time ")
    print(b - a)

    num_of_files =cluster_reads.cluster (outputfolder + "/fgsnored.faa",mink)
    of = open(outputfolder + "/assembled.fasta", "w")
    of.close()
    for i in range (1,num_of_files+1):

        assembler.PLOGAmain(outputfolder + "/fgsnored.faa"+ str(i) + ".fnn",mink, outputfolder + "/assembled.fasta" )

        cmd = "rm " + outputfolder + "/fgsnored.faa"+ str(i) + ".fnn"
        os.system(cmd)



    print ("is it finished")

    # print (cmd)






