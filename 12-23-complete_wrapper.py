def buildConcensus(kmers,looklength):
    ## takes a list (for each keyed anchor), and a sequence of bases up to looklength ; computes base composition of the next seq of bases and ouptuts
    # total num per base, and consensus fraction and base
    baseCount = [] # number at that position
    baseComp = '' # which is the most frequent base
    baseFrac = [] # what is its frequency
    print("CALLED")
    for j in range(1,(looklength-1)): #loop through each base ahead
        print(j)
        mycounter=[] #initialize dummy var
        ############################### CHECK 
        for eachseq in kmers: # loop through each sequence; 
            if len(eachseq)>(j-1): # test if the sequence is bigger than jth
                mycounter.append(eachseq[j:j+1])
            # record the ith base

        print(mycounter)
        
        vA = mycounter.count("A")
        vT = mycounter.count("T")
        vC = mycounter.count("C")
        vG = mycounter.count("G")
        n = vA + vT + vC + vG
        print(n)
        # add n to counter position j
        if n>0 :
            baseCount.append(n)
            if (vA == max(vA,vT,vC,vG)): # assign vT add a T to the string
                baseComp+='A'#.append("A")
                baseFrac.append(vA/n)
            if (vC == max(vA,vT,vC,vG)): # assign vT add a T to the string
                baseComp = baseComp + 'C' #.append("C")
                #baseComp.append("C")
                baseFrac.append(vC/n)
            if (vG == max(vA,vT,vC, vG)): # assign vT add a T to the string
                baseComp+='G'
                #baseComp.append("G")
                baseFrac.append(vG/n)
            if (vT == max(vA,vT,vC, vG)): # assign vT add a T to the string
                baseComp+='T'
                #baseComp.append("T")
                baseFrac.append(vT/n)
            
    return baseComp,baseFrac,baseCount

def recordNextKmers(anchorlist,looklength,myseqs,anchorlength,DNAdict):   #anchorlist is a list -- we will loopthrough the sequence myseq and check if any of the anchorlist kmers are defined
    # ,anchorlength is length of kmers in file
    import re
    ## LOOK AHEAD IN THE STRING
    
    for myseq in myseqs:
## loop through each kmer in the read; if it is part of the old dictionary, which gets passed in, record all of the mkers 
        for i in range(1,len(myseq)):
            # get substr at the ith position of length anchorlength to compare to anchorlist
            mystring= myseq[i:i+anchorlength]
             # check if it exists
            if mystring in anchorlist:
                ## find where the match is
                match = myseq.find(mystring)

                s = match
                e = match + len(mystring)
 #               print ("match endss")
 #               print(e)
        ## test if the stringlength is long enough to get the string 
                go=min(len(myseq), (e + looklength ))
                nextkmer = myseq [ e : e+looklength  ]
#               print("nextkmer")
#                    print(nextkmer)
                    ## keep dict small so less than 100 seqs:
                if len(DNAdict[mystring]) < 100 :
                    DNAdict[mystring].append(nextkmer)
    return DNAdict

################# GETTING REAL SEQS
def returnSeqs(fqfile,maxlines):
    import gzip
    print("Counting total number of reads...")
    myseqs=[]
    # Count reads                                                                                                                   
    tot_lines =  0
    with gzip.open(fqfile, "rt") as handle:# parse reads
        for read_seq in handle:
            # check we're in sequence line (remainder of 2)                                                                         
            tot_lines += 1
            if tot_lines%4!=2:
                continue
            # strip of new line character                                                                                               
            read_seq = read_seq.strip('\n')
            if len(myseqs)< maxlines:
                myseqs.append(read_seq)
    return(myseqs)

###############
################# GETTING REAL SEQS
def returnAnchors(infile):
    import gzip
#    import pandas
    print("Counting total number of reads...")
    anchors=[]
    # Count reads                                                                                                                       
    tot_lines =  0

    with gzip.open(infile, "rt") as handle:

        # parse read                                                                                                                    
        for line in handle:

            # check we're in sequence line (remainder of 2)                                                                             
            tot_lines += 1
            if tot_lines>1 :
  #              continue
            # strip of new line character                                                                                               
                qval = line.strip().split("\t")[5]
                seq = line.strip().split("\t")[6]
               

                if float(qval) < .001 : # lots of clusters 
#                if float(qval) < .001 and float(line.strip().split("\t")[2])>10: # lots of clusters 
                    #print("passed q val")
                    anchors.append(seq)
    return list(set(anchors))

## start running with fq file input and anchorInfile

## takes 2 files fqfile and infile and outputs 

## 3 files: concensus as a fasta file
## splits the fastq to give the prefix
## seq_id anchor \t consensusseq \t string of fraction values
## seq_id anchor \t consensusseq \t string of fraction values
 

#fqfile="/oak/stanford/groups/horence/julias/srr_outs/ERR2516403_2.fastq.gz"

looklength=200



#anchorlist=['TAC','GAG','CTC','TTA']

#anchorInfile = "/scratch/groups/horence/jordi/metagenomics/sarscov2/dgmfinder_v5/weinberger_long_8M_anchors_annot.txt.gz"
#anchorInfile = "/scratch/groups/horence/jordi/metagenomics/sarscov2/dgmfinder_v5/weinberger_s1_1M_anchors_annot.txt.gz"
#fqfile= "/scratch/groups/horence/jordi/metagenomics/sarscov2/fastq/weinberger_long_8M.fastq.gz"
#fqfile= "/scratch/groups/horence/jordi/metagenomics/sarscov2/fastq/Weinberger-SC-3397-01_S1_R1.fastq.gz"

def WriteSeqInfoFiles(anchorInfile, fqfile):

    # string processing                                                                                                                                                      
    fqnext=fqfile.split("fastq/")[1]
    prefix=fqnext.split(".fastq")[0]
    print(prefix)
    ## print outs                                                                                                                                                            
    fasta = prefix + "_concensus.fa"
    totalcountfile = prefix + "_counts.fa"
    fractions = prefix + "_fractions.tab"
    myseqs=returnSeqs(fqfile, maxlines=1000000)                                                                                                   
    anchorlist=returnAnchors(anchorInfile)                                                                                                                    

    #myfracs=['.2','.22']
    #myseqs = ['GATTACAGATTGAGATTACAGGATTACGATTACCAGATTACTAGATTACAAGATTACGATTACAATTACAATTACAGATACAGATTACA','GATTACAGATGGGTGAGATTACAGGATTACGATTACCAGATTACTAGATTAACTAACGATTACGATTACAACCCTTACAATTACAGATACAGATTACA']
    ################### RUN AND WRAP other work:                                                                                                                             
    ## read in the real files                                                                                                                                      #myseqs=returnSeqs(fqfile, maxlines=1000000)                                                                                                                   #anchorlist=returnAnchors(anchorInfile)                                                                                                                                  
    #anchorlist=['TA','AG']
    ## DNA dictionary stores the set of reads after each anchor in the angorlist                                                                                             
    DNAdict={}
    #initialize                                                                                                                                                              
    anchorlength =  len(anchorlist[0])
    for an in anchorlist:
        DNAdict[an]=[]
    ###### get all of the next kmers for the anchors!    in PREPARATION FOR BUILDING CONCENSUS                                                                               
    nextseqs = recordNextKmers(anchorlist,looklength,myseqs,anchorlength,DNAdict) # should be a list                                                                         
    for kk in nextseqs.keys():
    ## gets the value as an array?                                                                                                                                           
    # syntax for getting the values of a key                                                                                                                                 
        if len( nextseqs.get(kk) )>0 : # build concensus                                                                                                                 
            out=buildConcensus( nextseqs.get(kk) ,looklength)
                #print("GOT OUT for key ")
                #print (kk)
                #print(out[0]) ######### out0 is the sequence after the kmer                                                                                                  
            if len(out[1])>0:
                print(kk+"--->"+out[0])
               ######################################## write the fasta file                                                                                                 
                with open(fasta, "w") as external_file:
                    external_file.write('>')
                    external_file.write(kk)
                    external_file.write('\n')
                    external_file.write(out[0])
                    external_file.write('\n')
                    external_file.close()
                print("out[1]")
                print(out[1])
                with open(fractions, "w") as external_file:
                    external_file.write(kk)
                    external_file.write('\t')
                    external_file.write(out[0])
                    external_file.write('\t')
#                    external_file.write('\t'.join(out[1]) + '\n')
                    ns=str(out[1])
                    
                    external_file.write('\t'+ns + '\n')
                    external_file.close()
                print("out[2]")
                print(out[2])
                with open(totalcountfile, "w") as external_file:
                    external_file.write(kk)
                    external_file.write('\t')
                    external_file.write(out[0])
                    external_file.write('\t')

                    external_file.write('\t'+str(out[2]) + '\n')

                    external_file.close()

 


anchorInfile = "/scratch/groups/horence/jordi/metagenomics/e_coli/dgmfinder_v5/SRR8660932_pass_2_anchors_annot.txt.gz"
fqfile= "/scratch/groups/horence/jordi/metagenomics/e_coli/fastq/SRR8660932_pass_2.fastq.gz"
WriteSeqInfoFiles(anchorInfile, fqfile)
