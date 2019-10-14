'''
Created on Sep 16, 2019

@author: kyang
'''

def sumTPMs(input = "/Users/kyang/Dropbox/Rotation_2/RNASeq/JSL1/50primers_transcripts",
            kallisto_output="/Users/kyang/Dropbox/Rotation_2/RNASeq/JSL1/Unstim_TPM"):
    #Very simple function, given input text file, traverse the Kallisto output TPM file and find the sum of the TPMs.
    genelist,TPMlist,searchlist = [],[],[]
    with open(kallisto_output) as inA:
        for line in inA:
            if "Average" not in line:
                genelist.append(line.split("\t")[0])
                TPMlist.append(float(line.split("\t")[1][:-1]))
    with open(input) as inA:
        for line in inA:
            searchlist.append(line[:-1])
    count = 0
    total = 0
    for s in searchlist:
        found = False
        for i,g in enumerate(genelist):
            if s == g.split(".")[0]:
                count += 1
                total += TPMlist[i]
                found = True
        if found == False:
            print(str(s))
    print("Final TPM count = "+str(total))
    print("Found "+str(count)+" out of "+str(len(searchlist))+" genes")
    
def uniqueIDs(input="/Users/kyang/Dropbox/Rotation_2/RNASeq/JSL1/Gene_IDs"):
    templist = []
    with open(input) as inA:
        for line in inA:
            if line not in templist:
                templist.append(line)
    for r in templist:
        print(r[:-1])
    print(str(len(templist)))