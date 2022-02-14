#! /usr/bin/env python3
import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument("-i","--i",help="this is Input file")
parser.add_argument('-f', '--fold', type=int, default = 70)
args = parser.parse_args()
num = args.fold
if args.i:
        Input_sequence = str(args.i)
                  
def Fastq_to_Fasta(fn):
    flag=0
    with open(fn, 'r') as file_handler:
        final_output = []
        
        for line in file_handler:
            
            if (re.search(pattern = '\\+', string = line)is not None):
                
                final={}
            if (re.search(pattern = '\\?', string = line)is None) and(re.search(pattern = '@', string = line)is not None):
                final={}
                final['name']=line.replace("@",">").strip('\n')    
                
                
            elif (re.search(pattern = '\\?', string = line)is None) and(re.search(pattern = '[ATCGn]+', string = line)is not None):
                    final["seq"]= line.strip("\n")
                    if(len(re.findall(r'[^ATCGNatcgn]+', final["seq"]))!=0):
                        flag =1
                    final["len"]= len(line.strip("\n"))
                    final_output.append(final)
    return final_output,flag

def GBtoFasta(fn):
    flag=0
    with open(fn, 'r') as file_handler:
        lines = []
        final = {}
        a =0
        for line in file_handler:
            if a == 1:
                lines.append(line)
            if (line.strip() != ''):
                
                if (re.search(pattern = 'ACCESSION', string = line) is not None):      
                    final["name"]=">"+(line.split())[1].strip('\n')
                if (re.search(pattern = 'LOCUS', string = line) is not None): 
                    
                    final["len"]=line.split()[2].strip('\n')
                if (re.search(pattern = 'ORIGIN', string = line) is not None): 
                    a = 1
                if (re.search(pattern = '//', string = line) is not None): 
                    
                    seq='' 
                    for i in range(1, len(lines)-1):
                        seq += ''.join(lines[i].split()[1:])
                    final["seq"]=seq.upper().strip('\n')
                    if(len(re.findall(r'[^ATCGNatcgn]+', final["seq"]))!=0):
                        flag =1
    return [final],flag
def Mega_to_Fasta(fn):
    flag=0
    with open(fn, 'r') as file_handler:
        lines=[]
        header= []
        seq = []
        final_output =[]
        for line in file_handler:
            lines.append(line)
        content_list = ("".join(lines[2:])).split("#")
        for each in content_list:
            info = each.split("\n")
    
            header.append(info[0])
            seq.append("".join(info[1:]))
    header = header[1:]
    seq = seq[1:] 
    #print(header[0])
    if len(header)==0:
        final={}
        final["name"]=header[0].strip('\n')
        final["seq"]=seq[0].strip('\n')
        final["len"]=len(seq[0]).strip('\n')
        final_output.append(final)
    else:
        
        for i in range(0, len(header)) :
            
                final={}
                final["name"]=header[i].strip('\n')
                final["seq"]=seq[i].strip('\n')
                final["len"]=len(seq[i].strip('\n'))
                if(len(re.findall(r'[^ATCGNatcgn]+', seq[i]))!=0):
                    flag =1
                final_output.append(final)
    return final_output,flag


def vcfref(file):
    vcffile=open(file,"r")
    sample=[]
    result={}
    references = []
    for line in vcffile:
        ref=[]
        alt=[]
        genetype=[]
        if line.startswith("#CHROM"):
            sample=line.strip("\n").split("\t")[9:]
            #print(sample)  
            for i in sample:
                result[i]=""
        elif line.startswith("##"):
            pass
        else:
            genetype=line.split("\t")
            #print(genetype)
            realgenetype=[]
            genetype_dict={}
            ref=line.split("\t")[3].split(",")
            references.append(ref)
    return references

def vcf(file,colnum):
    vcffile=open(file,"r")
    sample=[]
    result={}
    rescon = []
    #refans = ""
    for line in vcffile:
        ref=[]
        alt=[]
        genetype=[]
        if line.startswith("#CHROM"):
            sample=line.strip("\n").split("\t")[9:]
            #print(sample)  
            for i in sample:
                result[i]=""
        elif line.startswith("##"):
            pass
        else:
            genetype=line.split("\t")
            #print(genetype)
            realgenetype=[]
            genetype_dict={}
            ref=line.split("\t")[3].split(",")
            alt=line.split("\t")[4].split(",")
            #print("Ref", ref)
            #print("Alt",alt)
            ninecol = line.split("\t")[colnum][0].split(":")
            #print(ninecol)
            if "1" in ninecol:
                rescon.append(alt)
            if "0" in ninecol:
                rescon.append(ref)
    
    print(">"+sample[colnum-9])
    return rescon

def printList():
    print(">NC_01234")
    refcon = vcfref(Input_sequence)
    #print(refcon)
    nac = ""
    for t in refcon:
        for s in t:
            nac = nac + s
    print(nac)

    for p in range(9,20):
        ans = vcf(Input_sequence,p)
        ohhoo = ""
        for i in ans:
            for j in i:
                ohhoo = ohhoo + j
        print(ohhoo)
printList()


def EMBL_to_Fasta(fn):
    flag=0
    with open(fn, 'r') as file_handler:
        lines = []
        final = {}
        a =0
        for line in file_handler:
            if a == 1:
                lines.append(line)
            if (line.strip() != ''):
                
                
                if (re.search(pattern = 'ID', string = line) is not None):      
                    final["name"]=">"+(line.split())[1].strip('\n')
                if (re.search(pattern = 'DE  ', string = line) is not None):               
                    final["descr"]=(line.strip('\n').replace('DE   ', 'descr='))
                if (re.search(pattern = 'SQ', string = line) is not None): 
                    a = 1
                    final["len"]=line.split()[2].strip('\n')
                    
                if (re.search(pattern = '//', string = line) is not None): 
                    
                    seq='' 
                    for i in range(1, len(lines)-1):
                        seq += ''.join(lines[i].split()[0:-1])
                    final["seq"]=seq.upper().strip('\n')
                    if(len(re.findall(r'[^ATCGNatcgn]+', final["seq"]))!=0):
                        flag =1    
    return [final],flag
    
with open(Input_sequence,"r") as fh1:
    first_line = fh1.readline()
    #print(first_line)
    match = re.search(pattern = '@', string = first_line)
    if (match is not None):
        final_output,flag = Fastq_to_Fasta(Input_sequence)
    match = re.search(pattern = 'ID', string = first_line)
    if (match is not None):
        final_output,flag = EMBL_to_Fasta(Input_sequence)
    match = re.search(pattern = 'LOCUS', string = first_line)
    if (match is not None):
        final_output,flag = GBtoFasta(Input_sequence)
    match = re.search(pattern = '#MEGA', string = first_line)
    if (match is not None):
        final_output,flag= Mega_to_Fasta(Input_sequence)
    
#print(final_output)
filename = Input_sequence.split(".")
if (flag == 0 ):
    appendix = "fna"
else:
    appendix = "faa"
if filename[-1]=="fna" or filename[-1]=="faa":
    output =".".join(filename)
    Output = open(output, "w")
elif len(filename)==1:
    ouput = filename[0]+"."+appendix
    Output = open(output, "w")
else:
    filename[-1]=appendix
    ouput = output =".".join(filename)
    Output = open(output, "w")

for final in final_output:
    Output.write(str(final['name'])+"|"+"|len="+str(final['len']))
    Output.write("\n")
    
    n = 0
    for c in list(str(final['seq'])):
        Output.write(str(c))
        n+=1
        if n > num :
            n=0
            Output.write("\n")
    Output.write("\n")




def Sam_to_Fasta(fn):
    with open (fn, 'r') as f:
        for line in f:
            if (re.search(pattern = '@', string = line) is None):
                cols=line.rstrip().split("\t")
                print(">"+ cols[0] + "\n" + cols[9])

with open(Input_sequence,"r") as fh2:
    first_line = fh2.readline()
    if first_line.startswith("@SQ"):
        Sam_to_Fasta(Input_sequence)




