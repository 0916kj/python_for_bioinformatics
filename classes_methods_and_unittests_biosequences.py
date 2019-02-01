### define nucleotide and protein sequence classes, 
### relevant methods to perform computation and analysis, 
### and unit tests to confirm methods

import random

class Sequence:
    ### def class to represent biological sequence,
    ### such as DNA, RNA, or protein
    def __init__(self, sequence):
        """special function used to create a class object
        In this example, an object is a specific sequence"""
        self.sequence=sequence
        """class is now created"""
    
    def get_sequence(self):
        return self.sequence
        
    def length(self):
        """Returns length of sequence"""
        return len(self.sequence)
        
    def calculate_mismatches(self,sequence2):
        """Given 2 sentences of equal length, count number of times
        the characters differ"""
        mismatches=0
        for i in range(0,len(self.sequence)):
            if self.sequence[i]!=sequence2[i]:
                mismatches+=1
            else:
                continue
        return mismatches
        
    def __getitem__(self,i):
        """Tells class what to do if someone tries to access the
        sequence using an index"""
        return self.sequence[i]
        
    def sequence_element_count(self):
        """Return a dictionary where the keys are a string in
        the sequence and the values are the numbers of time
        that string occurs in the sequence.
        General solution to nucleotide counting in sequence"""
        count_dic=dict()
        string_elements=set(self.sequence) #set returns a storage object that only contains unique items
        #can perform set operations on them (Union, intersection...)
        for element in string_elements:
            print(element)
            count_dic[element]=self.sequence.count(element)
        return count_dic
    def randomized_sequence(self):
        """
        Returns randomly shuffled version of sequence
        """
        return ''.join(random.sample(self.sequence,len(self.sequence)))
        
class Nucleotide(Sequence):
    def nucleotide_count(self):
        return self.sequence_element_count()
        
    def percent_gc(self):
        """Returns percentage of G or C bases in sequence"""
        return calculate_percent_gc(self.sequence)
    def mask(self,subseq):
        """
        Identifies subsequence in nucleotide sequence
        and converts the subsequence to lowercase
        """
        self.subseq=subseq
        masked_seq=self.sequence
        return masked_seq.replace(self.subseq,self.subseq.lower())

class DNA(Nucleotide):
    def transcribe_DNA(self):
        transcribed_seq=self.sequence.replace("T","U")
        return transcribed_seq
    def reverse_complement(self):
        complement=''
        for base in self.sequence:
            if base=='A':
                complement+='T'
            elif base=='T':
                complement+='A'
            elif base=='G':
                complement+='C'
            elif base=='C':
                complement+='G'
        return complement[::-1]
        
class RNA(Nucleotide):
    def translate_RNA(self):
        """Given an RNA sequence translate it into a protein sequence. 
        Don't include the stop codon"""

        protein_seq=''
        codon=''
    
        for base in range(0,len(self.sequence),3):
            codon=self.sequence[base:base+3]
            if codon in codon_table.keys():
                if codon=='UAA':
                    break
                elif codon=='UAG':
                    break
                elif codon=='UGA':
                    break    
                else:
                    protein_seq+=codon_table[codon]
            else:
                print('Not a codon')
                break
        return protein_seq
    def get_DNA_sequence(self):
        DNA_seq=self.sequence.replace("U","T")
        return DNA_seq

class protein(Sequence):
    def get_mass(self):
        mass=0
        for amino_acid in self.sequence:
            if amino_acid in amino_acid_mass_dic.keys():
                mass+=amino_acid_mass_dic[amino_acid]
        return mass

class FASTA_File:
    def __init__(self,filepath,seq_type):
        ##seq_types: DNA, RNA, protein
        self.filepath=filepath
        self.seq_type=seq_type
        self.fasta_dictionary=nested_dic(self.filepath)
        self.accessions=list(self.fasta_dictionary)
    def read_fasta_file_obj(self):
         return read_fasta_file(self.filepath)
    def parse_fasta_record_obj(self):
         return parse_fasta_record(self.filepath)
    def accession_random_sample(self,sample_size):
        return random.sample(self.accessions,sample_size)
    def __and__(self,other_filepath):
        self.other_filepath=other_filepath
        self.other=FASTA_File(other_filepath,self.seq_type)
        self.intersection=[]
        for acc in self.accessions:
            if acc in self.other.accessions:
                self.intersection.append(acc)
        return self.intersection
    def __or__(self,other_filepath):
        self.other_filepath=other_filepath
        self.other=FASTA_File(other_filepath,self.seq_type)
        self.union=self.accessions
        for acc in self.other.accessions:
            if acc not in self.union:
                self.union.append(acc)
        return self.union
    def subset_write_out(self,output_fasta):
        self.output_fasta=output_fasta
        #open output file path to write new fasta to
        with open(self.output_fasta, 'w') as file:
            for acc in self.accessions:
                file.write('>'+acc+'\n')
                #write out each accession no. as a line
    def in_metadata(self,search_term,case_sensitive):
        #accepts search term and boolean to trigger case sensitive
        ##or non-case sensitive search
        acc_list=[]
        if case_sensitive==True:
            for acc in self.accessions:
                if search_term in self.fasta_dictionary[acc]["Metadata"]:
                    acc_list.append(acc)
        elif case_sensitive==False:
            for acc in self.accessions:
                if search_term.lower() in self.fasta_dictionary[acc]["Metadata"].lower():
                    #convert both search string and metadata strings to lowercase
                    ##for non=case sensitive search
                    acc_list.append(acc)
        else:
            print("Please use boolean for case_sensitive")
        return acc_list

RNA_codon_table = {
# 			Second Base
#    U 	         C             A             G
# U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Urp', # UxG
# C
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}

singleletter = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 
                'Trp': 'W', 'Asn': 'N', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 
                'Ala': 'A', 'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H', 
                'Arg': 'R', 'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 
                '---': '*'}
#practice problem
codon_table = dict()
for codon in RNA_codon_table:
    amino_acid = RNA_codon_table[codon]
    if amino_acid in singleletter:#see if key in dictionary to prevent error
            single_letter_name = singleletter[amino_acid]
            codon_table = [single_letter_name]=amino_acid

amino_acid_mass_dic = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

def calculate_percent_gc(sequence):
    """Given a nucleotide sequence return the percentage of the sequence that
    is G or C. Round to the first decimal place"""
    count = 0
    percent = 0
    for i in range(0,len(sequence)):
        if sequence[i] == 'G':
            count += 1
        if sequence[i] == 'C':
            count += 1
    percent = (count/len(sequence))*100
    return round(percent,1)

def read_fasta_file(file_path):
    """
    Uses a generator to read large fasta file
    """
    with open(file_path,"r") as in_file:
        record_list=[] #make a list to store each record
        count=-1
        #new_record=None
        for line in in_file:
            if '>' in line: #checks if new record by looking for accession no.
                count+=1
                record_list.append(line)                
            elif line:
                record_list[count]+=line
    for item in record_list:
        if item: #don't return None
            yield item #yield each record as a string

def parse_fasta_record(file_path):
    """
    Return each record as a dictionary containing acc, meta, and seq
    Store sequence as seq obj
    """
    fasta_read=read_fasta_file(file_path)
    #create gen obj that yields each record
    for record in fasta_read:
        split_record=record.replace("\n","|").split("|") #turn record into list
        accession=''
        metadata=''
        sequence=''
        fasta_record_dic=dict()
        for item in split_record:
            if '>' in item:
                accession=item.replace('>','')
            else:
                if item.isupper(): #searches for items with all uppercase
                    ##ie sequence
                    sequence+=item
                else:
                    metadata+=item         
        fasta_record_dic["Accession"]=accession.strip()
        fasta_record_dic["Metadata"]=metadata.strip()
        fasta_record_dic["Sequence"]=Sequence(sequence)
        ##sequence stored in dictionary as sequence object;
        ##nucleotide or aa sequence can be accessed through 
        ##fasta_record_dic["Sequence"].sequence
        yield fasta_record_dic

                     
def nested_dic(file_path):
    gen_obj=parse_fasta_record(file_path)
    fasta_dic=dict()
    for item in gen_obj:
        dic_obj=item
        fasta_dic[dic_obj["Accession"]]={"Sequence":dic_obj["Sequence"],
                 "Metadata":dic_obj["Metadata"]}
    return fasta_dic


output_file_path = "output.fasta"  
        
fasta_file_obj1=FASTA_File("protein_fasta1.fasta","protein")
fasta_file_obj2=FASTA_File("protein_fasta2.fasta","protein")
DNA_obj=DNA("ATGGTCGTGAGATGCGTGAGTTCGGTATG")

"""
Create unittests to verify methods
"""

if __name__=="__main__":
    import unittest
    
    class Test_fasta_objs(unittest.TestCase):
        def test__and__(self):
            fasta_file_obj1=FASTA_File(fasta1_file_path,"protein")
            #fasta_file_obj2=FASTA_File(fasta2_file_path,"protein")
            answer=['PF3D7_1353700']
            self.assertEqual(fasta_file_obj1.__and__(fasta2_file_path)
            ,answer)
        def test__or__(self):
            fasta_file_obj1=FASTA_File(fasta1_file_path,"protein")
            #fasta_file_obj2=FASTA_File(fasta2_file_path,"protein")
            answer=['PF3D7_1353700','PF3D7_0626800','PF3D7_0903300',
                    'PF3D7_1467500','PF3D7_1033100','PF3D7_1330100',
                    'PF3D7_1143600','PF3D7_0221800','PF3D7_1448000',
                    'PF3D7_1457100','PF3D7_0626801']
            self.assertEqual(fasta_file_obj1.__or__(fasta2_file_path),
                             answer)
        def test_in_metadata(self):
            fasta_file_obj1=FASTA_File(fasta1_file_path,"protein")
            answer=['PF3D7_0626800']
            self.assertEqual(fasta_file_obj1.in_metadata("pyruvate",False),
                             answer)
        def test_mask(self):
            dna_obj=DNA("AGTCGGAGT")
            answer="agtCGGagt"
            self.assertEqual(dna_obj.mask("AGT"),
                             answer)
            
    unittest.main()
    
