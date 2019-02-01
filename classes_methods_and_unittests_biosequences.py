#####################In class Assignment########################################
###Create an DNA class that inherits from class Nucleotide and contains the methods:
            ### transcribe_dna,reverse_complement
###Create an RNA class that inherits from class Nucelotide and contains the methods:
    ### translate_rna and get_dna_sequence

###Create a Protein class that inherits from class Sequence and contains the methods:
    ###get_mass (use the amino_acid_mass_dic)
    
###Create unit tests to test your class methods (see nucleotide_sequence_functions.py)
    ###You should be able to figure out how to get unit-test working by
    ###using the tests in nucleotide_sequence_functions.py as a template
    
codon_table= {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
amino_acid_mass_dic = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }


def calculate_percent_gc(sequence):
    """Given a nucleotide sequence return the percentage of the sequence that
    is G or C. Round to the first decimal place"""
    count=0
    percent=0
    for i in range(0,len(sequence)):
        if sequence[i]=='G':
            print('G')
            count+=1
        if sequence[i]=='C':
            print('C')
            count+=1
    percent=(float(count)/len(sequence))*100
    return round(percent,1)

class Sequence:
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
        string_elements=set(self.sequence) ##set returns a storage object that 
        ## only contains unique items
        ## can perform set operations on them (Union, intersection...)
        for element in string_elements:
            print(element)
            count_dic[element]=self.sequence.count(element)
        return count_dic
        
class Nucleotide(Sequence):
    def nucleotide_count(self):
        return self.sequence_element_count()
        
    def percent_gc(self):
        """Returns percentage of G or C bases in sequence"""
        return calculate_percent_gc(self.sequence)

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

if __name__=="__main__":
    import unittest
    class TestMethods(unittest.TestCase):
        def test_transcribe_DNA(self):
            sequence=DNA("AGCT")
            answer="AGCU"
            self.assertEqual(sequence.transcribe_DNA(),answer)
            
        def test_reverse_complement(self):
            sequence=DNA("AGGCT")
            answer="AGCCT"
            self.assertEqual(sequence.reverse_complement(),answer)
    unittest.main()