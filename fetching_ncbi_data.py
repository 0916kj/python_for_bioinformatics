## Search for NCBI records by keyword (i.e., organism), 
## write results to out file (i.e., .FASTA or .gb) and read results

from Bio import Entrez
from Bio import SeqIO
from time import sleep ##pause between ncbi searches
## due to NCBI restrictions

Entrez.email = input("Email address: ") ##Entrez requires email
## to communicate with NCBI databases

def id_list(db, number_of_results, keyword):
    ##db is NCBI database being searched, i.e. "nucleotide"
    handle = Entrez.esearch(db, retmax=number_of_results,
                            term=keyword)
    record = Entrez.read(handle) #creates dictionary with search info
    handle.close()
    return record["IdList"] ##return list of NCBI id's, which
## can be used to fetch records

def write_ids_to_file(db, list_of_ids, retrieval_type,
                      record_type, file_type, outfile):
    ## retrieval_type determines format of retrieved data, 
    ## i.e. FASTA or gb (GenBank)
    ## record_type is type of record you wish to pull
    ## on specified id number, i.e. "genbank"
    ## file_type is type of file to write
    with open(outfile,"w") as f:
        f.write("")
    for id_num in list_of_ids:
    ##iterates through all id numbers in list, can use output of id_list(...)
        handle = Entrez.efetch(db, id=id_num, rettype=retrieval_type,
                             retmode="text")
        record = SeqIO.read(handle, record_type)
        handle.close()
        with open(outfile,"a") as out:
            SeqIO.write(record, out, file_type)
        sleep(1) ##one second delay between searches
            ## writes to file, formats as genbank file

def parse_output_file(outfile, record_type):
    with open(outfile, "r") as handle:
        for record in SeqIO.parse(handle, record_type) :
            print(repr(record)) #prints basic info about new genbank record
            print("\n")
    handle.close()

def search_term_to_parsed_output(db, number_of_results, keyword,
                                 retrieval_type, record_type, file_type, 
                                 outfile):
    ## combine functions to go directly from search to parsed output file
    ids = id_list(db, number_of_results, keyword)
    write_ids_to_file(db, ids, retrieval_type, record_type, file_type, outfile)
    ## write desired records to outfile, which can now be parsed
    parse_output_file(outfile, record_type)


search_term_to_parsed_output("nucleotide", 10, "Mus musculus",
                                 "fasta", "fasta", "fasta", 
                                 "mouse_fasta_record.fasta")
## parse results of first 10 mouse FASTA records on nucleotide database