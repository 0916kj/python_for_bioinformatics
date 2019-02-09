### Find relevant genes, look up proteins
### Write proteins to Fasta and search Blast
### Return DataFrame containing Blast results
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import SearchIO
import pandas as pd
from time import sleep ##pause between ncbi searches
## due to NCBI restrictions

Entrez.email=input("Email: ") ##give email to access NCBI servers

def gene_list(search_term,max_results):
    ##write list of genes associated with term to list
    handle=Entrez.esearch(db="gene",term=search_term,retmax=max_results)
    record=Entrez.read(handle)
    handle.close()
    ids=record["IdList"]
    title_list=[]
    for entry in ids:
        handle = Entrez.esummary(db="gene",id=entry)
        record = Entrez.read(handle)
        handle.close()
        title_list.append(record["DocumentSummarySet"]["DocumentSummary"]
        [0]["Name"])
        sleep(0.5)  ##delay between searches
    return title_list

def protein_search_dic(search_term,max_results):
    ##search protein db for each gene name, write results to list,
    ###store in dictionary where key = gene
    genes_of_interest = gene_list(search_term,max_results)
    gene_protein_dic={}
    for gene in genes_of_interest:
        protein_search_term=gene+" AND homo sapiens[Orgn]"
        handle=Entrez.esearch(db="protein",term=protein_search_term)
        record=Entrez.read(handle)
        handle.close()
        ids=record["IdList"]
        list_of_proteins=[]
        sleep(0.5)
        for entry in ids:
            handle = Entrez.esummary(db="protein",id=entry)
            record = Entrez.read(handle)
            handle.close()
            list_of_proteins.append(record[0]["Title"])
            sleep(0.5)
        gene_protein_dic[gene]=list_of_proteins
    return gene_protein_dic
        
def protein_search_fasta(search_term,fasta_file_out,max_results):
    ###search for proteins related to identified genes
    ####return FASTA file
    ids=[]
    genes_of_interest = gene_list(search_term,max_results)
    for gene in genes_of_interest:
        protein_search_term=gene+" AND homo sapiens[Orgn]"
        handle=Entrez.esearch(db="protein",term=protein_search_term,
                              retmax=5)
        record=Entrez.read(handle)
        handle.close()
        ids+=record["IdList"]
        sleep(0.5)
    records=[]
    for id_num in ids:
        handle=Entrez.efetch(db="protein",id=id_num,retmode='text',
                             rettype='fasta')
        record = SeqIO.read(handle,"fasta")
        handle.close()
        records.append(record)
        sleep(0.5)
    SeqIO.write(records, fasta_file_out, "fasta")

def blast_result_df(fasta_file_in,blast_type,db,
                    xml_file_out):
    ##run blast on list of proteins in fasta file
    ##parse blast results, return dictionary containing queries and hits
    with open(fasta_file_in) as f:
        fasta_string=f.read()
    with NCBIWWW.qblast(blast_type, db, fasta_string) as result_handle:
        with open(xml_file_out, 'w') as xml_file: 
            ##write blast results to xml
            xml_file.write(result_handle.read())
            ##parse xml with searchio to find hits
    blast_qresult = SearchIO.parse(xml_file_out, "blast-xml")
    id_list=[]
    desc_list=[]
    hits_list=[]
    for record in blast_qresult:
        if record.id in id_list:
            continue
        else:
            id_list.append(record.id)
            desc_list.append(record.description)
            hits_list.append(record.hits)
    result_dict={"Id":id_list,"Description":desc_list,"Hits":hits_list} 
    result_df=pd.DataFrame.from_dict(data=result_dict)
    return result_df
        
protein_search_fasta("tinnitus AND homo sapiens[Orgn]",
                     "tinnitus_proteins.fasta",5)
## create DataFrame of Blast results for search on human
## tinnitus genes/proteins
print(blast_result_df("tinnitus_proteins.fasta","blastp","swissprot",
                "protein_blast.xml"))