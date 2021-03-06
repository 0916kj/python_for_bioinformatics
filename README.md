# python_for_bioinformatics
Programs written for graduate Python for Bioinformatics course.
- Objectives: Create classes and methods for biological data (i.e., protein and nucleotide sequences), programs to parse, analyze, and write biological data in common bio research file formats, automate database searches and downloads from NCBI servers, and utilize various data structures and analysis tools to understand biological data
- Originally written for Python 2, edited for Python 3 when necessary.
 - classes_methods_and_unittests_biosequences.py: define nucleotide and protein sequence classes, relevant methods to perform computation and analysis, and unit tests to confirm methods
 - protein_fasta1 and protein_fasta2: fasta files used in unit tests in classes_methods_and_unittests_biosequences.py
 - fetching_ncbi_data: search for NCBI records by keyword, write to .fasta or .gb file, then parse file for name, sequence, and annotations summary
- search_genes_retrieve_protein.py: search proteins by relevant genes, print proteins to fasta, run BLAST alignment/search on proteins, then print BLAST results to DataFrame
- machine_learning_tutorial.py: ML Practice from https://medium.freecodecamp.org/the-hitchhikers-guide-to-machine-learning-algorithms-in-python-bfad66adb378
Practice on ML basics with Scikit-Learn: Regressions, Decision Trees, Support Vector Machines, K-Nearest Neighbors, Visualization
