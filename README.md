
#Leiara Rivera, Clemson University

#Description: 
This Python script is designed to download RNAseq data from the TCGA-BRCA GDC manifest, process the data, and create a Gene Expression Matrix (GEM). It assumes that the GDC client has already been installed with its latest version, and is used to download the necessary data. The script will estimate memory requirements, process the gene expression data from each sample, and save the results in a '.tsv' file. 

#Objective: 
The goal is to download and process RNAseq data from TCGA-BRCA, build a GEM that maps genes to expression values across samples, and the output is the full matrix in a '.tsv' file called 'genexp_optimize_matrix.tsv'. 

#Requirements 
- Python 3.x 
- Libraries: pandas, tqdm, psutil 

To install the required Libraries, use the following: 

pip install pandas tqdm psutil

#File Structure and Location to be used as example! 
- '/scratch/leiarar/tca-brca-trans/' is my directory containing RNAseq data samples 
- 'genexp_optimize_matrix.tsv' is my output file containing the Gene Expression Matrix 

#How to Run Script 
1. Make sure you have the required libraries installed. 
2. Place your downloaded RNAseq data under the directory '/scratch/leiarar/tca-brca-trans/'. 
3. Run the Python script: python genexp_optimize_matrix.py 

#Output 
Output Filename: 'genexp_optimize_matrix.tsv'
Matrix Dimensions: Rows represent genes and columns represent samples. 
Memory Usage: The script will display memory usage statistics when the matrix is saved. 

#Troubleshooting Tips 
1. Memory Issues: 
Problem: The script estimates that it needs more memory than is available on your system. 
Solution: Consider processing the data in small batches (like a test section), or use a system with more RAM (at least 16GB ideally). 
2. Missing TSV files: 
Problem: The script cannot find .tsv files in some sample directories. 
Solution: Ensure that the sample directories contain '.tsv' files with the necessary gene expression data. If a directory is missing '.tsv' files, it will be skipped with a warning. 

#License - this project is licensed under the MIT License. 


Contact 
Leiara Rivera 
Bioengineering Student, Clemson University 
Email: leiarar@clemson.edu 
