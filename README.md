
# Task: Build a Gene Expression Matrix with Manifest from TCGA (BRCA gene) 
### Leiara Rivera, Clemson University

## Description: 
This Python script is designed to efficiently download RNAseq data from the TCGA-BRCA GDC manifest, process the data, and create a Gene Expression Matrix (GEM). It assumes that the GDC client has already been installed with its latest version and is used to download the necessary data. The script will estimate memory requirements to know if storage is optimal for processing datasets. Then, process the gene expression data from each sample and save the results in a '.tsv' file. 

##  Objective: 
The goal is to download and process RNAseq data from TCGA-BRCA and build a GEM that maps genes to expression values across samples. The output is the full matrix in a '.tsv' file called 'genexp_optimize_matrix.tsv'. 

## Requirements 
- Python 3
- Install the following Libraries: 
```
import os
import glob
import pandas as pd
from tqdm import tqdm
import psutil
import gc
from datetime import datetime

```
### Set your working directory: 
```
def main():
    # Set the main directory
    main_dir = "/scratch/leiarar/tca-brca-trans/test"
    print(f"Working in directory: {main_dir}")
    
    # Change to the main directory
    os.chdir(main_dir)

```

## File Structure and Location to be used as example! 
- '/scratch/leiarar/tca-brca-trans/' is my directory containing RNAseq data samples 
- 'genexp_optimize_matrix.tsv' is my output file containing the Gene Expression Matrix 

### To estimate the memory required for downloading the datasets, use the following: 
```
 # Estimate memory requirements
    est_memory_gb = estimate_memory_requirements(sample_dirs)
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    
    print(f"Estimated memory requirement: {est_memory_gb:.2f} GB")
    print(f"Available memory: {available_memory_gb:.2f} GB")
    
    if est_memory_gb > available_memory_gb:
        print(f"WARNING: This process may require more memory than available.")
        print(f"Consider processing in batches or on a machine with at least {max(16, int(est_memory_gb*1.5))} GB RAM")

```
## Output 
* Output Filename: 'genexp_optimize_matrix.tsv'
* Matrix Dimensions: Rows represent genes and columns represent samples. 
* Memory Usage: The script will display memory usage statistics when the matrix is saved. 

## How to Run Script 
1. Make sure you have the required libraries installed. 
2. Place your downloaded RNAseq data under the directory '/scratch/leiarar/tca-brca-trans/'. 
3. Run the Python script: 'python genexp_optimize_matrix.py'

  ### Full Script to be processed: 
  ````
#!/usr/bin/env python3

import os
import glob
import pandas as pd
from tqdm import tqdm
import tempfile
import psutil
import gc
from datetime import datetime

def estimate_memory_requirements(sample_dirs, avg_genes=60000):
    """Estimate memory requirements based on number of samples and genes"""
    # Rough estimate: each value takes ~8 bytes (float64)
    # Plus overhead for pandas DataFrame
    estimated_bytes = avg_genes * len(sample_dirs) * 8 * 1.5  # 1.5 factor for overhead
    estimated_gb = estimated_bytes / (1024**3)
    return estimated_gb

def main():
    # Set the main directory
    main_dir = "/scratch/leiarar/tca-brca-trans/"
    print(f"Working in directory: {main_dir}")
    
    # Change to the main directory
    os.chdir(main_dir)
    
    # Get all sample directories
    sample_dirs = [d for d in os.listdir(main_dir) 
                  if os.path.isdir(d) and not d.startswith(".") and not d.startswith("temp_")]
    
    print(f"Found {len(sample_dirs)} sample directories")
    
    # Estimate memory requirements
    est_memory_gb = estimate_memory_requirements(sample_dirs)
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    
    print(f"Estimated memory requirement: {est_memory_gb:.2f} GB")
    print(f"Available memory: {available_memory_gb:.2f} GB")
    
    if est_memory_gb > available_memory_gb:
        print(f"WARNING: This process may require more memory than available.")
        print(f"Consider processing in batches or on a machine with at least {max(16, int(est_memory_gb*1.5))} GB RAM")
    
    # First, extract gene information from the first sample to get gene names
    print("Extracting gene information from first sample...")
    first_sample = sample_dirs[0]
    tsv_files = glob.glob(f"{first_sample}/*.tsv")
    
    if not tsv_files:
        print(f"Error: No TSV file found in first sample directory {first_sample}")
        return
    
    # Dictionary to store gene_name -> index mapping for faster lookups
    gene_names = []
    
    try:
        with open(tsv_files[0], 'r') as f:
            # Skip the first 6 lines (header and metadata)
            for _ in range(6):
                next(f)
            
            # Process each line to extract gene_name (column 1)
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) > 5:
                    gene_name = fields[1]  # gene_name column
                    gene_names.append(gene_name)
    except Exception as e:
        print(f"Error processing {tsv_files[0]}: {e}")
        return
    
    print(f"Extracted information for {len(gene_names)} genes")
    
    # Initialize the DataFrame with gene_name as index
    df = pd.DataFrame(index=gene_names)
    
    # Process each sample
    print("Processing samples...")
    for sample_dir in tqdm(sample_dirs):
        tsv_files = glob.glob(f"{sample_dir}/*.tsv")
        
        if not tsv_files:
            print(f"Warning: No TSV file found in {sample_dir}, skipping")
            continue
        
        try:
            # Create a temporary dictionary to hold this sample's data
            sample_data = {}
            
            with open(tsv_files[0], 'r') as f:
                # Skip the first 6 lines
                for _ in range(6):
                    next(f)
                
                # Read gene expression values
                for i, line in enumerate(f):
                    if i >= len(gene_names):
                        break
                        
                    fields = line.strip().split('\t')
                    if len(fields) > 5:
                        gene_name = fields[1]
                        expression_value = fields[5]  # stranded_second column
                        sample_data[gene_name] = expression_value
            
            # Add this sample's data as a new column
            df[sample_dir] = pd.Series(sample_data)
            
        except Exception as e:
            print(f"Error processing {sample_dir}: {e}")
            df[sample_dir] = pd.Series(["NA"] * len(gene_names), index=gene_names)
        
        # Force garbage collection to free memory
        gc.collect()
    
    # Fill NA values
    df.fillna("NA", inplace=True)
    
    # Set the output filename
    output_filename = "genexp_optimize_matrix.tsv"
    
    # Save the DataFrame to TSV
    print(f"Saving matrix to {output_filename}...")
    df.to_csv(output_filename, sep='\t')
    
    print(f"\nMatrix created: {output_filename}")
    print(f"Matrix dimensions: {df.shape[0]} genes × {df.shape[1]} samples")
    print(f"Memory usage: {df.memory_usage(deep=True).sum() / (1024**3):.2f} GB")

if __name__ == "__main__":
    start_time = datetime.now()
    main()
    end_time = datetime.now()
    print(f"Total execution time: {end_time - start_time}")
````

## Troubleshooting Tips 
1. Memory Issues: 
Problem: The script estimates that it needs more memory than is available on your system. 
Solution: Consider processing the data in small batches (like a test section), or use a system with more RAM (at least 16GB ideally). 
2. Missing TSV files: 
Problem: The script cannot find .tsv files in some sample directories. 
Solution: Ensure that the sample directories contain '.tsv' files with the necessary gene expression data. If a directory is missing '.tsv' files, it will be skipped with a warning.
3. Create a temporary file to run script. 
 

#### Libraries required: 
```
#!/usr/bin/env python3

import os
import glob
from tqdm import tqdm
import gc
import tempfile
import shutil
import datetime
```
## Use script to create a temporary file: 
```
  # Create a temporary directory
    temp_dir = tempfile.mkdtemp(prefix="gene_matrix_")
    print(f"Created temporary directory: {temp_dir}")
```

## Python Script for Creating temporary file:

```
#!/usr/bin/env python3

import os
import glob
from tqdm import tqdm
import gc
import tempfile
import shutil
import datetime

def main():
    # Set the main directory
    main_dir = "/scratch/leiarar/tca-brca-trans/test"
    print(f"Working in directory: {main_dir}")
    
    # Change to the main directory
    os.chdir(main_dir)
    
    # Create a temporary directory
    temp_dir = tempfile.mkdtemp(prefix="gene_matrix_")
    print(f"Created temporary directory: {temp_dir}")
    
    # Get all sample directories
    sample_dirs = [d for d in os.listdir(main_dir) 
                  if os.path.isdir(d) and not d.startswith(".") and not d.startswith("temp_")]
    
    print(f"Found {len(sample_dirs)} sample directories")
    
    # First, extract gene information from the first sample to get gene_id and gene_name
    print("Extracting gene information from first sample...")
    first_sample = sample_dirs[0]
    tsv_files = glob.glob(f"{first_sample}/*.tsv")
    
    if not tsv_files:
        print(f"Error: No TSV file found in first sample directory {first_sample}")
        return
        # Dictionary to store gene_id -> gene_name mapping
    gene_info = {}
    
    try:
        with open(tsv_files[0], 'r') as f:
            # Skip the first 6 lines (header and metadata)
            for _ in range(6):
                next(f)

            # Process each line to extract gene_id and gene_name
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) > 5:
                    gene_id = fields[0]  # gene_id column
                    gene_name = fields[1]  # gene_name column
                    gene_info[gene_id] = gene_name
    except Exception as e:
        print(f"Error processing {tsv_files[0]}: {e}")
        return
    
    print(f"Extracted information for {len(gene_info)} genes")
    
    # Set the output filename to 14th_try.tsv as requested
    output_filename = "14th_try.tsv"
    
     # Create the output file and write the header
    with open(output_filename, 'w') as out_f:
        # Write the header line
        header = ['gene_id', 'gene_name'] + sample_dirs
        out_f.write('\t'.join(header) + '\n')

        # Write data for each gene
        print("Writing gene expression data...")
        for gene_id, gene_name in tqdm(gene_info.items()):
            # Start a new line with gene_id and gene_name
            line_data = [gene_id, gene_name]

            # For each sample, extract the stranded_second value
            for sample_dir in sample_dirs:
                tsv_files = glob.glob(f"{sample_dir}/*.tsv")
                expression_value = "NA"  # Default if not found

                if tsv_files:
                    try:
                        with open(tsv_files[0], 'r') as f:
                            # Skip the first 6 lines
                            for _ in range(6):
                                next(f)

                            # Look for the gene in this sample
                            for line in f:
                                fields = line.strip().split('\t')
                                if len(fields) > 5 and fields[0] == gene_id:
                                    expression_value = fields[5]  # stranded_second column
                                    break
                    except Exception as e:
                        print(f"Error reading {tsv_files[0]}: {e}")

                line_data.append(expression_value)

            # Write the complete line to the output file
            out_f.write('\t'.join(line_data) + '\n')
            
  print(f"\nMatrix created: {output_filename}")
    print(f"Matrix dimensions: {len(gene_info)} genes × {len(sample_dirs)} samples")

if __name__ == "__main__":
    main()

```   


### License - this project is licensed under the MIT License. 


### Contact Information
Leiara Rivera 
Email: leiarar@clemson.edu 
