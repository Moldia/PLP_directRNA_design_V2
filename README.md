# Installation  
## Step 1: Install Docker
Before using the provided Dockerfile, ensure Docker is installed on your system. Follow the official instructions based on your operating system:

[Docker Installation Guide](https://docs.docker.com/engine/install/)  
Ensure Docker is properly installed by running the following command in your terminal:  
```
docker --version
```
## Step 2: Clone the GitHub Repository
Clone the repository containing the Dockerfile to your local machine. Replace <repository-url> with the actual URL of your repository:
```
git clone https://github.com/Moldia/PLP_directRNA_design_V2.git
cd PLP_directRNA_design_V2
```
## Step 3: Build the Docker Image
Build the Docker image using the docker build command. Replace <image-name> with a name of your choice for the image:
```
docker build -t plp_probe_design_v2 .
```
# Running the docker image   
After successfully building the image, create and run a container using the following command:
```
docker run -it -v $PWD:/app plp_probe_design_v2
```

You can test whether `cutadapt` is installed:  
```
cutadapt --version
```

# Running  
## Extract features
`python3 codes/extract_features.py --gtf data/tmp.gtf --genes Grik2 --identifier_type gene_name --gene_feature CDS --output extract_features_output.txt`

## Extract transcriptome
`python3 codes/extract_mrna.py --gtf data/tmp.gtf --fasta data/Mus.fa --output_file data/transcriptome_out.fa`  

_This step can be parallelized together with `extract_features.py` and the output is required by last step `find_target.py`_  

## Extract sequences
`python3 codes/extract_sequences.py --fasta data/Mus.fa --output_fasta extract_seqs_output.fa --identifier_type gene_name --plp_length 30 --gtf_output extract_features_output.txt`

## Find targets
`python3 codes/find_target.py --selected_features extract_features_output.txt --fasta_file extract_seqs_output.fa --output_file targets.txt --iupac_mismatches "5:R,10:G" --reference_fasta data/transcriptome_out.fa --max_errors 4 --Tm_min 58 --Tm_max 62 --lowest_percentile_Tm_score_cutoff 5 --min_dist_probes 8 --filter_ligation_junction`


## Run whole workflow  
```
python3 codes/run_probe_design.py \
--extract_features_gtf data/tmp.gtf \
--extract_features_genes Grik2 \
--extract_features_identifier_type gene_name \
--extract_features_gene_feature CDS \
--extract_features_output extract_features_output.txt \
--extract_transcriptome_gtf data/tmp.gtf \
--extract_transcriptome_fasta data/Mus.fa \
--extract_transcriptome_output_file data/transcriptome_out.fa \
--extract_sequences_fasta data/Mus.fa \
--extract_sequences_output_fasta extract_seqs_output.fa \
--extract_sequences_identifier_type gene_name \
--extract_sequences_plp_length 30 \
--extract_sequences_gtf_output extract_features_output.txt \
--find_target_selected_features extract_features_output.txt \
--find_target_fasta_file extract_seqs_output.fa \
--find_target_output_file targets.txt \
--find_target_iupac_mismatches 5:R,10:G \
--find_target_reference_fasta data/transcriptome_out.fa \
--find_target_max_errors 4 \
--find_target_Tm_min 58 \
--find_target_Tm_max 62 \
--find_target_lowest_percentile_Tm_score_cutoff 5 \
--find_target_min_dist_probes 8 \
--find_target_filter_ligation_junction \
--find_target_num_probes 15 
```mermaid
flowchart LR
  %% Inputs (Top Layer)
  gtf[gtf]
  genes[genes]
  identifier_type[identifier_type]
  gene_feature[gene_feature]
  fasta[fasta]
  min_coverage[min_coverage]
  gc_min[gc_min]
  gc_max[gc_max]
  num_probes[num_probes]
  iupac_mismatches[iupac_mismatches]
  max_errors[max_errors]
  check_specificity[check_specificity]
  plp_length[plp_length]
  Tm_min[Tm_min]
  Tm_max[Tm_max]
  lowest_percentile_Tm_score_cutoff[lowest_percentile_Tm_score_cutoff]
  min_dist_probes[min_dist_probes]
  filter_ligation_junction[filter_ligation_junction]

  %% Script
  features[extract_feature.py]
  transcriptome[extract_mrna.py]
  sequences[extract_sequences.py]
  probes[find_target.py]
  
  %% Output
  output1[extracted_features.txt]
  output2[transcriptome_output.fa]
  output3[extract_seqs_output.fa]
  probes.txt[probes.txt]

  %% Connections
  gtf --> features
  genes --> features
  identifier_type --> features
  gene_feature --> features
  features --> output1

  gtf --> transcriptome
  fasta --> transcriptome
  transcriptome --> output2

  fasta --> sequences
  identifier_type --> sequences
  plp_length --> sequences
  output1 --> sequences 
  sequences --> output3

  output1 --> probes
  output2 --> probes
  output3 --> probes
  iupac_mismatches --> probes
  fasta --> probes
  max_errors--> probes
  Tm_min--> probes
  Tm_max--> probes
  lowest_percentile_Tm_score_cutoff--> probes
  min_dist_probes--> probes
  filter_ligation_junction--> probes
  min_coverage--> probes
  gc_min--> probes
  gc_max--> probes
  num_probes--> probes
  check_specificity--> probes
  plp_length--> probes
  probes --> probes.txt
```
