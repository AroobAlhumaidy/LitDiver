#!/bin/bash
# This script aims to fetch the pubmed database and retrive the metadata of the papers. Based on a set of keywords provided in a text file 
# Log: this is prod version of the dev version: fetch_pubmed_lit.Edirect.v1.5.sh
# Usage: time bash fetch_pubmed_lit.prod.v1.0.sh mesh_keywords.txt
# By: Aroob ALhumaidy

if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"

if [ ! -f "$input_file" ]; then
  echo "Input file not found: $input_file"
  exit 1
fi

# temp files
command_file="commands_to_run.sh"
error_file="failed_keywords.txt"

> "$command_file"
> "$error_file"

# generate teh corresponding command to each keyword (loop) 
while IFS= read -r keyword || [ -n "$keyword" ]; do
  keyword=$(echo "$keyword" | sed 's/[“”]/"/g')

  clean_keyword=$(echo "$keyword" | sed 's/["()]/\\&/g')

  # Generate the PubMed search command and append it to the command file
  echo "esearch -db pubmed -query \"$clean_keyword\" | efetch -format medline > \"results_Pubmed_${clean_keyword}.txt\"" >> "$command_file"
done < "$input_file"

echo "Generated all commands. Now executing them..."

chmod +x "$command_file"

# Execute the entire command file at once and capture errors
bash "$command_file" || echo "Some commands failed. Check $error_file for details."

echo "Searches completed."
