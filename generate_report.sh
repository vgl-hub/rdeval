#!/bin/sh

# Initialize variables
RDEVAL=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
output_file=""
input_files=()

parsing_inputs=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output)
      output_file="$2"
      shift 2
      ;;
    -i|--input)
      parsing_inputs=true
      shift
      ;;
    -*)
      echo "Unknown option: $1"
      echo "Usage: $0 -o output_file.html -i input1 [input2 ...]"
      exit 1
      ;;
    *)
      if $parsing_inputs; then
        input_files+=("$1")
        shift
      else
        echo "Unexpected argument: $1"
        echo "Usage: $0 -o output_file.html -i input1 [input2 ...]"
        exit 1
      fi
      ;;
  esac
done

# Validate required arguments
if [[ -z "$output_file" ]] || [[ "${#input_files[@]}" -eq 0 ]]; then
  echo "Error: Both --output and at least one --input are required."
  echo "Usage: $0 -o output_file.html -i input1 [input2 ...]"
  exit 1
fi

r_vector="c("
for file in "${input_files[@]}"; do
  r_vector+="'$(readlink -f ${file})',"
done
r_vector="${r_vector%,})"  # Remove trailing comma and close

# Run the R command
R -e "rmarkdown::render('${RDEVAL}/figures.Rmd', output_file='$(pwd)/$output_file', params=list(input_files = $r_vector))"
