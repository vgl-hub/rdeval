#!/bin/sh

# Usage text
usage="Usage: $(basename "$0") -i input1.rd [input2.rd ...] -o output_file [--dynamic] [-h]

Arguments:
    -i, --input     Input Rmd file(s)
    -o, --output    Output filename (must end in .html or .pdf)
    --dynamic       Enable dynamic plotting, only supported in html (default: static)
    -h, --help      Show help
"

# Initialize variables
RDEVAL=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Initialize variables
output_file=""
input_files=()
dynamic="FALSE"
parsing_inputs=false

# Parse arguments
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
    -d|--dynamic)
      dynamic="TRUE"
      shift
      ;;
    -h|--help)
      echo "$usage"
      exit 0
      ;;
    -*)
      echo "Unknown option: $1"
      echo "$usage"
      exit 1
      ;;
    *)
      if $parsing_inputs; then
        input_files+=("$1")
        shift
      else
        echo "Unexpected argument: $1"
        echo "$usage"
        exit 1
      fi
      ;;
  esac
done

# Validate required arguments
if [[ -z "$output_file" ]] || [[ "${#input_files[@]}" -eq 0 ]]; then
  echo "Error: Both --output and at least one --input are required."
  echo "$usage"
  exit 1
fi

# Validate output extension and infer output format
if [[ "$output_file" == *.html ]]; then
  output_format="html"
elif [[ "$output_file" == *.pdf ]]; then
  output_format="pdf"
else
  echo "Error: Output file must end with .html or .pdf"
  exit 1
fi

# If dynamic is TRUE but format is PDF, override and warn
if [[ "$dynamic" == "TRUE" && "$output_format" == "pdf" ]]; then
  echo "Warning: Dynamic plotting is only supported with HTML. Switching to static plotting."
  dynamic="FALSE"
fi

# Construct input file vector for R
r_vector="c("
for file in "${input_files[@]}"; do
  r_vector+="'$(readlink -f "$file")',"
done
r_vector="${r_vector%,})"  # Remove trailing comma and close

# Run the R command with additional parameters
R -e "rmarkdown::render('${RDEVAL}/figures.Rmd', output_file='$(pwd)/$output_file', output_format='${output_format}_document', params=list(input_files = $r_vector, dynamic = $dynamic))"
