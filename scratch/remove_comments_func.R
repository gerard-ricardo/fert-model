input  = "./scratch/Fertmod_commented.R"
output = "./R_scripts/Backend.R"

#input  = './scratch/config_commented.R'
#output = "./R_scripts/config func.R"


process_script_v2 <- function(input_file, output_file) {
  content <- readLines(input_file)
  processed_content <- vector("character", length = length(content))
  for (i in seq_along(content)) {
    line <- content[i]
    # Directly keep lines starting with ##
    if (grepl("^\\s*##", line)) {
      processed_content[i] <- line
      # Also keep section headers used in RStudio (Ctrl + Shift + R)
    } else if (grepl("^\\s*#.*-{4,}", line)) {  # Look for # followed by at least four dashes
      processed_content[i] <- line
    } else {
      # Remove single # comments or inline comments, preserving code
      processed_content[i] <- sub("\\s*#.*$", "", line)
    }
  }
  # Optionally, filter out completely empty lines if desired
  processed_content <- processed_content[processed_content != ""]
  # Write the processed content to the output R script
  writeLines(processed_content, output_file)
  # Clear all lines in the input file by writing an empty string
  #writeLines("", input_file)
}
process_script_v2(input, output)

add_empty_lines_before_comments <- function(input_file) {
  content <- readLines(input_file)
  processed_content <- vector("character", length = length(content) * 2) # Allocate enough space
  j <- 1 # Index for processed_content
  for (i in seq_along(content)) {
    line <- content[i]
    # Check for lines starting with ## or section comments
    if (grepl("^\\s*##", line) || grepl("^\\s*#.*-{4,}", line)) {
      # Insert an empty line before this line if it's not the first line
      if (i > 1) {
        processed_content[j] <- "" # Add an empty line
        j <- j + 1
      }
    }
    processed_content[j] <- line
    j <- j + 1
  }
  # Trim the processed_content vector to remove unused elements
  processed_content <- processed_content[1:j-1]
  # Write the processed content back to the input file
  writeLines(processed_content, input_file)
}
add_empty_lines_before_comments(output)

clear_file_content <- function(input_file) {
  # Overwrite the file with an empty string, effectively clearing its content
  writeLines("", input_file)
}

# Example usage
clear_file_content(input)

