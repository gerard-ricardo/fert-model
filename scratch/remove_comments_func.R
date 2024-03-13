process_script_v2 <- function(input_file, output_file) {
  content <- readLines(input_file)
  processed_content <- vector("character", length = length(content))

  for (i in seq_along(content)) {
    line <- content[i]

    # Directly keep lines starting with ##
    if (grepl("^\\s*##", line)) {
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
  writeLines("", input_file)
}

process_script_v2("./scratch/Fertmod_commented.R", "./R/Fertmod.R")

