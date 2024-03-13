process_script_v2 <- function(input_file, output_file) {
  content <- readLines(input_file)
  processed_content <- vector("character", length = length(content))

  for (i in seq_along(content)) {
    line <- content[i]

    # Directly keep lines starting with ##
    if (grepl("^\\s*##", line)) {
      processed_content[i] <- line
    }
    # Remove whole-line comments starting with a single #
    # This regex checks if a line starts with optional whitespace followed by a single #
    # and removes the line if it matches.
    # else if (grepl("^\\s*#[^#]", line)) {
    #   processed_content[i] <- sub("^\\s*#.*$", "", line)
    # }
    # Remove inline comments, preserving code
    else {
      processed_content[i] <- sub("\\s*#.*$", "", line)
    }
  }

  # Optionally, filter out completely empty lines if desired
  processed_content <- processed_content[processed_content != ""]

  # Write the processed content to the output R script
  writeLines(processed_content, output_file)
}


process_script_v2("./scratch/Fertmod_comments.R", "./R/Fertmod.R")
