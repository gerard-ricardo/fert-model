# ============================================================================
# Front-End Test for Coral Model
# ============================================================================
# This script is designed to set up and run the coral model simulations.
# Users can select from predefined scenarios or customize settings via a Shiny app.

# ============================================================================
# Load Required Scripts
# ============================================================================
# Before running the model, source the necessary R scripts.
# These scripts define the model settings and the backend logic.

source("./R_scripts/settings.R")           # Model settings
source("./R_scripts/1_print_outputs.R") # Model outputs
source("./R_scripts/Backend.R") # Model backend logic

# Creating a task-specific folder inside 'out' directory to avoid overwriting results
dir_name <- "./out/test3"  #NAME THE FOLDER
dir.create(file.path(getwd(),  dir_name), showWarnings = T)
new_folder_path <- file.path(getwd(), dir_name)



# ============================================================================
# Choose a Scenario
# ============================================================================
# There are three predefined scenarios: 'default', 'test', and 'custom'.
# - The 'default' scenario uses a spacing of 2m between 9 total colonies, arranged in a square.
# - The 'test' scenario is similar but with reduced time parameters for quicker runs.
#   Note: The grid size is relative to no_timesteps, so fewer no_timesteps decreases computational time significantly.
# - The 'custom' settings can be adjusted through the Shiny app for more granular control.
#   To customize settings, run the Shiny app script within the 'inst' folder and use the UI to adjust your preferences.
#   Use patch_type = "patch_bot" to allow changes in the patch spacing and density.

# ============================================================================
# Add custom settings
# ============================================================================
(df = data.frame(
    colony_diam = c(20),
    depth_m = c(3),
    flow_rate = c(0.2)))

# ============================================================================
# Run the Model
# ============================================================================
# Execute the model using one of the scenarios. Replace "test" with your chosen scenario.
# Available options: 'default' or 'test'

for(i in 1:nrow(df)) {
  input1 <- paste(as.character(unname(df[i,])), collapse="_")
  input_df <- df[i,]
  run_model(scenario = "test", store = T, input1 = input1, input_df = input_df,
            colony_diam = df[i, 'colony_diam'], den_cell = df[i, 'den_cell'], int_col_spac = df[i, 'int_col_spac'], patch_type = df[i, 'patch_type'],
            depth_m = df[i, 'depth_m'], flow_rate = df[i, 'flow_rate'], no_height = df[i, 'no_height'], no_width = df[i, 'no_height'],
            tb = df[i, 'tb'], fe = df[i, 'fe'], E0 = df[i, 'E0'], patch_density  = df[i, 'patch_density'], spe_speed  = df[i, 'spe_speed'],
            long_const  = df[i, 'long_const'], bundlebreak = df[i, 'bundlebreak'], bundlebreak_sd = df[i, 'bundlebreak_sd'],
            polyp_den = df[i, 'polyp_den'], egg_dia = df[i, 'egg_dia'])  # Choose 'default', or 'test'.
}

# ============================================================================
# Note:
# After running the model, results and logs will be generated based on the scenario.
# For detailed insights, refer to the output files created by the model execution.
# ============================================================================



