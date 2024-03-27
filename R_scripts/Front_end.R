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
#   Alternatively, you can manually edit the 'settings_custom' within 'settings.R' and save.

# ============================================================================
# Run the Model
# ============================================================================
# Execute the model using one of the scenarios. Replace "test" with your chosen scenario.
# Available options: 'default', 'test', 'custom'.

run_model(scenario = "custom")  # Choose 'default', 'test', or 'custom'.

# ============================================================================
# Note:
# After running the model, results and logs will be generated based on the scenario.
# For detailed insights, refer to the output files created by the model execution.
# ============================================================================



