
# Fertmod

<!-- badges: start -->
<!-- badges: end -->

This project simulates coral fertilisation success. The core of Fertmod is built on R scripts designed for research and academic purposes, focusing on understanding and predicting coral fertilisation outcomes under various environmental and biological scenarios.

## Project Structure

The project comprises several R scripts, each serving a distinct role in the simulation process:
- **Settings**: Defines the model settings.
- **Print Outputs**: Handles output generation and printing.
- **Backend Logic**: Contains the core simulation logic.

Users can select from predefined scenarios or customize settings via a provided Shiny app interface.

### Running the Simulation

1. **Download the repository** to your local machine.

2. **Open Front_end.R script**

3. **Source the required scripts** in your R environment. Customized scenarios can be set via the Shiny app in the 'custom_settings' folder.

   ```r
   source("./R_scripts/settings.R")
   source("./R_scripts/1_print_outputs.R")
   source("./R_scripts/Backend.R")
   
   # Available options: 'default', 'test', 'custom'.
    run_model(scenario = "test")
   
   ```
4. **Outputs** will be in the 'out' folder


## Note

This is a 'cleaned script' with reduced functionality. 



