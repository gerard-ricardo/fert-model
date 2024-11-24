
# Fertmod

<!-- badges: start -->
<!-- badges: end -->

This project simulates coral fertilisation success. The core of Fertmod is built on R scripts designed for research and academic purposes, focusing on understanding and predicting coral fertilisation outcomes under various environmental and biological scenarios.

## Project Structure

The project comprises several R scripts, each serving a role in the simulation process:
- **Settings**: Defines the model initial settings.
- **Backend**: Contains the core model structure.
- **Print Outputs**: Handles output generation and printing.

Users can select from predefined scenarios or customize settings manually or via a provided Shiny app interface.

### Running the Model via the Frontend script

1. **Download the repository** to your local machine.

2. **Open Front_end.R script** in the R_scripts folder.

3. **Name the output folder** in 'dir_name'. Customized scenarios can be set via the Shiny app in the 'custom_settings' folder.

4. **Customize the input variables** if needed. This over-rides any initial settings.

5. **Set the scenario** in the run_model function

   ```r
   # Available options: 'default', or 'test'.
    run_model(scenario = "test",...)
   
   ```
6. **Run** the whole code.

7. **Outputs** Fertilisation success for each time step will be shown in the console. Further outputs and plots will be in the 'out' folder

### Running the Model via the Shiny-app frontend

1. **Open the app.R** 

2. **Click RunApp** 

3. **Customize the input variables** if needed. This over-rides any initial settings.

4. **Run model**

5. **Outputs** Fertilisation success for each time step will be shown in the console. Input setting will be saved as a .rds file in the 'out' folder. No plots will be saved using the shiny-app.


## Note

 - This is a 'cleaned script' with reduced functionality. 
 - Typical run times are ~2.5 minutes for a test run, and 3-5 hours for a default run. 



![Screenshot of Shiny App](https://raw.githubusercontent.com/gerard-ricardo/fert-model/main/scratch/Screenshot%202024-03-27%20092933.png
)


## System Information

### R Version and Platform
- **R Version:** 4.4.1 (2024-06-14 ucrt)  
- **Platform:** x86_64-w64-mingw32/x64  
- **Operating System:** Windows 11 x64 (build 22631)  

### Locale
- `LC_COLLATE=English_Australia.utf8`  
- `LC_CTYPE=English_Australia.utf8`  
- `LC_MONETARY=English_Australia.utf8`  
- `LC_NUMERIC=C`  
- `LC_TIME=English_Australia.utf8`  

### Attached Base Packages
- `stats`  
- `graphics`  
- `grDevices`  
- `utils`  
- `datasets`  
- `methods`  
- `base`  

### Other Attached Packages
- `truncnorm_1.0-9`  
- `RColorBrewer_1.1-3`  
- `dplyr_1.1.4`  
- `abind_1.4-5`  
- `gridExtra_2.3`  
- `lattice_0.22-6`  
- `ggplot2_3.5.1`  

### Loaded Namespaces
- `utf8_1.2.3`  
- `generics_0.1.3`  
- `stringi_1.7.12`  
- `digest_0.6.36`  
- `magrittr_2.0.3`  
- `grid_4.4.1`  
- `pkgload_1.3.4`  
- `fastmap_1.2.0`  
- `pkgbuild_1.4.4`  
- `sessioninfo_1.2.2`  
- `urlchecker_1.0.1`  
- `promises_1.2.0.1`  
- `purrr_1.0.1`  
- `fansi_1.0.4`  
- `scales_1.3.0`  
- `textshaping_0.4.0`  
- `cli_3.6.3`  
- `shiny_1.8.1.1`  
- `crayon_1.5.2`  
- `rlang_1.1.1`  
- `ellipsis_0.3.2`  
- `munsell_0.5.1`  
- `remotes_2.5.0`  
- `withr_3.0.0`  
- `cachem_1.1.0`  
- `devtools_2.4.5`  
- `tools_4.4.1`  
- `memoise_2.0.1`  
- `colorspace_2.1-0`  
- `httpuv_1.6.9`  
- `vctrs_0.6.5`  
- `R6_2.5.1`  
- `mime_0.12`  
- `lifecycle_1.0.4`  
- `stringr_1.5.0`  
- `fs_1.6.4`  
- `htmlwidgets_1.6.4`  
- `usethis_2.2.3`  
- `miniUI_0.1.1.1`  
- `ragg_1.3.2`  
- `pkgconfig_2.0.3`  
- `pillar_1.9.0`  
- `later_1.3.0`  
- `gtable_0.3.5`  
- `glue_1.6.2`  
- `profvis_0.3.8`  
- `Rcpp_1.0.12`  
- `systemfonts_1.1.0`  
- `tibble_3.2.1`  
- `tidyselect_1.2.1`  
- `rstudioapi_0.16.0`  
- `farver_2.1.1`  
- `xtable_1.8-4`  
- `htmltools_0.5.8.1`  
- `labeling_0.4.3`  
- `compiler_4.4.1`  
