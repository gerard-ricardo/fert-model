library(shiny)
library(withr)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-image: url('https://raw.githubusercontent.com/gerard-ricardo/fert-model/main/scratch/PA190103.jpg');
        background-size: cover;
        background-position: center center;
        background-attachment: fixed;
      }
      .time-parameters { background-color: #ADD8E6; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .patch-parameters { background-color: #90EE90; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .flow-parameters { background-color: #FF9999; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .kinetics-parameters { background-color: #9B59B6; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .parameter-heading { font-weight: bold; color: #FFFFFF; }
      .shiny-input-container { width: 100%; }
    "))
  ),
  titlePanel("Customize and Run Fertilisation Model"),
  sidebarLayout(
    sidebarPanel(
      class = "wide-sidebar",
      selectInput("patch_type", "Patch Type:",
                  choices = c("2022ahyac", "2021tenuis", "2021digit", "patch_bot"), selected = "patch_bot"),
      uiOutput("interColSpaceMessage"),
      div(class = "time-parameters",
          h4("Time Parameters", class = "parameter-heading"),
          numericInput("no_timestep", "Number of Timesteps:", value = 500),
          numericInput("time", "Time:", value = 1),
          numericInput("polarbody", "Polar Body:", value = 0),
          numericInput("bundlebreak", "Bundle Break:", value = 400),
          numericInput("bundlebreak_sd", "Bundle Break SD:", value = 100)
      ),

      div(class = "patch-parameters",
          h4("Patch Parameters", class = "parameter-heading"),
          numericInput("int_col_spac", "Inter-Colony Space:", value = 2),
          numericInput("egg_viab", "Egg Viability:", value = 0.92),
          numericInput("colony_diam", "Colony Diameter:", value = 27),
          numericInput("polyp_den", "Polyp Density:", value = 80.6),
          numericInput("fecun_zone", "Fecundity Zone:", value = 0.7),
          numericInput("bund_asc", "Bundle Ascent:", value = 0.008),
          numericInput("bund_asc_egg", "Bundle Ascent Egg:", value = 0.0027),
          numericInput("patch_density", "Patch Density:", value = 9)
      ),

      div(class = "flow-parameters",
          h4("Flow Parameters", class = "parameter-heading"),
          numericInput("Ks", "Ks:", value = 0.27),
          numericInput("flow_rate", "Flow Rate:", value = 0.10),
          numericInput("van_K_c", "Van K Constant:", value = 0.41),
          numericInput("long_const", "Longitudinal Constant:", value = 0.20),
          numericInput("trans_const", "Transverse Constant:", value = 0.135),
          numericInput("vert_const", "Vertical Constant:", value = 0.067),
          numericInput("depth_m", "Depth (m):", value = 4),
          numericInput("trans_dim", "Transverse Dimension:", value = 10)
      ),

      div(class = "kinetics-parameters",
          h4("Kinetics Parameters", class = "parameter-heading"),
          numericInput("tb", "tb:", value = 0.005138657),
          numericInput("fe", "fe:", value = 0.002391663),
          numericInput("E0", "E0:", value = 0.003764767),
          numericInput("spe_speed", "Spawning Speed:", value = 0.35),
          numericInput("spe_speed_sd", "Spawning Speed SD:", value = 0.049),
          numericInput("egg_dia", "Egg Diameter:", value = 0.517),
          numericInput("egg_dia_sd", "Egg Diameter SD:", value = 0.005),
          numericInput("eggperb", "Eggs per Bundle:", value = 6),
          numericInput("speperb", "Sperm per Bundle:", value = 2.1e6)
      ),

      actionButton("runModel", "Run Model"),
      width = 4
    ),
    mainPanel(
      width = 8  # Keep the width for layout consistency, but no content
    )
  )
)

# Server logic
server <- function(input, output, session) {
  print(paste("Current working directory is:", getwd()))

  #
  project_root <- getwd()
  print(paste("Project root is:", project_root))

  # Source required scripts
  required_scripts <- c("boundary_0_func.R", "1_print_outputs.R", "Backend.R")
  lapply(required_scripts, function(script) {
    file_path <- file.path(project_root, "R_scripts", script)
    print(paste("Attempting to source script:", file_path))
    if (!file.exists(file_path)) {
      stop(paste("Required script not found:", file_path))
    }
    source(file_path)
    print(paste("Successfully sourced script:", file_path))
  })


  # Reactive data frame to hold settings
  df <- reactive({
    data.frame(
      colony_diam = input$colony_diam,
      depth_m = input$depth_m,
      flow_rate = input$flow_rate,
      int_col_spac = input$int_col_spac,
      patch_type = input$patch_type,
      tb = input$tb,
      fe = input$fe,
      E0 = input$E0,
      patch_density = input$patch_density,
      spe_speed = input$spe_speed,
      long_const = input$long_const,
      bundlebreak = input$bundlebreak,
      bundlebreak_sd = input$bundlebreak_sd,
      polyp_den = input$polyp_den,
      egg_dia = input$egg_dia,
      stringsAsFactors = FALSE
    )
  })

  # Warning for certain patch types
  output$interColSpaceMessage <- renderUI({
    if (input$patch_type %in% c("2022ahyac", "2021tenuis", "2021digit")) {
      return(tags$div(style = "color: red;", "Inter-Colony Space or Patch Density cannot be altered for this patch configuration."))
    } else {
      return(NULL)  # No message for other patch types
    }
  })

  # Log of current settings
  output$currentSettings <- renderPrint({
    df()
  })

  # Run model and log output
  observeEvent(input$runModel, {
    log <- ""
    for (i in 1:nrow(df())) {
      input1 <- paste(as.character(unname(df()[i,])), collapse = "_")
      input_df <- df()[i,]

      # Simulated run_model function (replace with your actual function)
      log <- paste0(log, "Running model with settings: ", toString(input_df), "\n")

      run_model(
        scenario = "test",
        store = TRUE,
        input1 = input1,
        input_df = input_df,
        colony_diam = df()[i, 'colony_diam'],
        depth_m = df()[i, 'depth_m'],
        flow_rate = df()[i, 'flow_rate'],
        int_col_spac = df()[i, 'int_col_spac'],
        patch_type = df()[i, 'patch_type'],
        tb = df()[i, 'tb'],
        fe = df()[i, 'fe'],
        E0 = df()[i, 'E0'],
        patch_density = df()[i, 'patch_density'],
        spe_speed = df()[i, 'spe_speed'],
        long_const = df()[i, 'long_const'],
        bundlebreak = df()[i, 'bundlebreak'],
        bundlebreak_sd = df()[i, 'bundlebreak_sd'],
        polyp_den = df()[i, 'polyp_den'],
        egg_dia = df()[i, 'egg_dia']
      )
    }
    output$modelLog <- renderText(log)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
