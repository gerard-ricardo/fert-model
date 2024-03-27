library(shiny)

# Assume settings_custom initially has some default values
# settings_custom  <- list(
#   grid_size = 1,
#   sub_steps = 1,
#   no_timestep = 3000,
#   time = 1,
#   polarbody = 0,
#   bundlebreak = 400,
#   bundlebreak_sd = 100,
#
#   patch_type = "2021digit",
#   int_col_spac = 1,
#   egg_viab = 0.92,
#   colony_diam = 27,
#   polyp_den = 80.6,
#   fecun_zone = 0.7,
#   bund_asc = 0.008,
#   bund_asc_egg = 0.0027,
#   patch_density = 4,
#
#   Ks = 0.27,
#   flow_rate = 0.10,
#   van_K_c = 0.41,
#   long_const = 0.20,
#   trans_const = 0.135,
#   vert_const = 0.067,
#   depth_m = 4,
#   trans_dim = 10,
#
#   tb = 0.26,
#   fe = 0.0037,
#   E0 = 0.0364,
#   spe_speed = 0.35,
#   spe_speed_sd = 0.049,
#   egg_dia = 0.517,
#   egg_dia_sd = 0.005,
#   eggperb = 6,
#   speperb = 2.1e6
# )

# # UI definition
# ui <- fluidPage(
#   titlePanel("Customize Model Settings"),
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("patch_type", "Patch Type:",
#                   choices = c("2022ahyac", "2021tenuis", "2021digit"), selected = "2021digit"),
#       #numericInput("grid_size", "Grid Size:", value = 1),
#       #numericInput("sub_steps", "Sub Steps:", value = 1),
#       numericInput("no_timestep", "Number of Timesteps:", value = 3000),
#       numericInput("time", "Time:", value = 1),
#       numericInput("polarbody", "Polar Body:", value = 0),
#       numericInput("bundlebreak", "Bundle Break:", value = 400),
#       numericInput("bundlebreak_sd", "Bundle Break SD:", value = 100),
#
#       numericInput("int_col_spac", "Inter-Colony Space:", value = 1),
#       numericInput("egg_viab", "Egg Viability:", value = 0.92),
#       numericInput("colony_diam", "Colony Diameter:", value = 27),
#       numericInput("polyp_den", "Polyp Density:", value = 80.6),
#       numericInput("fecun_zone", "Fecundity Zone:", value = 0.7),
#       numericInput("bund_asc", "Bundle Ascent:", value = 0.008),
#       numericInput("bund_asc_egg", "Bundle Ascent Egg:", value = 0.0027),
#       numericInput("patch_density", "Patch Density:", value = 4),
#
#       numericInput("Ks", "Ks:", value = 0.27),
#       numericInput("flow_rate", "Flow Rate:", value = 0.10),
#       numericInput("van_K_c", "Van K Constant:", value = 0.41),
#       numericInput("long_const", "Longitudinal Constant:", value = 0.20),
#       numericInput("trans_const", "Transverse Constant:", value = 0.135),
#       numericInput("vert_const", "Vertical Constant:", value = 0.067),
#       numericInput("depth_m", "Depth (m):", value = 4),
#       numericInput("trans_dim", "Transverse Dimension:", value = 10),
#
#       numericInput("tb", "tb:", value = 0.26),
#       numericInput("fe", "fe:", value = 0.0037),
#       numericInput("E0", "E0:", value = 0.0364),
#       numericInput("spe_speed", "Spawning Speed:", value = 0.35),
#       numericInput("spe_speed_sd", "Spawning Speed SD:", value = 0.049),
#       numericInput("egg_dia", "Egg Diameter:", value = 0.517),
#       numericInput("egg_dia_sd", "Egg Diameter SD:", value = 0.005),
#       numericInput("eggperb", "Eggs per Bundle:", value = 6),
#       numericInput("speperb", "Sperm per Bundle:", value = 2.1e6),
#
#       actionButton("updateSettings", "Update Settings")
#     ),
#     mainPanel(
#       verbatimTextOutput("currentSettings")
#     )
#   )
# )

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-image: url('https://raw.githubusercontent.com/gerard-ricardo/fert-model/main/scratch/PA190103.jpg'); /* Set the background image */
        background-size: cover; /* Cover the entire page */
        background-position: center center; /* Center the background image */
        background-attachment: fixed; /* Fix the background image position */
      }
      .time-parameters { background-color: #ADD8E6; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .patch-parameters { background-color: #90EE90; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .flow-parameters { background-color: #FF9999; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .kinetics-parameters { background-color: #9B59B6; padding: 20px; border-radius: 5px; margin-bottom: 20px; opacity: 0.9; }
      .parameter-heading { font-weight: bold; color: #FFFFFF; }
      .shiny-input-container { width: 100%; }
    "))
  ),
  titlePanel("Customize Model Settings"),
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
          numericInput("tb", "tb:", value = 0.26),
          numericInput("fe", "fe:", value = 0.0037),
          numericInput("E0", "E0:", value = 0.0364),
          numericInput("spe_speed", "Spawning Speed:", value = 0.35),
          numericInput("spe_speed_sd", "Spawning Speed SD:", value = 0.049),
          numericInput("egg_dia", "Egg Diameter:", value = 0.517),
          numericInput("egg_dia_sd", "Egg Diameter SD:", value = 0.005),
          numericInput("eggperb", "Eggs per Bundle:", value = 6),
          numericInput("speperb", "Sperm per Bundle:", value = 2.1e6)
      ),

      actionButton("updateSettings", "Update Settings"),
      width = 4
    ),
    mainPanel(
      verbatimTextOutput("currentSettings"),
      width = 8
    )
  )
)


# Server logic
server <- function(input, output, session) {
  output$interColSpaceMessage <- renderUI({
    if (input$patch_type %in% c("2022ahyac", "2021tenuis", "2021digit")) {
      return(tags$div(style = "color: red;", "Inter-Colony Space or Patch Density cannot be altered for this patch configuration."))
    } else {
      return(NULL)  # No message for other patch types
    }
  })
  observeEvent(input$updateSettings, {
    # Update settings_custom with the input values
    settings_custom <<- list(
      #grid_size = input$grid_size,
      #sub_steps = input$sub_steps,
      no_timestep = input$no_timestep,
      time = input$time,
      polarbody = input$polarbody,
      bundlebreak = input$bundlebreak,
      bundlebreak_sd = input$bundlebreak_sd,

      patch_type = input$patch_type,
      int_col_spac = input$int_col_spac,
      egg_viab = input$egg_viab,
      colony_diam = input$colony_diam,
      polyp_den = input$polyp_den,
      fecun_zone = input$fecun_zone,
      bund_asc = input$bund_asc,
      bund_asc_egg = input$bund_asc_egg,
      patch_density = input$patch_density,

      Ks = input$Ks,
      flow_rate = input$flow_rate,
      van_K_c = input$van_K_c,
      long_const = input$long_const,
      trans_const = input$trans_const,
      vert_const = input$vert_const,
      depth_m = input$depth_m,
      trans_dim = input$trans_dim,

      tb = input$tb,
      fe = input$fe,
      E0 = input$E0,
      spe_speed = input$spe_speed,
      spe_speed_sd = input$spe_speed_sd,
      egg_dia = input$egg_dia,
      egg_dia_sd = input$egg_dia_sd,
      eggperb = input$eggperb,
      speperb = input$speperb
    )


    # Display the updated settings_custom list
    output$currentSettings <- renderPrint({ settings_custom })
  })
}

# Run the application
shinyApp(ui = ui, server = server)

