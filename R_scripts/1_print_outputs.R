library(ggplot2)
library(lattice)
library(gridExtra)
print_output_and_save_plots <- function(int_col_spac, no_timestep, smallmat, tot_spe_indiv, spemap3_3, eggmap3_4, surface, spemap_temp3_3,
                                        fertcounter, embmap_store, new_folder_path, execution_mode = "single", input_id = NULL, input_df = NULL,
                                        time_diff = NULL) {
  file_name_suffix <- if(execution_mode == "parallel" && !is.null(input_id)) {
    paste0("_", input_id)
  } else {
    ""
  }
  print(paste("Input ID:", input_id))
  print(paste("File name suffix:", file_name_suffix))
  print(paste("Time taken:", time_diff))
  save.image(file = file.path(new_folder_path, paste0("last_run", file_name_suffix, ".Rdata")))
  p0 = plot_fun1(spemap3_3, surface)
  p1 = plot_fun2(eggmap3_4, surface)
  coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  p2 = levelplot(log10(t(apply(spemap_temp3_3[[no_timestep/2]][,,2], 2, rev))/10^6), col.regions = coul, at = c(0, 1, 2, 3, 3.5, 4, 5, 6), labels = T, cuts  = 5,
                 contour = T, region = T,xlab = "Transverse", ylab = "Longitudinal", main = "Sperm conc. at mid-time \n(log10(sperm/mL))")
  p3 = levelplot(log10(t(apply(spemap_temp3_3[[no_timestep]][,,2], 2, rev))/10^6), col.regions = coul, at = c(0, 1, 2, 3, 3.5, 4, 5, 6), labels = T, cuts  = 5,
                 contour = T, region = T,xlab = "Transverse", ylab = "Longitudinal",  main = "Final sperm conc.\n(log10(sperm/mL))")
  coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  p5 <- levelplot(t(apply(smallmat, 2, rev)), col.regions = coul,
                  xlab = "Transverse", ylab = "Longitudinal",
                  main = list(label = "Conc. (cells/m^3)", cex = 1),
                  colorkey = FALSE,
                  panel = function(x, y, ...) {
                    panel.levelplot(x, y, ...)
                    panel.grid(h = -1, v = -1, col = "black")
                  })
  lay <- rbind(c(0,0,1,1,2,2,3,3),
               c(0,0,1,1,2,2,3,3),
               c(0,0,1,1,2,2,3,3),
               c(0,0,1,1,2,2,3,3),
               c(0,0,1,1,2,2,3,3),
               c(0,0,1,1,2,2,3,3)
  )
  lay <- rbind(lay,
               c(4, 4, 4, 4, 4, 4, 4, 4),
               c(4, 4, 4, 4, 4, 4, 4, 4),
               c(4, 4, 4, 4, 4, 4, 4, 4),
               c(4, 4, 4, 4, 4, 4, 4, 4),
               c(4, 4, 4, 4, 4, 4, 4, 4))
  gs = list(p0, p1, p2, p3, p5)
  plots1 = grid.arrange(grobs = gs, layout_matrix = lay)

  ####fert counter####
  fert_df = data.frame(time = 1:(no_timestep), fert = fertcounter)
  tail(fert_df)
  max_fert = max(fert_df$fert, na.rm = T)
  input_df$max_fert = max_fert
  save(input_df, file = paste0(new_folder_path, "/input", file_name_suffix, ".Rdata"))
  filename <- file.path(new_folder_path, paste0("fertcounter_", ".RData"))
  source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
  p4 <- ggplot() +
    geom_point(fert_df,
               mapping = aes(x = time, y = fert), position = position_jitter(width = .00, height = .00),
               alpha = 0.50, size = 3) +
    theme_sleek2()
  p4 <- p4 + annotate("text",
                      x = max(fert_df$time, na.rm = T) * 0.3, y = max(fert_df$fert, na.rm = T) * 0.5,
                      label = paste("Max fert:", round(max_fert * 100, 2), "%"),
                      hjust = 1, vjust = -1)
  p4 <- p4 + annotate("text",
                      x = max(fert_df$time, na.rm = T) * 0.3, y = max(fert_df$fert, na.rm = T) * 0.4,
                      label = paste("Total embryos:", round(sum(embmap_store[[no_timestep - 1]][, , 2], na.rm = T), 0)),
                      hjust = 1, vjust = -1)
  p4 <- p4 + labs(
    x = expression(Time ~ (s)),
    y = expression(Cumulative ~ fert. ~ success ~ (prop.)))
  ggsave(file = paste0("clouds", file_name_suffix, ".pdf"), plot  = plots1, width = 10, height = 8, path = new_folder_path)
  ggsave(file = paste0("fert_counter", file_name_suffix, ".pdf"), plot = p4, width = 10, height = 8, path = new_folder_path)
}
