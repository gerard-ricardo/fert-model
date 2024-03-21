print_output_and_save_plots <- function(int_col_spac, no_timestep, spemap3_3, eggmap3_4, surface, spemap_temp3_2, fertcounter, embmap_store, new_folder_path) {

  save.image(file = file.path(path = new_folder_path, paste0("last_run_",  ".Rdata")))
  #pdf(paste0('spemap_', 'int_col', int_col_spac,".pdf"));plot_fun1(spemap3_3, surface);dev.off() #last sperm map
  #pdf(paste0('eggmap_', 'int_col', int_col_spac,".pdf"));plot_fun2(eggmap3_4, surface);dev.off() #last sperm map
  p0 = plot_fun1(spemap3_3, surface)
  p1 = plot_fun2(eggmap3_4, surface)

  #sperm/mL plot
  #pdf(paste0('sperm_mL_', 'int_col', int_col_spac,".pdf"), width = 10, height = 8)
  coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  p2 = levelplot(log10(t(apply(spemap_temp3_2[[no_timestep]][,,2], 2, rev))/10^6), col.regions = coul, at = c(0, 1, 2, 3, 3.5, 4, 5, 6), labels = T, cuts  = 5,
                 contour = T, region = T,xlab = "Transverse", ylab = "Longitudinal", main = "Final sperm concentration (log10(sperm/mL))")

  #dev.off()
  #ggsave("sperm_mL.png", height = 8, width = 12)  #check png via Files tab
  #ggsave(file = paste0("sperm_mL", 'int_col', int_col_spac,".pdf"), width = 10, height = 8)

  lay <- rbind(c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2)
  )
  gs = list(p0, p1, p2)
  plots1 = grid.arrange(grobs = gs, layout_matrix = lay)


  ####fert counter####
  fert_df = data.frame(time = 1:(no_timestep), fert = fertcounter)
  tail(fert_df)
  max_fert = max(fert_df$fert, na.rm = T)
  #df_name <- paste0("fert_df", vec_x[i])
  #assign(df_name, fert_df)
  filename <- file.path(new_folder_path, paste0("fertcounter_", ".RData")) #set name and path to save df
  #save(get(df_name), file = filename)  #save df
  #save(list = df_name, file = filename)
  #plot(fert_df$time, fert_df$cumfert)
  source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code
  p4 <- ggplot() +
    geom_point(fert_df,
               mapping = aes(x = time, y = fert), position = position_jitter(width = .00, height = .00),
               alpha = 0.50, size = 3) +
    theme_sleek2()
  # p4 = p4 + scale_y_continuous( limits = c(0, 1))
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
  # p1
  # ggsave("fert_counter.png", height = 8, width = 12)  #check png via Files tab

  # saved files
  ggsave(file = paste0("clouds_", ".pdf"), plots1, width = 10, height = 8, path = new_folder_path)
  ggsave(p4, file = paste0("fert_counter_", ".pdf"), width = 10, height = 8, path = new_folder_path)
  # savehistory(file.path(new_folder_path, "my_history.R"))

}
