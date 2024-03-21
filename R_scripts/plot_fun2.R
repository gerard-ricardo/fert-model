plot_fun2 = function(x, i=2){
  coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
  levelplot(t(apply(x[,,i], 2, rev)), col.regions = coul2, xlab = "Transverse",
            ylab = "Longitudinal", main = "Conc. (cells/m^3)")}
