plot_fun3 = function(x, i=2, map){
  if(map == 'sperm')
  {
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    levelplot(t(apply(x[,,i], 2, rev)), col.regions = coul, xlab = "Transverse",
              ylab = "Longitudinal", main = "Concentration (cells/m^3)")
  }
  
  if(map == 'egg')
  {
    coul <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    levelplot(t(apply(x[,,i], 2, rev)), col.regions = coul, xlab = "Transverse",
              ylab = "Longitudinal", main = "Concentration (cells/m^3)")
  }
}