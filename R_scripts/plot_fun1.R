plot_fun1 = function(x, i=2){

    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    levelplot(t(apply(x[,,i], 2, rev)), col.regions = coul, xlab = "Transverse",
            ylab = "Longitudinal", main = "Conc. (cells/m^3)")

}
