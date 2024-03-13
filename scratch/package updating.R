
#copy code to Fertmod_comments and then run remove code




# package updating and checking -------------------------------------------


## After developing functions, document the package
devtools::document()  #update functions

devtools::test()

devtools::check()


devtools::install("C:/Users/gerar/OneDrive/1_Work/R/ECx", upgrade = 'never')
