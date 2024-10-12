latex = function( metricas, names, multrow.name, align, digits ){
  
  dados = data.frame( Priori = c(paste0( "\\multirow{", 
                                         ncol( metricas ), 
                                         "}{*}{",
                                         multrow.name,
                                         "}"
                                         ), 
                                 "", 
                                 "" ),
                      Metricas = colnames( metricas ),
                      t( metricas ) 
  )
  row.names( dados ) = NULL
  
  # kbl
  latex = kableExtra::kbl( dados, 
                           format = 'latex',
                           booktabs = TRUE,
                           align = align,
                           col.names = names,
                           escape = FALSE,
                           digits = digits)
  return( latex )
}
