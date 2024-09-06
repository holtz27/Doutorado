# Função para estimar a densidade em x0 usando a função density do R
estimador_densidade <- function(x, x0, bw = "nrd0", kernel = "gaussian") {
  # Calcula a densidade usando a função density do R
  densidade <- density(x, bw = bw, kernel = kernel)
  
  # Aproxima a densidade no ponto x0
  densidade_x0 <- approx(densidade$x, densidade$y, xout = x0)$y
  
  return(densidade_x0)
}
