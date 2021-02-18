library(mcclust)
library(mcclust.ext)
library(partitions)

parts = t(partitions::setparts(5))
# mcclust::vi.dist(parts[4,], parts[1,])

vi_mat = matrix(0, nrow(parts), nrow(parts))
for(i in 2:nrow(parts))
{
  for(j in 1:(i-1))
  {
    vi_mat[i,j] = mcclust::vi.dist(parts[i,], parts[j,])
    vi_mat[j,i] = vi_mat[i,j] 
  }
}


### CASO ASIMMETRICO 
p_c = function(c, c0, psi_inf_eq, psi_sup)
{
  psi = psi_inf_eq
  id_c0 = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c0)==1 ) )
  id_c = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c)==1 ) )
  if( length(unique(c)) > length(unique(c0)) ) psi = psi_sup
  exp(- psi * vi_mat[id_c, id_c0] )
}

psi_seq = seq(0, 3, by = 0.02)
length(psi_seq)



plot_partitions = function(row_c0, psi_inf_eq_seq, psi_sup_seq)
{
  if(length(psi_inf_eq_seq)!=length(psi_sup_seq)) return("le sequenze devono avere uguale lunghezza \n")
  
  res = matrix(NA, nrow(parts), length(psi_inf_eq_seq))
  for(i in 1:length(psi_inf_eq_seq))
  {
    tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[row_c0,], psi_inf_eq = psi_inf_eq_seq[i], psi_sup =  psi_sup_seq[i]) )
    res[,i] = tmp/sum(tmp)
  }
  
  par(mar = c(4, 4.7, 1, 4) + 0.1, xpd = NA)
  plot(0, xlim = range(psi_inf_eq_seq), ylim = c(0,1), type = "n", axes = F,  
       xaxs = "i", yaxs = "i", xlab = expression(psi[1]), ylab = "")
  axis(1, tck=-0.01)
  axis(4, tck=-0.01, las = 1)
  mtext("Cumulative probability", side = 4, line = 3, cex = par("cex.lab"))
  res2 = res
  for(i in 1:length(psi_inf_eq_seq))
  {
    res2[,i] = cumsum(res[,i])
  }
  
  cbind(res2[,1], parts, apply(parts, 1, function(x) length(unique(x))))
  axis(2, at = res2[,1][c(1,16,41,51)], labels = c("","","",""), tick = F)
  print(res2[,1][c(1,16,41,51)])
  axis(2, at = c(0.013,0.15,0.55,0.88,0.993), labels = c("1 block", "2 blocks", "3 blocks", "4 blocks", "5 blocks"), 
       tick = F, las = 1, line = 0.5)
  
  polygon(x = c(psi_inf_eq_seq, rev(psi_inf_eq_seq)),
          y = c( res2[row_c0-1,], rev(res2[row_c0,]) ),
          col = "#BDEBF2", border = NA)
  for(i in 1:52)
  {
    lines( smooth.spline(x = psi_inf_eq_seq, y = res2[i,]) )
  }
  box()
  segments(0, c(0, res2[,1][c(1,16,41,51)], 1), -0.1, c(0, res2[,1][c(1,16,41,51)], 1), lty = 3)
}


### c0 = (1,2) (3,4,5)
# riga 16
plot_partitions(16, psi_seq, psi_seq)
plot_partitions(16, psi_inf_eq_seq = psi_seq, psi_sup_seq = 1.5*psi_seq)


### c0 = (1,2) (3,4) (5)
# riga 32
plot_partitions(32, psi_seq, psi_seq)
plot_partitions(32, psi_inf_eq_seq = psi_seq, psi_sup_seq = 1.5*psi_seq)
