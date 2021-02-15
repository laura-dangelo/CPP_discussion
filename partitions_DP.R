library(mcclust)
# library(mcclust.ext)
library(partitions)

parts = t(partitions::setparts(5))
# mcclust::vi.dist(parts[4,], parts[1,])


DP_prob = function(c, alpha)
{
  N = length(c)
  K = length(unique(c))
  lambdas = sapply(1:K, function(x) sum(c==x))
  out = K * log(alpha) - lfactorial(alpha + N - 1) + lfactorial(alpha - 1) + sum( lfactorial(lambdas - 1) )
  return(exp(out))
}
prob_parts = apply(parts, 1, function(x) DP_prob(x, alpha = 1))
sum(prob_parts)

vi_mat = matrix(0, nrow(parts), nrow(parts))
for(i in 2:nrow(parts))
{
  for(j in 1:(i-1))
  {
    vi_mat[i,j] = mcclust::vi.dist(parts[i,], parts[j,])
    vi_mat[j,i] = vi_mat[i,j] 
  }
}


p_c_DP = function(c, c0, psi_inf_eq, psi_sup, alpha)
{
  psi = psi_inf_eq
  id_c0 = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c0)==1 ) )
  id_c = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c)==1 ) )
  if( length(unique(c)) > length(unique(c0)) ) psi = psi_sup
  exp(- psi * vi_mat[id_c, id_c0] ) * DP_prob(c, alpha)
}

# p_c_DP(parts[1,], parts[3,], psi_inf_eq = 1, psi_sup = 1.4,  alpha = 1)

psi_seq = seq(0, 3, by = 0.02)
length(psi_seq)



plot_partitions = function(row_c0, psi_inf_eq_seq, psi_sup_seq)
{
  if(length(psi_inf_eq_seq)!=length(psi_sup_seq)) return("le sequenze devono avere uguale lunghezza \n")
  
  res = matrix(NA, nrow(parts), length(psi_inf_eq_seq))
  for(i in 1:length(psi_inf_eq_seq))
  {
    tmp = apply( parts, 1, function(x) p_c_DP(c = x, c0 = parts[row_c0,], psi_inf_eq = psi_inf_eq_seq[i], psi_sup =  psi_sup_seq[i], alpha = 1) )
    res[,i] = tmp/sum(tmp)
  }
  
  par(mar = c(5, 2, 4, 3) + 0.1)
  plot(0, xlim = range(psi_inf_eq_seq), ylim = c(0,1), type = "n", axes = F,  xaxs = "i", yaxs = "i", xlab = "", ylab = "")
  axis(1, tck=-0.01)
  axis(4, tck=-0.01, las = 1)
  res2 = res
  for(i in 1:length(psi_inf_eq_seq))
  {
    res2[,i] = cumsum(res[,i])
  }
  
  cbind(res2[,1], parts, apply(parts, 1, function(x) length(unique(x))))
  axis(2, at = res2[,1][c(1,16,41,51)], labels = c("","","",""))
  axis(2, at = c(0.1,0.42,0.77,0.95,0.999), labels = c("1", "2", "3", "4", "5"), tick = F, las = 1)
  
  
  polygon(x = c(psi_inf_eq_seq, rev(psi_inf_eq_seq)),
          y = c( res2[row_c0-1,], rev(res2[row_c0,]) ),
          col = "#BDEBF2", border = NA)
  for(i in 1:52)
  {
    lines( smooth.spline(x = psi_inf_eq_seq, y = res2[i,]) )
  }
  box()
  
}



### c0 = (1,2) (3,4,5)
# riga 16
par(mar = c(5, 2, 4, 3) + 0.1)
plot_partitions(16, psi_seq, psi_seq)

##### asimmetrico
plot_partitions(16, psi_seq, psi_seq*2) # penalizzo di più quelle con più clusters
plot_partitions(16, 2*psi_seq, psi_seq) # penalizzo di più quelle con meno (o uguale??? numero di clusters) 



### c0 = (1,2) (3,4) (5)
# riga 32
par(mar = c(5, 2, 4, 3) + 0.1)
plot_partitions(32, psi_seq, psi_seq)

##### asimmetrico
plot_partitions(32, psi_seq, psi_seq*2) # penalizzo di più quelle con più clusters
plot_partitions(32, 2*psi_seq, psi_seq) # penalizzo di più quelle con meno (o uguale??? numero di clusters) 
