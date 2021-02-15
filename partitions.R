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

p_c = function(c, c0, psi)
{
  id_c0 = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c0)==1 ) )
  id_c = which(sapply(1:nrow(parts), function(x) prod(parts[x,] == c)==1 ) )
  exp(- psi * vi_mat[id_c, id_c0] )
}

psi_seq = seq(0, 3, by = 0.02)
length(psi_seq)

### c0 = 1 blocco
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[1,], psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F)
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}


### c0 = (1,2) (3,4,5)
# riga 16
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[16,], psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F)
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}


### c0 = (1,2) (3,4) (5)
# riga 32
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[32,], psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F)
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}


### c0 = (1) (2) (3) (4,5)
# riga 51
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[51,], psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F)
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}




################################################################################
############### CASO ASIMMETRICO 



vi_mat = matrix(0, nrow(parts), nrow(parts))
for(i in 2:nrow(parts))
{
  for(j in 1:(i-1))
  {
    vi_mat[i,j] = mcclust::vi.dist(parts[i,], parts[j,])
    vi_mat[j,i] = vi_mat[i,j] 
  }
}

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

### c0 = (1,2) (3,4,5)
# riga 16
par(mfrow=c(1,3))
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[16,], psi_inf_eq = psi_seq[i], psi_sup = psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"="~psi[sup]) )
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}
polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[15,], rev(res2[16,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)


res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[16,], psi_inf_eq = psi_seq[i], psi_sup = 2 * psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"= 0.5 "~psi[sup]))
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[15,], rev(res2[16,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)


res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[16,], psi_inf_eq = psi_seq[i], psi_sup = .5 * psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"= 2 "~psi[sup]))
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}
polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[15,], rev(res2[16,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)








### c0 = (1,2) (3,4) (5)
# riga 32
par(mfrow=c(1,3))
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[32,], psi_inf_eq = psi_seq[i], psi_sup = psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"="~psi[sup]) )
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}
polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[31,], rev(res2[32,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)


res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[32,], psi_inf_eq = psi_seq[i], psi_sup = 2 * psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"= 0.5 "~psi[sup]) )
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}
polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[31,], rev(res2[32,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)


res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[32,], psi_inf_eq = psi_seq[i], psi_sup = .5 * psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"= 2 "~psi[sup]) )
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}
polygon(x = c(psi_seq, rev(psi_seq)),
        y = c( res2[31,], rev(res2[32,]) ),
        col = "#BDEBF2", border = NA)
for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)




### c0 = (1) (2) (3) (4,5)
# riga 51
res = matrix(NA, nrow(parts), length(psi_seq))
for(i in 1:length(psi_seq))
{
  tmp = apply( parts, 1, function(x) p_c(c = x, c0 = parts[51,], psi_seq[i]) )
  res[,i] = tmp/sum(tmp)
}
str(res)
head(res)
plot(0, xlim = range(psi_seq), ylim = c(0,1), type = "n", axes = F, main=expression(psi[inf]~"="~psi[sup]) )
axis(1)
axis(2, at = seq(0,1,length.out = 52)[c(1,2,17,42,51)], labels = c("","","","",""))
axis(2, at = seq(0,1,length.out = 52)[c(1,10,30,46,52)], labels = c("1", "2", "3", "4", "5"), tick = F)
axis(4)
box()
res2 = res
for(i in 1:length(psi_seq))
{
  res2[,i] = cumsum(res[,i])
}

for(i in 1:52)
{
  lines( smooth.spline(x = psi_seq, y = res2[i,]) )
}
segments(0,0,3,0)
