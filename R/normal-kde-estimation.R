kdeF <- function(x){

  bw <- stats::bw.nrd0(x)
  lims <- range(x) + c(-1, 1) * 3 * bw

  x_seq <- seq(from = lims[1], to = lims[2], length = 512)
  Fhat <- rowMeans(stats::pnorm(q = outer(x_seq, x, "-"), mean = 0, sd = bw))

  return(list("x" = x_seq,
              "Fhat" = Fhat,
              "bw" = bw))
}

sd_floor <- function(x){
  s <- stats::sd(x)
  if(s == 0){
    return(1e-10)
  }else{
    return(s)
  }
}
