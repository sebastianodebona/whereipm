#' Test Vital rate function
#' 
#' Creates dummy (default) parameters and visually test vital rate functions
#' 
#' @export
#' 
#' @param type Defines which functions to be tested: S, survival; R, reproduction; G, movement.
#' 

test.vital.fun <- function(type){
  
  kp <- create.kernel.params()
  vp <- mock.vital.params()
  
  if(type=='S'){
    
    S <- diag(sx(density=c(10, rep(0, 100-1))/1, 
            s.int=vp$s.int, 
            s.slo=vp$s.slo))
    
    image(t(S))
    return(invisible(S))
    
  } else if (type=='R'){
    
    R <- diag(rx(density=c(10, rep(0, 100-1))/1,
                 r.int=vp$r.int,
                 r.slo=vp$r.slo))
    
    image(t(R))
    return(invisible(R))
    
  } else if (type=='G'){
    
    b <- kp$L+c(0:kp$n.bins)*(kp$U-kp$L)/kp$n.bins # meshpoints
    y <- 0.5*(b[1:kp$n.bins]+b[2:(kp$n.bins+1)])
    
    G <- outer(y,y, FUN=gx_rw, p.stat=vp$p.stat, 
               s.stat=vp$s.stat, 
               s.mob=vp$s.mob)
    image(t(G))
    return(invisible(G))
    
  }
    
}