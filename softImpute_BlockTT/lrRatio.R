lrRatio <- function(u_old, d_old, v_old, u, d, v)
{
  #Compute ratio for stopping criterion
  rat = 1
  
  rat = 
    sqrt(sum((u - u_old)^2)) / max( sqrt(sum(u_old^2)) , sqrt(sum(u^2)) ) 
    
    
    # 0.5 * (
    # sqrt(sum((u - u_old)^2)) / max( sqrt(sum(u_old^2)) , sqrt(sum(u^2)) ) + 
    # sqrt(sum((d_old - d)^2)) / max( sqrt(sum(d_old^2)) , sqrt(sum(d^2)) ) )
  
 
  return(rat) 
}