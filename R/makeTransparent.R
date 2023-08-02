makeTransparent<-function(someColor, alpha=100)
{
    stopifnot(exprs = {
      is.numeric(alpha)
    })
  
    newColor <- grDevices::col2rgb(someColor)
    return(apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)}))
}
