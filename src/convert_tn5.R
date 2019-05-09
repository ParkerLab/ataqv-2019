# for converting between numeric Tn5 and Factor Tn5
# Example usage:
# as.factor.tn5('2')
# as.numeric.tn5(as.factor.tn5('2'))

tn5_conversions <- data.frame(num=c(0.2, 0.5, 0.66, 1, 1.5, 2, 5))
tn5_conversions$fac <- factor(paste(tn5_conversions$num, 'X', sep=''), paste(tn5_conversions$num, 'X', sep = ''), ordered = T)
rownames(tn5_conversions) <- as.character(tn5_conversions$fac)

as.numeric.tn5 <- function(tn5) {
  return(as.numeric(as.character(gsub('X', '', tn5))))
}

as.factor.tn5 <- function(tn5) {
  if (is.character(tn5)) {
    return(tn5_conversions[tn5,'fac'])
  } else if (is.numeric(tn5)) {
    return(tn5_conversions[paste(tn5, 'X', sep=''),'fac'])
  }
}