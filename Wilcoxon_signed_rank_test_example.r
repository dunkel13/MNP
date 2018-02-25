ttest.wilcox.examples <- function(x,y,z,k)
{
  subjects <- x
  mean1 <- y
  mean2 <- z
  standarddev <- k
  print( c("Number of measurements:   ", x))
  print( c("Mean of group 1: ", y))
  print( c( "Mean of group 2: ",z))
  print( c("Standard deviation:  ", k))
  group1 <- rnorm(x, y, k)
  group2 <- rnorm(x, z, k)
  framedata <- cbind(group1, group2)
  print(framedata)
  print( list (t.test(group1, group2, var.equal = T), wilcox.test(group1, group2)))
}

set.seed(570)
ttest.wilcox.examples(13,90,100,10)
# <https://github.com/cyrilobolonsky/hypothesis-tests-in-R/blob/master/hypothesis_tests_in_R.R>
