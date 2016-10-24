
### A simple 2-dimensional example from R 

x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
cl <- kmeans(x, 2, 20)
plot(x, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8)


### Determining the number of clusters $K$

x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 2, sd = 0.3), ncol = 2))

cl2 <- kmeans(x, 2, 20)
cl3 <- kmeans(x,3,20)
cl4 <- kmeans(x,4,20)
cl5 <- kmeans(x,5,20)
cl6 <- kmeans(x,6,20)

wc2 <- sum(cl2$withinss)
wc3 <- sum(cl3$withinss)
wc4 <- sum(cl4$withinss)
wc5 <- sum(cl5$withinss)
wc6 <- sum(cl6$withinss)

plot(c(2:6),c(wc2,wc3,wc4,wc5,wc6),xlab="Number of clusters",
   ylab="Sum of squares",type="o")
lines(c(2:6),c(wc2,wc3,wc4,wc5,wc6))



