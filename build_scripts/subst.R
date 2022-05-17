library(terra)
packageVersion("terra")

r <- rast(ncols=5, nrows=5, xmin=0, xmax=1, ymin=0, ymax=1, crs="")
r <- init(r, 1:6)
x <- subst(r, 3, 7); x11(); plot(r)
x <- subst(r, 2:3, NA); x11(); plot(x)
x <- subst(x, NA, 10); x11(); plot(x)
x <- subst(r, c(5,3), c(3,5)); x11(); plot(x)

R <- as.matrix( r, wide=TRUE )
new <- c(3,5)
old <- c(5,3)
x <- as.vector(R)
y <- c(new,x)[match(x, c(old,x))]
