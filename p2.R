rm(list=ls())

# function needed for visualization purposes

sigmoid = function(params, x) {

  params[1] / (1 + exp(-params[2] * (x - params[3])))

}



x = 1:53

y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.18, 0.18, 0.18, 0.33, 0.33, 0.33, 0.33, 0.41,   0.41, 0.41, 0.41, 0.41, 0.41, 0.5, 0.5, 0.5, 0.5, 0.68, 0.58, 0.58, 0.68, 0.83, 0.83, 0.83, 0.74, 0.74, 0.74, 0.83, 0.83, 0.9, 0.9, 0.9, 1, 1, 1, 1, 1, 1, 1)

jpeg(file = "p2.jpeg")

plot(x,y)


# fitting code

fitmodel <- nls(y ~ a/(1 + exp(-b * (x-c))), start=list(a=1,b=.5,c=30))

# visualization code

# get the coefficients using the coef function

params=coef(fitmodel)
print(params)
var = paste(round(params,3), collapse = ', ')

val = paste(' valores a, b, c: ', var)
text(2, .8, val, pos=4)

y2 <- sigmoid(params,x)

points(y2,type="l")

dev.off()

#points(y)
