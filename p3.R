rm(list=ls())

sigmoid = function(params, x) {

  params[1] / (1 + exp(params[2] * (x - params[3])))

}

data = read.csv("datos.csv", header=T)
print(data)
attach(data)

for (i in 3:3){
	j = 2*i-1
	x = data[,j]
	y = data[,j+1]
y = y-min(y)

jpeg(file = "p3.jpeg")

plot(x, y)

#fitmodel <- nls(y ~ a/(1 + exp(b * (x-c))), start=list(a=100,b=.5,c=55))
#222222222222
#fitmodel <- nls(y ~ a/(1 + exp(b * (x-c))), start=list(a=100,b=.5,c=10))

#3333333333
fitmodel <- nls(y ~ a/(1 + exp(b * (x-c))), start=list(a=100,b=.5,c=10))


# visualization code

# get the coefficients using the coef function

params=coef(fitmodel)
print(params)
var = paste(round(params,3), collapse = ', ')
val = paste(' valores a, b , c: ', var)

print(val)
text(30, 80, val, pos=4)

xx = seq(0,70,.1)
y2 <- sigmoid(params,xx)

points(xx, y2, type="l")
dev.off()
}