rm(list=ls())

source('Graf_Ajuste.R')

data = read.csv("datos.csv", header=T)
print(data)
attach(data)

x = as.double(B90_5)
y = as.double(Area_5)
plot(x,y)

ii = x < 50
x = x[ii]
y = y[ii]
plot(x,y)

plot(x,y)

m.s <- nls(y ~ b + ((a-b)/(1 + exp(-c * (x - d)))), start = list(a = min(y), 
    b = max(y), c = 1, d = round(median(x))), trace = TRUE)
    
print(m.s)

lines(x, fitted(m.s), lty = 2, lwd = 2, col = "red")



