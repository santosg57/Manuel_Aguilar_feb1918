rm(list=ls())

source('fun.R')
source('EncuentraPunto.R')

parametros_1 = c(.9506,         4.3094,   .2511,   33.01) # BIEN
parametros_2 = c(6.7207,   100.4621,   .5987,    8.7781) # BIEN
parametros_3 = c(10.5005,   102.2501, 0.7652, 14.4703) # BIEN
parametros_4 = c(12.8236, 112.3567, 0.6993, 20.7787) # BIEN
parametros_5 = c(5.768,       96.282,   1.452,    20.718) # BIEN
parametros_6 = c(9.380559, 91.335671,  1.914156, 34.432734) #
parametros_7 = c(10.0252, 89.8495, 0.7803, 38.3616) # BIEN
parametros_8 = c(16.537239, 104.978845,   2.244498,  41.034993) #
parametros_9 = c(12.6247, 96.3998, 1.1135, 53.1809) # BIEN

pp = matrix(c(parametros_1, parametros_2, parametros_3, parametros_4,
parametros_5, parametros_6, parametros_7, parametros_8, parametros_9), byrow = T, ncol=4)

#print(pp)

#data = read.csv("datos.csv", header=T)
#print(data)
#attach(data)

#x = as.double(B90_9)
#y = as.double(Area_9)

#ii = y < 10
#x = x[ii]
#y = y[ii]

#y2 = max(y)
#y = y*100/y2

ii <- as.integer(readline("Curva 1-9 o 0[TERMINA], Mete un numero= "))

if (1 <= ii & ii <= 9){
	xx= seq(0,100,.1)
	p = pp[ii,]
	a = p[1]
	b = p[2]
	c = p[3]
	d = p[4]
	yy = fun(xx, a, b, c, d)
	yy2 = max(yy)
	yy1 = min(yy)
	yn = yy*100/yy2
	ya = yy1*100/yy2
	plot(xx,yn, type='l', xlim=c(-5,100))
	points(c(xx[1], xx[length(xx)]), c(ya, ya), type='l', col='red')
	text(-3,ya, round(ya,1), col='red')
	y1 <- as.double(readline("Introduce el valor de Y0 para calcular X0= "))
    y0 = y1*yy2/100
    x0 = EncuentraPunto(p, y0)
    print(paste(' Valor de X0: ', round(x0,2)))
    points(c(0,x0), c(y1, y1), type='l')
    points(c(x0,x0), c(0, y1), type='l')
}


