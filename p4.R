rm(list=ls())

source('Graf_Ajuste.R')

data = read.csv("datos.csv", header=T)
print(data)
attach(data)

#jpeg(file = "p4.jpeg")
par(mfcol=c(3,3))
Graf_Ajuste(x=B90_1, y=Area_1, a0=100, b0=.5, c0=50, let=1)
Graf_Ajuste(x=B90_2, y=Area_2, a0=100, b0=.5, c0=10, let=2)
Graf_Ajuste(x=B90_3, y=Area_3, a0=100, b0=.5, c0=20, let=3)
Graf_Ajuste(x=B90_4, y=Area_4, a0=100, b0=.5, c0=30, let=4)
Graf_Ajuste(x=B90_5, y=Area_5, a0=100, b0=.5, c0=30, let=5)
Graf_Ajuste(x=B90_6, y=Area_6, a0=84, b0=.7, c0=35, let='6NO', ajus=1)
Graf_Ajuste(x=B90_7, y=Area_7, a0=80, b0=.5, c0=20, let=7)
Graf_Ajuste(x=B90_8, y=Area_8, a0=80, b0=.5, c0=20, let=8)
Graf_Ajuste(x=B90_9, y=Area_9, a0=80, b0=.5, c0=50, let=9)

#dev.off()


