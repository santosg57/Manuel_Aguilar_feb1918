rm(list=ls())

data = read.csv("datos.csv", header=T)
print(data)
attach(data)

jpeg(file = "p1.peg")
par(mfcol=c(3,3))
plot(B90_1, Area_1)
plot(B90_2, Area_2)
plot(B90_3, Area_3)
plot(B90_4, Area_4)
plot(B90_5, Area_5)
plot(B90_6, Area_6)
plot(B90_7, Area_7)
plot(B90_8, Area_8)
plot(B90_9, Area_9)

dev.off()


