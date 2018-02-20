Graf_Ajuste <- function(x=0, y=0, a0=0, b0=0, c0=0, let=0, ajus=1){

   sigmoid = function(params, x) {
      params[1] / (1 + exp(params[2] * (x - params[3])))
   }

y0 = min(y)
print(y0)
plot(x, y, xlim=c(0, 80), main=paste(let,' - ecuacion: a/(1 + exp(b (x-c)))'))

cat('ajus=', ajus)
if (ajus == 1){
	y = y-y0
	fitmodel <- nls(y ~ a/(1 + exp(b * (x-c))), start=list(a=a0,b=b0,c=c0))
    params=round(coef(fitmodel),3)
    print(params)
    xx = seq(0,80,.1)
    y2 <- sigmoid(params,xx)
   y2 = y2 + y0
   points(xx, y2, type="l")
   sa = paste('a: ', params[1], sep='')
   sb = paste('b: ', params[2], sep='')
   sc = paste('c: ', params[3], sep='')

   text(40, 85, sa, pos=4)
   text(40, 75, sb, pos=4)
   text(40, 65, sc, pos=4)
   #text(40, 55, num, pos=4)
   }
   else {
   	x = seq(0,80,.1)
   	y = a0/(1 + exp(b0 * (x-c0))) + y0
   	points(x,y, type='l')
   }
 }
