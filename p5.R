#http://rstudio-pubs-static.s3.amazonaws.com/7812_5327615eb0044cf29420b955ddaa6173.html

sigmoid <- function(x, lower_asymptote, carrying_capacity, growth_rate, time_max) {
    return(lower_asymptote + ((carrying_capacity - lower_asymptote)/(1 + exp(-growth_rate * 
        (x - time_max)))))
}
x <- 1:100
y <- sigmoid(1:100, 1, 50, 0.1, 50) + rnorm(100, 0, 20)

plot(x,y)

m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = min(y), 
    b = max(y), c = 1, d = round(median(x))), trace = TRUE)
    
print(m.s)

lines(x, fitted(m.s), lty = 2, lwd = 2, col = "red")

