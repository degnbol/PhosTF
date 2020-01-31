
# purpose is visualizing a loss function that can be alternative to MSE, 
# since we might care more about explaining the presence of a perturbation effect rather than be worried if we overshoot the effect

linex = function() {
    x = pred - true
    signs = sign(true)
    signs[signs == 0] = sample(c(-1,1), sum(signs == 0), replace=TRUE)
    a = -signs
    linex = exp(a*x) - a*x - 1 
    linex
}


mse = function(true) {(true - pred)^2}
quadquad = function(true, a) {
    if(missing(a)) a = abs(true)/(abs(true)+.5)
    (1 - a * (pred*true > true*true)) * (true - pred)^2
}

# try different values of 'true'
true = 1
pred = seq(true-10,true+10,.1)
plot(pred, mse(true), type="l", col="red")
lines(pred, quadquad(true), col="blue")













