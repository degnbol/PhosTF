
# purpose is visualizing a loss function that can be alternative to MSE, 
# since we might care more about explaining the presence of a perturbation effect rather than be worried if we overshoot the effect

# try different values of true:
true = 1.2
pred = seq(true-3,true+3,.01)

mse = (pred - true)^2
x = pred - true
a = -sign(true)
linex = exp(a*x) - a*x - 1 


plot(pred, mse, type="l")
lines(pred, linex)
