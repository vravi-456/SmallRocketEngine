tankRadius = 3;
fun = @(x)sqrt(tankRadius^2-(x-tankRadius).^2);
toMinimize = @(x)integral(fun,0,x) - V/(2*pi);
h = fzero(toMinimize,[0,tankRadius]) % m