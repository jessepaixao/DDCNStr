function I=trapeze(signal,dt,cdtinit)
n=length(signal);
I(1,1) =cdtinit;
for i=2:n
I(1,i) = I(1,i-1)+(dt/2)*(signal(1,i)+signal(1,i-1));
 %I(1,i+1) = I(1,i-1)+(dt)*(0.3854*signal(1,i+1)+1.2832*signal(1,i)+0.3854*signal(1,i-1));
end
end