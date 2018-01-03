function eta = ema(S)
% targeted sampling rate
del_t = 0.5;
[ro,co] = size(S);
tm = ro;
% maximum duration of time series for one spectra
tmax = 3600; 
fnyq = 1/(2*del_t);
del_f = 1/tmax;
% freuency bins of randomized interval
%r = ((0.1*rand(3600,1)-0.05)*del_f)+del_f;
%sum = 0;
%for i=1:length(r)
%   sum = sum + r(i);
%   f(i,:) = sum;
%end
f = del_f:del_f:fnyq;
pi2 = 2*pi();

a=(2*S*del_f).^0.5;
disp(a(1,1:10));

t = 0:del_t:tmax;
for i=1:tm
    for j=1:length(t)
        eta(i,j)=0;
        for k=1:length(f)
           eta(i,j)=eta(i,j)+(a(i,k)*cos((pi2*f(k))*t(j)+(rand(1)*pi2))); 
        end
    end
end
