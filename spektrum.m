function S = spektrum(Hs,Ts)

    %addpath('D:\Paper EMA\Matlab EMA');
    %File = xlsread('hasil.xls', 3);
    %Hs=File(1,1); %significant wave height
    %Ts=File(1,2); %significant wave period
    
    tmax = 3600; % maximum duration of time series for one spectra
    del_t = 0.5; % targeted sampling rate
    fnyq = 1/(2*del_t);
    del_f = 1/tmax;
    f = del_f:del_f:fnyq;
    
    % spectrum generator
    gamma_1 = 3.3; % Assume that gamma_1 value is mean value
    beta_i = 0.0624*(1.094 - 0.01915*log(gamma_1))/(0.230 + 0.0336*gamma_1 - (0.185/(1.9 + gamma_1)));
    divtstp = 1 - 0.132*(gamma_1 + 0.2)^(-0.559);
    divltfp = 2*0.07^2;
    divgtfp = 2*0.09^2;
    Tp = Ts/divtstp;
    Hssq = Hs^2;
    Tppow4 = Tp^(-4);

    n = length(f);
    fp = 1/Tp;
    for i=1:n
        %two condition of sepctra depending on the value of fp
        if f(i)<=fp
           S(i) = beta_i*Hssq*Tppow4*f(i)^(-5)*exp(-1.25*(f(i)*Tp)^(-4))*3.3^(exp(-((f(i)*Tp - 1)^2)/divltfp));
        else 
           S(i) = beta_i*Hssq*Tppow4*f(i)^(-5)*exp(-1.25*(f(i)*Tp)^(-4))*3.3^(exp(-((f(i)*Tp - 1)^2)/divgtfp));
        end
    end
