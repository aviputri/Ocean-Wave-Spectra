function [hmean, hmax, hs, h10, hsp, heta] = sintetis(file_input)
    pdata = xlsread(file_input,3);
    pdata(:,1) = pdata(:,1)*0.01;   % convert cm -> m
    len = length(pdata);

    hs = [0,0];
    hmean = [0,0];
    hmax = [0,0];
    h10 = [0,0];

    % Calculate JONSWAP constants before entering loop
    gamma_1 = 3.3; % Assume that gamma_1 value is mean value
    beta_i = 0.0624*(1.094 - 0.01915*log(gamma_1))/(0.230 + 0.0336*gamma_1 - (0.185/(1.9 + gamma_1)));
    divtstp = 1 - 0.132*(gamma_1 + 0.2)^(-0.559);
    divltfp = 2*0.07^2;
    divgtfp = 2*0.09^2;

    for i=1:len
        disp(i);

        % JONSWAP spectrum generation
        S = spektrum2(pdata(i,1),pdata(i,2),beta_i,divtstp,divltfp,divgtfp);
        %hsp(i,:) = S;
        disp('spektrum selesai');

        % Inverse Discrete Fourier Transform to obtain sea level elevation
        E = ema(S);
        %heta(i,:) = E;
        disp('IDFT selesai');

        % Zerodown crossing
        [hmean(i,:), hmax(i,:), hs(i,:), h10(i,:)] = zerodown(E(1,:),0.5);
        disp('zerodown selesai');
        % Output Hs,Ts of synthetic and measured data, respectively
        disp(hs(i,:));
        disp(pdata(i,:));
    end

    [success,message] = xlswrite('sintetis.xls', hmean, 1);
    [success,message] = xlswrite('sintetis.xls', hmax, 2);
    [success,message] = xlswrite('sintetis.xls', hs, 3);
    [success,message] = xlswrite('sintetis.xls', h10, 4);
    
