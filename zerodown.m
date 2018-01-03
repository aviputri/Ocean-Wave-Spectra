function [M_mean, M_max, M_significant, M_10, M, M_sorted] = zerodown(E, sr)
M = [0 0];		% matrix to contain (unsorted) H and T:
% first column contains H, second column contains T

n = length(E);	% number of data points
iz = 1;			% index of current zero down crossing
tz = 1.0;		% interpolated time of current zero down crossing
tz_prev = 1.0;	% interpolated time of previous zero down crossing
i_prev = 1;

% determine mean sea level for this function call
msl = mean(E);

% normalize elevation
for i = 1:n
    E(i) = E(i) - msl;
end

% determine height and length of each zero down crossing
for i = 1:n-1
    % if elevation is positive and then becomes negative
    if (E(i) >=0 && E(i+1) < 0)
        % period of current wave is equal to the (interpolated) time of zero down crossing subtracted by time of previous zero down crossing
        tz = lagrpol(E(i:i+1), i:i+1, 2, 0.0);
        T = sr * (tz - tz_prev);
        M(iz,2) = T;

        % height of the wave equals maximum elevation during this period subtracted by minimum elevation during period
        % SHOULD USE INTERPOLATION
        Z = E(i_prev:i);
        M(iz,1) = max(Z) - min(Z);

        iz = iz + 1;
        tz_prev = tz;
        i_prev = i;
    end
end

% select column 1 (wave amplitude) and sort in descending order
[H,order] = sort(M(:,1), 1, 'descend');
% 
M_sorted = M(order,:);

% calculate Hs, Ts, H10, T10, and Hmean and Tmean in one pass
n = length(M);
ns = floor(n/3);
n10 = floor(n/10);
Havg = 0;
Tavg = 0;
for i = 1:n
    Havg = Havg + (M_sorted(i,1) - Havg) / i;
    Tavg = Tavg + (M_sorted(i,2) - Tavg) / i;
    if (i == n10)
        M_10 = [Havg, Tavg];
    end
    if (i == ns)
        M_significant = [Havg, Tavg];
    end
end

M_mean = [Havg, Tavg];
M_max = M_sorted(1,:);
