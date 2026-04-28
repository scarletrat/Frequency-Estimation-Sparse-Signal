function [X, Xe_out, Xo_out] = ditfft_stop(x, N0, stop_level)
if nargin == 1
    N0 = length(x);
    stop_level = log2(N0);
elseif nargin == 2
    stop_level = log2(N0);
end

N = length(x);

% --- Base case: always go all the way to N=1 ---
if N == 1
    X = x;
    Xe_out = [];
    Xo_out = [];
    return;
end

% Always recurse fully down regardless of stop_level
[Xe, ~, ~] = ditfft_stop(x(1:2:end), N0, stop_level);
[Xo, ~, ~] = ditfft_stop(x(2:2:end), N0, stop_level);

% current_level: how many butterfly stages have been completed coming back up
% e.g. N=2 means 1 level up from bottom, N=8 means 3 levels up
levels_from_bottom = log2(N);

% Stop combining if we've gone up more stages than allowed 
if levels_from_bottom > stop_level
    X = [Xe, Xo];      
    Xe_out = Xe;
    Xo_out = Xo;
    return;
end

k = 0:N/2-1;
W = exp(-1j*2*pi*k/N);

if N == N0
    Xe_out = Xe;
    Xo_out = Xo;
else
    Xe_out = [];
    Xo_out = [];
end

X = [Xe + W.*Xo, Xe - W.*Xo];
end
