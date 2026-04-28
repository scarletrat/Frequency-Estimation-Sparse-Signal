clc;
clear;
delta = 1; % Time shift

% Open file for writing
fid = fopen('Log1.txt', 'w');
%% Sweep over increasing N
        SNR_dB = [0 5 20 40];                 % Desired Signal-to-Noise Ratio in dB
for snr = SNR_dB
        fprintf(fid,'\n=== SNR_db %d ===\n', snr);
        fprintf('\n=== SNR_db %d ===\n', snr);
for exp_n = 12
    N = 2^exp_n;
    %% Sweep k_true size from 1% to 10% of N
    k_ratios = [1];
    fprintf(fid,'\n=== N = 2^%d = %d ===\n', exp_n, N);

    for r = 1
        num_tones = r;
        N_temp = N/2;
        %k_true    = randperm(N_temp, num_tones);
        k_true = randperm(N_temp,1);
        %% Build signal (sum of tones)
        n = (0:N-1);
        x = zeros(1, N);
        for i = 1:length(k_true)
            x = x + exp(1j * 2*pi * k_true(i)/N * n);
        end
        x = awgn(x, snr, 'measured'); % Add AWGN        fprintf('\n  --- num_tones = %d---\n', num_tones);

        %% Inner loop: stages k = 1 .. log2(N)
        for k = 1
            [~, X, ~]         = ditfft_stop(x, N, k);
            x_shifted         = circshift(x, delta);
            [~, X_shifted, ~] = ditfft_stop(x_shifted, length(x_shifted), k);
            %figure;
            %plot(1:N/2,abs(X));
            %% Find strongest aliased bins (2x oversampled)
            [~, sorted_idx] = sort(abs(X), 'descend');
            idx = sorted_idx(1 : num_tones * 2);

            %% Recover true bins from phase difference
            k_est = zeros(1, num_tones * 2);
            for i = 1:length(k_est)
                dphi     = angle(X_shifted(idx(i))) - angle(X(idx(i)));
                k_est(i) = mod(-round(dphi * N / (2*pi * delta)), N);
            end

            %% Deduplicate
            if length(unique(k_est)) ~= num_tones
                k_unique = k_est(1:num_tones);
            else
                k_unique = unique(k_est);
            end

            %% Report
            match = all(sort(k_unique) == sort(k_true));
           
              fprintf(fid,'    Stage %2d | N = %7d | tones = %6d | Match: %s\n', ...
                    k, N, num_tones, string(match));
fprintf(fid,'True bins   : %d\n', k_true);
fprintf(fid,'Two Strongest Aliased bins: %d %d\n', idx-1);
fprintf(fid,'Estimated   : %d %d\n', (k_est));
        end
    end
end
end
% Close file
fclose(fid);
