clc;
clear;
delta = 1; % Time shift

%% Sweep over increasing N
for exp_n = 14:15
    N = 2^exp_n;

    %% Sweep k_true size from 1% to 10% of N
    k_ratios = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01]; % 1%, 2%, 5%, 10%

    fprintf('\n=== N = 2^%d = %d ===\n', exp_n, N);

    for r = k_ratios
        num_tones = max(1, round(N/2 * r));
        N_temp = N/2;
        k_true    = randperm(N_temp, num_tones);

        %% Build signal (sum of tones)
        n = (0:N-1);
        x = zeros(1, N);
        for i = 1:length(k_true)
            x = x + exp(1j * 2*pi * k_true(i)/N * n);
        end

        fprintf('\n  --- num_tones = %d (%.2f%% of N) ---\n', num_tones, r*100);

        %% Inner loop: stages k = 1 .. log2(N)
        for k = 1:log2(N)
            [~, X, ~]         = ditfft_stop(x, N, k);
            x_shifted         = circshift(x, delta);
            [~, X_shifted, ~] = ditfft_stop(x_shifted, length(x_shifted), k);

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
            fprintf('    Stage %2d | N = %7d | tones = %6d | Match: %s\n', ...
                    k, N, num_tones, string(match));
        end
    end
end
