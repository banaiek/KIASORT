function adjusted_thr = adjust_correlation_threshold(N, M, original_thr)
    df_N = N - 2;
    t = original_thr * sqrt(df_N) / sqrt(1 - original_thr^2); % T-statistic
        df_M = M - 2;
    adjusted_thr = t / sqrt(t^2 + df_M); % New threshold for M
end

