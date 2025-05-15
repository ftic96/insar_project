function spi = calc_spi(precipitation, num_days, overlap) 

arguments
    precipitation (:,:) timetable
    num_days (1,1) double {mustBePositive} = 30
    overlap (1,1) double {mustBeNonnegative} = num_days - 1
end
addpath pearsons
% overlap must be less than num_days! 

% Pre-allocate Data
length_spi = ceil((height(precipitation) - num_days + 1) / (num_days - overlap));
rolling_window = zeros(length_spi, width(precipitation));
spi_mat = zeros(length_spi, width(precipitation));

% Set Dates
name_date = precipitation.Properties.DimensionNames{1};
i_SPI = num_days:(num_days-overlap):height(precipitation);
dates = precipitation.(name_date)(i_SPI); 

for ii = 1:width(precipitation)
    current_var = precipitation.(ii);
    for n = 1:length_spi
        idx = i_SPI(n) - num_days + 1: i_SPI(n);
        rolling_window = sum(current_var(idx));
    end

    zero_days = nnz(rolling_window == 0);
    chance_of_zero = zero_days / length(rolling_window);

    rw_dist = rolling_window;
    rw_dist(rolling_window == 0) = [];

    [alpha,beta,xi,Gamma] = pearson3_fit(rw_dist);
    cumulative_dist = pearson3_cdf(rolling_window, alpha, beta, xi, Gamma);
    cumulative_dist(isnan(cumulative_dist)) = 0;
    cumulative_dist = cumulative_dist + chance_of_zero .* (1 - cumulative_dist);

    spi_mat(:,ii) = norminv(cumulative_dist,0,1);
end

spi = array2timetable(spi_mat,...
    "VariableNames",precipitation.Properties.VariableNames, ...
    "RowTimes", dates);
end
