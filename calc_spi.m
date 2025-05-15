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

% for ii = 1:width(precipitation)
current_var = (precipitation);%.(ii);
if ~isempty(current_var)
    for n = 1:length_spi
        daterange = timerange(dates(n), dates(n) + num_days);
        rolling_window(n,:) = sum(current_var{daterange,:},1);
    end

    zero_days = sum(rolling_window == 0);
    chance_of_zero = zero_days / length(rolling_window);
    
    for ii = 1:width(precipitation)
        rw_dist = rolling_window(:,ii);
        rw_dist(rolling_window(:,ii) == 0) = [];

        [alpha,beta,xi,Gamma] = pearson3_fit(rw_dist);
        cumulative_dist = pearson3_cdf(rolling_window(:,ii), alpha, beta, xi, Gamma);
        cumulative_dist(isnan(cumulative_dist)) = 0;
        cumulative_dist = cumulative_dist + chance_of_zero(ii) .* (1 - cumulative_dist);

        spi_mat(:,ii) = norminv(cumulative_dist,0,1);
    end
end

spi = array2timetable(spi_mat,...
    "VariableNames",precipitation.Properties.VariableNames, ...
    "RowTimes", dates);
end
