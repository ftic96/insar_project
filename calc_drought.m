function [period, intensity] = calc_drought(spi, min_length)
% calc_drought is a tool for assigning periods of drought based on the
% standardized precipitation index (SPI). As defined by the World
% Meteorological Organization, a drought is a period where the SPI is
% consistently lower than 0, and at some point reaches -1. 

arguments
    spi (:,1) timetable
    min_length (1,1) double {mustBePositive} = 30
end

dry = spi < 0;
very_dry = spi < -1;

d = diff(dry);

%  Find the name of the timetable component
time_name = spi.Properties.DimensionNames{1};

% Finds starting and ending dates for a dry period
dry_starts = 1 + find(d.(1) == 1);
if dry.(1)(1) == 1
    dry_starts = [1;dry_starts];
end

dry_ends = 1 + find(d.(1) == -1);
if dry.(1)(height(dry)) == 1
    dry_ends = [dry_ends;height(dry)];
end

% Checks that the dry period reaches "very-dry" status
severe_enough = zeros(length(dry_starts),1);
intensity = zeros(length(dry_starts),1);
for jj = 1:length(dry_starts)
    severe_enough(jj) = any(very_dry.(1)(dry_starts(jj):dry_ends(jj)));
    intensity(jj) = sum(spi.(1)(dry_starts(jj):dry_ends(jj)));
end
dry_starts(severe_enough == 0) = [];
dry_ends(severe_enough == 0) = [];
intensity(severe_enough == 0) = [];

% Check that the drought is long enough
long_enough = dry_ends - dry_starts >= min_length;
dry_starts(long_enough == 0) = [];
dry_ends(long_enough == 0) = [];

start_dates = spi.(time_name)(dry_starts);
end_dates = spi.(time_name)(dry_ends);

period = table(start_dates, end_dates);
intensity(long_enough == 0) = [];

end