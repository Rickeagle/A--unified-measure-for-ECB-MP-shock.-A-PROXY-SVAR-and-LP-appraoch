% Interpolate Quarterly Data to Monthly Frequency

quarterlyData = readtable('QuarterlyData.xlsx');


dates = datetime(quarterlyData.DATE, 'InputFormat', 'yyyy-MM-dd');
values = quarterlyData.LaborCostIndex;


monthlyDates = dates(1):calmonths(1):dates(end);

% Perform linear interpolation
monthlyValues = interp1(dates, values, monthlyDates, 'linear');
monthlyData = table(monthlyDates', monthlyValues', 'VariableNames', {'Date', 'Labor_Cost_Index'});

writetable(monthlyData, 'MonthlyData.xlsx');
