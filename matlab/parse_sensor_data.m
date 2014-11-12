function X = parse_sensor_data(sensor_data)

% Dataset from http://www.select.cs.cmu.edu/data/labapp3/index.html

% In this case, epoch is a monotonically increasing sequence number from each
% mote. Two readings from the same epoch number were produced from different
% motes at the same time. There are some missing epochs in this data set. 
% Moteids range from 1-54; data from some motes may be missing or truncated.
% Temperature is in degrees Celsius. Humidity is temperature corrected relative
% humidity, ranging from 0-100%. Light is in Lux (a value of 1 Lux corresponds 
% to moonlight, 400 Lux to a bright office, and 100,000 Lux to full sunlight.) 
% Voltage is expressed in volts, ranging from 2-3; the batteries in this case 
% were lithium ion cells which maintain a fairly constant voltage over their 
% lifetime; note that variations in voltage are highly correlated with 
% temperature.
% 
% We also have a Matlab-friendly version of the dataset that is sorted 
% sequentially by the time stamp in seconds. First measurement with time stamp 0
% is at 2004-02-28 00:58:15.315133. The file is 16MB gzipped, 88MB uncompressed.
% Each line contains the following entries:
% time:real 	moteid:int 	temperature:real 	humidity:real 	light:real 	voltage:real
% 
% Finally, we have aggregate connectivity data averaged over all time. This data
% consists of sender id, receiver id, and probability of a message from a sender
% successfully reaching a receiver. Note that this is not a symmetric
% relationship -- sensor A may hear B better than B hears A. 

n_day = 2;
% t = 30:30:n_day*24*60*60; % time is seconds for n_day days, in 30 seconds interval
% t = 30:30:n_day*24*60*1; % One hour for debug
t = 30:30:3233050;
ids = 1:54;
X = NaN * ones(length(t),length(ids));
for i = 1:length(t),
    d = sensor_data(sensor_data(:,1) == t(i), :);
    for j = 1:size(d,1),
        id = d(j,2);
        if find(ids == id),
            X(t(i)/30, id) = d(j,3);
        end
    end
end

% Remove columns that contain no data (full of NaN)
X(:,sum(isnan(X)) / size(X,1) == 1) = [];

% Interpolate to fill the missing data
% Approach taken from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/156746

for i = 1:size(X,2),
    d = X(:,i);
    d_x =find(~isnan(d));
    d_y = d(~isnan(d));
    X(:,i)=interp1(d_x,d_y,1:length(d));
end
