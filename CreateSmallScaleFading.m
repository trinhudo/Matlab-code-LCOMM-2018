function [Channel] = CreateSmallScaleFading( NoSources, ...
    NoDestinations, AntennaPerSource, AntennaPerDestination, dist)
% This function is to generate the channel
% NoSources : Number of source objects
% NoDestinations: Number of destination objects
% AntennaPerSource: Number of antennas per source (scalar or vector)
% AntennaPerDestination: Number of antennas per destination (scalar or vector)
%
if (nargin < nargin('CreateSmallScaleFading'))
    dist = {[]};
end
%
if (length(AntennaPerSource)>1)
    if (length(AntennaPerSource)<NoSources)
        disp('Error: Length of vector providing the number of anntennas per source is not enough');
        exit(0)
    end
else
    AntennaPerSource = repmat(AntennaPerSource,1,NoSources);
end
%
if (length(AntennaPerDestination)>1)
    if (length(AntennaPerDestination)<NoDestinations)
        disp('Error: Length of vector providing the number of anntennas per destination is not enough');
        exit(0)
    end
else
    AntennaPerDestination = repmat(AntennaPerDestination,1,NoDestinations);
end
%
for i = 1:1:NoSources
    for j=1:1:NoDestinations
        if (strcmp(dist{1},'Rice'))
            v = dist{2}*ones(AntennaPerSource(i),AntennaPerDestination(j));
            s = dist{3}*ones(AntennaPerSource(i),AntennaPerDestination(j));
            Channel(:,:,j,i) = sqrt(1/2)*(ricernd(v,s) + 1i*ricernd(v,s));
        else
            Channel(:,:,j,i) = sqrt(1/2)...
                *(randn(AntennaPerSource(i),AntennaPerDestination(j)) ...
                + 1i*randn(AntennaPerSource(i),AntennaPerDestination(j)));
        end
    end
end
end