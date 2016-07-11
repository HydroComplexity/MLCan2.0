function [day,month] = doy2date(yearV,doyV)
% The DOY2DATE function converts the day of year to a Gregorian date. 
%
% Syntax:
%
%   [DAY,MONTH] = DOY2DATE(DOY,YEAR)
%
% Inputs:
%
%   DOY - Day of year, enumerated in [1,366]
%   YEAR - Year
%
% Outputs:
%   
%   DAY - Day of month, enumerated in [1,31]
%   MONTH - Month, enumerated in [1,12]
%
% Deprecated Syntax (do not use in new code):
%
%   [X] = DOY2DATE(DOY,YEAR)
%
%   X(:,1) - DAY
%   X(:,2) - MONTH
%
% DOY2DATE was written by:
%
% Ethan Miller (esmiller@uiuc.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:  2009 February 3 / Version 1.1 (vectorized)
% Original:  before 2007 December 2 / Version 1.0

    if (length(doyV)~=length(yearV))
        error('doy2date: inputs must have the same dimensions');
        return;
    end

    day = zeros(size(doyV));
    month = zeros(size(doyV));
    
    if (nargout~=2)
        warning('doy2date: deprecated syntax, use separate DAY and MONTH vectors');
        day = zeros(length(doyV),2);
    end
    
    for k = 1:length(doyV)
        year = yearV(k);
        doy = doyV(k);

        daysmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
        if (((mod(year,4)==0) & (mod(year,100)~=0)) | (mod(year,400)==0))
            daysmonth(2)=daysmonth(2)+1;
        end
        for i = 1:12
            if (doy <= daysmonth(i))
                if (nargout==2)
                    day(k) = doy;
                    month(k) = i;
                elseif (nargout==1)
                    day(k,1) = doy;
                    day(k,2) = i;
                end
                break;
            else
                doy = doy - daysmonth(i);
            end
        end
    end
    

% function [dd,mm] = doy2date(year,doy)
% 
% if mod(year,4) ~= 0 
%     dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
% else
%     if mod(theYear,100) == 0
%         dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
%         if mod(theYear,1000) == 0
%             dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
%         end
%     else
%         dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
%     end
% end
%     
% mm = ones(size(doy))*NaN;
% dd = ones(size(doy))*NaN;
% 
% for im = 1:12
%    I = find(doy > sum(dayspermonth(1:im)) & doy <= sum(dayspermonth(1:im+1)));
%    mm(I) = ones(size(I)).*im;
%    dd(I) = doy(I) - ones(size(I))*sum(dayspermonth(1:im));
% end
% 
% 
% return;
% 
% if mod(year,4) ~= 0 
%     dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
% else
%     if mod(theYear,100) == 0
%         dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
%         if mod(theYear,1000) == 0
%             dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
%         end
%     else
%         dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
%     end
% end
% 
% im = 0;
% found = 0;
% while ~found
%    im = im+1;
%    if doy > sum(dayspermonth(1:im)) & doy <= sum(dayspermonth(1:im+1))
%       found = 1;
%       mm = im;
%       dd = doy - sum(dayspermonth(1:im));
%    end
% end
