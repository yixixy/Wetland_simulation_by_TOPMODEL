
%%% Fitting the parameters of the function from STocker et al. (2014)
% by Shushi Peng and Yi Xi update date: 2019/12/08
% contact yixi@pku.edu.cn

filein = 'G:\DATA\TOPMODEL\ga2\ga2.nc'; % CTI data
% Please replace it with your directory The CTI data can be available from
% Marthews et al. (2015)
% https://catalogue.ceh.ac.uk/documents/6b0c4358-2bf3-4924-aa8f-793d468b92be
Paradir = 'F:\1InundatedArea_180219\data\Parameter\';
% Please replace it with your directory

ul_lat = 86.09999996880516; % the upper boundery of the  CTI data
ul_lon = -180; % the left boundery of the CTI data
resol = 1/240; % resolution of the CTI data

lat = ncread(filein, 'lat');
lon = ncread(filein, 'lon');

ul_lat = ul_lat-1/240*length(lat);

resos = [0.1, 0.25, 0.5, 1, 2];
% you can only choose one for you foci

for M = 1:15
    
    % water table depth
    WT = -2:0.05:1;
    
    % function to fit the relationship between fsat and WT
    F = @(x, xdata)(1+x(1).*exp(-x(2).*(xdata-x(3)))).^(-1/x(1));
    x0 = [4, 10, 0];
    
    for rr = 1:5
        reso = resos(rr);
        vp = NaN(360/reso, 180/reso);
        kp = NaN(360/reso, 180/reso);
        qp = NaN(360/reso, 180/reso);

        for ii = 1:360/reso
            parfor jj = 1:180/reso
                
                latp_up = -90+reso*jj;
                latp_dn = -90+reso*(jj-1);
                latp = -90+reso*(jj-1/2);
                lonp_up = -180+reso*(ii-1);
                lonp_dn = -180+reso*ii;
                lonp = -180+reso*(ii-1/2);
                
                if latp_up <= max(lat) && latp_dn >= min(lat) &&...
                        lonp_up <= max(lon) && lonp_dn >= min(lon)
                    indlat = floor((latp_dn-ul_lat)/resol);
                    indlon = floor((lonp_up-ul_lon)/resol);
                    
                    % read CTI matrix
                    data_temp = ncread(filein, 'Band1', [indlon+1, indlat+1],...
                        [min(reso/resol, length(lon)-indlon), min(reso/resol, length(lat)-indlat)]);
                    data_temp(data_temp==-1) = nan;
                    data_temp = reshape(data_temp, [size(data_temp, 1)*size(data_temp, 2), 1]);
                    
                    if ~isempty(data_temp)
                        CTIam = nanmean(data_temp(:));
                        fsat = zeros(length(WT), 1);
                        
                        for pp = 1:length(WT)
                            % CTI threshold for inundation
                            CTIx = CTIam-M*WT(pp);
                            cnt = length(find(data_temp >= CTIx));
                            % calculate inundated fraction
                            fsat(pp) = cnt/length(data_temp);
                        end
                        
                        if max(fsat)>0
                            ind=find(fsat==max(fsat));
                            % remove the flat parts of the curve
                            [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0, WT(1:ind)', fsat(1:ind));
                            if ~isreal(x)
                                x(1:3) = NaN;
                            end
                        else
                            x(1:3) = NaN;
                        end
                        
                        vp(ii, jj) = x(1);
                        kp(ii, jj) = x(2);
                        qp(ii, jj) = x(3);
                    end
                end
            end
        end
        
        if ~exist([Paradir, 'reso', num2str(reso), '\'], 'dir')
            mkdir([Paradir, 'reso', num2str(reso), '\'])
        end
        save([[Paradir, 'reso', num2str(reso), '\'], ...
            'gridPara_reso', num2str(reso), '_M', num2str(M), '.mat'], 'vp', 'kp', 'qp');
        
    end
    disp(M)
    
end
