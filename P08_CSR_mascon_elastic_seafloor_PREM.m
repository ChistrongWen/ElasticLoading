function OUT = P08_CSR_mascon_elastic_seafloor_PREM()
% P08 | CSR mascon trend -> PREM elastic loading displacement of seafloor.
% NorthPacific_SeaLevel; 2005-01 through 2024-12.
% Input CSR rg MUST be equivalent water thickness (EWH), NOT a sea-level
% trend already fitted, and use the full GLOBAL load, not just regional load.
% Output: U (up), E (east), N (north) in mm/yr.
%
% Implements the user's Farrell/PREM convention:
% G_U = ur2/(1e12 * Re * angular_distance_rad), source mass = h*area*rho.
% The numerical cell-centre approximation is regularized for near-field
% distances. Compare several nearAngleDeg/targetStepDeg before publication.
%
% No source data are modified. This does NOT automatically remove GIA,
% land leakage, or ocean/land mass ambiguities of GRACE mascon solutions.

%% ---------------------------- SETTINGS -----------------------------
root = '/Volumes/PolyUSSD/HongKong/p-2026-NorthPacific_SeaLevel';
srcFile = fullfile(root,'data','processed','GRACE','GRACE_mascon_4_20260707.mat');
boundFile = fullfile(root,'config','regions','NorthPacific_160E_130W_35N.txt');
greenFile = fullfile(root,'config','greens','cf_PREM.txt');
outDir = fullfile(root,'results','solid_earth','CSR_elastic_loading_2005_2024');
gmtDir = fullfile(root,'gmt','input','fig_CSR_elastic_loading');

% Verify from the original CSR data generator; never infer EWH units by size.
sourceEWHUnit = 'cm';  % 'cm', 'mm' or 'm'; change if necessary
rho = 1000;           % kg/m^3; EWH-to-mass conversion
Re = 6371e3;          % m, consistent with the user's mass_load_hv
startMonth = datetime(2005,1,1);
endMonth   = datetime(2024,12,1);
months = (startMonth:calmonths(1):endMonth)';
minTrendMonths = 100;   % minimum valid months to fit spatial trend
minSourceAbsGt = 0;     % 0 => do not threshold source mass trends
% Grid spacing of target seafloor sites. 2 degrees for a first test;
% use 1 degree for the publication-quality spatial/area mean test.
targetStepDeg = 2;
% Critical sensitivity setting: minimum Green-function angular distance.
% Avoids 1/theta singularity of point-centre integration. See caveats.
nearAngleDeg = 0.25;
% Number of GLOBAL source cells to process per batch (memory control).
sourceChunk = 30000;
% Set to false if original PREM table theta is not in degrees.
greenThetaInDegrees = true;

assert(isfile(srcFile),'Missing CSR mascon file: %s',srcFile);
assert(isfile(boundFile),'Missing study polygon: %s',boundFile);
assert(isfile(greenFile),['Missing cf_PREM.txt at %s. ' ...
    'Set greenFile to the actual PREM table path.'],greenFile);
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(gmtDir), mkdir(gmtDir); end
assert(numel(months)==240,'Analysis calendar should contain 240 months.');

%% 1. Load CSR (the three GRACE solutions are in one MAT file)
fprintf('Loading CSR mascon: %s\n',srcFile);
S = load(srcFile,'mascon_csr');
assert(isfield(S,'mascon_csr'),'Missing mascon_csr structure.');
D = S.mascon_csr; clear S
assert(all(isfield(D,{'lon','lat','rg','tt'})), ...
    'Expected mascon_csr.lon, lat, rg, tt.');
t0 = parseCSRtime(D.tt);  % [decimal_year, year, month]
[nr,nc,nt] = size(D.rg);
assert(nt==numel(t0),'CSR rg time dimension mismatches tt.');
[LON,LAT] = gridCoords(D.lon,D.lat,nr,nc);
assert(max(abs(LAT(:)))<=90+1e-6,'Latitude values invalid.');

% Grid-cell area, spherical and longitude-cyclic, inferred from axes.
lonVec = mod(double(LON(1,:)),360);
latVec = double(LAT(:,1));
assert(numel(lonVec)>1 && numel(latVec)>1);
dlon=median(abs(diff(unwrap(deg2rad(lonVec)))))*180/pi;
dlat=median(abs(diff(latVec)));
assert(dlon>0 && dlat>0 && dlon<3 && dlat<3, ...
    'Unexpected grid spacing, check longitude/latitude arrangement.');
latLo=max(-90,latVec-dlat/2);latHi=min(90,latVec+dlat/2);
areaRow=(Re^2)*deg2rad(dlon).*(sind(latHi)-sind(latLo));
areaM2=repmat(abs(areaRow),1,nc);

switch lower(sourceEWHUnit)
    case 'mm', scaleToMetres=1e-3;
    case 'cm', scaleToMetres=1e-2;
    case 'm',  scaleToMetres=1;
    otherwise, error('Invalid sourceEWHUnit, use mm/cm/m.');
end

%% 2. Explicit monthly aggregation: duplicate same-month samples are averaged
% Spatial OLS trend on unsmoothed monthly EWH grids. No interpolation.
yr=year(months)+(month(months)-0.5)/12;
td=yr-mean(yr);
keys0=year(t0)*12+month(t0);
keys=year(months)*12+month(months);
s0=zeros(nr,nc);
sx=zeros(nr,nc); sxx=zeros(nr,nc); sy=zeros(nr,nc); sxy=zeros(nr,nc);
validMonths=zeros(numel(months),1); duplicates=0;
for im=1:numel(months)
    idx=find(keys0==keys(im));
    if isempty(idx), continue; end
    duplicates=duplicates+max(numel(idx)-1,0);
    monthSum=zeros(nr,nc); monthN=zeros(nr,nc);
    for jj=1:numel(idx)
        z=double(D.rg(:,:,idx(jj)));
        v=isfinite(z)&abs(z)<1e10;
        monthSum(v)=monthSum(v)+z(v);
        monthN(v)=monthN(v)+1;
    end
    v=monthN>0;
    if ~any(v(:)),continue;end
    z=monthSum(v)./monthN(v); % units are still source EWH units
    x=td(im);
    s0(v)=s0(v)+1;
    sx(v)=sx(v)+x;
    sxx(v)=sxx(v)+x*x;
    sy(v)=sy(v)+z;
    sxy(v)=sxy(v)+x.*z;
    validMonths(im)=nnz(v);
end
clear D monthSum monthN
fprintf('CSR duplicate extra epochs within analysis calendar: %d\n',duplicates);
den=sxx-sx.^2./max(s0,1);
trendEWH=nan(nr,nc);ok=s0>=minTrendMonths & den>0;
trendEWH(ok)=(sxy(ok)-sx(ok).*sy(ok)./s0(ok))./den(ok); % source EWH / yr
trendEWH_mm=trendEWH * scaleToMetres * 1000;
clear sx sxx sy sxy den
fprintf('Valid CSR source trend cells: %d / %d\n',nnz(ok),numel(ok));

%% 3. Global source mass-change trend kg/year per mascon grid cell
load_kgyr=trendEWH(:) .* scaleToMetres .* areaM2(:) * rho;
sv=isfinite(load_kgyr)&isfinite(LAT(:))&isfinite(LON(:));
if minSourceAbsGt>0
    sv=sv & abs(load_kgyr)>=minSourceAbsGt*1e12;
end
srcLon=mod(double(LON(sv)),360);
srcLat=double(LAT(sv));
srcMass=load_kgyr(sv);
assert(~isempty(srcMass),'No valid global CSR mass load trend cells.');
clear load_kgyr
fprintf('Global source cells: %d; grid %.4f x %.4f deg\n',numel(srcMass),dlat,dlon);

%% 4. Read normalized PREM Green functions from the user's table
T=readmatrix(greenFile);
assert(size(T,2)>=5,'PREM table must have >=5 columns.');
q=T(:,1);ur=T(:,4);vh=T(:,5);
good=isfinite(q)&isfinite(ur)&isfinite(vh)&q>=0;
q=q(good);ur=ur(good);vh=vh(good);
if ~greenThetaInDegrees, q=rad2deg(q); end
[q,idx]=unique(q,'stable');ur=ur(idx);vh=vh(idx);
[q,ix]=sort(q);ur=ur(ix);vh=vh(ix);
assert(q(1)<=nearAngleDeg && q(end)>=180-1e-4, ...
    'Green table should cover near-field and ~180 degrees. Check theta unit/range.');

%% 5. Sample seafloor points on a REGULAR regional grid
parts=readBoundaryParts(boundFile);
Tlon=0.5*targetStepDeg:targetStepDeg:(360-0.5*targetStepDeg);
Tlat=(-90+0.5*targetStepDeg):targetStepDeg:(90-0.5*targetStepDeg);
[TLO,TLA]=meshgrid(Tlon,Tlat);
mask=makeRegionMask(TLO,TLA,parts);
targetLon=TLO(mask);targetLat=TLA(mask);
assert(~isempty(targetLon),'No grid points inside regional polygon.');
% Regular cell areas; cos(latitude) proportional to area for fixed grid.
weights=cosd(targetLat);
N=numel(targetLon);
U=nan(N,1);E=nan(N,1);North=nan(N,1);
fprintf('Computing elastic PREM displacement at %d regional seafloor sites...\n',N);

% Interpolation angular degrees. source -> observation azimuth (degrees).
% To limit memory, process sources in chunks for every observation site.
for it=1:N
    lo=targetLon(it);la=targetLat(it);
    tsum=[0 0 0];
    for st=1:sourceChunk:numel(srcMass)
        en=min(st+sourceChunk-1,numel(srcMass)); ix=st:en;
        [arc,az]=sphericalDistanceAzimuth(srcLon(ix),srcLat(ix),lo,la);
        arcSafe=max(arc,nearAngleDeg); % regularized cell-centre integral
        gU=interp1(q,ur,arcSafe,'linear',NaN) ./ (1e12*Re*deg2rad(arcSafe));
        gH=interp1(q,vh,arcSafe,'linear',NaN) ./ (1e12*Re*deg2rad(arcSafe));
        assert(all(isfinite(gU))&&all(isfinite(gH)), ...
            'Green interpolation outside table or missing values.');
        masses=srcMass(ix);
        tsum(1)=tsum(1)+sum(gU.*masses);
        tsum(2)=tsum(2)+sum(-sind(az).*gH.*masses);
        tsum(3)=tsum(3)+sum(-cosd(az).*gH.*masses);
    end
    U(it)=1000*tsum(1);E(it)=1000*tsum(2);North(it)=1000*tsum(3);
    if mod(it,20)==0 || it==N
        fprintf('  %d/%d sites\n',it,N);
    end
end

%% 6. Regional area-weighted mean and outputs
uMean=sum(U.*weights)/sum(weights);
eMean=sum(E.*weights)/sum(weights);
nMean=sum(North.*weights)/sum(weights);
OUT=struct('source','CSR','source_file',srcFile,'Green_file',greenFile, ...
    'region_file',boundFile,'time_start',startMonth,'time_end',endMonth, ...
    'EWH_input_unit',sourceEWHUnit,'density_kg_m3',rho, ...
    'source_lon',LON,'source_lat',LAT,'source_trend_EWH_mm_per_yr',trendEWH_mm, ...
    'target_lon',targetLon,'target_lat',targetLat,'target_weight',weights, ...
    'up_mm_per_yr',U,'east_mm_per_yr',E,'north_mm_per_yr',North, ...
    'regional_up_mm_per_yr',uMean,'regional_east_mm_per_yr',eMean, ...
    'regional_north_mm_per_yr',nMean,'target_step_degree',targetStepDeg, ...
    'nearfield_cutoff_degree',nearAngleDeg, ...
    'source_valid_months',s0,'monthly_source_cells_count',validMonths, ...
    'note',['Elastic displacement caused by global CSR EWH trend using ' ...
    'point-centre PREM Green functions with nearfield angular regularization. ' ...
    'Sensitivity to target step and cutoff must be tested. ' ...
    'Not automatically a correction to altimetry or full sea-level budget.']);
save(fullfile(outDir,'P08_CSR_elastic_PREM_results.mat'),'OUT','-v7.3');
writetable(table(targetLon,targetLat,U,E,North, ...
    'VariableNames',{'lon','lat','up_mm_per_yr','east_mm_per_yr','north_mm_per_yr'}), ...
    fullfile(outDir,'P08_CSR_elastic_displacement_sites.csv'));
writematrix([targetLon,targetLat,U], ...
    fullfile(gmtDir,'CSR_elastic_up_trend_xyz.txt'),'Delimiter','tab');
writematrix([mod(LON(ok),360),LAT(ok),trendEWH_mm(ok)], ...
    fullfile(gmtDir,'CSR_EWH_trend_global_xyz.txt'),'Delimiter','tab');

f=figure('Color','w','Position',[100 100 950 470]);
scatter(targetLon,targetLat,22,U,'filled');axis tight;grid on;
cb=colorbar;cb.Label.String='Elastic vertical displacement (mm/yr)';
xlabel('Longitude (deg E)');ylabel('Latitude (deg N)');
title(sprintf('CSR/PREM elastic seabed motion 2005-2024 | mean %+.3f mm/yr',uMean));
exportgraphics(f,fullfile(outDir,'P08_CSR_elastic_up_map.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'P08_CSR_elastic_up_map.pdf'),'ContentType','vector');
close(f);

fprintf('\n===== CSR global mass trend => regional elastic seafloor rate =====\n');
fprintf('Period: 2005-01 to 2024-12 | %d valid source cells\n',nnz(ok));
fprintf('Study-region sample points: %d | step: %.2f deg\n',N,targetStepDeg);
fprintf('PREM near-field cutoff: %.3f deg\n',nearAngleDeg);
fprintf('Regional mean UP   : %+.6f mm/yr (positive = uplift)\n',uMean);
fprintf('Regional mean EAST : %+.6f mm/yr\n',eMean);
fprintf('Regional mean NORTH: %+.6f mm/yr\n',nMean);
fprintf('Saved outputs: %s\n',outDir);
end

function t=parseCSRtime(a)
if isdatetime(a)
    t=datetime(year(a),month(a),1);return;
end
a=double(a);
assert(ismatrix(a)&&size(a,2)>=3,'Expected CSR tt Nx3: [decimalYear,year,month].');
y=a(:,2); m=a(:,3);
assert(all(isfinite(y)&isfinite(m)&y>=1990&y<=2100&m>=1&m<=12& ...
    abs(y-round(y))<1e-6&abs(m-round(m))<1e-6), ...
    'CSR tt columns 2 and 3 are not year/month.');
t=datetime(y,m,ones(size(y)));
end

function [LON,LAT]=gridCoords(lon,lat,nr,nc)
lon=double(lon);lat=double(lat);
if isvector(lon)&&isvector(lat)
    assert(numel(lon)==nc&&numel(lat)==nr,'Grid coordinates inconsistent with rg.');
    [LON,LAT]=meshgrid(lon(:)',lat(:));
elseif isequal(size(lon),[nr nc]) && isequal(size(lat),[nr nc])
    LON=lon;LAT=lat;
else
    error('CSR lon/lat must be vectors or matrices with rg spatial size.');
end
end

function parts=readBoundaryParts(file)
B=readmatrix(file);assert(size(B,2)>=2,'Boundary requires lon and lat.');
x=B(:,1);y=B(:,2);v=isfinite(x)&isfinite(y)&abs(y)<=90;
cuts=[0;find(~v);numel(x)+1];parts={};
for j=1:numel(cuts)-1
    ii=(cuts(j)+1):(cuts(j+1)-1);if numel(ii)<3,continue;end
    xx=rad2deg(unwrap(deg2rad(mod(x(ii),360)))); yy=y(ii);
    keep=[true;abs(diff(xx))+abs(diff(yy))>1e-10];
    xx=xx(keep);yy=yy(keep);
    if numel(xx)>=3,parts{end+1}=[xx(:),yy(:)];end %#ok<AGROW>
end
assert(~isempty(parts),'No valid polygon rings.');
end

function mask=makeRegionMask(lon,lat,parts)
lon=mod(lon,360);mask=false(size(lon));
for j=1:numel(parts)
    p=parts{j};
    for shift=-2:2
        [in,on]=inpolygon(lon+360*shift,lat,p(:,1),p(:,2));
        mask=mask|in|on;
    end
end
mask=mask&isfinite(lon)&isfinite(lat)&abs(lat)<=90;
end

function [arc,az]=sphericalDistanceAzimuth(lon1,lat1,lon2,lat2)
% Source (lon1,lat1) to receiver (lon2,lat2); az clockwise from north.
ph1=deg2rad(lat1);ph2=deg2rad(lat2);
dl=deg2rad(mod(lon2-lon1+180,360)-180);
c=sin(ph1).*sin(ph2)+cos(ph1).*cos(ph2).*cos(dl);
arc=rad2deg(acos(min(1,max(-1,c))));
y=sin(dl).*cos(ph2);
x=cos(ph1).*sin(ph2)-sin(ph1).*cos(ph2).*cos(dl);
az=mod(rad2deg(atan2(y,x)),360);
end
