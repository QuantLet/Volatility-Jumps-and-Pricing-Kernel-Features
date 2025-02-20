
function [sampleestimate returns_axis] = CDI_estimator(realizedKhRet, realizedQdenRet, splinefit, moments, bases, step ) %splinefit to be called as @splinefit
% INPUTS:   realizedKhRet = cell with each entry a vector representing support for returns (x-axis) only up to realized returns for the given  month
%           realizedQdenRet = cell with each entry a vector representing Q-density values corresponding to the realizedKhRet returns axis 
%           splinefit=user-suplied function used to fit the spline call with @splinefit (@OptSDF)
%           moments=number of moments used at each step in GMM 
%           bases = vector w number of bases from each model
%           step = spacing between values of estimated SDF
%a
%
% OUTPUTS:  sampleestimate = estimated sdf over entire returns space.  Will need to be truncated to reasonable return values where data is reasonable
%           returns_axis = axis of returns corresponding to sampleestimate sdf values
% ================================================================================================================         

n=length(realizedKhRet);
for h=1:length(realizedKhRet);
    mx(h)=max(realizedKhRet{h});
    mn(h)=min(realizedKhRet{h});
end
maxRet=max(mx);
minRet=min(mn);
returns_axis=minRet:step:maxRet;
r=length(returns_axis);
%numbmodels=length(moments);

%--------------------------------------------------------------------------
%               Create spline bases 
%--------------------------------------------------------------------------
interiorknots = minRet:(maxRet-minRet)/(bases-3):maxRet;
knots = augknt(interiorknots,4); %augmented knot sequence which contains the first and last knot
% with exact multiplicity k (k=4 for cubic-b-splines)
numbases=size(interiorknots,2)+2; %number of basis functions (note interiorknots actually contains 2 non-interior points)
Mspline=spmak(knots,eye(numbases));

tempsplineint=zeros(length(realizedKhRet),numbases);
for i=1:length(realizedKhRet);%i indexes day
    realizeKh=realizedKhRet{i};
    realizeQden=realizedQdenRet{i};
    splineday=(spval(Mspline,realizeKh))';%gives basis values for all return values x below the day's realized returns
    for j=1:numbases                      %where each column of splineday is a basis function over returns
        y=splineday(:,j);
        tempsplineint(i,j)=trapz(realizeKh,y'.*realizeQden);
    end
end
splineint=tempsplineint;%number of days-by-number of bases



%--------------------------------------------------------------------------
%     Sample estimate  (linear function of 'splinint' by matching moments)
%--------------------------------------------------------------------------
optimal_splinefit=@(x) splinefit(x,moments,interiorknots); %splinefit as function of only data, with moments and interiorknots fixed
sampleestimate = optimal_splinefit(splineint);
