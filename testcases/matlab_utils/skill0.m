function [ s ] = skill0( measured, computed, varargin )
   %UNTITLED3 Summary of this function goes here
   %   Detailed explanation goes here

   OPT.estimate=[];
   OPT.nfirst = 0; % skip nfirst
   OPT.nlast = 0; % skip nlast
   OPT=setproperty(OPT,varargin{:});

   x       = unique([measured(:,1)]);
   x       = x(~isnan(x));

   zmt     = interp1(measured(:,1), measured(:,2), x); zmt=zmt(1+OPT.nfirst:end-OPT.nlast);
   zct     = interp1(computed(:,1), computed(:,2), x); zct=zct(1+OPT.nfirst:end-OPT.nlast);

   if ~isempty(OPT.estimate)
      zet = OPT.estimate;
   else
      zet=nan(size(zmt));
   end


   %% remove nanned data
   nd=ndims(zmt);

   N1=isnan(zmt);
   N2=isnan(zct);
   n=0;
   for i=1:nd
      n=n+sum(N1+N2);
   end

   if n>0
      disp('NaN''s found in data. Nan''s will not be used in calculation of skill');
      zmt(N1+N2>0)=[];
      zct(N1+N2>0)=[];
      zet(N1+N2>0)=[];
   end

   %% compute skills
   s.r2     = mean((zct-mean(zct)).*(zmt-mean(zmt)))/(std(zmt)*std(zct));
   rmse    = sqrt(mean((zct-zmt).^2));
   rmsm    = sqrt(mean(zmt.^2));
   s.sci     = rmse/max(rmsm,abs(mean(zmt)));

   s.relbias = mean(zct-zmt)/max(rmsm,abs(mean(zmt)));
   s.bias    = mean(zct-zmt);

   s.skill   = 1- sum((zct-zmt).^2) / sum((zmt).^2);
   s.rmse    = rmse;
   if ~isempty(OPT.estimate)
      zct(isnan(zet))=[];
      zmt(isnan(zet))=[];
      zet(isnan(zet))=[];
      s.bss     = 1-(sum((zct-zmt).^2))/(sum((zet-zmt).^2));
   end
end

