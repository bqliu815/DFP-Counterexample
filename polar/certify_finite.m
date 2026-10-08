function cert=certify_finite(stem)
% SPDX-License-Identifier: MIT
if nargin<1,stem='finite';end
load([stem,'.mat'],'objective','finiteOrbit');
D=[objective.points,objective.corrections,objective.radii];
S=sym(D,'f');
[~,ix]=sort(D(:,2));
gap=S(ix(2:end),2)-S(ix(1:end-1),2)-S(ix(2:end),5)-S(ix(1:end-1),5);
assert(all(isAlways(gap>0)),'Projection separation certificate failed');
R=sym(ceil(objective.ratio*1e8))/sym(100000000);
margin=R^2*S(:,5).^2-S(:,3).^2-S(:,4).^2;
assert(all(isAlways(margin>=0)),'Correction bound certificate failed');
s3=sym(1732050808)/sym(1000000000);assert(isAlways(s3^2>3));
p=(sym(135)/16+15*s3/2)*R;
pUpper=sym(ceil(double(p)*1e6))/sym(1000000);
assert(isAlways(pUpper>=p)&&isAlways(pUpper<1));
cert=struct('status','exact_dyadic_checks_passed','points',size(D,1), ...
    'ratio_upper',double(R),'perturbation_upper',double(pUpper),'hessian_band',double([1-pUpper,1+pUpper]), ...
    'sqrt3_rational_upper','1732050808/1000000000', ...
    'separation','disjoint projections on the second coordinate', ...
    'data_interpretation','stored binary64 values treated as exact dyadic rationals', ...
    'minimum_projection_margin',double(min(gap)),'matlab',version);
if strcmp(stem,'finite'),path='certificate.json';else,path=[stem,'_certificate.json'];end
PolarDFP.writeJSON(path,cert);disp(cert);
end
