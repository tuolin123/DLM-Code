function dist2 = distance_p2(ref, pval)

if any(size(ref) ~= size(pval))
    error('Error. Please make sure the input have same dimension')
end

ind = ref <= 0.05;
sqdiff = (ref - pval).^2;
sqdiff005 = sqdiff(ind);
%diff = abs(ref - pval);

dist2 = sqrt(mean(sqdiff005(~isnan(sqdiff005))));
end