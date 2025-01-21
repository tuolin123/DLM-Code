function dist = distance_p(ref, pval)

if any(size(ref) ~= size(pval))
    error('Error. Please make sure the input have same dimension')
end

ind = ref <= 0.05;
ratiodiff = pval./ref;
ratiodiff005 = ratiodiff(ind);
%diff = abs(ref - pval);

dist = mean(ratiodiff005(~isnan(ratiodiff005)));
end

