function [Gamma]=pseudoxcorr(v1,v2,n,M)

Gamma=zeros(1,M);

for s=1:M
    for l=n-M+1:n-s+1
        Gamma(s)=Gamma(s)+v1(l)*v2(l+s-1);
    end
end

end