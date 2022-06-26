function [annualEms]=emissionDetector(surudekiEnIyiBireyinYeri,eff,emissionPrm,plantType,time)


for i=1:time
    annualEms(1,i)=0;
    for j=1:plantType
        annualEms(1,i)=annualEms(1,i)+emissionPrm(1,j)*surudekiEnIyiBireyinYeri(j,i)/eff(1,j);
    end
end

end