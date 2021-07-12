function [L_new] = relabelL(L)

for ii = 1:size(L,3)
    temp = imcomplement(L(:,:,ii));
    temp = watershed(temp,8);
    L_new(:,:,ii) = bwlabel(temp);
end

end