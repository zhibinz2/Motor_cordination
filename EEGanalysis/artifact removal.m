imgaesc(trialdatamax)
imagesc(trialdatamax);colorbar
imagesc(log10(trialdatamax));colorbar
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > 100) = 1;
imagesc(binarythreshhold);
chansum = sum(binarythreshhold,2);
epochsum = sum(binarythreshhold,1);
bar(chansum)
badchanedit = find(chansum > 100);
badchanedit
goodchanedit = setdiff([1:128],badchanedit);
goodchanedit
bar(epochsum)
badepochedit = find(epochsum > 30);
goodepochedit = setdiff([1:350],badepochedit);
chansum2 = sum(binarythreshhold(goodchanedit,goodepochedit),2);
epochsum2 = sum(binarythreshhold(goodchanedit,goodepochedit),1);
whos
bar(goodchanedit,chansum2)
bar(goodepochedit,epochsum2)
bar(goodchanedit,chansum2)
bar(goodepochedit,epochsum2)
newbadepochedit =  goodepochedit(find(epochsum2 > 15))
badepochedit = [badepochedit newbadepochedit];
goodepochedit_final = setdiff([1:350],badepochedit);
chansum2 = sum(binarythreshhold(goodchanedit,goodepochedit),2);
bar(goodchanedit,chansum2)
newbadchanedit = goodchanedit(find(chansum2 > 50))
badchanedit = [badchanedit; newbadchanedit']
goodchanedit_final = setdiff([1:128],badchanedit)
goodepochedit_final
imagesc(trialdatamax(goodchan_final,goodepoch_final))
imagesc(trialdatamax(goodchanedit_final,goodepochedit_final))
imagesc(log10(trialdatamax(goodchanedit_final,goodepochedit_final)))