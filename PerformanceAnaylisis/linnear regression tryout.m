load carsmall
Year = categorical(Model_Year);
tbl = table(MPG,Weight,Year);
mdl = fitlm(tbl,'MPG ~ Year + Weight^2');

plot(mdl)

%%
x=[1 2 3 4 5 6 7 8 9 10]

y=[4 5 2 7 2 8 10 2 1 5]

tbl = table(x',y')

mdl = fitlm(tbl,'linear')
plot(mdl)