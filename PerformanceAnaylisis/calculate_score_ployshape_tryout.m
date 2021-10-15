poly1 = polyshape([0 0 1 1],[1 0 0 1]);
poly2 = polyshape([0.75 1.25 1.25 0.75],[0.25 0.25 0.75 0.75]);
plot(poly1)
hold on
plot(poly2)
hold off

polyout = intersect(poly1,poly2)
polyout = union(poly1,poly2)

plot(polyout)

area(polyout)