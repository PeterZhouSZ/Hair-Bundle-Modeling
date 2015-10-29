% Hair Bundle Geometry Script
clf

L1 = 4.7;   % length of shortest stereocilia
L2 = 8.3;   % length of longest stereocilia
theta1 = 10;% angle of shortest stereocilia to the vertical
theta2 = 11;% angle of longest stereocilia to the vertical
base = 4;   % length of base of the hair bundle

line([0 base],[0 0])
line([0 L1*sind(theta1)], [0 L1*cosd(theta1)])
line([base base-L2*sind(theta2)], [0 L2*cosd(theta2)])
line([L1*sind(theta1) base-L2*sind(theta2)],[L1*cosd(theta1) L2*cosd(theta2)])

A = L1*sind(theta1);
B = L1*cosd(theta1);
C = base-L2*sind(theta2);
D = L2*cosd(theta2);
hypotenuse = sqrt(((C-A)^2) + ((D-B)^2));
theta3 = asind((C-A)/hypotenuse);
x = (hypotenuse/2)*sind(theta3);
y = (hypotenuse/2)*cosd(theta3);
midpoint_x = A + x;
midpoint_y = B + y;

xx = cosd(theta3);
yy = sind(theta3);

ionto_x = midpoint_x - xx;
ionto_y = midpoint_y + yy;

axis([0 10 0 10])
hold on
plot(midpoint_x,midpoint_y,'.')
plot(L1*sind(theta1),L1*cosd(theta1),'.')
plot(base-L2*sind(theta2),L2*cosd(theta2),'.')
plot(ionto_x, ionto_y,'.')
axis square
set(gca, 'Visible','off')

r1 = distance(A,B,ionto_x, ionto_y);
r2 = distance(midpoint_x, midpoint_y, ionto_x, ionto_y);
r3 = distance(C,D,ionto_x, ionto_y);

