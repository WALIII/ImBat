function[outlier_flight_indexes] = find_outlier_flights(P,F)

% Function that identifies the outlier flights, in order to change the
% start time or exclude that flight from the clustering.

xv = [-2.4,2.0,2.0,-2.4,-2.4];
yv = [2.2,2.2,-1.55,-1.55,2.2];

xq = [F(1,:)];
yq = [F(2,:)];

[in,on] = inpolygon(xq,yq,xv,yv);

figure(); hold on
plot(xv,yv) % polygon
plot(xq(in&~on),yq(in&~on),'r+') % points strictly inside
plot(xq(on),yq(on),'k*') % points on edge
plot(xq(~in),yq(~in),'bo') % points outside
hold off

outlier_flight_indexes=in;

disp(strcat(num2str(sum(in))," flights identified as not starting on the wall or feeder"));

end