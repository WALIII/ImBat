function mbdis=get_3d_dist(micpos,batpos)
%output is sound delay in ms dependent on mic to bat distance in 3D
if size(batpos,2)>1
mbdis=sqrt((micpos(1) - batpos(1,:)).^2 + (micpos(2) - batpos(2,:)).^2 ...
+ (micpos(3) - batpos(3,:)).^2);
else
  mbdis=sqrt((micpos(1) - batpos(1)).^2 + (micpos(2) - batpos(2)).^2 ...
+ (micpos(3) - batpos(3)).^2);  
end
mbdis=((mbdis./1e3)./340);