
FileTif='DSampled_190528_fly1.tif';
InfoImage=imfinfo(FileTif);
Ysiz = [InfoImage(1).Height InfoImage(1).Width length(InfoImage)];
%mImage=InfoImage(1).Width;
%nImage=InfoImage(1).Height;
%NumberImages=length(InfoImage);
Y=zeros(Ysiz(1),Ysiz(2),Ysiz(3),'single');
 
TifLink = Tiff(FileTif, 'r');
for i=1:Ysiz(3)
   TifLink.setDirectory(i);
   Y(:,:,i)=TifLink.read();
end
TifLink.close();