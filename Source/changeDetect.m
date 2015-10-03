function [rgb_mean,rgb_sd]=changeDetect
warning off all
%--------------------Initial Setup---------------------------%
vid = videoinput('winvideo',1,'YUY2_640X480');
triggerconfig(vid,'manual');
set(vid,'FramesPerTrigger',1);
set(vid,'TriggerRepeat', Inf);
start(vid)
%------------------------------------------------------------%

%------------------------------------------------------------%
f1=10;

for i=1:f1
    trigger(vid)
    temp=getdata(vid);
    rgb(:,:,:,i)=double(temp(50:400,150:480,:));  
    
    norm=rgb(:,:,1,i).^2+rgb(:,:,2,i).^2+rgb(:,:,3,i).^2;
    norm=(norm).^(.5);
    rgb(:,:,1,i)=rgb(:,:,1,i)./norm;
end

rgb_mean=sum(rgb,4)/f1;
rgb_sd=std(rgb,0,4);
%------------------------------------------------------------%

stop(vid)
end