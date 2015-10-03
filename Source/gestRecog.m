function gestRecog(ref)
warning off all

%--------------------Initial Setup---------------------------%
vid = videoinput('winvideo',1,'YUY2_640X480');
triggerconfig(vid,'manual');
set(vid,'FramesPerTrigger',1);
set(vid,'TriggerRepeat', Inf);
%------------------------------------------------------------%


%------------------Aquire Initial Frame----------------------%
%------------To Determine Valid Number of Clusters-----------%
'Detemining Valid Number of Clusters'                                       %#ok<NOPRT>
temp=input('Keep Your Hand Infront of Camera');                           %#ok<NASGU>
CF=getsnapshot(vid);                                                          %Aquiring Image from Device
'Frame Aquired'                                                           %#ok<NOPRT>
CF=CF(1:441,1:411,:);                                                     %Image Cropping

CF(:,:,2)=medfilt2(CF(:,:,2),[5,5],'symmetric');                          %Color Smoothening
CF(:,:,3)=medfilt2(CF(:,:,3),[5,5],'symmetric');                          %

temp(:,:,1)=imresize(CF(:,:,2),[90 ,120]); 
temp(:,:,2)=imresize(CF(:,:,2),[90 ,120]);                                  %Resizing Image 
temp(:,:,3)=imresize(CF(:,:,3),[90 ,120]);                                  %
CF=temp;
ab = double(CF(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

for i=2:3
    nColors =i;
    [cluster_idx cluster_center]= kmeans(ab,nColors,'distance','cosine','emptyaction','singleton');
    [silh3,h] = silhouette(ab,cluster_idx,'cosine');
    result(i)=mean(silh3);
end

nColors=find(result==max(result))                                         %Valid Number of Clusters
%------------------------------------------------------------%
    

'Number of Valid Clusters : '                                             %#ok<NOPRT>
nColors                                                               %#ok<NOPRT>
[mean_frame sd]=changeDetect;
sd=3*sd;
up=mean_frame+sd;
down=mean_frame-sd;

start(vid) 
for idx=1:75
    trigger(vid)
    CF=getdata(vid);                                                      %Aquiring Frame
    CF=CF(50:400,150:480,:);                                              %Image Cropping
    CF(:,:,2)=medfilt2(CF(:,:,2),[5,5],'symmetric');                      %Color Smoothening
    CF(:,:,3)=medfilt2(CF(:,:,3),[5,5],'symmetric');                      %
    
    %------------------------Change Detection Mask------------------------%
    temp=ycbcr2rgb(CF);
    temp=double(temp);
    norm=temp(:,:,1).^2+temp(:,:,2).^2+temp(:,:,3).^2;
    norm=(norm).^(.5);
    temp(:,:,1)=temp(:,:,1)./norm;
    temp(:,:,2)=temp(:,:,2)./norm;
    temp(:,:,3)=temp(:,:,3)./norm;
    
    test=(temp>down)&(temp<up);
    CDM=test(:,:,1)&test(:,:,2)&test(:,:,3);
    CDM=~CDM;
    
    %---------------------------------------------------------------------%
   
    
    %----------------------------Clustering-------------------------------%
    nrows = size(CF,1);
    ncols = size(CF,2);
    ab = double(CF(:,:,2:3));
    ab = reshape(ab,nrows*ncols,2);
    [cluster_idx cluster_center]= kmeans(ab,nColors,'distance','cosine','emptyaction','singleton','replicates',1);
    pixel_labels = reshape(cluster_idx,nrows,ncols);
    
    %-----Detemining The Cluster Closest To Skin Color Distribution-------%
    
    skin_cluster(1)=150;   skin_cluster(2)=110;                           %Centre Of Skin Color Distribution In CB-CR Color Space
    cluster_center(:,1)=cluster_center(:,1)-skin_cluster(1);
    cluster_center(:,2)=cluster_center(:,2)-skin_cluster(2);
    [temp,cluster_dist]=cart2pol(cluster_center(:,2),cluster_center(:,1));
    SDM=(pixel_labels==find(min(cluster_dist)==cluster_dist));
    SDM=~SDM;
    
    %---------------------------------------------------------------------%
    
    %---------------------------------------------------------------------%
    SDM=SDM&CDM;
    SDM=imfill(SDM,'holes');
    SDM=bwlabel(SDM);
    s=regionprops(SDM,'BoundingBox','Area','Orientation');
    n=size(s,1);
    max_area=0;  
    k=1;
    for i=1:n
        if max_area<s(i,1).Area
            max_area=s(i,1).Area;k=i;
        end
    end
    
    if size(s)~=0
        x1=s(k,1).BoundingBox(1);   x2=s(k,1).BoundingBox(1)+s(k,1).BoundingBox(3);
        y1=s(k,1).BoundingBox(2);   y2=s(k,1).BoundingBox(2)+s(k,1).BoundingBox(4);
    else continue;
    end
    angle=s(k,1).Orientation;
    CF=ycbcr2rgb(CF);
    CF(y1:y2,x1:x1,:)=0;  CF(y1:y2,x2:x2,:)=0;  
    CF(y1:y1,x1:x2,:)=0;  CF(y2:y2,x1:x2,:)=0;
    SDM=SDM(y1:y2-30,x1:x2-1);                                            %Cropping Image to Fit the Boundig Box
    if size(SDM)>0
        SDM=imresize(SDM,[280,240]);                                          %Resizing Image to  a Standard Size for Matching
    else continue
    end
    imshow(CF)
    [y,x]=find(SDM~=0);                                                   %Determining the Centroid of the Cropped Image
    cx=sum(x)/sum(sum(SDM~=0));                                           %
    cy=sum(y)/sum(sum(SDM~=0));                                           %
    %----------------Extracting Boundary from the SDMary Image------------%
    bound=bwboundaries(SDM,'noholes');                                                
    n=size(bound,1);
    max_bound=0;k=1;
    for i=1:n
        s=size(bound{i,1},1);
        if max_bound<s
            max_bound=s;k=i;
        end
    end
    bound=bound{k,1};
    bound(:,1)=bound(:,1)-cy;
    bound(:,2)=bound(:,2)-cx;
    z=frdescp(bound);
    z30_bound=ifrdescp(z,30);
    z30_im=bound2im(z30_bound);
    %---------------------------------------------------------------------%
    
    %--------------------Obtaining Signature of Boundry-------------------%
    [theta,rho]=cart2pol(z30_bound(:,2),z30_bound(:,1));
    rho=imresize(rho,[1000,1]);
   
    co=find(rho==min(rho));                                               %Transformations for Rotation Invariance
    rho2(1:1000-co+1)=rho(co:1000);                                       %
    rho2(1000-co+2:1000)=rho(1:co-1);                                     %
    rho2=rho2-min(rho2);
    rho2=rho2/(max(rho2)/100);                                            %Transformation for Size Invariance  
    %---------------------------------------------------------------------%
   
    %---------Determinig Eucledian Distance Between Gestures--------------%
    ed=zeros(1,83);
    for i=1:83
        ed(i)=sum((rho2(:)-ref(:,i)).^2);
    end
    ed=sqrt(ed);
    result=find(ed==min(ed));
    %---------------------------------------------------------------------%
     
    %----------------------Matching the Gesture --------------------------%
    if ed(result)<500
        switch result
            case {1 ,26 ,27 ,54}
                text(cx,cy,'FIVE'); 
            case {2 ,28 ,29 ,55, 56}
                 text(cx,cy,'FOUR'); 
            case {3}
                text(cx,cy,'L shape'); 
            case {4}
                text(cx,cy,'THUMB'); 
            case {5, 30 ,31}
                text(cx,cy,'CLOSED HAND'); 
            case {6 , 60, 61}
                text(cx,cy,'ONE'); 
            case {7 ,32 ,62 ,63}
                text(cx,cy,'TWO'); 
            case {8, 57,58 ,59}
                text(cx,cy,'THREE'); 
            case {9}
                text(cx,cy,'SEVEN'); 
            case {10 ,33}
                text(cx,cy,'EIGHT'); 
            case {11, 34 ,35}
                text(cx,cy,'NINE'); 
            case {12}
                text(cx,cy,'TEN'); 
            case {13, 64, 51, 52 ,53}
                text(cx,cy,'ELEVEN'); 
            case {14 ,65, 66}
                text(cx,cy,'TWELVE'); 
            case {15, 36,37}
                text(cx,cy,'THIRTEEN'); 
            case {16, 49 ,50}
                text(cx,cy,'FOURTEEN'); 
            case {17}
                text(cx,cy,'FIFTEEN'); 
            case {18}
                text(cx,cy,'Gun'); 
            case {19 ,40 ,41}
                text(cx,cy,'FIST'); 
            case {20}
                text(cx,cy,'PINKY'); 
            case {21}
                text(cx,cy,'SIXTEEN'); 
            case {22, 67, 68 ,69}
                text(cx,cy,'SEVENTEEN'); 
            case {23, 47, 48}
                text(cx,cy,'EIGHTEEN'); 
            case {24 ,44 ,45, 46}
                text(cx,cy,'SIX'); 
            case {25 ,38, 39, 42 ,43}
                text(cx,cy,'PINKY'); 
            case {70, 71, 72}
                if abs(angle)>=0 && abs(angle)<=19
                    text(cx,cy,'FIVE ROTATED'); 
                end
            case{73, 74}
                if abs(angle)>=17 && abs(angle)<=30
                    text(cx,cy,'ONE ROTATED'); 
                end
            case{75 ,76 ,77}
                if abs(angle)>=0 && abs(angle)<=20
                    text(cx,cy,'TWO ROTATED'); 
                end
            case{78, 79, 80}
                if abs(angle)>=0 && abs(angle)<=20
                    text(cx,cy,'THREE ROTATED'); 
                end
            case{81 ,82 ,83}
                if abs(angle)>=0 && abs(angle)<=20
                    text(cx,cy,'FOUR ROTATED'); 
                end
        end
    end 
    %---------------------------------------------------------------------%
end

stop(vid)
end
