%% every iteration,the motile cells that are on the edge of the biofilm swim away 
% finds the number of cells, ratio, and proportion over time 

clear all; close all; tic;
ratio_sum = [];
Rsq_sum = [];
area_sum = [];
mottot_sum = [];

for N=1:100
global mfill bio_matrix area_m area_s row_num col_num neighbors perim edge_perimeter perimeter rad_r theta_r theta radius area init_row init_col mottot

sd = 5;
%number of rows and columns in grid
row_num=100;
col_num=100;

init_row = 50;
init_col = 50;

maxiter= 16000; %1000; %max iterations 
piter=10; %plots per x iterations 
% piter2=10000; %plots per x iterations 
% piter3=250; %plots per x iterations 
precision = 2;

%initial Radii 
R1 = 7;
R2 = 5;
R0 = 6;
eps_layer = 2;

%setting up biofilm and growth matrices 
bio_matrix = zeros(row_num, col_num);
neighbors = zeros(row_num, col_num);
perimeter = zeros(row_num, col_num);
bio_matrix0 = zeros(row_num, col_num);
grow_time = zeros(row_num, col_num);
eps = zeros(row_num, col_num);
rad = zeros(row_num, col_num);
radius = zeros(row_num, col_num);
theta = zeros(row_num, col_num);



%number of cells that can be pushed in order for a surrounded cell to grow 
%max_push = 8; 
max_eps_push = 18;

%for cell growth tracking
row_fill=[];
col_fill=[];
time_keeper=[];
bio_height=[];
bio_width=[];
radial_avg = [];
radii = [];
Rtime = [];
ratio_time = [];
num_motile = [];
num_matrix = [];
rad_r = [];
theta_r = [];
pro_m = [];
pro_s = [];
motsave = [];
areasave = [];
numfactors = [];


%setting up radius matrix 
for i=1:row_num
    for j=1:col_num
        y = ((i-double(init_row)));
        x = ((j-double(init_col)));
        rad(i,j) = sqrt((y)^2 + (x)^2);
        theta(i,j) = atand(y/x);
        if i>init_row && j<init_col
            theta(i,j) = theta(i,j) + 180;
        end 
        if i<=init_row && j<init_col
            theta(i,j) = theta(i,j) + 180;
        end 
        if i<init_row && j>=init_col
            theta(i,j) = theta(i,j) + 360;
        end     
    end 
end  

for i=1:row_num
    for j=1:col_num
        for n=.5:360.5 %(round(sqrt(2*(double(init_row)^2)))+.5)
            if rad(i,j)>n && rad(i,j)<(n+1)
                radius(i,j)=n+.5;
            end  
        end 
    end 
end 

% figure(1) 
% pcolor(radius)
% colorbar
% figure(2)
% pcolor(theta)
% colorbar 


%setting power spectrum plot matrices 
power = double.empty;
phi = double.empty;


%initializing cell types in initial radii 
for i=1:row_num
    for j=1:col_num
        num = rand;
        if radius(i,j)<=R1 && radius(i,j)>R2
                bio_matrix(i,j) = 2;
                bio_matrix0(i,j) = 2;
                eps(i,j) = ((randn/sd) + eps_layer)+rand;
        end 
        if radius(i,j)<=R2
            bio_matrix(i,j) = 1;
            bio_matrix0(i,j) = 1;
            eps(i,j) = ((randn/sd) + 1)+rand;
        end 
    end
end

% if you want to make it initially random 
% for i=1:row_num
%     for j=1:col_num
%         num = rand;
%         if radius(i,j)<=R0
%            num = rand; 
%            if rand<.5  
%               bio_matrix(i,j) = 1;
%               bio_matrix0(i,j) = 1;
%               eps(i,j) = ((randn/sd) + 1)+rand;
%            else 
%               bio_matrix(i,j) = 2;
%               bio_matrix0(i,j) = 2;
%               eps(i,j) = ((randn/sd) + eps_layer)+rand;
%            end
%         end 
%         
%     end
% end

%setting initial growth times
for i=1:row_num
    for j=1:col_num
        if bio_matrix0(i,j)~=0
           grow_time(i,j)=((randn/10) + 1)+rand; %making initial cell growth times out of phase 
        end 
    end 
end  



% find_area()
time_keeper = [time_keeper, 0];

% num_motile = [num_motile, area_s];
% num_matrix = [num_matrix, area_m];
% pro_m = [pro_m, area_m/(area_s + area_m)];
% pro_s = [pro_s, area_s/(area_s + area_m)];

perim_fill()
num_neighbors_fill()
find_perimeter_fill()

perimi = [];
perimj = [];
for i=1:row_num
    for j=1:col_num
        if perimeter(i,j)==1
            perimi = [perimi; i];
            perimj = [perimj; j];
        end 
    end 
end 

minj = min(perimj);
maxj = max(perimj);
mini = min(perimi);
maxi = max(perimi);
diffj = maxj-minj+1;
diffi = maxi-mini+1;
if rem(diffj,2)~=0
    diffj=diffj+1;
end 
if rem(diffi,2)~=0
    diffi=diffi+1;
end 

if diffj>diffi
    diff=diffj;
    mini = mini-(diffj-diffi)/2;
elseif diffi>diffj 
    diff=diffi;
    minj = minj-(diffi-diffj)/2;
else 
    diff=diffj;
end 


fac=[];
for l=1:diff
    if rem(diff,l)==0 %&& rem(diffi,l)==0
        fac=[fac;l];
    end 
end 

for ll=length(fac):-1:1
    if fac(ll)>(0.5*diff)
        fac(ll) = [];
    end 
end 

factors = fac;

numfactors = [numfactors; length(factors)];

fact2 = diff./factors;
fractdim=[];
for kk=1:length(factors)
    keep = [];
    for n=1:fact2(kk)
        for m=1:fact2(kk)
        dude=0;
        %n
        for ii=mini-1+factors(kk)*n-(factors(kk)-1):mini-1+factors(kk)*n
            for jj=minj-1+(factors(kk)*m-(factors(kk)-1)):minj-1+(factors(kk)*m)
                %fprintf('%d %d ',ii,jj)
                if perimeter(ii,jj)==1
                   dude=dude+1;
                end
            end 
        end 
        if dude~=0
            dude=1;
        end 
        keep = [keep; dude];
        dude=0;
        end 
    end 
    fractdim(kk)=sum(keep);
end 

[slope,S] = polyfit(log(factors),log(fractdim),1);
Rsq = 1 - (S.normr/norm((log(fractdim)) - mean((log(fractdim)))))^2;
fract = -slope(1);

Rtime = [Rtime; Rsq];     
ratio_time = [ratio_time; fract];

motsave = [motsave; mottot];
areasave = [areasave; area];
% fract
%%
for k=0:maxiter

%allowing motile cells to leave if they are on the edge 
for u=1:row_num
    for p=1:col_num      
        if bio_matrix0(u,p) == 1 && (bio_matrix0(u+1,p) == 0 || bio_matrix0(u-1,p) == 0 || bio_matrix0(u,p+1) == 0 || bio_matrix0(u,p-1) == 0)
           bio_matrix(u,p) = 2;
           eps(u,p) = ((randn/sd) + eps_layer);
        end 
    end 
end


% while medge>0
%     for i=1:row_num
%         for j=1:col_num
%             
% num_neighbors_edge()
% find_perimeter_edge()
% % 
% % figure(2)
% % pcolor(edge_perimeter)
% 
% %allowing motile cells to leave if they are on the edge 
% for u=1:row_num
%     for p=1:col_num      
%         if bio_matrix(u,p) == 1 && edge_perimeter(u,p)==1
%              bio_matrix(u,p) = 0;
%         end 
%     end 
% end
% end 
        [i,j] = find(grow_time==min(grow_time(grow_time> 0)));
        eps_up = 0;
        eps_down = 0;
        eps_left = 0;
        eps_right = 0;
 %------------------------------------------------------------------------%       
        if i>1 && i<row_num && j>1 && j<col_num 
%                    i
%                    j
                  %-------------------------------------------------%
                           %finds first zero above 
                           a=i;
                           while bio_matrix0(a,j)~=0 && a<row_num
                                 a=a+1; 
                           end  
                           max_row=a-1;
                           length_above = max_row-i;
                              %amount of eps above 
                              for n = max_row:-1:i+1
                                  eps_up = eps_up + eps(n,j);
                              end
                           %finds first zero below
                           b=i;
                           while bio_matrix0(b,j)~=0 && b>1
                                 b=b-1; 
                           end  
                           min_row=b+1;
                           length_below = i-min_row;
                              %amount of eps below 
                              for m = min_row:i-1
                                  eps_down = eps_down + eps(m,j);
                              end
                           %finds first zero to the right 
                           c=j;
                           while bio_matrix0(i,c)~=0 && c<col_num
                                 c=c+1;
                           end  
                           max_col=c-1;
                           length_right = max_col-j;
                              %amount of eps to right 
                              for l = max_col:-1:j+1
                                  eps_right = eps_right + eps(i,l);
                              end
                           %finds first zero to the left 
                           d=j;
                           while bio_matrix0(i,d)~=0 && d>1
                                 d=d-1; 
                           end  
                           min_col=d+1;
                           length_left = j-min_col;
                              %amount of eps to left
                              for p = min_col:j-1
                                  eps_left = eps_left + eps(i,p);
                              end
                   %-------------------------------------------------%
                   num = rand;
                   if bio_matrix0(i,j)==1
                       %goes left 
                       if num < .25 %left_prob
                           if eps_left<=max_eps_push 
                               for left_shift=min_col:j-1
                                   bio_matrix(i,left_shift-1) = bio_matrix0(i,left_shift);
                                   grow_time(i,left_shift-1) = grow_time(i,left_shift);
                                   eps(i,left_shift-1) = eps(i,left_shift);
                               end
                               bio_matrix(i,j-1)=1;
                               eps(i,j-1)=((randn/sd) + 1);
                               grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes right 
                       if num >= .25 && num < .50 
                           if eps_right<=max_eps_push  
                               for right_shift=max_col:-1:j+1
                                   bio_matrix(i,right_shift+1) = bio_matrix0(i,right_shift);
                                   grow_time(i,right_shift+1) = grow_time(i,right_shift);
                                   eps(i,right_shift+1) = eps(i,right_shift);
                               end
                               bio_matrix(i,j+1)=1;
                               eps(i,j+1)=((randn/sd) + 1);
                               grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes down
                       if num >= .50 && num < .75
                           if eps_down<=max_eps_push 
                               for down_shift=min_row:i-1
                                   bio_matrix(down_shift-1,j) = bio_matrix0(down_shift,j);
                                   grow_time(down_shift-1,j) = grow_time(down_shift,j);
                                   eps(down_shift-1,j) = eps(down_shift,j);
                               end
                               bio_matrix(i-1,j)=1;
                               eps(i-1,j)=((randn/sd) + 1);
                               grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes up 
                       if num > .75 
                           if eps_up<=max_eps_push  
                               for up_shift=max_row:-1:i+1
                                   bio_matrix(up_shift+1,j) = bio_matrix0(up_shift,j);
                                   grow_time(up_shift+1,j) = grow_time(up_shift,j);
                                   eps(up_shift+1,j) = eps(up_shift,j);
                               end
                               bio_matrix(i+1,j)=1;
                               eps(i+1,j)=((randn/sd) + 1);
                               grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                   end 
                  %______________________________________________________________________%
                  if bio_matrix0(i,j)==2
                       %goes left 
                       if num < .25 
                           if bio_matrix0(i,j-1)==0
                              bio_matrix(i,j-1)=2;
                              eps(i,j-1)=((randn/sd) + eps_layer);
                              grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_left<=max_eps_push  
                               for left_shift=min_col:j-1
                                   bio_matrix(i,left_shift-1) = bio_matrix0(i,left_shift);
                                   grow_time(i,left_shift-1) = grow_time(i,left_shift);
                                   eps(i,left_shift-1) = eps(i,left_shift);
                               end
                               bio_matrix(i,j-1)=2;
                               eps(i,j-1)=((randn/sd) + eps_layer);
                               grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes right 
                       if num >= .25 && num < .50 
                           if bio_matrix0(i,j+1)==0
                              %disp("right")
                              bio_matrix(i,j+1)=2;
                              eps(i,j+1)=((randn/sd) + eps_layer);
                              grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_right<=max_eps_push  
                               for right_shift=max_col:-1:j+1
                                   bio_matrix(i,right_shift+1) = bio_matrix0(i,right_shift);
                                   grow_time(i,right_shift+1) = grow_time(i,right_shift);
                                   eps(i,right_shift+1) = eps(i,right_shift);
                               end
                               bio_matrix(i,j+1)=2;
                               eps(i,j+1)=((randn/sd) + eps_layer);
                               grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes down
                       if num >= .50 && num < .75 
                           if bio_matrix0(i-1,j)==0
                              %disp("down")
                              bio_matrix(i-1,j)=2;
                              eps(i-1,j)=((randn/sd) + eps_layer);
                              grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_down<=max_eps_push 
                               for down_shift=min_row:i-1
                                   bio_matrix(down_shift-1,j) = bio_matrix0(down_shift,j);
                                   grow_time(down_shift-1,j) = grow_time(down_shift,j);
                                   eps(down_shift-1,j) = eps(down_shift,j);
                               end
                               bio_matrix(i-1,j)=2;
                               eps(i-1,j)=((randn/sd) + eps_layer);
                               grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes up 
                       if num > .75 
                           if bio_matrix0(i+1,j)==0
                              %disp("up")
                              bio_matrix(i+1,j)=2;
                              eps(i+1,j)=((randn/sd) + eps_layer);
                              grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_up<=max_eps_push  
                               for up_shift=max_row:-1:i+1
                                   bio_matrix(up_shift+1,j) = bio_matrix0(up_shift,j);
                                   grow_time(up_shift+1,j) = grow_time(up_shift,j);
                                   eps(up_shift+1,j) = eps(up_shift,j);
                               end
                               bio_matrix(i+1,j)=2;
                               eps(i+1,j)=((randn/sd) + eps_layer);
                               grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                   end
                  %______________________________________________________________________%
                   time = grow_time(i,j);
                   grow_time(i,j)= grow_time(i,j)+((randn/10) + 1);                  
        end 
                
bio_matrix0 = bio_matrix;


if rem(k,piter)==0


time_keeper = [time_keeper, time];
find_area()
num_motile = [num_motile, area_s];
num_matrix = [num_matrix, area_m];
pro_m = [pro_m, area_m/(area_s + area_m)];
pro_s = [pro_s, area_s/(area_s + area_m)];

% num_neighbors_edge()
% find_perimeter_edge()
% fractal_unwrap()
%ratio = perim/sqrt(area);
perim_fill()
num_neighbors_fill()
find_perimeter_fill()

perimi = [];
perimj = [];
for i=1:row_num
    for j=1:col_num
        if perimeter(i,j)==1
            perimi = [perimi; i];
            perimj = [perimj; j];
        end 
    end 
end 

minj = min(perimj);
maxj = max(perimj);
mini = min(perimi);
maxi = max(perimi);
diffj = maxj-minj+1;
diffi = maxi-mini+1;
if rem(diffj,2)~=0
    diffj=diffj+1;
end 
if rem(diffi,2)~=0
    diffi=diffi+1;
end 

if diffj>diffi
    diff=diffj;
    mini = mini-(diffj-diffi)/2;
elseif diffi>diffj 
    diff=diffi;
    minj = minj-(diffi-diffj)/2;
else 
    diff=diffj;
end 


fac=[];
for l=1:diff
    if rem(diff,l)==0 %&& rem(diffi,l)==0
        fac=[fac;l];
    end 
end 

for ll=length(fac):-1:1
    if fac(ll)>(0.5*diff)
        fac(ll) = [];
    end 
end 

factors = fac;

numfactors = [numfactors; length(factors)];
% dim

fact2 = diff./factors;
fractdim=[];
for kk=1:length(factors)
    keep = [];
    for n=1:fact2(kk)
        for m=1:fact2(kk)
        dude=0;
        %n
        for ii=mini-1+factors(kk)*n-(factors(kk)-1):mini-1+factors(kk)*n
            for jj=minj-1+(factors(kk)*m-(factors(kk)-1)):minj-1+(factors(kk)*m)
                %fprintf('%d %d ',ii,jj)
                if perimeter(ii,jj)==1
                   dude=dude+1;
                end
            end 
        end 
        if dude~=0
            dude=1;
        end 
        keep = [keep; dude];
        dude=0;
        end 
    end 
    fractdim(kk)=sum(keep);
end 

[slope,S] = polyfit(log(factors),log(fractdim),1);
Rsq = 1 - (S.normr/norm((log(fractdim)) - mean((log(fractdim)))))^2;
fract = -slope(1);
% length(factors)
% figure(k+1)
x=0:0.1:3;
y=slope(1)*x+slope(2);
% plot(x,y,'LineWidth',2);
% hold on 
% plot(log(factors),log(fractdim),'o','LineWidth',2,'MarkerSize', 10);
% xlabel('log(box side length)')
% ylabel('log(number of occupied boxes)')
% set(gca,'FontSize',14)
% legend('slope = ',num2str(slope(1)))
ratio_time = [ratio_time; fract];
Rtime = [Rtime; Rsq];
% fract
motsave = [motsave; mottot];
areasave = [areasave; area];


% if rem(k,piter2)==0
% % mymap = [1 1 1; 0 .75 0; 0.75 0 0.75];
% % figure(k+1)
% % pcolor(bio_matrix)
% % title(['fractal dimension= ',num2str(slope(1)),' time= ',num2str(time)])
% % colormap (mymap)
% % drawnow
% 
% % figure(k+1)
% % pcolor(perimeter)
% end 
% 

% mymap = [1 1 1; 0 .75 0; 0.75 0 0.75];
% figure(k+1)
% pcolor(bio_matrix)
% title(['fractal dimension= ',num2str(slope(1)),' time= ',num2str(time)])
% colormap (mymap)
% drawnow

end

end 


%%       
figure(4)
plot(time_keeper, ratio_time)
hold on
title('fractal dimension vs. time')
xlabel('time')
ylabel('fractal dimension')
grid on        

N



time_keeper = round(time_keeper,precision);
time_fill = 0:1*10^(-precision):13;
ratio_fill = zeros(length(time_fill),1);
Rsq_fill = zeros(length(time_fill),1);
area_fill = zeros(length(time_fill),1);
mottot_fill = zeros(length(time_fill),1);

for t=1:length(time_keeper)
    for p=1:length(time_fill)
        if time_keeper(:,t) == time_fill(:,p)
            ratio_fill(p:length(ratio_fill),:) =  ratio_time(t,:);
            Rsq_fill(p:length(Rsq_fill),:) =  Rtime(t,:);
            area_fill(p:length(area_fill),:) =  areasave(t,:);
            mottot_fill(p:length(mottot_fill),:) =  motsave(t,:);
        end 
    end 
end 


% if N == 1
%     ratio_sum = zeros(length(ratio_fill),1);
% end 
area_sum = [area_sum, area_fill];
mottot_sum = [mottot_sum, mottot_fill];
ratio_sum = [ratio_sum, ratio_fill];
Rsq_sum = [Rsq_sum, Rsq_fill];
% ratio_sum = ratio_sum + ratio_fill;



% ratio_final = ratio_sum./N;

%%
end 

toc;


ratio_avg = zeros(length(ratio_sum(:,1)),1);
Rsq_avg = zeros(length(ratio_sum(:,1)),1);
area_avg = zeros(length(ratio_sum(:,1)),1);
mottot_avg = zeros(length(ratio_sum(:,1)),1);
for kp=1:length(ratio_sum(:,1))
    Rsq_avg(kp,1) = mean(Rsq_sum(kp,:));
    area_avg(kp,1) = mean(area_sum(kp,:));
    mottot_avg(kp,1) = mean(mottot_sum(kp,:));
    ratio_avg(kp,1) = mean(ratio_sum(kp,:));
    ratio_min(kp,1) = min(ratio_sum(kp,:));
    ratio_max(kp,1) = max(ratio_sum(kp,:));
end 
figure(3)
plot(time_fill,ratio_avg) 
% plot(time_fill,ratio_avg./ratio_avg(length(ratio_avg))) 
% hold on 
% plot(time_fill,Rsq_avg) 
% 
% plot(time_keeper, (numfactors./numfactors(length(numfactors))./4)+0.75)



% figure(4)
% plot(time_fill,area_avg./mottot_avg(length(mottot_avg)))
% hold on 
% plot(time_fill,mottot_avg./mottot_avg(length(mottot_avg)))

%save shape_fill_200_push14 time_fill ratio_avg
%save fractdim20 time_fill ratio_avg

%fractsave = ratio_avg;

%%
% perimeter = zeros(row_num, col_num);
% for i=1:row_num
%     for j=1:col_num
%         if radius(i,j)<=10
%             perimeter(i,j)=1;
%         end 
% %         if i==50 && j>=15 && j<=86
% %             perimeter(i,j)=1;
% %         end 
%     end 
% end
% 
% % num_neighbors()
% % find_perimeter()
% 
% perimi = [];
% perimj = [];
% for i=1:row_num
%     for j=1:col_num
%         if perimeter(i,j)==1
%             perimi = [perimi; i];
%             perimj = [perimj; j];
%         end 
%     end 
% end 
% 
% minj = min(perimj);
% maxj = max(perimj);
% mini = min(perimi);
% maxi = max(perimi);
% diffj = maxj-minj+1;
% diffi = maxi-mini+1;
% if rem(diffj,2)~=0
%     diffj=diffj+1;
% end 
% if rem(diffi,2)~=0
%     diffi=diffi+1;
% end 
% 
% if diffj>diffi
%     diff=diffj;
% else 
%     diff=diffi;
% end 
% 
% fac=[];
% for l=1:diff
%     if rem(diffj,l)==0 && rem(diffi,l)==0
%         fac=[fac;l];
%     end 
% end 
% 
% for ll=length(fac):-1:1
%     if fac(ll)>10
%         fac(ll) = [];
%     end 
% end 
% factors = fac;
% 
% 
% fact2 = diff./factors;
% fractdim=[];
% for kk=1:length(factors)
%     keep = [];
%     for n=1:fact2(kk)
%         for m=1:fact2(kk)
%         dude=0;
%         %n
%         for ii=mini-1+factors(kk)*n-(factors(kk)-1):mini-1+factors(kk)*n
%             for jj=minj-1+(factors(kk)*m-(factors(kk)-1)):minj-1+(factors(kk)*m)
%                 %fprintf('%d %d ',ii,jj)
%                 if perimeter(ii,jj)==1
%                    dude=dude+1;
%                 end
%             end 
%         end 
%         if dude~=0
%             dude=1;
%         end 
%         keep = [keep; dude];
%         dude=0;
%         end 
%     end 
%     fractdim(kk)=sum(keep);
% end 
% 
% slope = polyfit(log(factors),log(fractdim),1);
% fract = -slope(1);
% plot(log(factors),log(fractdim),'o')
% hold on 
% x=0:0.1:3;
% y=slope(1)*x+slope(2);
% plot(x,y)
% xlabel('log(box side length)')
% ylabel('log(number of occupied boxes)')
% set(gca,'FontSize',12)
% 
% figure(2)
% plot(time_fill,ratio_avg)

%save fractmet3_4 time_fill ratio_avg
% fractsave = ratio_avg; 

%1101

%pcolor(bio_matrix)

% mymap = [1 1 1; 0 .75 0; 0.75 0 0.75];
% figure(N)
% pcolor(bio_matrix)
% colormap (mymap)
% drawnow
% 
% 
% 
% figure(N+100)
% mymap = [1 1 1; 0 0 0];
% pcolor(perimeter)
% colormap (mymap)
% drawnow

% save FDevoluSC18 time_fill ratio_avg