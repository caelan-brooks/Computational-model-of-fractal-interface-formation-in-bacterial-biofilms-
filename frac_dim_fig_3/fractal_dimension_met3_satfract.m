%% every iteration,the motile cells that are on the edge of the biofilm swim away 
% finds the number of cells, ratio, and proportion over time 

clear all; close all; tic;
ratio_sum = [];
Rsq_sum = [];

for M=10:2:30
max_eps_push = M;
for N=1:20
global mfill bio_matrix area_m area_s row_num col_num neighbors perim perimeter rad_r theta_r theta radius area init_row init_col
sd = 5;
%number of rows and columns in grid
row_num=100;
col_num=100;

init_row = 50;
init_col = 50;

maxiter= 16000; %1000; %max iterations 
piter=25000; %plots per x iterations 
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


%setting initial growth times
for i=1:row_num
    for j=1:col_num
        if bio_matrix0(i,j)~=0
           grow_time(i,j)=((randn/10) + 1)+rand; %making initial cell growth times out of phase 
        end 
    end 
end  

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

end 


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

rsaves(N,1) = Rsq;
fractdimsave(N,1) = fract;

end 

r2 = mean(rsaves);
fractdim = mean(fractdimsave);
fracttrack(M,1) = fractdim; 
r2track(M,1) = r2;
M

end 
 
% save fractsat10to12 fracttrack