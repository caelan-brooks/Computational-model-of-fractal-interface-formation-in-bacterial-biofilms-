%% every iteration,the motile cells that are on the edge of the biofilm swim away 
% finds the number of cells, ratio, and proportion over time 

clear all; close all; tic;
for M=1:14
ratio_sum = [];
for N=1:10
global bio_matrix area_m area_s row_num col_num neighbors perim perimeter rad_r theta_r theta radius area

%number of rows and columns in grid
row_num=100;
col_num=100;

init_row = 50;
init_col = 50;

maxiter= 4000; %1000; %max iterations 
piter=10; %plots per x iterations 
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
max_eps_push = M;

%for cell growth tracking
row_fill=[];
col_fill=[];
time_keeper=[];
bio_height=[];
bio_width=[];
radial_avg = [];
radii = [];
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
                eps(i,j) = ((randn/5) + eps_layer)+rand;
        end 
        if radius(i,j)<=R2
            bio_matrix(i,j) = 1;
            bio_matrix0(i,j) = 1;
            eps(i,j) = ((randn/5) + 1)+rand;
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
%               eps(i,j) = 1;
%            else 
%               bio_matrix(i,j) = 2;
%               bio_matrix0(i,j) = 2;
%               eps(i,j) = 2;
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



find_area()
time_keeper = [time_keeper, 0];
% num_motile = [num_motile, area_s];
% num_matrix = [num_matrix, area_m];
% pro_m = [pro_m, area_m/(area_s + area_m)];
% pro_s = [pro_s, area_s/(area_s + area_m)];

num_neighbors_edge()
find_perimeter_edge()
edgetrack = [];
for i=1:row_num
    for j=1:col_num 
        if perimeter(i,j)==1
            edgetrack = [edgetrack ; rad(i,j)];
        end 
    end
end 
        
ratio_time = [ratio_time; mean(edgetrack)];

%%
for k=0:maxiter

%allowing motile cells to leave if they are on the edge 
for u=1:row_num
    for p=1:col_num      
        if bio_matrix0(u,p) == 1 && (bio_matrix0(u+1,p) == 0 || bio_matrix0(u-1,p) == 0 || bio_matrix0(u,p+1) == 0 || bio_matrix0(u,p-1) == 0)
           bio_matrix(u,p) = 2;
           eps(u,p) = ((randn/5) + eps_layer);
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
                               eps(i,j-1)=((randn/10) + 1)+rand;
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
                               eps(i,j+1)= ((randn/5) + 1);
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
                               eps(i-1,j)= ((randn/5) + 1);
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
                               eps(i+1,j)= ((randn/5) + 1);
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
                              eps(i,j-1)=((randn/5) + eps_layer);
                              grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_left<=max_eps_push  
                               for left_shift=min_col:j-1
                                   bio_matrix(i,left_shift-1) = bio_matrix0(i,left_shift);
                                   grow_time(i,left_shift-1) = grow_time(i,left_shift);
                                   eps(i,left_shift-1) = eps(i,left_shift);
                               end
                               bio_matrix(i,j-1)=2;
                               eps(i,j-1)=((randn/5) + eps_layer);
                               grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes right 
                       if num >= .25 && num < .50 
                           if bio_matrix0(i,j+1)==0
                              %disp("right")
                              bio_matrix(i,j+1)=2;
                              eps(i,j+1)=((randn/5) + eps_layer);
                              grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_right<=max_eps_push  
                               for right_shift=max_col:-1:j+1
                                   bio_matrix(i,right_shift+1) = bio_matrix0(i,right_shift);
                                   grow_time(i,right_shift+1) = grow_time(i,right_shift);
                                   eps(i,right_shift+1) = eps(i,right_shift);
                               end
                               bio_matrix(i,j+1)=2;
                               eps(i,j+1)=((randn/5) + eps_layer);
                               grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes down
                       if num >= .50 && num < .75 
                           if bio_matrix0(i-1,j)==0
                              %disp("down")
                              bio_matrix(i-1,j)=2;
                              eps(i-1,j)=((randn/5) + eps_layer);
                              grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_down<=max_eps_push 
                               for down_shift=min_row:i-1
                                   bio_matrix(down_shift-1,j) = bio_matrix0(down_shift,j);
                                   grow_time(down_shift-1,j) = grow_time(down_shift,j);
                                   eps(down_shift-1,j) = eps(down_shift,j);
                               end
                               bio_matrix(i-1,j)=2;
                               eps(i-1,j)=((randn/5) + eps_layer);
                               grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes up 
                       if num > .75 
                           if bio_matrix0(i+1,j)==0
                              %disp("up")
                              bio_matrix(i+1,j)=2;
                              eps(i+1,j)=((randn/5) + eps_layer);
                              grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_up<=max_eps_push  
                               for up_shift=max_row:-1:i+1
                                   bio_matrix(up_shift+1,j) = bio_matrix0(up_shift,j);
                                   grow_time(up_shift+1,j) = grow_time(up_shift,j);
                                   eps(up_shift+1,j) = eps(up_shift,j);
                               end
                               bio_matrix(i+1,j)=2;
                               eps(i+1,j)=((randn/5) + eps_layer);
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
% figure(1)
% pcolor(bio_matrix)
% title(['time= ',num2str(time)])
% colorbar
% drawnow 


time_keeper = [time_keeper, time];
find_area()
num_motile = [num_motile, area_s];
num_matrix = [num_matrix, area_m];
pro_m = [pro_m, area_m/(area_s + area_m)];
pro_s = [pro_s, area_s/(area_s + area_m)];

num_neighbors_edge()
find_perimeter_edge()
% fractal_unwrap()
%ratio = perim/sqrt(area);
edgetrack = [];
for i=1:row_num
    for j=1:col_num 
        if perimeter(i,j)==1
            edgetrack = [edgetrack ; rad(i,j)];
        end 
    end
end 
        
ratio_time = [ratio_time; mean(edgetrack)];

end
end 


 
%%

% sat = ratio_time(length(ratio_time),:);
% ratio_time = ratio_time./sat;
%%

% k = length(ratio_time);
% while ratio_time(k,:) == 1
%     sat_spot = k; 
%     k = k-1;
% end


% rat_round = round(ratio_time,1);
% time_r = round(sat/2,1);
% for ind = 1:length(ratio_time)
%     if rat_round(ind) == time_r
%         time_spot = ind;
%     elseif rat_round(ind) ==round(time_r+0.1,1)
%         time_spot = ind;
%     elseif rat_round(ind) == round(time_r-0.1,1)
%         time_spot = ind;
%     end 
% end

%sat_time = time_keeper(:,sat_spot);
%sat_time = time_keeper(:,ind);
time_keeper = time_keeper; %./sat_time;
%%       
        

% figure(M)
% plot(time_keeper, ratio_time,'.')
% hold on
% title('time vs perimeter/sqrt(area)')
% xlabel('time')
% ylabel('ratio')
% grid on 


time_keeper = round(time_keeper,precision);
time_fill = 0:1*10^(-precision):6;
ratio_fill = zeros(length(time_fill),1);

for t=1:length(time_keeper)
    for p=1:length(time_fill)
        if time_keeper(:,t) == time_fill(:,p)
            ratio_fill(p:length(ratio_fill),:) =  ratio_time(t,:);
        end 
    end 
end 


% if N == 1
%     ratio_sum = zeros(length(ratio_fill),1);
% end 

ratio_sum = [ratio_sum, ratio_fill];
% ratio_sum = ratio_sum + ratio_fill;


end 
M
% ratio_final = ratio_sum./N;

%%


toc;

% 
% save curve_14 time_fill ratio_final


%%

ratio_avg = zeros(length(ratio_sum(:,1)),1);
for k=1:length(ratio_sum(:,1))
    ratio_avg(k,1) = mean(ratio_sum(k,:));
    ratio_min(k,1) = min(ratio_sum(k,:));
    ratio_max(k,1) = max(ratio_sum(k,:));
end 
plot(time_fill,ratio_avg,'o')
timesave(:,M) = time_fill;
linslope(:,M) = ratio_avg;
hold on 
 [slope,S] = polyfit(time_fill(200:500),ratio_avg(200:500),1);
% 
grate(M,1) = slope(1);

end 
% 
% time_f = [];
% ratio_a = [];
% err_f = [];
% err = ratio_max - ratio_min;
% for ps=1:200:length(ratio_sum(:,1))
% time_f = [time_f, time_fill(1,ps)];
% ratio_a = [ratio_a, ratio_avg(ps,1)];
% err_f = [err_f, err(ps,1)];
% end
% %%
% figure(5)
% errorbar(time_f,ratio_a,err_f)
% title('time vs perimeter/sqrt(area) with error')
% xlabel('time')
% ylabel('ratio')
% grid on 
% 
% 
% 
% save eps_layer_2 
%save edgepush16 time_fill ratio_avg

x=2:0.1:5;
y=slope(1)*x+slope(2);
plot(x,y,'LineWidth',2);


figure(10)
for i=2:2:14 
    plot(timesave(1:20:601,i),linslope(1:20:601,i),'o-','LineWidth',2,'MarkerSize', 6)
    hold on 
end 

plot(grate)
%save linear_edge_growth_data timesave linslope
%save growth_rate_SC grate