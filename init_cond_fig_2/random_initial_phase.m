clear all; close all; tic;
global bio_matrix area_m area_s row_num col_num neighbors perim perimeter rad_r theta_r theta radius area min_d kilt edge
%% every iteration,the motile cells that are on the edge of the biofilm swim away 
% finds the number of cells, ratio, and proportion over time 
% can find intensity as radius increases as well as fractal outline 

for R1=4:4 %8
for ii=1:300
%number of rows and columns in grid
row_num=100;
col_num=100;

%initial incoculation postion in grid 
init_row = 50;
init_col = 50;

maxiter= 10000; %1000; %max iterations 
piter=10000; %plots per x iterations 

%initial Radii 
%R1 = 7;
R2 = R1-2;
R0 = 6;

%how much harder it is to push matrix cells 
eps_layer = 2;

%amount of eps that can be pushed in order for a surrounded cell to grow 
max_eps_push = 12;


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


% % if you want to start the initial biofilm with random mixture 
for i=1:row_num
    for j=1:col_num
        num = rand;
        if radius(i,j)<=R1
           num = rand; 
           if rand<.5  
              bio_matrix(i,j) = 1;
              bio_matrix0(i,j) = 1;
              eps(i,j) = 1;
           else 
              bio_matrix(i,j) = 2;
              bio_matrix0(i,j) = 2;
              eps(i,j) = 2;
           end
        end 
        
    end
end

%initalizing biofilm with initial radii 
% for i=1:row_num
%     for j=1:col_num
%         num = rand;
%         if radius(i,j)<=R1 && radius(i,j)>R2
%                 bio_matrix(i,j) = 2;
%                 bio_matrix0(i,j) = 2;
%                 eps(i,j) = eps_layer;
%         end 
%         if radius(i,j)<=R2
%             bio_matrix(i,j) = 1;
%             bio_matrix0(i,j) = 1;
%             eps(i,j) = 1;
%         end 
%     end
% end


% %concentric circle start 
% for i=1:row_num
%     for j=1:col_num
%         for k=R1:-1:0
%             if mod(R1,2)==0
%                 if mod(k,2)==0 && radius(i,j)==k
%                         bio_matrix(i,j) = 2;
%                         bio_matrix0(i,j) = 2;
%                         eps(i,j) = eps_layer;
%                 end 
%                 if mod(k,2)~=0 && radius(i,j)==k
%                     bio_matrix(i,j) = 1;
%                     bio_matrix0(i,j) = 1;
%                     eps(i,j) = 1;
%                 end 
%             end 
%             
%             
%             
%              if mod(R1,2)~=0
%                 if mod(k,2)~=0 && radius(i,j)==k
%                         bio_matrix(i,j) = 2;
%                         bio_matrix0(i,j) = 2;
%                         eps(i,j) = eps_layer;
%                 end 
%                 if mod(k,2)==0 && radius(i,j)==k
%                     bio_matrix(i,j) = 1;
%                     bio_matrix0(i,j) = 1;
%                     eps(i,j) = 1;
%                 end 
%              end 
%             
%              
%              
%         end  
%     end
% end

%setting initial growth times
for i=1:row_num
    for j=1:col_num
        if bio_matrix0(i,j)~=0
           grow_time(i,j)=((randn/10) + 1)+rand; 
        end 
    end 
end  

mat0=0;
mot0=0;

for i=1:row_num
    for j=1:col_num
        if bio_matrix(i,j)==2
            mat0 = mat0+1;
        end 
        if bio_matrix(i,j)==1
            mot0 = mot0+1;
        end 
    end 
end 

%% initial conditions of the biofilm are tracked 
find_area()
time_keeper = [time_keeper, 0];
num_motile = [num_motile, area_s];
num_matrix = [num_matrix, area_m];
pro_m = [pro_m, area_m/(area_s + area_m)];
pro_s = [pro_s, area_s/(area_s + area_m)];

num_neighbors()
find_perimeter()
% fractal_unwrap()
ratio = perim/sqrt(area);
ratio_time = [ratio_time; ratio];

%%
for k=0:maxiter

%letting motile cells leave every iteration if they are on the edge 
for u=1:row_num
    for p=1:col_num      
        if bio_matrix0(u,p) == 1 && (bio_matrix0(u+1,p) == 0 || bio_matrix0(u-1,p) == 0 || bio_matrix0(u,p+1) == 0 || bio_matrix0(u,p-1) == 0)
           bio_matrix(u,p) = 2;
           eps(u,p) = eps_layer;
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
                               eps(i,j-1)=1;
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
                               eps(i,j+1)=1;
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
                               eps(i-1,j)=1;
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
                               eps(i+1,j)=1;
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
                              eps(i,j-1)=eps_layer;
                              grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_left<=max_eps_push  
                               for left_shift=min_col:j-1
                                   bio_matrix(i,left_shift-1) = bio_matrix0(i,left_shift);
                                   grow_time(i,left_shift-1) = grow_time(i,left_shift);
                                   eps(i,left_shift-1) = eps(i,left_shift);
                               end
                               bio_matrix(i,j-1)=2;
                               eps(i,j-1)=eps_layer;
                               grow_time(i,j-1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes right 
                       if num >= .25 && num < .50 
                           if bio_matrix0(i,j+1)==0
                              %disp("right")
                              bio_matrix(i,j+1)=2;
                              eps(i,j+1)=eps_layer;
                              grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_right<=max_eps_push  
                               for right_shift=max_col:-1:j+1
                                   bio_matrix(i,right_shift+1) = bio_matrix0(i,right_shift);
                                   grow_time(i,right_shift+1) = grow_time(i,right_shift);
                                   eps(i,right_shift+1) = eps(i,right_shift);
                               end
                               bio_matrix(i,j+1)=2;
                               eps(i,j+1)=eps_layer;
                               grow_time(i,j+1)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes down
                       if num >= .50 && num < .75 
                           if bio_matrix0(i-1,j)==0
                              %disp("down")
                              bio_matrix(i-1,j)=2;
                              eps(i-1,j)=eps_layer;
                              grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_down<=max_eps_push 
                               for down_shift=min_row:i-1
                                   bio_matrix(down_shift-1,j) = bio_matrix0(down_shift,j);
                                   grow_time(down_shift-1,j) = grow_time(down_shift,j);
                                   eps(down_shift-1,j) = eps(down_shift,j);
                               end
                               bio_matrix(i-1,j)=2;
                               eps(i-1,j)=eps_layer;
                               grow_time(i-1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                       %goes up 
                       if num > .75 
                           if bio_matrix0(i+1,j)==0
                              %disp("up")
                              bio_matrix(i+1,j)=2;
                              eps(i+1,j)=eps_layer;
                              grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           elseif eps_up<=max_eps_push  
                               for up_shift=max_row:-1:i+1
                                   bio_matrix(up_shift+1,j) = bio_matrix0(up_shift,j);
                                   grow_time(up_shift+1,j) = grow_time(up_shift,j);
                                   eps(up_shift+1,j) = eps(up_shift,j);
                               end
                               bio_matrix(i+1,j)=2;
                               eps(i+1,j)=eps_layer;
                               grow_time(i+1,j)= grow_time(i,j)+((randn/10) + 1);
                           end 
                       end 
                   end
                  %______________________________________________________________________%
                   time = grow_time(i,j);
                   grow_time(i,j)= grow_time(i,j)+((randn/10) + 1);                  
        end               
bio_matrix0 = bio_matrix;

mymap = [.1 .1 .1; 0 .75 0; 0.75 0 0.75];
if rem(k,piter)==0
%% plots biofilm every certain amount of interations (piter) 
% figure(k+1)
% pcolor(bio_matrix)
% colormap (mymap)
% title(['time= ',num2str(time)])
% %colorbar
% drawnow 
%%
%         figure(60)
%       h = pcolor(bio_matrix);
%       set(h, 'EdgeColor', 'none');
%       colormap (mymap)
%       title(['time= ',num2str(time)])
%%

% time_keeper = [time_keeper, time];
% find_area()
% num_motile = [num_motile, area_s];
% num_matrix = [num_matrix, area_m];
% pro_m = [pro_m, area_m/(area_s + area_m)];
% pro_s = [pro_s, area_s/(area_s + area_m)];

% num_neighbors()
% find_perimeter()
% % fractal_unwrap()
% ratio = perim/sqrt(area);
% ratio_time = [ratio_time; ratio];

end
%time
%     if time>=(8-0.05) && time<=(8+0.05)
%         for radval=1:40
%             innit = [];
%             motprop = 0;
%             for i=1:row_num
%                 for j=1:col_num
%                     if radius(i,j)==radval
%                         innit = [innit; bio_matrix(i,j)];
%                         if bio_matrix(i,j)==1
%                             motprop = motprop+1;
%                         end 
%                     end
%                 end 
%             end 
%   
%         propsave(radval,1) = motprop/length(innit);
%         end 
%     
% end 



end 

% proportion of motile cells at different radius values 
for radval=1:30
    innit = [];
    motprop = 0;
    matprop = 0;
    for i=1:row_num
        for j=1:col_num
            if radius(i,j)==radval
                innit = [innit; bio_matrix(i,j)];
                if bio_matrix(i,j)==1
                    motprop = motprop+1;
                end 
                if bio_matrix(i,j)==2
                    matprop = matprop+1;
                end 
            end
        end 
    end 

propsave(radval,1) = motprop/length(innit);
propsavemat(radval,1) = matprop/length(innit);
end 


proptot(:,ii) = propsave;
proptotmat(:,ii) = propsavemat;
% ratsave(ii) = ratio; 
end 
%rats(iii+1) = mean(ratsave); 




%% average intensity around a certain radii and plotting radii around the fractal pattern 
%figure (2) and (3)
% for n=0:1:30
%     for i=1:row_num
%         for j=1:col_num
%             if radius(i,j)==n
%                 power = [power; bio_matrix(i,j)];
%                 phi = [phi; theta(i,j)];
% %                 bio_matrix(i,j)=5;
%             end 
%         end 
%     end 
%     compile = cat(2,phi,power);
%     compile2 = sortrows(compile); 
%     x1=compile2(:,1);
%     y1=compile2(:,2);
%     radial_avg = [radial_avg ; sum(y1)/length(y1)];
%     radii = [radii ; n];
%     power = double.empty;
%     phi = double.empty;
% end 

%figure(2)    
%plot(radii, radial_avg)
% xlabel('Radius')
% ylabel('Average Intensity')
% grid on
% title('Average Intensity vs Radius')


%figure(3)
%plot(r_fract(:,1),r_fract(:,2))



%% cell quantity over time plot 
% figure(4) 
% plot(time_keeper, num_matrix) %,'-o'
% hold on 
% plot(time_keeper, num_motile)
% hold on 
% title('number of cell types')
% xlabel('time')
% ylabel('cell quantity')
% legend('matrix cells','motile cells')
% grid on 

%% proportion of cells over time plot 
% figure(5)
% plot(time_keeper, pro_m)
% hold on 
% plot(time_keeper, pro_s)
% hold on
% title('proportion of cell types')
% xlabel('time')
% ylabel('proportion')
% legend('proportion of matrix cells','proportion of motile cells')
% grid on 

%% perimeter/sqrt(area) plot 
% figure(6)
% plot(time_keeper, ratio_time)
% hold on
% title('time vs perimeter/sqrt(area)')
% xlabel('time')
% ylabel('ratio')
% grid on 



%% closest matrix cell to a motile cell 
% matrix_motile_relation()
% 
% for k=1:length(min_d)
%     for n=.5:360.5 %(round(sqrt(2*(double(init_row)^2)))+.5)
%         if min_d(k)>n && min_d(k)<(n+1)
%             min_d2(k)=n+.5;
%         end  
%     end 
% end 
% 
% min_d2 = transpose(min_d2);

% figure(7)
% edges = 0:10;
% histogram(min_d) 
% min_d = sortrows(min_d);
% min_d = round(min_d,4);
%%
% num_map=[];
% num_map(1)=1;
% for ik=2:length(min_d)
%     if min_d(ik)>min_d(ik-1)
%         num_map=[num_map; min_d(ik)];
%     end 
% end 
%%

% for w=1:length(num_map)
%     num_hist(w)=sum(min_d(:)==num_map(w));
% end 
% 
% % figure(8)
% % plot(num_map,num_hist)
% 
% 
% num_bar = [];
% for t=1:10
%    num_bar = [num_bar; sum(min_d2(:)==t)];
% end 
% x = 1:10;
% y1 = num_bar;
% toc;

% Save data_y2 y2
% y3 = cat(2,y1,y2);
% p0_box=0;
% for i=1:row_num
%     for j=1:col_num
%         if perimeter(i,j)==1
%             p0_box = p0_box+1;
%         end 
%     end 
% end 
% %%
% split = row_num/20;
% dimen = [1 split (2*split) (3*split) (4*split) (5*split) (6*split) (7*split) (8*split) (9*split) (10*split) (11*split) (12*split) (13*split) (14*split) (15*split) (16*split) (17*split) (18*split) (19*split) (20*split)];
% %%
% pl_box=0;
% p_box=0;


% for h=1:length(dimen)
%     for l=1:length(dimen)    
%         for i=dimen(h):
%             for j=1:split
%                 if perimeter(i,j)==1
%                 p_box = p_box+1;
%                 end 
%             end 
%         end 
%     end 
% end 
% if p_box>0
%     pl_box=plbox+1;
% end 
    

%%

% frac_grid = zeros(row_num,col_num);
% for i=1:row_num
%     for j=1:col_num
%         if i==dimen(1) || i==dimen(2)  || i==dimen(3) || i==dimen(4) || i==dimen(5) || i==dimen(6) || i==dimen(7) || i==dimen(8)  || i==dimen(9) || i==dimen(10) || i==dimen(11) || i==dimen(12)  || i==dimen(13) || i==dimen(14) || i==dimen(15) || i==dimen(16) || i==dimen(17) || i==dimen(18)  || i==dimen(19) || i==dimen(20) || i==dimen(21) 
%             frac_grid(i,j)=1;
%         elseif j==dimen(1) || j==dimen(2)  || j==dimen(3) || j==dimen(4) || j==dimen(5) || j==dimen(6) || j==dimen(7)  || j==dimen(8) || j==dimen(9) || j==dimen(10) || j==dimen(11) || j==dimen(12)  || j==dimen(13) || j==dimen(14) || j==dimen(15) || j==dimen(16) || j==dimen(17)  || j==dimen(18) || j==dimen(19) || j==dimen(20) || j==dimen(21)
%              frac_grid(i,j)=1;
%         end 
%     end 
% end 
        
%%       
% for i=split+1:(2*split)
%     for j=split+1:(2*split)
%         if perimeter(i,j)==1
%            p_box = p_box+1;
%         end 
%     end 
% end 
% 
% if p_box>0
%     pl_box=pl_box+1;
% end 

% for i=1:row_num
%     for j=1:col_num
%         if frac_grid(i,j) == 1
%             perimeter(i,j) = 3;
%         end 
%     end 
% end 

%save bfilm_data3 bio_matrix rad perimeter

%save pushtrack20 rats 
% 
% pushtracks = [];
% pushtracks = [pushtracks, rats];
%epsmount = [1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3];
for u=1:30
    maptot(u,1)=mean(proptotmat(u,:));
    ptot(u,1)=mean(proptot(u,:));
end 

ptotmat = ptot./mat0;
ptotmot = ptot./mot0;
maptotmat = maptot./mat0;
maptotmot = maptot./mot0;
plot(ptot,'o')
xlabel('radius')
ylabel('proportion of motile cells')

radtracksave = [];
for radval=1:30
    for i=1:row_num
        for j=1:col_num
            if radius(i,j)==radval
                if bio_matrix(i,j)~=0
                    radtracksave = [radtracksave; radval];
                end 
            end
        end 
    end 
end 

ptotalbul(:,(R1-3))=ptot;
ptotalmat(:,(R1-3))=ptotmat;
ptotalmot(:,(R1-3))=ptotmot;
maptotalmat(:,(R1-3))=maptotmat;
maptotalmot(:,(R1-3))=maptotmot;
end 
figure(4)
for pl=1:5
    %plot(ptotalbul(:,pl),'o')
    %plot(ptotalmat(:,pl),'o')
    plot(ptotalmot(:,pl),'-o','LineWidth',2,'MarkerSize', 8,'color', [0.9290    0.6940    0.1250])
    hold on 
end
% legend('R2 = 4', ' R2 = 5', ' R2 = 6', ' R2 = 7', ' R2 = 8', 'fontsize', 12)
% save psavetot ptotal 
% 
% for pl=1:3
%     plot(ptotal(:,pl),'o')
%     hold on 
% end 

%save R1mult_conc_push14 ptotalbul
%

%save RAND200 ptotalmot ptotalmat 

ptotsdude = ptotalmot;