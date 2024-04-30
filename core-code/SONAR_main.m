clearvars -except thepath h;
tic

% Construct full paths for the input files
y_file = fullfile(thepath, 'y.txt');
u_file = fullfile(thepath, 'u.txt');
N_file = fullfile(thepath, 'N.txt');
label_file = fullfile(thepath, 'label.txt');
coord_file = fullfile(thepath, 'coord.txt');

% Import data using the constructed paths
y = importdata(y_file, ",");
y = y.data; y = y';
u = importdata(u_file, ",");
u = u.data; u = u'; temp = repelem(1, size(u, 2)); u = [temp; u];
N = importdata(N_file, ",");
N = N.data;
label = importdata(label_file, ",");
label = label.data;
coord = importdata(coord_file, ",");
coord = coord.data;

D=pdist(coord);
D=squareform(D);
[r,c]=size(D);
w=zeros(length(N),length(N));
%radius
%h=1.2;
for w1=1:r
    for w2=1:c
        if (D(w1,w2)<h)&&(label(w1)==label(w2))
            w(w1,w2)=(1-(D(w1,w2)/h)^2)^2;
        else
            w(w1,w2)=0;
        end
    end
end
norm_y=y;
for y_i=1:size(y,1)
    norm_y(y_i,:)=norm_y(y_i,:)/sum(norm_y(y_i,:));
end
cos_dist=squareform(pdist(norm_y,'cosine'));
%auto
w=w.^(100*cos_dist);
G=size(u,2);
x0=rand(1,size(u,1));
x0=[x0,10];
lb=zeros(size(u,1)+1,1)';
ub=[];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon =[];
options = optimoptions('fmincon','UseParallel',true,'SpecifyObjectiveGradient',true);
JIE=zeros(length(N),(size(u,1)+1));
parfor i = 1:length(N)
    JIE(i,:)=fmincon(@Grad_of_my_GWNBR_likeli,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,u,N,w,G,y,i);
end
JIE(:,1)=[];
JIE(:,end)=[];
for guiyi_i = 1:size(JIE,1)
    JIE(guiyi_i,:)=JIE(guiyi_i,:)/sum(JIE(guiyi_i,:));
end

% Save results using fullfile
results_file = fullfile(thepath, 'SONAR_results.mat');
save(results_file, 'JIE');

toc
fprintf("Complete\n")
