%% AUTHORS : SHIVA KARTHIK REDDY KOYYA , RAMESH NAIR
%% The project is intended for study purpose and is trying to 
%% understand the Human Brain in terms of the stress induced when 
%% view certain images.


clear
clc
%% imorting data
load('nsclassdata.mat')
load('sClaasData.mat')
% importing data for robustness test
 load('matlab.mat')
 r_ns=rob_ns;
 r_s=rob_s;
SData=stress;
NSData=class1Data;
figure(1)
plot(SData(1:256*5,1),'r');
hold on
plot(NSData(1:256*5,1),'b');
grid on
title('plot of rawdata for both classes channel 1 ');
legend('stress','nonstress');
hold off
%% call to function to perform filtering, artifact removal,feature genration and feature reduction
[W,Xm,Label]=allMajor(SData,NSData);
%% robustness check data
[W_r,Xm_r,Label_r]=allMajor(r_s,r_ns);

%% data for classifier
X=Xm*W;
[IDX,Z]=rankfeatures(X',Label','Criterion','bhattacharyya');
X_S=X(1:56,IDX(1:10,1));
X_NS=X(57:end,IDX(1:10,1));
L_S=Label(1:56,:);
L_NS=Label(57:end,:);
% feature matrix for robost data
X_r=Xm_r*W_r;
[m,n]=size(X_r);
for i=1:n
    X_r(:,i)=(X_r(:,i)-min(X_r(:,i)))/(max(X_r(:,i))-min(X_r(:,i)));
end  
Y_r=Label_r;

[IDX_r,Z_r]=rankfeatures(X_r',Label_r','Criterion','bhattacharyya');

X_rs=X_r(1:56,IDX_r(1:10,1));
X_rns=X_r(57:end,IDX_r(1:10,1));
L_rs=Y_r(1:56,:);
L_rns=Y_r(57:end,:);






   %% Diving the data into TRAINING , CROSSVALIDATION , PERFORMANCE .
P=randperm(56);
Xt=[X_S(P(1:28),:);X_NS(P(1:28),:)];
Yt=[L_S(P(1:28),:);L_NS(P(1:28),:)];

Xv=[X_S(P(29:28+14),:);X_NS(P(29:28+14),:)];
Yv=[L_S(P(29:28+14),:);L_NS(P(29:28+14),:)];

Xp=[X_S(P(43:end),:);X_NS(P(43:end),:)];
Yp=[L_S(P(43:end),:);L_NS(P(43:end),:)];

Xr=[X_rs(P(43:end),:);X_rns(P(43:end),:)];
Yr=[L_rs(P(43:end),:);L_rns(P(43:end),:)];

%% implementation of SVM with Gaussion RBF

fprintf('\nstart of implementation of svm \n')
model = svmtrain(Xt, Yt, 'kernel_function', 'rbf', 'rbf_sigma', 1 );
p = svmclassify(model, Xt);
fprintf('Training Accuracy: %f\n', mean(double(p == Yt)) * 100);

fprintf('performing computation for optimal value Sigma \n')
error= 999999;
s =  1:0.2:10;


for j=1:0.2:10
    sigma_c = j;
    model = svmtrain(Xt, Yt, 'kernel_function', 'rbf', 'rbf_sigma', j );
    d = svmclassify(model, Xt);
    fprintf('\n Training Accuracy: %f\n', mean(double(d == Yt)) * 100);
    p = svmclassify(model, Xv);
    fprintf('Validation set Accuracy: %f\n', mean(double(p == Yv)) * 100);
    fprintf('value of sigma %f\n',j)
    pred = svmclassify(model, Xv);
    new_error = mean(double(pred~=Yv));
    if ( new_error < error &&(mean(double(d == Yt)) * 100>=mean(double(p == Yv)) * 100))
        error = new_error;
        sigma = sigma_c;
        
    end
end



fprintf('optimal value of sigma is  sigma= %f\n',sigma)
fprintf('traning for optimal values of c and sigma \n')

model = svmtrain(Xt, Yt, 'kernel_function', 'rbf', 'rbf_sigma', sigma );
%% final set of acctracy for the system

fprintf('\ndisplaying the final set of accuracies for the model with all types of data sets\n')
%% Test system on Trainign Dataset
%Xt, Yt
p = svmclassify(model, Xt);
fprintf('Training Accuracy SVM: %f\n', mean(double(p == Yt)) * 100);
%%  Cross Validating system  using
% Xv, Yv
p = svmclassify(model, Xv);
fprintf('Validation set Accuracy SVM: %f\n', mean(double(p == Yv)) * 100);

%% Computing  performance of the system using
% Xp , Yp

p = svmclassify(model, Xp);
fprintf('performance set Accuracy SVM: %f\n', mean(double(p == Yp))* 100);
 
%% Computing  robustness of the system using
% Xr , Yr

p = svmclassify(model, Xr);
fprintf('robustness of SVM: %f\n', mean(double(p == Yr))* 100);




    %% implementation of neural network for classification 
 %% Setup the parameters you will use for this exercise
    
fprintf('\n start of implementation of neural network classifier' )   
input_layer_size  = size(Xt,2);      
hidden_layer_size = 25;              
num_labels = 2;                       
                         
%%
fprintf('Initializing Neural Network Parameters ...\n')

initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);

% Unroll parameters
initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];

%% 

fprintf('Training Neural Network... \n')

%  After you have completed the assignment, change the MaxIter to a larger
%  value to see how more training helps.
options = optimset('MaxIter', 50);

%  You should also try different values of lambda
lambda = 1;
Yt(Yt==0)=2;
% Create "short hand" for the cost function to be minimized
costFunction = @(initial_nn_params) nnCostFunction(initial_nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, Xt, Yt, lambda);

% Now, costFunction is a function that takes in only one argument (the
% neural network parameters)
[nn_params, cost] = fmincg(costFunction, initial_nn_params, options);

% Obtain Theta1 and Theta2 back from nn_params
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));


   
    %% Test system on Trainign Dataset
    pred = predict(Theta1, Theta2, Xt);
    fprintf('Training Set Accuracy NEURAL: %f\n', mean(double(pred == Yt)) * 100);

    %%  Cross VAlidati system with using
    % Xv, Yv
    Yv(Yv==0)=2;
    pred = predict(Theta1, Theta2, Xv);
    fprintf('validation Set Accuracy NEURAL: %f\n', mean(double(pred == Yv)) * 100);

    %% Compute performance of the system 
    % Xr , Yr
   Yp(Yp==0)=2;
    pred = predict(Theta1, Theta2, Xp);
    fprintf('performance Set Accuracy NEURAL: %f\n', mean(double(pred == Yp)) * 100);

 
  % compute robustness of the model using neural network
    Yr(Yr==0)=2;
    pred = predict(Theta1, Theta2, Xr);
    fprintf('Robustness of NEURAL: %f\n', mean(double(pred == Yr)) * 100);


















