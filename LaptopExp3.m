clear all
close all
clc
tic
%% Introducing Elementary Variables
cam             = webcam; %Introduce the webcam. Since the laptop has only one cam, it will select the default. 
cam.Resolution  = '320x240'; %Choose a resolution that supports. We select the lowest since we are interested in SPC.

n               = 8; %The square length of the desired object scene which is nxn. 
m               = 50; %The measurement amount
N               = n^2; %The total amount of data. 

%% Fixing a Reference Frame. (Not practical, but for learning purposes, we need to see the original image in the end.)
orgImg          = snapshot(cam); % Take a picture with the laptop cam.
imshow(orgImg)
orgImg          = rgb2gray(orgImg); %Transform it into grayscale image. 
orgImg          = imresize(orgImg, [n n]); %Resize it nxn so that we can compare our nxn result. 
orgImgDouble    = double(orgImg); %Convert to double.
% imshow(orgImg)

%% Initiate Zero Matrices to Enhance Speed. (MATLAB works slower if a variable changes dimensions every iteration.)
y               = zeros(m,1) ;% Measurement vector of mx1. ( y = Phi * x)
Theta           = zeros(m,N); %Identifying matrix (y = Phi * x >> y = Phi * Psi * s>> y = Theta * s >>>> Theta = Phi * Psi)
x1              = zeros(N,1); % The object vector for reconstruction with l1-magic.
x2              = zeros(N,1); % The object vector for reconstruction with least squares.
PatternCellArray = cell(1,m); %Array for keeping the illumination patterns.
Phi             = zeros(m,N); %The collective pattern matrix. (y = Phi * x)
CellArray       = cell(1,m);
%% Constructing the Random Pattern Array
rng(1)



 for i=1:m % Select the corresponding rows of the Hadamard matrix, and put them in the array by resizing to nxn
    A           = (sign(randn(n,n))+ones(n,n))/2

    % A = load("Atest2_32x32Binary.mat");

    PatternCellArray{i}=A;
    CellArray {i} = A(:);
 end
 PatternCellArray(cellfun('isempty', PatternCellArray)) = []; %Remove the empty cells from the array.
 
 for i=1:m %Insert the random patterns into Phi
dummyVariable   = PatternCellArray{i}; 
dummyVariable   = dummyVariable(:).';
Phi(i,:)        = dummyVariable;
 end

 %% Display the Patterns and Take Measurements
figure('Toolbar', 'none', 'Menu', 'none', 'Position',[1920 0 1920 1080], 'WindowState', 'fullscreen');  %Fullscreen with no bars, menu etc.
set(gca, 'Unit', 'normalized', 'Position', [0 0 1 1]);
set(gcf, 'color', 'k'); %Background black
imshow(zeros(n)) % leave an extra layer of black at the background,

for i=1:m
    figure('Toolbar', 'none', 'Menu', 'none', 'Position',[1920 0 1920 1080], 'WindowState', 'fullscreen');  %Fullscreen with no bars, menu etc.
    set(gcf, 'color', 'k'); %Background black

    set(gca, 'Unit', 'normalized', 'Position', [0 0 1 1]);
    imshow(PatternCellArray{i}); %Show the i th pattern in array. 
    pause(2)
    img         = snapshot(cam); % Take a screenshot
    pause(8  ) % Pause for a sec to observe
    img         = rgb2gray(img); % Convert to grayscale
    img         = imresize(img, [n n]); %Resize
    imgd        = double(img); % Convert to double 
    imgd        = imgd(:); % Turn it to a double vector
    y(i)        = mean(imgd,'all');
    %y(i) = imresize(img, [1,1]);
    close

end
%% Extension Code for Fixed Experiment
%{
y = [
    


    ]'

%}
%% Kutz Style
%{
A_mat       = zeros(ny,nx);

for j=1:m
    A_mat(r1k(j))=1;
    Adel        = reshape(idct2(A_mat),nx*ny,1);
    Adelta(j,:) = Adel;
    A_mat(r1k(j)) = 0;
end
%}
AA = zeros(n);

%for j=1:m
   % AA


%% Construct Theta matrix ( y = Theta * s >>> Theta = Phi * Psi)
% Essentially, construct 1D DCT Basis vectors, 
 for i=1:N
     i;
     B          = zeros(1,N); % This is practically for the vectors from an eye matrix, or Kronecker matrix
     B(i)       = 1;  % For the elementary basis vectors . 
     Psi        = idct(B)'; % Now, D <-> K are complementary. Meaning, it takes infinite basis K vectors to form a sinusoid, and vice versa
     Theta(:,i) = Phi*Psi; % Form the reconstruction matrix Theta sequentially, from Basis vector Psi and the sensing matrix Phi
     toc
 end
 %% l2-norm Solution
 s2             = pinv(Theta)*y; % least squares
 %% l1-norm Solution
 % l1-magic
 s1             = l1eq_pd(s2*0,Theta,Theta',y,5e-3,20); % first input is the starting point... 
 % ...either pick up the vector from LS, or start from zero 

 %x = l1eq_pd(y,A,A',b,5e-3,32);
 %% Image Reconstructions for L1
for j = 1:N
    j;  % counter
    B           = zeros(1,N); % Same as above, however this time construct x1 sequentially
    B(j)        = 1;
    Psi         = idct(B)';
    x1          = x1+Psi*s1(j);
    toc
end
  %% Image Reconstructions for L2 
 for j = 1:N
    j;
    B           = zeros(1,N);
    B(j)        = 1;
    Psi         = idct(B)';
    x2          = x2+Psi*s2(j);
    toc
 end
%%
%{
A_mat = zeros(n,n);

for j=1:N
    A_mat(j)=1;
    Psi9 = reshape(idct2(A_mat),N,1);
    Adelta(j,:) = Psi9;
    A_mat(j) = 0;
    
end
%}


%% CVX L1 Reconstruction
cvx_begin;
    variable x3(size(Theta,2));
    minimize ( norm(x3,1) );
%    minimize ( norm( Theta*x3-y,1 ) ); 

    subject to 
        Theta*x3 == y; % constraint function
cvx_end
%% CVX L2 Reconstruction

cvx_begin;
    variable x4(size(Theta,2));
    minimize ( norm(Theta*x4-y,2) );
%   minimize ( norm(x4,2) );
cvx_end


 %% Display the Reconstructed Images
 figure('name','Compressive sensing image reconstructions')
subplot(2,3,1), imagesc(reshape(img,n,n)), xlabel('original'), axis image
subplot(2,3,2), imagesc(reshape(x1,n,n)), xlabel('l1-magic'), axis image
subplot(2,3,3), imagesc(reshape(x2,n,n)), xlabel('l2-norm'), axis image
subplot(2,3,4), imagesc(reshape(x3,n,n)), xlabel('l1-CVX'), axis image
subplot(2,3,5), imagesc(reshape(x4,n,n)), xlabel('l2-CVX'), axis image

colormap gray
%%
toc
     
     
 
 
 
