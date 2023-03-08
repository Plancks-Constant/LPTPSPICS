clear all
close all
clc
tic
%% Introducing Elementary Variables
cam             = webcam; %Introduce the webcam. Since the laptop has only one cam, it will select the default. 
cam.Resolution  = '320x240'; %Choose a resolution that supports. We select the lowest since we are interested in SPC.

n               = 10; %The square length of the desired object scene which is nxn. 
m               = 95; %The measurement amount
N               = n^2;
p               = 0.5; %The total amount of data. 

%% Fixing a Reference Frame. (Not practical, but for learning purposes, we need to see the original image in the end.)
orgImg          = snapshot(cam); % Take a picture with the laptop cam.
%imshow(orgImg)
orgImg          = rgb2gray(orgImg); %Transform it into grayscale image. 
orgImg          = imresize(orgImg, [n n]); %Resize it nxn so that we can compare our nxn result. 
orgImgDouble    = double(orgImg); %Convert to double.
%imshow(orgImg)

%% Initiate Zero Matrices to Enhance Speed. (MATLAB works slower if a variable changes dimensions every iteration.)
y               = zeros(m,1) ;% Measurement vector of mx1. ( y = Phi * x)
Theta           = zeros(m,N); %Identifying matrix (y = Phi * x >> y = Phi * Psi * s>> y = Theta * s >>>> Theta = Phi * Psi)
x1              = zeros(N,1); % The object vector for reconstruction with l1-magic.
x2              = zeros(N,1); % The object vector for reconstruction with least squares.
PatternCellArray = cell(1,m); %Array for keeping the illumination patterns.
Phi             = zeros(m,N); %The collective pattern matrix. (y = Phi * x)
CellArray       = cell(1,m);
ObservePhotosCellArray = cell(1,m);
%%

A               = rand(n);
A               = (A<p)


%% Constructing the Hadamard Pattern matrix 

H               = hadamard(N); %Construct NxN Hadamard matrix. 
H(find(H==-1))  = 0; %Exchange -1 to 0.

%H(1,:) = zeros(1,N)

%% Or constructing random pattern matrix
%{
 rng(2)
for i=1:m % Random matrix
     

    A = (sign(randn(n,n))+ones(n,n))/2
    PatternCellArray{i}=A;
CellArray {i} = A(:);
PatternCellArray(cellfun('isempty', PatternCellArray)) = []; %Remove the empty cells from the array.

 end
%}
%%


u = 0; % how many times 1 come up counter


randRows        = zeros(1,m); % 
randRows        = randperm(N,m); % Select m hadamard patterns out of N, index vector

while ismember(1, randRows) == true % Do not include the first pattern, which is the full bright pattern.
    u = u+1;
    randRows        = randperm(N,m); % Choose m unique numbers from N.
    
end
%{
for i = randRows % Select the corresponding rows of the Hadamard matrix, and put them in the array by resizing to nxn
 A              = H(i,:); 

 PatternCellArray{i}=reshape(A, n, n);
end

%}
% PatternCellArray(cellfun('isempty', PatternCellArray)) = []; %Remove the empty cells from the array.
 
 %% Extension Code for fixed pattern
%PatternCellArray = {[1,1,1,0,1,0,0,0;1,1,1,0,1,0,1,0;0,0,1,0,0,1,0,0;1,0,1,1,1,0,0,0;0,0,1,1,0,0,1,0;0,1,1,1,0,0,0,1;0,1,0,1,0,1,0,0;0,0,1,0,0,1,0,1],};
 for i=1:m %Insert the random patterns into Phi
dummyVariable   = PatternCellArray{i}; 
dummyVariable   = dummyVariable(:).';
Phi(i,:)        = dummyVariable;
 end
%% Extension COde for Black background
figure('Toolbar', 'none', 'Menu', 'none', 'Position',[1920 0 1920 1080], 'WindowState', 'fullscreen');  %Fullscreen with no bars, menu etc.
set(gca, 'Unit', 'normalized', 'Position', [0 0 1 1]);
set(gcf, 'color', 'k'); %Background black
imshow(zeros(n)); %Show the first pattern in array to give a black background, while switching patterns


%{
%% Extension Code for pure white screen
 figure('Toolbar', 'none', 'Menu', 'none', 'Position',[1920 0 1920 1080], 'WindowState', 'fullscreen');  %Fullscreen with no bars, menu etc.
set(gca, 'Unit', 'normalized', 'Position', [0 0 1 1]);
set(gcf, 'color', 'k'); %Background black
imshow(ones(n)); %Show the i th pattern in array. 
pause(0.5)
 close
%}

 %% Display the Patterns and Take Measurements
for i=1:m
figure('Toolbar', 'none', 'Menu', 'none', 'Position',[1920 0 1920 1080], 'WindowState', 'fullscreen');  %Fullscreen with no bars, menu etc.
set(gca, 'Unit', 'normalized', 'Position', [0 0 1 1]);
set(gcf, 'color', 'k'); %Background black

%displayPattern = (PatternCellArray{i}); % direct illumination projection

displayPattern = flip(PatternCellArray{i},2); % Laptop screen projection does not correspond to the same pixel as the matrix projection

imshow(displayPattern); %Show the i th pattern in array.4

pause(1) % Pause for a sec to observe
img = snapshot(cam); % Take a screenshot


pause(0.2) % Pause for a sec to observe
img = rgb2gray(img); % Convert to grayscale
ObservePhotosCellArray{i} = img; % Just for checking the images later, not really necessary. 
img = imresize(img, [n n]); %Resize

%imgd = double(img); % Convert to double 
img = img(:); % Turn it to a double vector
%imgd = imgd(:); % Turn it to a double vector

y(i) = sum((img),'all')/(n^2); % Make it single pixel 

%y(i) = sum((imgd),'all')/(n^2);
%y(i) = imresize(img, [1,1]);
close

end
%%Extension code for manual y input


%% Construct Theta matrix ( y = Theta * s >>> Theta = Phi * Psi)
 for i=1:N
     i; % counter to keep an eye at the command window
     B = zeros(1,N);
     B(i) = 1; % Kronecker matrix, basis vectors sequentially 
     Psi = idct(B)'; % 1D DCT basis vectors
     Theta(:,i) = Phi*Psi; % Sequentially build up the reconstruction matrix
 end



 %% l2-norm Solution
     s2 = pinv(Theta)*y; % LS 
 %% l1-norm Solution (Basis Pursuit)
 s1 = l1eq_pd(s2,Theta,Theta',y,1e-3,50);
 %x = l1eq_pd(y,A,A',b,5e-3,32);
  %% Image Reconstructions for L2 
 for j = 1:N
    j;
    B = zeros(1,N);
    B(j) = 1;
    Psi = idct(B)';
    x2 = x2+Psi*s2(j); % Sequentially build up signal
 end
%% Image Reconstructions for L1
for j = 1:N
    j;
    B = zeros(1,N);
    B(j) = 1;
    Psi = idct(B)';
    x1 = x1+Psi*s1(j);
end
 %% Display the Reconstructed Images
 figure('name','Compressive sensing image reconstructions')
subplot(1,3,1), imagesc(reshape(img,n,n)), xlabel('original'), axis image
subplot(1,3,2), imagesc(reshape(x1,n,n)), xlabel('l1-norm'), axis image
subplot(1,3,3), imagesc(reshape(x2,n,n)), xlabel('l2-norm'), axis image
colormap gray
%%
toc
     
     
 
 
 
