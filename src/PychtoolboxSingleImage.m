% Clear the workspace
close all;
clear all;
sca;

%% Data Structure initialization variables
%dataSample = struct('name',{},'rawdata',{},'label',{},'StartTime',{},'StimulusTime',{},'StopTime',{},'Trials',{});
global  channelNumbers collectionInterval countNew isCollecting indexFlip collectFlag inNam dataSample imgArray Hd stimTime startTime stopTime rawWindow trials labeled  mode  repeat_trials numDevices samplingRate xx yy T IC;
dataSample = struct('name',{},'rawdata',[],'label',{},'StartTime',{},'StimulusTime',{},'StopTime',{},'Trials',{});
trials=0;
countNew =0;
collectFlag=0;
indexFlip=1;
inNam = [];
imgArray={};
stimTime ={};
startTime ={};
stopTime={};
labeled ={};
%% Data Collection Parameters
% Mode of Operation mode =1 -> Testing, mode =0 -> Training
mode =0;
% Flag controlling the DataCollection in the Buffer isCollecting=1 ->
% Buffer Data else =0 stop buffering
isCollecting =1;
% Variable regarding Data Collection System
repeat_trials=15; % trials of visual stimulus to be repeated
numDevices=1; % num of Bio-radio s 150 Connected, and collecting EEG
samplingRate =256;
collectionInterval =samplingRate * 21; % Raw Data buffer size collectionInterval = sampleingRate * timeTo buffer in seconds
channelNumbers = [1 2 3 4 5]; % Number of Channels Used
rawWindow = zeros( length(channelNumbers),collectionInterval); % Raw Data Buffer Allocation

%% Setting up a Timer for a background process for Collecting Bio-Radio Data (Control)
t = timer;
t.StartFcn = @bioRadio_Start_fcn;
t.Period = 0.5;%What time period?
t.StartDelay=2;
t.TimerFcn=@bioRadio_Execute_fcn4;
t.TasksToExecute =inf;
t.ExecutionMode='fixedSpacing';
t.StopFcn =@bioRadio_Stop_fcn;

%% Reading the Images

[fileImages, pathImages, wtv] = uigetfile('*.jpg','*.png','Load Images','MultiSelect', 'on');

if ischar(fileImages)
    nImages = 1;
else
    [wtv nImages] = size(fileImages);
end

for i=1:nImages
    
    if ischar(fileImages)
        img_name = strcat(pathImages,fileImages);
    else
        img_name = strcat(pathImages,fileImages{i});
    end
    
    im = imread(img_name);
    imGray = mat2gray(im);
    
%     figure, imshow(imGray), title('Original Image');
end
%% Filter Design
freq_range = [1 30];
order =6;
filter_type = 'butter';
d = fdesign.bandpass('N,F3dB1,F3dB2',order,freq_range(1),freq_range(2),samplingRate);
Hd = design(d,filter_type);
%% Image Display Code
% Here we call some default settings for setting up Psychtoolbox
KbName('UnifyKeyNames');
deviceIndex =[];

PsychDefaultSetup(2);

% Get the screen numbers
try
    
 escapeKey = KbName('ESCAPE');
    spaceKey = KbName('space');
    keysOfInterest=zeros(1,256);
    keysOfInterest(escapeKey)=1;
    keysOfInterest(spaceKey)=1;
    
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);
W=windowRect(RectRight); % screen width
H=windowRect(RectBottom); % screen height
% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Here we load in an image from file. This one is a image of rabbits that
% is included with PTB
% theImageLocation = [PsychtoolboxRoot 'PsychDemos' filesep...
%     'AlphaImageDemo' filesep 'konijntjes1024x768.jpg'];
% theImage = imread('Hopetoun_falls.jpg');
% theimage = imresize(theImage,[512 512]);
% % Get the size of the image
% [s1, s2, s3] = size(theimage);

% Here we check if the image is too big to fit on the screen and abort if
% it is. See ImageRescaleDemo to see how to rescale an image.
% if s1 > screenYpixels || s2 > screenYpixels
%     disp('ERROR! Image is too big to fit on the screen');
%     sca;
%     return;
% end
% 
% % Make the image into a texture
% imageTexture = Screen('MakeTexture', window, theImage);
% 
% % Draw the image to the screen, unless otherwise specified PTB will draw
% % the texture full size in the center of the screen. We first draw the
% % image in its correct orientation.
% Screen('DrawTexture', window, imageTexture, [], [], 0);
% 
% % Flip to the screen
% Screen('Flip', window);
% 
% % Wait for two seconds
% WaitSecs(2);
% 
% % Now fill the screen green
% Screen('FillRect', window, [0 1 0]);
% 
% % Flip to the screen
% Screen('Flip', window);
% 
% % Wait for two seconds
% WaitSecs(2);
% 
% % Draw the image to the screen for a second time this time upside down and
% % drawn onto our updated blue background
% Screen('DrawTexture', window, imageTexture, [], [], 180);
% 
% % Flip to the screen
% Screen('Flip', window);
% 
% % Wait for one second
% WaitSecs(2);

% 
% image1 = imread('7001.jpg');
% image2 = imread('7002.jpg');
% image3 = imread('7003.jpg');
% image4 = imread('7004.jpg');
% image5 = imread('7006.jpg');
% image6 = imread('7009.jpg');
% image7 = imread('7010.jpg');
% image8 = imread('7017.jpg');
% image9 = imread('3010.jpg');
% image10 = imread('3015.jpg');
% image11 = imread('3019.jpg');
% image12 = imread('3051.jpg');
% image13 = imread('3053.jpg');
% image14 = imread('3053.jpg');
% image15 = imread('3053.jpg');
% image16 = imread('3071.jpg');
% imgcell={image1,image2,image3,image4,image5,image6,image7,image8, image9,image10,image11,image12,image13,image14,image15,image16};

baseRect = [0 0 400 400];
crossSize = 20;
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);

rectColor = [1 1 1];
order=2;

 collectFlag=1;
        trials=trials+1;
        startTime= [startTime;clock()]; % Start Event
        tic;
%  barLength = 16; % in pixels
%     barWidth = 2; % in pixels
%     barColor = 0.5; % number from 0 (black) to 1 (white) 
    
    
    
  
% Now fill the screen green
Screen('FillRect', window, [0 0 0]);

% Flip to the screen
Screen('Flip', window);

% Wait for two seconds
WaitSecs(3);
start(t); %Start timer here
index = 0;
for i = 1:1
   
    
   Screen('FillRect', window, rectColor, centeredRect);
   
%    Screen('SetMouseHelper', window, x, y , 2);
%     Screen('Flip', window);
   
    fixationDuration = 0.5; % Length of fixation in seconds
      %Screen = drawCross(window,W,H);
    barLength = 32; % in pixels
    barWidth = 4; % in pixels
    barColor = 0; % number from 0 (black) to 1 (white) 
    Screen('FillRect', window, barColor,[ (W-barLength)/2 (H-barWidth)/2 (W+barLength)/2 (H+barWidth)/2]);
    Screen('FillRect', window, barColor ,[ (W-barWidth)/2 (H-barLength)/2 (W+barWidth)/2 (H+barLength)/2]);
 Screen('Flip', window)
 
%       tFixation = Screen('Flip', window);
% function drawCross(window,W,H)

% Show fixation cross

     WaitSecs(5);
     c=clock;
%      [mstart sstart]= c(:,5:end);
    for j=1:nImages
         
        % index of image
        index = index + 1;
        im=strcat(pathImages,fileImages{j});
        img=imread(im); 
        imageTexture = Screen('MakeTexture', window, img);

        Screen('DrawTexture', window, imageTexture, [], [], 0);

        % Flip to the screen
        Screen('Flip', window);

        % Wait for two seconds
        WaitSecs(3);
        
        if j ==8
%     [mend send]= c(:,5:end);        
             % Now fill the screen green
%              extractDataSample% Calling the epoching function
Screen('FillRect', window, [1 1 1]);

% Flip to the screen
Screen('Flip', window);

% Wait for two seconds
WaitSecs(10);
    end
    end
%     timeStimulus = [timeStimulus;toc]; % Stimulus Event
%     stimTime = [stimTime;timeStimulus];
% stopTime= [stopTime;clock()];
%  extractMyDataSample(order);

end 

stop(t);
% % Now fill the screen green
% Screen('FillRect', window, [1 1 0]);
% 
% % Flip to the screen
% Screen('Flip', window);
% 
% % Wait for two seconds
% WaitSecs(2);
%   KbQueueCreate(deviceIndex, keysOfInterest);
%     KbQueueStart(deviceIndex);
%         KbQueueRelease(deviceIndex);
Screen('CloseAll');
% Clear the screen
% sca;
catch
    stop(t);
    KbQueueRelease(deviceIndex);
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(pyschlasterror);
end


 