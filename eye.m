function varargout = eye(varargin)
% EYE MATLAB code for eye.fig
%      EYE, by itself, creates a new EYE or raises the existing
%      singleton*.
%
%      H = EYE returns the handle to a new EYE or the handle to
%      the existing singleton*.
%
%      EYE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EYE.M with the given input arguments.
%
%      EYE('Property','Value',...) creates a new EYE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eye_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eye_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eye

% Last Modified by GUIDE v2.5 28-Sep-2019 01:15:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eye_OpeningFcn, ...
                   'gui_OutputFcn',  @eye_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before eye is made visible.
function eye_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eye (see VARARGIN)

% Choose default command line output for eye
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes eye wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = eye_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input_image h
figure
RH = imhist(input_image(:,:,1), 256);
GH = imhist(input_image(:,:,2), 256);
BH = imhist(input_image(:,:,3), 256);

h(1) = stem(1:256, RH);
h(2) = stem(1:256 + 1/3, GH);
h(3) = stem(1:256 + 2/3, BH);title('Output Image Histogram');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input_image
global Z B out h
B = imresize(input_image, [584 565]);
% Read image
im = im2double(B);
% Convert RGB to Gray via PCA
lab = rgb2lab(im);
f = 0;
wlab = reshape(bsxfun(@times,cat(3,1-f,f/2,f/2),lab),[],3);
[C,S] = pca(wlab);
S = reshape(S,size(lab));
S = S(:,:,1);
gray = (S-min(S(:)))./(max(S(:))-min(S(:)));
% Contrast Enhancment of gray image using CLAHE
J = adapthisteq(gray,'numTiles',[8 8],'nBins',128);
%Background Exclusion
% Apply Average Filter
h = fspecial('average', [9 9]);
JF = imfilter(J, h);
figure, imshow(JF)
% Take the difference between the gray image and Average Filter
Z = imsubtract(JF, J);
figure, imshow(Z)
% Threshold using the IsoData Method
level= isodata(Z); % this is our threshold level
%level = graythresh(Z)
% Convert to Binary
BW = im2bw(Z, level-.008);
% Remove small pixels
BW2 = bwareaopen(BW, 100);
% Overlay
BW2 = imcomplement(BW2);
%figure;imshow(BW2)

out = imoverlay(B, BW2, [0 0 0]);
%figure, imshow(out);
h=im2bw(out,0.1)
figure;imshow(h)

I=double(imresize(input_image,[416,486]));
mI=double(imresize(out,[416,486]));
save M.mat I mI
% b = bwboundaries(V2);
% figure;set(imshow(input_image));
% hold on
% for k = 1:numel(b)
%     plot(b{k}(:,2), b{k}(:,1), 'b', 'Linewidth', 1)
% end
% 
% global IMG
% IMG=double(V2);
% CLM = 2;
% [ IN, CEN, NOF ] = SUB_FUN8( IMG, CLM );
% lb = zeros(2,1);
% figure;
% imshow(IN(:,:,1),[]);title('segmented image 1')
% figure
% imshow(IN(:,:,2),[]);title('segmented image 2')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input_image
global I
I = rgb2gray(input_image);
axes(handles.axes2);
imshow(I,[]);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I V2
Ip = single(I);
thr = prctile(Ip(Ip(:)>0),1) * 0.9;
Ip(Ip<=thr) = thr;
Ip = Ip - min(Ip(:));
Ip = Ip ./ max(Ip(:));    

% compute enhancement for two different tau values
V1 = vesselness2D(Ip, 0.5:0.5:2.5, [1;1], 1, false);
V2 = vesselness2D(Ip, 0.5:0.5:2.5, [1;1], 0.5, false);
figure; 
imshow(V1,[])
title('Enhanced Image for tau=1')

figure; 
imshow(V2,[])
axes(handles.axes2);
imshow(V2,[])
title('Enhanced Image for tau=0.5')

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input_image
[filename, pathname]=uigetfile('file selector');
a=strcat([pathname filename]);
input_image = double(imread(a))/255;
axes(handles.axes1);
imshow(input_image,[])


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global V2 h
% % dim = ndims(V2);
% % if(dim == 3)
% %     %Input is a color image
% %     input_image = rgb2gray(V2);
% % end
% 
% %Extract Blood Vessels
% Threshold = 1;
% bloodVessels = VesselExtract(V2, Threshold);
% 
% %Output Blood Vessels image
% 
% figure;imshow(V2);title('Input Image');
% figure;imshow(bloodVessels,[]);title('Extracted Blood Vessels');
binary=im2bw(V2,0.7)
%figure;imshow(binary);title('binary')
CC = bwconncomp(binary, 8);
S = regionprops(CC, 'Area');
L = labelmatrix(CC);
BW2 = ismember(L, find([S.Area] >= 50));
%figure;imshow(BW2);title('Region with boundary')
labelled = bwlabel(BW2); % binary image is the one on the right, which I'm assuming is binary
 object_number = 1; % 1 is the object which is closest to the left hand side edge, and if there are two, then it choses the one nearest the top of the image. 
 labelled(labelled == object_number) = 0;
 binary_image_corrected = labelled > 0;
 figure;imshow(binary_image_corrected);title('Vessel extracted')
  BW2=imresize(binary_image_corrected,[640,480]);
   h=imresize(h,[640,480]);
 r=h-BW2;
  BW = im2bw(r,0.05);
 BW= bwareaopen(BW,10);
 figure;imshow(BW,[]);
 save N.mat binary_image_corrected
