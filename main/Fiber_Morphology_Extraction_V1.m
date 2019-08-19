function varargout = Fiber_Morphology_Extraction_V1(varargin)
% FIBER_MORPHOLOGY_EXTRACTION_V1 MATLAB code for Fiber_Morphology_Extraction_V1.fig
%      FIBER_MORPHOLOGY_EXTRACTION_V1, by itself, creates a new FIBER_MORPHOLOGY_EXTRACTION_V1 or raises the existing
%      singleton*.
%
%      H = FIBER_MORPHOLOGY_EXTRACTION_V1 returns the handle to a new FIBER_MORPHOLOGY_EXTRACTION_V1 or the handle to
%      the existing singleton*.
%
%      FIBER_MORPHOLOGY_EXTRACTION_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIBER_MORPHOLOGY_EXTRACTION_V1.M with the given input arguments.
%
%      FIBER_MORPHOLOGY_EXTRACTION_V1('Property','Value',...) creates a new FIBER_MORPHOLOGY_EXTRACTION_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fiber_Morphology_Extraction_V1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fiber_Morphology_Extraction_V1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fiber_Morphology_Extraction_V1

% Last Modified by GUIDE v2.5 07-Aug-2019 16:58:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fiber_Morphology_Extraction_V1_OpeningFcn, ...
                   'gui_OutputFcn',  @Fiber_Morphology_Extraction_V1_OutputFcn, ...
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


% --- Executes just before Fiber_Morphology_Extraction_V1 is made visible.
function Fiber_Morphology_Extraction_V1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fiber_Morphology_Extraction_V1 (see VARARGIN)

% Choose default command line output for Fiber_Morphology_Extraction_V1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fiber_Morphology_Extraction_V1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fiber_Morphology_Extraction_V1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

    %set(handles.Sigma,'String','1');
    set(handles.rhoRes,'String','1');
    set(handles.Peaks,'String','100');
    set(handles.Theta,'String','1');
    set(handles.Gap,'String','40');
    set(handles.Length,'String','40');
    set(handles.Scale,'String','50');
    set(handles.Thresh,'String','0.5');
    set(handles.SE,'String','80');
    set(handles.DT,'String','2');
    set(handles.CT,'String','50');
    set(handles.TT,'String','4');


% --- Executes on button press in SHT.
function SHT_Callback(hObject, eventdata, handles)
% hObject    handle to SHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SHT


% --- Executes on button press in OHT.
function OHT_Callback(hObject, eventdata, handles)
% hObject    handle to OHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OHT


% --- Executes on button press in PHT.
function PHT_Callback(hObject, eventdata, handles)
% hObject    handle to PHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PHT

% --- Executes on button press in Gradient.
function Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gradient

% --- Executes on button press in BRM.
function BRM_Callback(hObject, eventdata, handles)
% hObject    handle to BRM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BRM

% --- Executes on button press in detectFiber.
function detectFiber_Callback(hObject, eventdata, handles)
    handles.output = hObject;
    axes(handles.axes2); % Chage -> axes4 to axes2
    hold off;
    cla reset;
    
    % Image preprocessing step
    global I;
    S = 2;
    global gaussian_im;
    gaussian_im = imgaussfilt(I,S);
    global gray_im;
    gray_im = mat2gray(gaussian_im,[0 255]);
    level = graythresh(gray_im);
    global bw_im;
    bw_im = imbinarize(gray_im,level);
    global centroildAll;
    global observedOri;
    %===============================%
    
    method_selection = get(handles.method,'selectedObject');
    Method = get(method_selection,'String');
    %global I;
    %global bw_im;
    %global gray_im;
    axes(handles.axes2); % change -> axes4 to axes2
    imshow(gray_im);
    hold on;
    rho = get(handles.rhoResolution,'value');
    theta = get(handles.thetaResolution,'value');
    thresh = get(handles.threshCof,'value');
    p = get(handles.noPeak,'value');
    scale = get(handles.suppressScale,'value');
    gap = get(handles.fillGap,'value');
    mlength = get(handles.minLength,'value');

    switch Method
        case 'Simple HT'
            nfiber = 0;
            SK = bwmorph(bw_im,'skel',inf);
            [H,Theta,Rho] = hough(SK,'RhoResolution',rho,'Theta',-90:theta:89);
            hood = floor(size(H)/scale);
            sizeN = hood - (1 - mod(hood,2));
            peaks = houghpeaks(H,p,'Threshold', thresh*max(H(:)),'NHoodSize',sizeN);
            lines = houghlines(SK,Theta,Rho,peaks,'FillGap',gap,'MinLength',mlength);
            centroildAll = [];
            observedOri = [];
            for k = 1:length(lines)
                xy = [lines(k).point1; lines(k).point2];
                plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

                % Plot beginnings and ends of lines
                plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
                plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
                centroid = (xy(1,:)+xy(2,:))/2;
                angl = -1*atand((xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2)));
                centroildAll = [centroildAll;centroid];
                observedOri = [observedOri; angl];
            end
            nfiber = nfiber + length(lines);
            val = num2str(nfiber);
            set(handles.fiberNum,'string',val);

        case 'Opening HT'
            SE_L = get(handles.lengthSE,'value');
            DT = get(handles.deltaTheta,'value');
            CT = get(handles.centerTol,'value');
            TT = get(handles.thetaTol,'value');
            nfiber = 0;
            [all_boundary,center,false_boundary,boundary,fiber_ori,fib_orientation] = opening_method(bw_im,SE_L,DT,CT,TT);
            nlines = apply_hough(gray_im,boundary,rho,theta,thresh,p,scale,gap,mlength);
            centroildAll = [];
            for k = 1:length(nlines)
                %if ne(lines(k).theta,angle)
                xy = [nlines(k,1), nlines(k,2); nlines(k,3), nlines(k,4)];
                plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
                % Plot beginnings and ends of lines
                plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
                plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
                centroid = (xy(1,:)+xy(2,:))/2;
                angl = -1*atand((xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2)));
                centroildAll = [centroildAll;centroid];
                observedOri = [observedOri; angl];
            end
            nfiber = nfiber + length(nlines);
            val = num2str(nfiber);
            set(handles.fiberNum,'string',val);

        case 'Partitioned HT'
            SK = bwmorph(bw_im,'skel',inf);
            Fibers = bwconncomp(SK);
            fib = cell(Fibers.NumObjects,1);
            nfiber = 0;
            N_line = 0;
            TF = 0;
            Line_index = [];
            centroildAll = [];
            observedOri = [];
            for j = 1:Fibers.NumObjects
                temp = zeros(length(Fibers.PixelIdxList{j}),2);
                ImageTemp = zeros(size(gray_im));
                for k = 1:length(Fibers.PixelIdxList{j})
                    [temp(k,1),temp(k,2)] = ind2sub(size(gray_im), Fibers.PixelIdxList{j}(k));
                    ImageTemp(temp(k,1),temp(k,2)) = 1;
                end
                fib{j}=temp;
                [nlines,lines_index]=findFiber(ImageTemp,rho,theta,scale,p,thresh,gap,mlength);
                N_line=N_line+1;
                Line_index=[Line_index;lines_index];
                for k = 1:size(lines_index,1)
                    xy=lines_index(k,:);
                    plot(xy([1,3]),xy([2,4]),'LineWidth',2,'Color','green');
                    plot(xy(1),xy(2),'x','LineWidth',2,'Color','yellow');
                    plot(xy(3),xy(4),'x','LineWidth',2,'Color','red');
                    centroid = (xy([1,2])+xy([3,4]))/2;
                    angl = -1*atand((xy(3)-xy(1))/(xy(4)-xy(2)));
                    centroildAll = [centroildAll;centroid];
                    observedOri = [observedOri; angl];
                    nfiber = nfiber+1;
                end
            end
            TF = TF + nfiber;
            val = num2str(TF);
            set(handles.fiberNum,'string',val);
            
        case 'Gradient HT'
            SK = bwmorph(bw_im,'skel',inf);
            Fibers = bwconncomp(SK);
            fib = cell(Fibers.NumObjects,1);
            nfiber = 0;
            N_line = 0;
            TF = 0;
            Line_index = [];
            centroildAll = [];
            observedOri = [];
            for j = 1:Fibers.NumObjects
                temp = zeros(length(Fibers.PixelIdxList{j}),2);
                ImageTemp = zeros(size(gray_im));
                for k = 1:length(Fibers.PixelIdxList{j})
                    [temp(k,1),temp(k,2)] = ind2sub(size(gray_im), Fibers.PixelIdxList{j}(k));
                    ImageTemp(temp(k,1),temp(k,2)) = 1;
                end
                fib{j}=temp;
                [nlines,lines_index]=findFiber(ImageTemp,rho,theta,scale,p,thresh,gap,mlength);
                N_line=N_line+1;
                Line_index=[Line_index;lines_index];
                for k = 1:size(lines_index,1)
                    xy=lines_index(k,:);
                    plot(xy([1,3]),xy([2,4]),'LineWidth',2,'Color','green');
                    plot(xy(1),xy(2),'x','LineWidth',2,'Color','yellow');
                    plot(xy(3),xy(4),'x','LineWidth',2,'Color','red');
                    centroid = (xy([1,2])+xy([3,4]))/2;
                    angl = -1*atand((xy(3)-xy(1))/(xy(4)-xy(2)));
                    centroildAll = [centroildAll;centroid];
                    observedOri = [observedOri; angl];
                    nfiber = nfiber+1;
                end
            end
            TF = TF + nfiber;
            val = num2str(TF);
            set(handles.fiberNum,'string',val);
            
        case 'Break Merge '
            XY = {};
            tol1 = 8;
            tol2 = 5;
            th_image = bwmorph(bw_im, 'thin', inf);
            Fibers = bwconncomp(th_image);
            fib = cell(Fibers.NumObjects,1);
            for i = 1:Fibers.NumObjects
                temp = zeros(length(Fibers.PixelIdxList{i}),2);
                ImageTemp = zeros(size(I));
                for k = 1:length(Fibers.PixelIdxList{i})
                    [temp(k,1),temp(k,2)] = ind2sub(size(I), Fibers.PixelIdxList{i}(k));
                    ImageTemp(temp(k,1),temp(k,2)) = 1;
                end
                fib{i} = temp;

                [y, x] = find(ImageTemp);
                data = [x, y];

                [IDX, noise, D]=DBSCAN(data,1.5,4);

                K = find(IDX == 0);
                NOISE = data(K,:); % Junction points
                [IDXN, noiseN, DN]=DBSCAN(NOISE,1.5,3);

                SEDGES = {}; % Seperated Edges
                A = []; % Array for angels 
                for j = 1:max(IDXN)
                    L = find(IDXN == j);
                    EP = NOISE(L,:); % EP for edges pixels (edges are seperated at juntions points)
                    if size(EP,1) >= 10 % Skip the small edges
                        X = FurthestPixel(EP); % Find the pixels at max distance of a cluster
                        A = [A, atand((X(1,2)-X(2,2))/(X(1,1)-X(2,1)))];
                        SEDGES = [SEDGES;EP];
                    end
                end

                B = A;
                C = SEDGES;
                MEDGES = {}; % Merged Edges

                while isempty(B) == 0    
                    [IN,dir] = ClosestAngel(B,tol1);
                    Merge = []; % Array of edges to merge
                    Merge = [Merge;cell2mat(C(IN(1)))];
                    DI = 1; % Index of cluster to be deleted
                    P1 = cell2mat(C(IN(1)));
                    for k = 2:length(IN)
                        P2 = cell2mat(C(IN(k)));
                        MP = [P1(round(size(P1,1)/2),:);P2(round(size(P2,1)/2),:)];
                        % MP for Middle Point of two edges
                        tempA = atand((MP(1,2)-MP(2,2))/(MP(1,1)-MP(2,1)));
                        if tempA >= dir-tol2 && tempA <= dir+tol2
                            Merge = [Merge;cell2mat(C(IN(k)))];
                            DI = [DI;IN(k)];
                        end
                    end
                    MEDGES = [MEDGES;Merge];
                    B(DI) = [];
                    C(DI) = [];
                end
                XY = [XY;MEDGES];
            end

            global SE
            SE = {}; % Start and End point of each Merged edges
            for i = 1:size(XY,1)
                SE = [SE;FurthestPixel(XY{i})];
            end

            %Img_gaus = imgaussfilt(Image, sigma);
            %gray_image = mat2gray(Img_gaus, [0 255]);
            %figure, imshow(gray_image)
            %hold on
            TF = 0;
            centroildAll = [];
            observedOri = [];
            for i = 1:size(SE,1)
                cor = SE{i};
                if pdist(cor) > 20
                    plot(cor(:,1), cor(:,2), 'color','g','LineWidth',1.5)
                    plot(cor(1,1), cor(1,2), 'x','color','r','LineWidth',1.5)
                    plot(cor(2,1), cor(2,2), 'x','color','y','LineWidth',1.5)
                    centroid = (cor(1,:)+cor(2,:))/2;
                    angl = -1*atand((cor(2,1)-cor(1,1))/(cor(2,2)-cor(1,2)));
                    centroildAll = [centroildAll;centroid];
                    observedOri = [observedOri; angl];
                    TF = TF+1;
                end
            end
            %hold off
            val = num2str(TF);
            set(handles.fiberNum,'string',val);
    end
    hold off;
    guidata(hObject,handles);



% --- Executes on slider movement.
function rhoResolution_Callback(hObject, eventdata, handles)
s2val = get(handles.rhoResolution,'value');
s2val = num2str(s2val);
set(handles.rhoRes,'string',s2val);


% --- Executes during object creation, after setting all properties.
function rhoResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhoResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function thetaResolution_Callback(hObject, eventdata, handles)
s3val = get(handles.thetaResolution,'value');
s3val = num2str(s3val);
set(handles.Theta,'string',s3val);


% --- Executes during object creation, after setting all properties.
function thetaResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetaResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function noPeak_Callback(hObject, eventdata, handles)
s4val = get(handles.noPeak,'value');
s4val = num2str(s4val);
set(handles.Peaks,'string',s4val);


% --- Executes during object creation, after setting all properties.
function noPeak_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noPeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function suppressScale_Callback(hObject, eventdata, handles)
s5val = get(handles.suppressScale,'value');
s5val = num2str(s5val);
set(handles.Scale,'string',s5val);


% --- Executes during object creation, after setting all properties.
function suppressScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to suppressScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function fillGap_Callback(hObject, eventdata, handles)
s7val = get(handles.fillGap,'value');
s7val = num2str(s7val);
set(handles.Gap,'string',s7val);


% --- Executes during object creation, after setting all properties.
function fillGap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fillGap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function minLength_Callback(hObject, eventdata, handles)
s8val = get(handles.minLength,'value');
s8val = num2str(s8val);
set(handles.Length,'string',s8val);

% --- Executes during object creation, after setting all properties.
function minLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function deltaTheta_Callback(hObject, eventdata, handles)
s10val = get(handles.deltaTheta,'value');
s10val = num2str(s10val);
set(handles.DT,'string',s10val);


% --- Executes during object creation, after setting all properties.
function deltaTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in resetAll.
function resetAll_Callback(hObject, eventdata, handles)
    axes(handles.axes1);
    hold off;
    cla reset;

    axes(handles.axes2);
    hold off;
    cla reset;

    axes(handles.axes3);
    hold off;
    cla reset;
    
    axes(handles.axes4);
    hold off;
    cla reset;
    
    set(handles.fiberNum,'String','');
    set(handles.MSEo,'String','');
    set(handles.MSEc,'String','');
    %set(handles.slider1,'Value',1);
    set(handles.rhoRes,'String','1');
    set(handles.rhoResolution,'Value',1);
    set(handles.Peaks,'String','100');
    set(handles.noPeak,'Value',100);
    set(handles.Theta,'String','1');
    set(handles.thetaResolution,'Value',1);
    set(handles.Gap,'String','40');
    set(handles.fillGap,'Value',40);
    set(handles.Length,'String','40');
    set(handles.minLength,'Value',40);
    set(handles.Scale,'String','50');
    set(handles.suppressScale,'Value',50);
    set(handles.Thresh,'String','0.5');
    set(handles.threshCof,'Value',0.5);
    set(handles.SE,'String','80');
    set(handles.lengthSE,'Value',80);
    set(handles.DT,'String','2');
    set(handles.deltaTheta,'Value',2);
    set(handles.CT,'String','50');
    set(handles.centerTol,'Value',50);
    set(handles.TT,'String','4');
    set(handles.thetaTol,'Value',4);
    
    set(handles.method,'selectedObject',handles.SHT);

% --- Executes on button press in exitGUI.
function exitGUI_Callback(hObject, eventdata, handles)
    close all;


% --- Executes on slider movement.
function lengthSE_Callback(hObject, eventdata, handles)
s9val = get(handles.lengthSE,'value');
s9val = num2str(s9val);
set(handles.SE,'string',s9val);

% --- Executes during object creation, after setting all properties.
function lengthSE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lengthSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function centerTol_Callback(hObject, eventdata, handles)
s11val = get(handles.centerTol,'value');
s11val = num2str(s11val);
set(handles.CT,'string',s11val);


% --- Executes during object creation, after setting all properties.
function centerTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centerTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function thetaTol_Callback(hObject, eventdata, handles)
s12val = get(handles.thetaTol,'value');
s12val = num2str(s12val);
set(handles.TT,'string',s12val);


% --- Executes during object creation, after setting all properties.
function thetaTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetaTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function threshCof_Callback(hObject, eventdata, handles)
s13val = get(handles.threshCof,'value');
s13val = num2str(s13val);
set(handles.Thresh,'string',s13val);



% --- Executes during object creation, after setting all properties.
function threshCof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshCof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in method.
function method_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in method 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function rhoRes_Callback(hObject, eventdata, handles)
    rho1 = get(hObject,'String');
    set(handles.rhoResolution,'value',str2num(rho1));
    guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rhoRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhoRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Peaks_Callback(hObject, eventdata, handles)
    peak1 = get(hObject,'String');
    set(handles.noPeak,'value',str2num(peak1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Peaks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gap_Callback(hObject, eventdata, handles)
    gap1 = get(hObject,'String');
    set(handles.fillGap,'value',str2num(gap1));
    guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function Gap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Thresh_Callback(hObject, eventdata, handles)
    thresh1 = get(hObject,'String');
    set(handles.threshCof,'value',str2num(thresh1));
    guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function Thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_Callback(hObject, eventdata, handles)
    theta1 = get(hObject,'String');
    set(handles.thetaResolution,'value',str2num(theta1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SE_Callback(hObject, eventdata, handles)
    se1 = get(hObject,'String');
    set(handles.lengthSE,'value',str2num(se1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function SE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DT_Callback(hObject, eventdata, handles)
    dt1 = get(hObject,'String');
    set(handles.deltaTheta,'value',str2num(dt1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function DT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CT_Callback(hObject, eventdata, handles)
    ct1 = get(hObject,'String');
    set(handles.centerTol,'value',str2num(ct1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TT_Callback(hObject, eventdata, handles)
    tt1 = get(hObject,'String');
    set(handles.thetaTol,'value',str2num(tt1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function TT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_Callback(hObject, eventdata, handles)
    scale1 = get(hObject,'String');
    set(handles.suppressScale,'value',str2num(scale1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_Callback(hObject, eventdata, handles)
    mlength1 = get(hObject,'String');
    set(handles.minLength,'value',str2num(mlength1));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NOF_Callback(hObject, eventdata, handles)
% hObject    handle to NOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NOF as text
%        str2double(get(hObject,'String')) returns contents of NOF as a double


% --- Executes during object creation, after setting all properties.
function NOF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on selection change in popupmenu2.
% function popupmenu2_Callback(hObject, eventdata, handles)
% % hObject    handle to popupmenu2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from popupmenu2
% 
% 
% % --- Executes during object creation, after setting all properties.
% function popupmenu2_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to popupmenu2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on selection change in Orientation.
function Orientation_Callback(hObject, eventdata, handles)
% hObject    handle to Orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Orientation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Orientation


% --- Executes during object creation, after setting all properties.
function Orientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulateImage.
function simulateImage_Callback(hObject, eventdata, handles)
% hObject    handle to simulateImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axes(handles.axes1);
    hold off;
    cla reset;
    
    Nx = 2400;
    Ny = 1800;
    %rhores = 0.00005;
    %lambda = Nx*Ny*rhores;
    
    global N;
    N = str2double(get(handles.NOF,'string'));
    ori_type =get(handles.Orientation,'value');
    
    alpha = zeros(N,1);
    global L
    L = zeros(N,1);
    W = zeros(N,1);
    Density = zeros(N,1);
    S = 2;
    
    for i = 1:N
        alpha(i) = orientation(ori_type);
        L(i) = normrnd(100,20);
        W(i) = 2*S;
        Density(i) = rtnorm(0,255,192,32);
    end
    global fiber_angle
    fiber_angle = alpha;
    
    X = randi([1,Nx],[N,1]);
    Y = randi([1,Ny],[N,1]);
    
    global I
    I = zeros(Ny,Nx);
    [XX,YY] = meshgrid(1:Nx,1:Ny);
    
    global actualCentroid;
    actualCentroid = [];
    for i = 1:N
        p1 = [-W(i)/2;-L(i)/2];
        p2 = [W(i)/2;-L(i)/2];
        p3 = [W(i)/2;L(i)/2];
        p4 = [-W(i)/2;L(i)/2];
        
        Theta = alpha(i);
        R = [cos(Theta) -sin(Theta);sin(Theta) cos(Theta)];
        center = [X(i);Y(i)];
        actualCentroid = [actualCentroid; center'];
        
        q1 = R*p1 + center;
        q2 = R*p2 + center;
        q3 = R*p3 + center;
        q4 = R*p4 + center;
        
        V_XY = [q1,q2,q3,q4];
        Wire=inpolygon(XX,YY,V_XY(1,:),V_XY(2,:))*Density(i);
        I = max(I, Wire);
    end
    axes(handles.axes1);
    imshow(I);   

% --- Executes on button press in extractMorphology.
function extractMorphology_Callback(hObject, eventdata, handles)
% hObject    handle to extractMorphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % ===============Orientation Distribution======================%
    handles.output = hObject;
    axes(handles.axes3);
    hold off;
    cla reset;
    
    global SE
    global fiber_angle
    deltatheta = 5;
    num = zeros(180/deltatheta,1);
    num(:,1) = -90+deltatheta:deltatheta:90;
    distributions = flipud(num);
    
    Ori = []; % Emperical Orientation
    EL = []; % Emperical Length 
    for i = 1:size(SE,1)
        d = pdist(SE{i});
        if d > 20
            len = sqrt((SE{i}(2,1)-SE{i}(1,1))^2 + (SE{i}(2,2)-SE{i}(1,2))^2);
            EL = [EL; len];
            angl = -1*atan((SE{i}(2,1)-SE{i}(1,1))/(SE{i}(2,2)-SE{i}(1,2)));
            Ori = [Ori;angl];
        end
    end
    
    PD1 = findpb(num,deltatheta,fiber_angle);
    PD2 = findpb(num,deltatheta,Ori);
    
    distribution1 = [distributions PD1];
    distribution2 = [distributions PD2];
    
    %figure;
    axes(handles.axes3);
    hold on;
    bar(distribution1(:,1),distribution1(:,2),'g');
    hold on
    scatter(distribution2(:,1),distribution2(:,2),'r *');
    hold on
    grid off
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    ylim([0,0.25])
    xlim([-90,90])
    xlabel('\theta')
    ylabel('P(\theta)')
    % legend('Location','northwest','True data','Emoerical distribution');
    legend('True data','Emperical distribution');
    title('Orientation Distribution')
    %hold off
    
    hold off;
    guidata(hObject,handles);
    
%     ori_type =get(handles.Orientation,'value');
%     if ori_type == 1 
%         S = 'U(-90,90)';
%         set(handles.MSEo,'String',S);
%     elseif ori_type == 2
%         S = 'N(50,15)';
%         set(handles.MSEo,'String',S);
%     elseif ori_type == 3
%         S = 'N(-30,5)';
%         set(handles.MSEo,'String',S);
%     elseif ori_type == 4
%         S = 'N(30,10)&N(-60,10)';
%         set(handles.MSEo,'String',S);
%     elseif ori_type == 5
%         S = '(-70,-30,30)';
%         set(handles.MSEo,'String',S);
%     elseif ori_type == 6
%         S = '(45)';
%         set(handles.MSEo,'String',S);    
%     end
        
    %====================================================================%
    
    %=========================Length Distribution========================%
%     handles.output = hObject;
%     axes(handles.axes4);
%     hold off;
%     cla reset;
%     
%     global L
%     deltaL = 4;
%     %maxL = round(max(EL))+10;
%     %minL = round(min(EL))-10;
%     Length = zeros(160/deltaL,1);
%     Length(:,1) = 20+deltaL:4:180;
%     %Length = zeros((round((maxL-minL)/deltaL)*deltaL)/deltaL,1);
%     %Length(:,1) = minL+deltaL:deltaL:maxL;
%     
%     LD1 = lengthD(Length,deltaL,L);
%     LD2 = lengthD(Length,deltaL,EL);
%     
%     distributionL1 = [Length LD1];
%     distributionL2 = [Length LD2];
%     
%     bar(distributionL1(:,1), distributionL1(:,2), 'g')
%     hold on
%     scatter(distributionL2(:,1),distributionL2(:,2),'r *');
% 
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
%     ylim([0,0.5])
%     xlim([20,180])
%     xlabel('Fiber length')
%     ylabel('P(length)')
%     % legend('Location','northwest','True data','Emoerical distribution');
%     legend('True data','Emperical distribution');
%     title('Length Distribution')
%     
%     hold off;
%     guidata(hObject,handles);
%     
%     S = 'N(100,20)';
%     set(handles.MSEc,'String',S);
    
    %====================================================================%
    
    %============================Plot Centroid===========================%
    handles.output = hObject;
    axes(handles.axes4);
    hold off;
    cla reset;
   
    global centroildAll;
    axes(handles.axes4);
    hold on;
    axis ij
    plot(centroildAll(:,1), centroildAll(:,2), 'o', 'color', 'b', 'LineWidth', 0.5)
    title('Spatial Distribution')
    hold off
    
    %==================MSE of Orientation and Centroid====================%
    global actualCentroid;
    global observedOri;
    ActualMorphology = [actualCentroid, fiber_angle*180/pi];
    ObservedMorphology = [centroildAll, observedOri];

    M = [];
    for i = 1:length(ObservedMorphology)
        A = ObservedMorphology(i,:);
        distances = sqrt(sum(bsxfun(@minus, ActualMorphology, A).^2,2));
        if min(distances) <= 40
            closest = ActualMorphology(find(distances==min(distances)),:);
            M = [M; A, closest];
        end
    end

    MSEc = mean((M(:,1)-M(:,4)).^2 + (M(:,2)-M(:,5)).^2);
    MSEo = mean((M(:,3)-M(:,6)).^2);
    So = num2str(MSEo);
    Sc = num2str(MSEc);
    set(handles.MSEc,'String',Sc);
    set(handles.MSEo,'String',So);
