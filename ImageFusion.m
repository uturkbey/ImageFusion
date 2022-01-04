classdef ImageFusion < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SaveImage_Button                matlab.ui.control.Button
        ImageFusionMethodsButtonGroup   matlab.ui.container.ButtonGroup
        Average_Button                  matlab.ui.control.RadioButton
        Min_Button                      matlab.ui.control.RadioButton
        Max_Button                      matlab.ui.control.RadioButton
        PCA_Button                      matlab.ui.control.RadioButton
        StationaryWavleletTransformL1_Button  matlab.ui.control.RadioButton
        DiscreteCosineTransform_Button  matlab.ui.control.RadioButton
        DiscreteWaveletTransform_Button  matlab.ui.control.RadioButton
        wfusingBuiltInMATLABFunction_Button  matlab.ui.control.RadioButton
        StationaryWavleletTransformL2_Button  matlab.ui.control.RadioButton
        ImagePanel_A                    matlab.ui.control.UIAxes
        ImagePanel_Fused                matlab.ui.control.UIAxes
        ImagePanel_B                    matlab.ui.control.UIAxes
        LoadImageB_Button               matlab.ui.control.Button
        LoadImageA_Button               matlab.ui.control.Button
        FUSE_Button                     matlab.ui.control.Button
        ImageFusionLabel                matlab.ui.control.Label
        Image                           matlab.ui.control.Image
        MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label  matlab.ui.control.Label
        Image2                          matlab.ui.control.Image
        MessagePromptTextAreaLabel      matlab.ui.control.Label
        MessagePromptTextArea           matlab.ui.control.TextArea
        withVariousSignalProcessingMethodsLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        ImageA = 0; %Variable to hold the Image A data
        ImageB = 0; %Variable to hold the Image B data
        ImageFused = 0; %Variable to hold the Image Data which is calculated by the fusion of Image A and Image B using one of the methods provided.
        method = "x"; %To keep track of last used method for saved images
    end
    
    methods (Access = private)
        
        function imfuse_avg(app)
            app.ImageFused = (app.ImageA + app.ImageB) / 2;
        end
        
        function imfuse_min(app)
            mm = app.ImageA < app.ImageB;
            app.ImageFused  = (mm.*app.ImageA) + ((~mm).*app.ImageB);
        end
        
        function imfuse_max(app)
            mm = app.ImageA > app.ImageB;
            app.ImageFused  = (mm.*app.ImageA) + ((~mm).*app.ImageB);
        end
        
        function imfuse_pca(app)
            % getting Image dimensions
            size1=size(app.ImageA);
            size2=size(app.ImageB);
            
            %vectorizing Images
            vec1=app.ImageA(:);
            vec2=app.ImageB(:);
            
            %covariance of vectors
            cov_mat=cov(vec1,vec2);
            
            %eigenvalues and eigenvectors 
            [V,D] = eig(cov_mat);
            
            sum1=sum(V,1);
            
            % Calculating PCA
            if D(1,1) >= D(2,2)
                pca = V(:,1)/sum1(1);
            else
                pca = V(:,2)/sum1(2);
            end
            
            app.ImageFused = pca(1)*app.ImageA + pca(2)*app.ImageB;
        end
        
        function imfuse_dctma(app)
            block_size=8;
            [m,n] = size(app.ImageA);
            
            for i=1:block_size:m
                for j=1:block_size:n
                    
                    %dividing image to blocks
                    block1 = app.ImageA(i:i+block_size-1,j:j+block_size-1);
                    block2 = app.ImageB(i:i+block_size-1,j:j+block_size-1);
                    
                    block1_dct = dct2(block1);
                    block2_dct = dct2(block2);
                    
                    %choosing the largest magnitude AC coefficients
                    dl = abs(block1_dct)-abs(block2_dct)>=0;
                    CBF = dl.*block1_dct+(~dl).*block2_dct;
                    
                    CBF(1,1)=(block1_dct(1,1)+block2_dct(1,1))/2; %averaging dc 
                    
                    app.ImageFused(i:i+block_size-1,j:j+block_size-1)=idct2(CBF); %inverse discrete cosine transform
                end
            end
        end
        
        function imfuse_dwt(app)
            %This function implements image fusion with Discrete Wavelet Transform(DWT)
            [a1,b1,c1,d1]=dwt2(app.ImageA,'db2');
            [a2,b2,c2,d2]=dwt2(app.ImageB,'db2');
            [k1,k2]=size(a1);
 
            %% Average Rule
            for i=1:k1
                for j=1:k2
                  a3(i,j)=(a1(i,j)+a2(i,j))/2;
               end
            end
            
            %% Max Rule
            for i=1:k1
              for j=1:k2
                  b3(i,j)=max(b1(i,j),b2(i,j));
                  c3(i,j)=max(c1(i,j),c2(i,j));
                  d3(i,j)=max(d1(i,j),d2(i,j));
              end
            end
            
            %% Inverse Wavelet Transform 
            app.ImageFused=idwt2(a3,b3,c3,d3,'db2');
        end
        
        function imfuse_swtL1(app)
            %This function implements image fusion with one level Discrete Stationary Wavelet Transform (SWT)

            %When appliying one level SWT, image dimensions should be even.  
            
            %Image decomposition using SWT
            [A1L1,H1L1,V1L1,D1L1] = swt2(app.ImageA,1,'sym2');
            [A2L1,H2L1,V2L1,D2L1] = swt2(app.ImageB,1,'sym2');
            
            %Fusion starts - Level 1
            AfL1 = 0.5*(A1L1+A2L1);
            D = (abs(H1L1)-abs(H2L1))>=0;
            HfL1 = D.*H1L1 + (~D).*H2L1;
            D = (abs(V1L1)-abs(V2L1))>=0;
            VfL1 = D.*V1L1 + (~D).*V2L1;
            D = (abs(D1L1)-abs(D2L1))>=0;
            DfL1 = D.*D1L1 + (~D).*D2L1;
            
            %Fused image formation
            app.ImageFused = iswt2(AfL1,HfL1,VfL1,DfL1,'sym2');
        end
        
        function imfuse_swtL2(app)
            %This function implements image fusion with two level Discerete Stationary Wavelet Transform (SWT)

            %When appliying SWT, 2^Level has to divide the size of the image. 
            %In this case, it should be at least even. Be careful: level is always one
            %in swt2 function calls and both A1L1 and A2L1 are at the same size with 
            %img1 and img2. Therefore even the level of SWT is 2, being even is enough.
            
            %Image decomposition using SWT
            [A1L1,H1L1,V1L1,D1L1] = swt2(app.ImageA,1,'sym2');
            [A2L1,H2L1,V2L1,D2L1] = swt2(app.ImageB,1,'sym2');
            [A1L2,H1L2,V1L2,D1L2] = swt2(A1L1,1,'sym2');
            [A2L2,H2L2,V2L2,D2L2] = swt2(A2L1,1,'sym2');
            %Fusion at Level 2
            AfL2 = 0.5*(A1L2+A2L2);
            D = (abs(H1L2)-abs(H2L2))>=0;
            HfL2 = D.*H1L2 + (~D).*H2L2;
            D = (abs(V1L2)-abs(V2L2))>=0;
            VfL2 = D.*V1L2 + (~D).*V2L2;
            D = (abs(D1L2)-abs(D2L2))>=0;
            DfL2 = D.*D1L2 + (~D).*D2L2;
            
            %Fusion at Level 1
            D = (abs(H1L1)-abs(H2L1))>=0;
            HfL1 = D.*H1L1 + (~D).*H2L1;
            D = (abs(V1L1)-abs(V2L1))>=0;
            VfL1 = D.*V1L1 + (~D).*V2L1;
            D = (abs(D1L1)-abs(D2L1))>=0;
            DfL1 = D.*D1L1 + (~D).*D2L1;
            
            %Fused image formation
            AfL1 = iswt2(AfL2,HfL2,VfL2,DfL2,'sym2');
            app.ImageFused = iswt2(AfL1,HfL1,VfL1,DfL1,'sym2');
        end
        
        function imfuse_wfusimg(app)
            %This function applies image fusing using a built in Matlab function called wfusimg
            %For more detail about wfusimg function in Wavelet Toolbox, see https://www.mathworks.com/help/wavelet/ref/wfusimg.html

            wv = "db2";
            lv = 5;
            app.ImageFused = wfusimg(app.ImageA, app.ImageB, wv, lv, "min", "max");
        end
        
        function gray_image = looseColors(app, image)
            %This function first checks if the input image is RGB.
            %Then converts it to double valued grayscale if RGB 
            %or just makes it double if already grayscale.
            
            if size(image, 3) == 3
                gray_image = double(rgb2gray(image));
            else
                gray_image = double(image);
            end
            
            %Finally checks if the input image size is even
            %if not pads zeros to make image sizes even
            sizes=size(gray_image);
            if mod(sizes(1),2)~=0 
                gray_image(sizes(1)+1,:)=0;
            end
                
            if mod(sizes(2),2)~=0 
                gray_image(:,sizes(2)+1)=0;
            end    
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %Show default images -- Later
        end

        % Button pushed function: LoadImageA_Button
        function LoadImageA_ButtonPushed(app, event)
            [pathA,~] = imgetfile(); %Ask user for an input image to be fused
            try
                app.ImageA = looseColors(app, imread(pathA)); %Read the choosen image and turn it into grayscale
                imshow(app.ImageA, [], 'Parent', app.ImagePanel_A); %Prompt the image to the related panel
                app.MessagePromptTextArea.Value = "Image A loaded successfully.";
            catch
                if size(app.ImageA) == 1
                    app.MessagePromptTextArea.Value = "Invalid path to Image A file! Try choosing a different path.";
                end
            end
        end

        % Button pushed function: LoadImageB_Button
        function LoadImageB_ButtonPushed(app, event)
            [pathB,~] = imgetfile(); %Ask user for an input image to be fused
            try
                app.ImageB = looseColors(app, imread(pathB)); %Read the choosen image and turn it into grayscale
                imshow(app.ImageB, [], 'Parent', app.ImagePanel_B); %Prompt the image to the related panel
                app.MessagePromptTextArea.Value = "Image B loaded successfully.";
            catch
                if size(app.ImageB) == 1
                    app.MessagePromptTextArea.Value = "Invalid path to Image B file! Try choosing a different path.";
                end
            end
        end

        % Button pushed function: FUSE_Button
        function FUSE_ButtonPushed(app, event)
          if size(app.ImageA) ~= 1 
              if size(app.ImageB) ~= 1 
                  if size(app.ImageA) == size(app.ImageB)
                    try
                        %Apply choosen fusion method and keep the method name as a
                        %string to use in save function
                        if app.Average_Button.Value == true
                            imfuse_avg(app);
                            app.method = "avg";
                        elseif app.Min_Button.Value == true
                            imfuse_min(app);
                            app.method = "min";
                        elseif app.Max_Button.Value == true
                            imfuse_max(app);
                            app.method = "max";
                        elseif app.PCA_Button.Value == true
                            imfuse_pca(app);
                            app.method = "pca";
                        elseif app.DiscreteWaveletTransform_Button.Value == true
                            imfuse_dwt(app);
                            app.method = "dwt";
                        elseif app.StationaryWavleletTransformL1_Button.Value == true
                            imfuse_swtL1(app);
                            app.method = "swtL1";
                        elseif app.StationaryWavleletTransformL2_Button.Value == true
                            imfuse_swtL2(app);
                            app.method = "swtL2";
                        elseif app.DiscreteCosineTransform_Button.Value == true
                            imfuse_dctma(app);
                            app.method = "dctma";
                        elseif app.wfusingBuiltInMATLABFunction_Button.Value == true
                            imfuse_wfusimg(app);
                            app.method = "wfusimg";
                        end
                        %Normalize the Fused Image
                        app.ImageFused = (app.ImageFused - min(app.ImageFused(:))) ./ (max(app.ImageFused(:)) - min(app.ImageFused(:)));
                        imshow(app.ImageFused, 'Parent', app.ImagePanel_Fused); %Prompt the fused image to the related panel
                        app.MessagePromptTextArea.Value = strcat("Image Fused successfully with ", app.method, " method.");
                    catch
                        app.MessagePromptTextArea.Value = "An unknown error occured during Image Fusion while applying the choosen fusion method!";
                    end
                  else
                    app.MessagePromptTextArea.Value = "Images to be fused are not of the same size! Please choose appropriate images.";
                  end
              else
                  app.MessagePromptTextArea.Value = "Image B is missing! Please load Image B.";
              end
          else
              app.MessagePromptTextArea.Value = "Image A is missing! Please load Image A.";
          end
        end

        % Button pushed function: SaveImage_Button
        function SaveImage_ButtonPushed(app, event)
            c = clock; %To save the fused images with date and hour in adddition to the fusion method
            filename = strcat('fused_image_', app.method, '_', string(c(1)), '-', ...
                              string(c(2)), '-', string(c(3)), '_', string(c(4)), '-', string(c(5)), '-', string(uint8(c(6))), '.jpg');
            if size(app.ImageFused) == 1
                app.MessagePromptTextArea.Value = "There is no Fused Image to be saved, yet! Try to upload some images and clicking the FUSE!!! button :)";
            else
                try
                    imwrite(app.ImageFused, filename);
                    app.MessagePromptTextArea.Value = strcat("Fused Image saved to the current working directory with the filename: ", filename);
                catch
                    app.MessagePromptTextArea.Value = "An unknown error occured during saving the Fused Image.";
                end            
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.8 0.8 0.8];
            app.UIFigure.Position = [100 100 849 653];
            app.UIFigure.Name = 'MATLAB App';

            % Create SaveImage_Button
            app.SaveImage_Button = uibutton(app.UIFigure, 'push');
            app.SaveImage_Button.ButtonPushedFcn = createCallbackFcn(app, @SaveImage_ButtonPushed, true);
            app.SaveImage_Button.BackgroundColor = [0.7176 0.2745 1];
            app.SaveImage_Button.FontName = 'Comic Sans MS';
            app.SaveImage_Button.FontSize = 20;
            app.SaveImage_Button.FontWeight = 'bold';
            app.SaveImage_Button.Position = [389 150 100 74];
            app.SaveImage_Button.Text = {'SAVE'; 'IMAGE'};

            % Create ImageFusionMethodsButtonGroup
            app.ImageFusionMethodsButtonGroup = uibuttongroup(app.UIFigure);
            app.ImageFusionMethodsButtonGroup.TitlePosition = 'centertop';
            app.ImageFusionMethodsButtonGroup.Title = 'Image Fusion Methods';
            app.ImageFusionMethodsButtonGroup.BackgroundColor = [0.6314 0.902 0.6627];
            app.ImageFusionMethodsButtonGroup.FontName = 'Comic Sans MS';
            app.ImageFusionMethodsButtonGroup.FontWeight = 'bold';
            app.ImageFusionMethodsButtonGroup.FontSize = 15;
            app.ImageFusionMethodsButtonGroup.Position = [394 362 424 153];

            % Create Average_Button
            app.Average_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.Average_Button.Text = 'Average';
            app.Average_Button.FontName = 'Comic Sans MS';
            app.Average_Button.Position = [11 100 68 22];
            app.Average_Button.Value = true;

            % Create Min_Button
            app.Min_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.Min_Button.Text = 'Min';
            app.Min_Button.FontName = 'Comic Sans MS';
            app.Min_Button.Position = [11 78 65 22];

            % Create Max_Button
            app.Max_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.Max_Button.Text = 'Max';
            app.Max_Button.FontName = 'Comic Sans MS';
            app.Max_Button.Position = [11 56 65 22];

            % Create PCA_Button
            app.PCA_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.PCA_Button.Text = 'PCA';
            app.PCA_Button.FontName = 'Comic Sans MS';
            app.PCA_Button.Position = [11 35 44 22];

            % Create StationaryWavleletTransformL1_Button
            app.StationaryWavleletTransformL1_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.StationaryWavleletTransformL1_Button.Text = '1 Level Stationary Wavelet Transform (SWT1L)';
            app.StationaryWavleletTransformL1_Button.FontName = 'Comic Sans MS';
            app.StationaryWavleletTransformL1_Button.Position = [95 78 294 22];

            % Create DiscreteCosineTransform_Button
            app.DiscreteCosineTransform_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.DiscreteCosineTransform_Button.Text = 'Discrete Cosine Transform with Max Coeffs. (DCTma)';
            app.DiscreteCosineTransform_Button.FontName = 'Comic Sans MS';
            app.DiscreteCosineTransform_Button.Position = [95 35 330 22];

            % Create DiscreteWaveletTransform_Button
            app.DiscreteWaveletTransform_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.DiscreteWaveletTransform_Button.Text = 'Discrete Wavelet Transform (DWT)';
            app.DiscreteWaveletTransform_Button.FontName = 'Comic Sans MS';
            app.DiscreteWaveletTransform_Button.Position = [95 100 228 22];

            % Create wfusingBuiltInMATLABFunction_Button
            app.wfusingBuiltInMATLABFunction_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.wfusingBuiltInMATLABFunction_Button.Text = 'wfusimg (Built-In MATLAB Function)';
            app.wfusingBuiltInMATLABFunction_Button.FontName = 'Comic Sans MS';
            app.wfusingBuiltInMATLABFunction_Button.Position = [95 14 231 22];

            % Create StationaryWavleletTransformL2_Button
            app.StationaryWavleletTransformL2_Button = uiradiobutton(app.ImageFusionMethodsButtonGroup);
            app.StationaryWavleletTransformL2_Button.Text = '2 Level Stationary Wavelet Transform (SWT2L)';
            app.StationaryWavleletTransformL2_Button.FontName = 'Comic Sans MS';
            app.StationaryWavleletTransformL2_Button.Position = [95 56 298 22];

            % Create ImagePanel_A
            app.ImagePanel_A = uiaxes(app.UIFigure);
            title(app.ImagePanel_A, 'IMAGE A')
            xlabel(app.ImagePanel_A, '')
            ylabel(app.ImagePanel_A, '')
            app.ImagePanel_A.FontName = 'Comic Sans MS';
            app.ImagePanel_A.FontSize = 15;
            app.ImagePanel_A.FontWeight = 'bold';
            app.ImagePanel_A.TickLabelInterpreter = 'none';
            app.ImagePanel_A.Box = 'on';
            app.ImagePanel_A.XTick = [];
            app.ImagePanel_A.XTickLabel = '';
            app.ImagePanel_A.YTick = [];
            app.ImagePanel_A.YTickLabel = '';
            app.ImagePanel_A.ZTick = [];
            app.ImagePanel_A.Color = [0.6706 0.8118 0.9804];
            app.ImagePanel_A.BackgroundColor = [0.8 0.8 0.8];
            app.ImagePanel_A.Position = [14 318 321 181];

            % Create ImagePanel_Fused
            app.ImagePanel_Fused = uiaxes(app.UIFigure);
            title(app.ImagePanel_Fused, 'FUSED IMAGE')
            xlabel(app.ImagePanel_Fused, '')
            ylabel(app.ImagePanel_Fused, '')
            app.ImagePanel_Fused.FontName = 'Comic Sans MS';
            app.ImagePanel_Fused.FontSize = 15;
            app.ImagePanel_Fused.FontWeight = 'bold';
            app.ImagePanel_Fused.TickLabelInterpreter = 'none';
            app.ImagePanel_Fused.Box = 'on';
            app.ImagePanel_Fused.XTick = [];
            app.ImagePanel_Fused.XTickLabel = '';
            app.ImagePanel_Fused.YTick = [];
            app.ImagePanel_Fused.YTickLabel = '';
            app.ImagePanel_Fused.ZTick = [];
            app.ImagePanel_Fused.Color = [0.8902 0.4431 0.4431];
            app.ImagePanel_Fused.BackgroundColor = [0.8 0.8 0.8];
            app.ImagePanel_Fused.Position = [500 97 307 250];

            % Create ImagePanel_B
            app.ImagePanel_B = uiaxes(app.UIFigure);
            title(app.ImagePanel_B, 'IMAGE B')
            xlabel(app.ImagePanel_B, '')
            ylabel(app.ImagePanel_B, '')
            app.ImagePanel_B.FontName = 'Comic Sans MS';
            app.ImagePanel_B.FontSize = 15;
            app.ImagePanel_B.FontWeight = 'bold';
            app.ImagePanel_B.TickLabelInterpreter = 'none';
            app.ImagePanel_B.Box = 'on';
            app.ImagePanel_B.XTick = [];
            app.ImagePanel_B.XTickLabel = '';
            app.ImagePanel_B.YTick = [];
            app.ImagePanel_B.YTickLabel = '';
            app.ImagePanel_B.ZTick = [];
            app.ImagePanel_B.Color = [0.6706 0.8118 0.9804];
            app.ImagePanel_B.BackgroundColor = [0.8 0.8 0.8];
            app.ImagePanel_B.Position = [8 68 319.777202072539 180.714640198511];

            % Create LoadImageB_Button
            app.LoadImageB_Button = uibutton(app.UIFigure, 'push');
            app.LoadImageB_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadImageB_ButtonPushed, true);
            app.LoadImageB_Button.Icon = 'insert-image-icon.jpg';
            app.LoadImageB_Button.IconAlignment = 'right';
            app.LoadImageB_Button.FontName = 'Comic Sans MS';
            app.LoadImageB_Button.FontSize = 15;
            app.LoadImageB_Button.FontWeight = 'bold';
            app.LoadImageB_Button.Position = [116 21 119 48];
            app.LoadImageB_Button.Text = {'Load'; 'Image B'};

            % Create LoadImageA_Button
            app.LoadImageA_Button = uibutton(app.UIFigure, 'push');
            app.LoadImageA_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadImageA_ButtonPushed, true);
            app.LoadImageA_Button.Icon = 'insert-image-icon.jpg';
            app.LoadImageA_Button.IconAlignment = 'right';
            app.LoadImageA_Button.FontName = 'Comic Sans MS';
            app.LoadImageA_Button.FontSize = 15;
            app.LoadImageA_Button.FontWeight = 'bold';
            app.LoadImageA_Button.Position = [116 271 120 48];
            app.LoadImageA_Button.Text = {'Load'; 'Image A'};

            % Create FUSE_Button
            app.FUSE_Button = uibutton(app.UIFigure, 'push');
            app.FUSE_Button.ButtonPushedFcn = createCallbackFcn(app, @FUSE_ButtonPushed, true);
            app.FUSE_Button.BackgroundColor = [0 1 0];
            app.FUSE_Button.FontName = 'Comic Sans MS';
            app.FUSE_Button.FontSize = 25;
            app.FUSE_Button.FontWeight = 'bold';
            app.FUSE_Button.FontColor = [1 0 0];
            app.FUSE_Button.Position = [390 248 100 71];
            app.FUSE_Button.Text = 'FUSE!!!';

            % Create ImageFusionLabel
            app.ImageFusionLabel = uilabel(app.UIFigure);
            app.ImageFusionLabel.BackgroundColor = [0.8 0.8 0.8];
            app.ImageFusionLabel.FontName = 'Comic Sans MS';
            app.ImageFusionLabel.FontSize = 45;
            app.ImageFusionLabel.FontWeight = 'bold';
            app.ImageFusionLabel.FontAngle = 'italic';
            app.ImageFusionLabel.FontColor = [0.0706 0.3412 0.0314];
            app.ImageFusionLabel.Position = [249 575 295 69];
            app.ImageFusionLabel.Text = 'Image Fusion';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [15 514 190 153];
            app.Image.ImageSource = 'Matrix_Screens.jpg';

            % Create MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label = uilabel(app.UIFigure);
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.HorizontalAlignment = 'center';
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.FontName = 'Comic Sans MS';
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.FontSize = 15;
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.FontColor = [0.0706 0.3412 0.0314];
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.Position = [561 544 135 83];
            app.MehmetYiitAvcUtkuTrkbeyBOUNEEV10Label.Text = {'Mehmet Yiğit Avcı'; 'Utku Türkbey'; 'BOUN EE'; 'V1.0'};

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Position = [719 528 109 115];
            app.Image2.ImageSource = 'Boun_Logo.png';

            % Create MessagePromptTextAreaLabel
            app.MessagePromptTextAreaLabel = uilabel(app.UIFigure);
            app.MessagePromptTextAreaLabel.HorizontalAlignment = 'center';
            app.MessagePromptTextAreaLabel.FontName = 'Comic Sans MS';
            app.MessagePromptTextAreaLabel.FontSize = 14;
            app.MessagePromptTextAreaLabel.FontWeight = 'bold';
            app.MessagePromptTextAreaLabel.Position = [327 49 63 38];
            app.MessagePromptTextAreaLabel.Text = {'Message'; 'Prompt'};

            % Create MessagePromptTextArea
            app.MessagePromptTextArea = uitextarea(app.UIFigure);
            app.MessagePromptTextArea.HorizontalAlignment = 'center';
            app.MessagePromptTextArea.FontName = 'Comic Sans MS';
            app.MessagePromptTextArea.FontWeight = 'bold';
            app.MessagePromptTextArea.Position = [400 13 417 85];
            app.MessagePromptTextArea.Value = {'Welcome to the Image Fusion App :)  Starter guide below!'; ''; 'You may start by uploading two images of the same size.'; 'Later, you can chooce any method you want from the table.'; 'Then you can push the FUSE!!! button to see some magic(?).'; 'Finally, you can save your cool fused image and show it to your'; 'friends, so that they will believe the magic(?) you experienced :)'};

            % Create withVariousSignalProcessingMethodsLabel
            app.withVariousSignalProcessingMethodsLabel = uilabel(app.UIFigure);
            app.withVariousSignalProcessingMethodsLabel.HorizontalAlignment = 'right';
            app.withVariousSignalProcessingMethodsLabel.FontName = 'Comic Sans MS';
            app.withVariousSignalProcessingMethodsLabel.FontSize = 17;
            app.withVariousSignalProcessingMethodsLabel.FontColor = [0.0706 0.3412 0.0314];
            app.withVariousSignalProcessingMethodsLabel.Position = [222 550 329 26];
            app.withVariousSignalProcessingMethodsLabel.Text = ' with Various Signal Processing Methods';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ImageFusion

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end