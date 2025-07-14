classdef Project_of_Dsp_App_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        Project_Dsp_2024UIFigure        matlab.ui.Figure
        Panel_13                        matlab.ui.container.Panel
        Transfer1Button_6               matlab.ui.control.Button
        CopyButton_6                    matlab.ui.control.Button
        Clear1Button_6                  matlab.ui.control.Button
        SaveInputAButton_2              matlab.ui.control.Button
        Transfer1Button_5               matlab.ui.control.Button
        CopyButton_5                    matlab.ui.control.Button
        Clear1Button_5                  matlab.ui.control.Button
        SaveInputBButton_3              matlab.ui.control.Button
        Clear1Button_4                  matlab.ui.control.Button
        SaveOutputButton_2              matlab.ui.control.Button
        UIAxes_6                        matlab.ui.control.UIAxes
        UIAxes_5                        matlab.ui.control.UIAxes
        UIAxes_4                        matlab.ui.control.UIAxes
        Panel_12                        matlab.ui.container.Panel
        Clear1Button_2                  matlab.ui.control.Button
        SaveOutputButton                matlab.ui.control.Button
        UIAxes_3                        matlab.ui.control.UIAxes
        Panel_11                        matlab.ui.container.Panel
        Transfer1Button_3               matlab.ui.control.Button
        CopyButton_3                    matlab.ui.control.Button
        Clear1Button_3                  matlab.ui.control.Button
        SaveInputBButton                matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
        Panel_8                         matlab.ui.container.Panel
        ExcelButton                     matlab.ui.control.Button
        Button_28                       matlab.ui.control.Button
        MonoStereoSignalsPanel_3        matlab.ui.container.Panel
        UndersampleOversampleButton_3   matlab.ui.control.Button
        MergeMonoStereoButton_3         matlab.ui.control.Button
        StereoMonoButton_3              matlab.ui.control.Button
        MonoStereoButton_3              matlab.ui.control.Button
        NumberofSamplesEditField        matlab.ui.control.NumericEditField
        NumberofSamplesEditFieldLabel   matlab.ui.control.Label
        RecordingParamatersPanel        matlab.ui.container.Panel
        ChannelsSwitch                  matlab.ui.control.Switch
        ChannelsSwitchLabel             matlab.ui.control.Label
        PreprocessingNormalizationPanel  matlab.ui.container.Panel
        StandardNormalizationButton_2   matlab.ui.control.StateButton
        Button_20                       matlab.ui.control.StateButton
        Button_19                       matlab.ui.control.StateButton
        Button_15                       matlab.ui.control.Button
        DAQButton_2                     matlab.ui.control.Button
        Button_13                       matlab.ui.control.Button
        AudioFileButton                 matlab.ui.control.Button
        Button_12                       matlab.ui.control.Button
        StarRecordingButton             matlab.ui.control.Button
        Panel_9                         matlab.ui.container.Panel
        FIRIIRFiltersPanel_2            matlab.ui.container.Panel
        Filters_2                       matlab.ui.control.DropDown
        Button_27                       matlab.ui.control.Button
        ConvolutionButton_3             matlab.ui.control.Button
        Button_25                       matlab.ui.control.Button
        PotentiationButton_2            matlab.ui.control.Button
        Button_24                       matlab.ui.control.Button
        DivisionButton                  matlab.ui.control.Button
        Button_23                       matlab.ui.control.Button
        ProductButton_2                 matlab.ui.control.Button
        Button_22                       matlab.ui.control.Button
        SubstractButton_2               matlab.ui.control.Button
        OperationsBetweenSignalsButton  matlab.ui.control.Button
        OperationsOnSignalsPanel_2      matlab.ui.container.Panel
        TimeShiftingButton              matlab.ui.control.Button
        TemporalReflectionButton_3      matlab.ui.control.Button
        Button_21                       matlab.ui.control.Button
        AdditionButton                  matlab.ui.control.Button
        Panel_4                         matlab.ui.container.Panel
        Button_29                       matlab.ui.control.Button
        DataSerieButton                 matlab.ui.control.Button
        Button_8                        matlab.ui.control.Button
        TriangularButton_2              matlab.ui.control.Button
        Button_7                        matlab.ui.control.Button
        StepButton                      matlab.ui.control.Button
        Button_6                        matlab.ui.control.Button
        RampButton_2                    matlab.ui.control.Button
        Button_5                        matlab.ui.control.Button
        SawtoothButton_2                matlab.ui.control.Button
        Button_4                        matlab.ui.control.Button
        CosButton_2                     matlab.ui.control.Button
        Button_3                        matlab.ui.control.Button
        sinButton                       matlab.ui.control.Button
        Button_2                        matlab.ui.control.Button
        Button                          matlab.ui.control.Button
        SincButton_2                    matlab.ui.control.Button
        ChirpButton                     matlab.ui.control.Button
        Panel_7                         matlab.ui.container.Panel
        Transfer1Button                 matlab.ui.control.Button
        CopyButton                      matlab.ui.control.Button
        Clear1Button                    matlab.ui.control.Button
        SaveInputAButton                matlab.ui.control.Button
        UIAxes_2                        matlab.ui.control.UIAxes
        Panel_2                         matlab.ui.container.Panel
        ShowOutputPanel                 matlab.ui.container.Panel
        ButtonGroup_5                   matlab.ui.container.ButtonGroup
        OUTPUTCheckBox                  matlab.ui.control.CheckBox
        DataInputPanel                  matlab.ui.container.Panel
        ButtonGroup_6                   matlab.ui.container.ButtonGroup
        INPUTBButton                    matlab.ui.control.RadioButton
        INPUTAButton                    matlab.ui.control.RadioButton
        AudioPlaybackPanel_2            matlab.ui.container.Panel
        Image5_6                        matlab.ui.control.Image
        Image5_5                        matlab.ui.control.Image
        Image5_4                        matlab.ui.control.Image
        Image5_3                        matlab.ui.control.Image
        EchoButton_2                    matlab.ui.control.StateButton
        ReverseButton_2                 matlab.ui.control.StateButton
        StopButton_2                    matlab.ui.control.Button
        PlayButton_3                    matlab.ui.control.Button
        VolumeControlPanel              matlab.ui.container.Panel
        ButtonGroup_2                   matlab.ui.container.ButtonGroup
        Image5_2                        matlab.ui.control.Image
        Image5                          matlab.ui.control.Image
        ConstantVolumeButton            matlab.ui.control.RadioButton
        DecreasingVolumeButton_2        matlab.ui.control.RadioButton
        IncreasingVolumeButton_2        matlab.ui.control.RadioButton
        ConstantVolumeSlider            matlab.ui.control.Slider
        MenuButton                      matlab.ui.control.Button
        GraphsPanel                     matlab.ui.container.Panel
        ButtonGroup_7                   matlab.ui.container.ButtonGroup
        AlltheGrahsButton               matlab.ui.control.RadioButton
        GraphOutputButton               matlab.ui.control.RadioButton
        GraphInputBButton               matlab.ui.control.RadioButton
        GraphInputAButton               matlab.ui.control.RadioButton
        GraphsButton                    matlab.ui.control.Button
        DropDown_6                      matlab.ui.control.DropDown
        AudioButton                     matlab.ui.control.Button
        DataButton                      matlab.ui.control.Button
        DropDown_5                      matlab.ui.control.DropDown
        DropDown_4                      matlab.ui.control.DropDown
        DropDown_2                      matlab.ui.control.DropDown
        UniversidadDelValleLabel        matlab.ui.control.Label
        Image2                          matlab.ui.control.Image
        PresentProjectButton            matlab.ui.control.Button
        Panel_3                         matlab.ui.container.Panel
        JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel  matlab.ui.control.Label
        DIGITALSIGNALSPROCESSINGPROJECT1Label  matlab.ui.control.Label
        Image4                          matlab.ui.control.Image
    end

    
    properties (Access = private)
       volumen;
       input_A;input_B;output;%Señales de entrada
       timeA;timeB;timeout;%vectores de tiempo
       freqA=0;freqB=0;freqout=0;%frecuencias de muestreo
       fA=0;fB=0;fout=0;
       h;ak_block;bk_block;h_filter;
       ak;bk;control;salidaconv;an;bn;salidaedf;
        
    end
    
    methods (Access = private)

       
        function f = calculateDominantFrequency(~,signal, fs)
            % Calcula la frecuencia dominante de una señal utilizando FFT
            n = length(signal);  % Número de muestras
            Y = fft(signal);  % FFT de la señal
            P2 = abs(Y/n);  % Magnitud de la FFT
            P1 = P2(1:floor(n/2)+1);  % Solo considerar la mitad positiva del espectro
            P1(2:end-1) = 2*P1(2:end-1);  % Ajustar la amplitud

            % Crear vector de frecuencias
            fVec = fs*(0:(n/2))/n;

            % Encontrar la frecuencia dominante (máxima magnitud)
            [~, idx] = max(P1);
            f = fVec(idx);  % La frecuencia dominante en Hz
        end


        function Y=ecuaciondif(~,x,dato1,dato2)   %%función ecuación  diferencia
                   ak1=dato1;bk1=dato2;inp1=x;
                   inp_f=zeros(1,(length(inp1)));
                   inp_r=zeros(1,(length(inp1)));
                   temp1=zeros(1,(length(bk1)));
                   temp2=zeros(1,(length(bk1)-1));
                   lx=(length(inp1)-1);    
                   lbk=(length(bk1)-1);
                   % Respuesta en estado vector (Forzada)
                   for j=0:lx
                       for k=0:lbk
                           ind1= j-k;
                           if ind1<0
                               xn=0;
                           else
                               xn=inp1(ind1+1);
                           end
                           temp1(k+1)=bk1(k+1)*xn;
                       end   
                       y0=sum(temp1);
                       inp_f(j+1)=y0;
                       lak=(length(ak1)-1);
                       for m=1 : lak  %Respuesta entrada cero (natural)
                           ind2= j-m;
                           if ind2<0
                               xn2=0;
                           else
                               xn2=inp_r(ind2+1);
                           end        
                           temp2(m)=ak1(m+1)*xn2;
                       end
                       y1=inp_f(j+1)-sum(temp2);
                       inp_r(j+1)=y1;
                   end
                   Y=inp_r';
        end
       function [y, ni_y,nf_y]=Convolucion(~,x, ni_x, h, ni_h )

            % Cálculo de la convolución : y(n) = suma_k [x(k) h(n-k)]
            lx=length(x); % longitud de x(n)
            lh=length(h); % longitud de h(n)  
            ly=lh+lx-1;   % longitud de y(n)

            nf_x=lx+ni_x-1;      %instante final de x(n)
            nf_h=lh+ni_h-1;      %instante final de h(n)

            % Datos de salida y(n)
            ni_y=ni_x+ ni_h;  % Instante inicio de y(n)
            nf_y=nf_x+ nf_h;  % Instante final de y(n)

            y=zeros(1,ly);
            i=0;
            for  n=ni_y:nf_y
              i=i+1;  
              temp=0; 
              for k= ni_x: nf_x
                kx=k-ni_x+1; %ajuste del indice de x(n) para evitar salir rango de Matlab
                ind= n-k;
                if  ind >= ni_h && ind <= nf_h
                  n_k= ind-ni_h+1;  %ajuste del indice de h(n-k) para evitar salir rango de Matlab
                  temp=temp+x(kx)*h(n_k);
               end
              end
          y(i)=temp;
            end
        end 
        %% new
        function [y, t_y] = cambiarFrecuenciaMuestreo(~,x, t, k, tipo)
    % x: vector de datos discretos
    % t: vector de instantes correspondientes a los datos
    % k: factor de submuestreo o sobremuestreo (entero)
    % tipo: 'submuestreo' o 'sobremuestreo'

     if strcmp(tipo,'1')
        % Submuestreo
        y = x(1:k:end); % Tomar cada k-ésimo valor
        t_y = t(1:k:end); % Vector de instantes correspondientes

    elseif strcmp(tipo,'2')
        % Sobremuestreo
        % Generar nuevos instantes
        t_y = linspace(t(1), t(end), length(x) * k); 
        % Interpolación para estimar nuevos valores
        y = interp1(t, x, t_y, 'linear');
    else
        error('Tipo no válido. Usa "1" o "2".');
    end
end

        %%% nueva funcion

        function [xf1, xf2, tf1, tf2, ffs1, ffs2] = resampleSignals(~, x1, x2, t1, t2, fs1, fs2)
            % Verificar si las frecuencias son iguales
            if fs1 == fs2
                xf1 = x1;
                xf2 = x2;
                tf1 = t1;
                tf2 = t2;
                ffs1 = fs1;
                ffs2 = fs2;
            else
                % Resamplear ambas señales a la frecuencia más alta
                if fs1 > fs2
                    xf1 = x1;
                    xf2 = resample(x2, fs1, fs2);  % Resamplear x2 a fs1
                    tf1 = t1;
                    tf2 = (0:length(xf2)-1) / fs1;  % Ajustar tiempo de x2 al nuevo número de muestras
                    ffs1 = fs1;
                    ffs2 = fs1;
                else
                    xf1 = resample(x1, fs2, fs1);  % Resamplear x1 a fs2
                    xf2 = x2;
                    tf1 = (0:length(xf1)-1) / fs2;  % Ajustar tiempo de x1 al nuevo número de muestras
                    tf2 = t2;
                    ffs1 = fs2;
                    ffs2 = fs2;
                end
            end

            % Asegurarse de que ambas señales tienen el mismo número de canales
            if size(xf1, 1) ~= size(xf2, 1)
                if size(xf1, 1) == 1
                    xf1 = [xf1; xf1];  % Convertir mono a estéreo duplicando el canal
                elseif size(xf2, 1) == 1
                    xf2 = [xf2; xf2];  % Convertir mono a estéreo duplicando el canal
                end
            end

            % Asegurarse de que las longitudes de las señales coinciden
            min_len = min(length(xf1), length(xf2));
            xf1 = xf1(1:min_len);
            xf2 = xf2(1:min_len);
            tf1 = tf1(1:min_len);  % Ajustar tiempo de xf1
            tf2 = tf2(1:min_len);  % Ajustar tiempo de xf2
        end


        function [x1f, x2f, t1f, t2f] = pad(~, x1, x2, t1, t2, fs)
            % Asegurarse de que las señales tengan el mismo número de muestras
            len_x1 = size(x1, 2);
            len_x2 = size(x2, 2);

            if len_x1 > len_x2
                x2f = [x2, zeros(size(x2, 1), len_x1 - len_x2)];  % Pad de x2
                x1f = x1;
                t2f = [t2, t2(end) + (1/fs)*(1:(len_x1 - len_x2))];  % Ajustar el tiempo de t2
            elseif len_x2 > len_x1
                x1f = [x1, zeros(size(x1, 1), len_x2 - len_x1)];  % Pad de x1
                x2f = x2;
                t1f = [t1, t1(end) + (1/fs)*(1:(len_x2 - len_x1))];  % Ajustar el tiempo de t1
            else
                x1f = x1;
                x2f = x2;
                t1f = t1;
                t2f = t2;
            end

            % Ajustar el tiempo
            t1f = t1f(1:size(x1f, 2));  % Asegurarse de que t1f tiene la misma longitud que x1f
            t2f = t2f(1:size(x2f, 2));  % Asegurarse de que t2f tiene la misma longitud que x2f
        end

        %%learning funcion

           %funtion of Plot Graphs
           function plot_Graphs(app,axe,x1,t1,fs)
          
            if app.DropDown_6.Value == "Time"
               %% plot o stem
                plot(axe,t1,x1);

            elseif app.DropDown_6.Value == "Frequency"
                fft_signal = fft(x1);
                N = length(x1);
                magnitude = abs(fft_signal/N);
                magnitude = magnitude(1:floor(N/2)+1);
                if N > 1
                magnitude(2:end) = magnitude(2:end)*2;
                end
                f = (0:N-1)*(fs/N);
                f = f(1:floor(N/2)+1);
                plot(axe,f,magnitude);  
            elseif app.DropDown_6.Value  == "Spectogram"
            windowLength = 256;    
            noverlap = 128;        
            nfft = 512;  
            if size(x1,1)==2
                x1 = (x1(1,:)+x1(2,:))/2;
                [~, F, T, P] = spectrogram(x1, windowLength, noverlap, [], fs, 'yaxis');
                ylim(app.UIAxes_2, [0 300]);
                surf(axe,T, F, 10*log10(abs(P)), 'edgecolor', 'none');view(axe,0, 90);colormap(axe,'jet');
            else
                [~, F, T, P] = spectrogram(x1, windowLength, noverlap, [], fs, 'yaxis');
                ylim(app.UIAxes_2, [0 300]);
                surf(axe,T, F, 10*log10(P), 'edgecolor', 'none');view(axe,0, 90);colormap(axe,'jet');
            end
            end
            end
    end
% function Input1Button_2Pushed(app, event) % es un boton 
% 
% window_length = 256; % Longitud de la ventana (en muestras)
% overlap = 128; % Solapamiento entre ventanas (en muestras)
% 
% % Calcular el espectrograma
% [s, f, t_spectrogram] = spectrogram(app.senal, window_length, overlap, [], app.fs);
% 
% spectrogram_dB = 10*log10(abs(s)); % Convertir la magnitud del espectrograma a dB (escala logarítmica)
% 
% imagesc(app.spectogram, t_spectrogram, f, spectrogram_dB);
% axis(app.spectogram, 'xy'); % Invertir el eje

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            % muestra de paneles completo
            app.Panel_3.Visible = 'on';
            app.Panel_4.Visible = 'off'; 
            app.Panel_7.Visible = 'off'; 
            app.Panel_8.Visible = 'off'; 
            app.Panel_9.Visible = 'off'; 
            app.VolumeControlPanel.Visible = 'off'; 
            app.AudioPlaybackPanel_2.Visible ="off";
            app.DataInputPanel.Visible ="off";
            app.ShowOutputPanel.Visible="off";
            app.GraphsPanel.Visible = 'off';
            %grafiaca 1 input a
            app.Panel_7.Visible ="off";
            %grafica 2 input b
            app.Panel_11.Visible ="off";
            %grafica 3 output
            app.Panel_12.Visible ="off";
            %todas las graficas
            app.Panel_13.Visible ="off";
        end

        % Button pushed function: PresentProjectButton
        function PresentProjectButtonPushed(app, event)

            % muestra de paneles completo
            app.Panel_3.Visible = 'on';
            app.Panel_4.Visible = 'off'; 
            app.Panel_7.Visible = 'off'; 
            app.Panel_8.Visible = 'off'; 
            app.Panel_9.Visible = 'off'; 
           
        end

        % Button pushed function: MenuButton
        function MenuButtonPushed(app, event)
             
              if strcmp(app.Panel_4.Visible, 'on')  % Si el panel está visible
                  
                  app.Panel_3.Visible = 'off';
                  app.Panel_4.Visible = 'off'; 
                  app.Panel_7.Visible = 'off'; 
                  app.Panel_8.Visible = 'off'; 
                  app.Panel_9.Visible = 'off'; 
                 
            else                                % Si el panel está oculto
                  app.Panel_3.Visible = 'off';
                  app.Panel_4.Visible = 'off'; 
                  app.Panel_7.Visible = 'on'; 
                  app.Panel_8.Visible = 'on'; 
                  app.Panel_9.Visible = 'off';
             end
          
        end

        % Value changed function: DropDown_4
        function DropDown_4ValueChanged(app, event)
            value = app.DropDown_4.Value;
             if strcmp(value,'Playback')
                   app.AudioPlaybackPanel_2.Visible = 'on'; 
                   app.VolumeControlPanel.Visible = 'off'; 
             end
             if strcmp(value,'Volume')
                    app.AudioPlaybackPanel_2.Visible = 'off'; 
                   app.VolumeControlPanel.Visible = 'on'; 
             end
             
        end

        % Button pushed function: AudioButton
        function AudioButtonPushed(app, event)
 
            if strcmp(app.AudioPlaybackPanel_2.Visible, 'on')  % Si el panel está visible
                    app.AudioPlaybackPanel_2.Visible = 'off'; 
                    app.VolumeControlPanel.Visible = 'off'; 
            else                                % Si el panel está oculto
                    app.AudioPlaybackPanel_2.Visible = 'on'; 
             end
        end

        % Value changed function: DropDown_2
        function DropDown_2ValueChanged(app, event)
            value = app.DropDown_2.Value;
         % Si se selecciona "Input/Output", muestra los paneles
    if strcmp(value,'Input')
        app.DataInputPanel.Visible = 'on';
        app.ShowOutputPanel.Visible = 'on';
        
    % Si se selecciona "Seleccionar", oculta los paneles
    elseif strcmp(value,'Output')
        app.DataInputPanel.Visible = 'off';
        app.ShowOutputPanel.Visible = 'off';
    end
        end

        % Button pushed function: DataButton
        function DataButtonPushed(app, event)
        
            if strcmp(app.DataInputPanel.Visible, 'on')  % Si el panel está visible
                   app.DataInputPanel.Visible = 'off';
                   app.ShowOutputPanel.Visible = 'off';
                   app.AudioPlaybackPanel_2.Visible ="off";

            else                                % Si el panel está oculto
                   app.DataInputPanel.Visible = 'on';
                   app.ShowOutputPanel.Visible = 'on';
             end
        end

        % Value changed function: DropDown_6
        function DropDown_6ValueChanged(app, event)
            %value = app.DropDown_6.Value;
           % if strcmp(value,'Time')
            %app.GraphsPanel.Visible = 'on';
            %end
            %if strcmp(value,'Frequency')
            %app.GraphsPanel.Visible = 'on';
            %end
            %if strcmp(value,'Spectogram')
            %app.GraphsPanel.Visible = 'on';
            %end
        end

        % Button pushed function: GraphsButton
        function GraphsButtonPushed(app, event)
           
              if strcmp(app.GraphsPanel.Visible, 'on')  % Si el panel está visible
                   app.GraphsPanel.Visible = 'off';
            else                                % Si el panel está oculto
                 app.GraphsPanel.Visible = 'on';
             end
        end

        % Value changed function: DropDown_5
        function DropDown_5ValueChanged(app, event)
            value = app.DropDown_5.Value;
              if strcmp(value,'Signal Inputs')
                   app.Panel_4.Visible = 'off'; 
                   app.Panel_8.Visible = 'on'; 
                   app.Panel_9.Visible = 'off'; 
             end
             if strcmp(value,'Synthetic Signals')
                   app.Panel_4.Visible = 'on'; 
                   app.Panel_8.Visible = 'off'; 
                   app.Panel_9.Visible = 'off'; 
             end
             if strcmp(value,'Processing Functions')
                   app.Panel_4.Visible = 'off'; 
                   app.Panel_8.Visible = 'off'; 
                   app.Panel_9.Visible = 'on'; 
             end
             if strcmp(value,'Others')
                   app.Panel_4.Visible = 'off'; 
                   app.Panel_8.Visible = 'off'; 
                   app.Panel_9.Visible = 'off'; 
             end
        end

        % Selection changed function: ButtonGroup_7
        function ButtonGroup_7SelectionChanged(app, event)
            selectedButton = app.ButtonGroup_7.SelectedObject;
            if strcmp(selectedButton.Text, 'Graph Input A')
                %grafiaca 1 input a
                app.UIAxes_2.Visible ="on";
                app.Panel_7.Visible ="on";
                %grafica 2 input b
                app.UIAxes.Visible ="off";
                app.Panel_11.Visible ="off";
                %grafica 3 output
                app.UIAxes_3.Visible ="off";
                app.Panel_12.Visible ="off";
                 %Todas las graficas 
                app.UIAxes_4.Visible ="off";
                app.UIAxes_5.Visible ="off";
                app.UIAxes_6.Visible ="off";
                app.Panel_13.Visible ="off";
            end

            if strcmp(selectedButton.Text, 'Graph Input B')
                %grafiaca 1 input a
                app.UIAxes_2.Visible ="off";
                app.Panel_7.Visible ="off";
                %grafica 2 input b
                app.UIAxes.Visible ="on";
                app.Panel_11.Visible ="on";
                %grafica 3 output
                app.UIAxes_3.Visible ="off";
                app.Panel_12.Visible ="off";
                 %Todas las graficas 
                app.UIAxes_4.Visible ="off";
                app.UIAxes_5.Visible ="off";
                app.UIAxes_6.Visible ="off";
                app.Panel_13.Visible ="off";
            end

            if strcmp(selectedButton.Text, 'Graph Output ')
               %grafiaca 1 input a
                app.UIAxes_2.Visible ="off";
                app.Panel_7.Visible ="off";
                %grafica 2 input b
                app.UIAxes.Visible ="off";
                app.Panel_11.Visible ="off";
                %grafica 3 output
                app.UIAxes_3.Visible ="on";
                app.Panel_12.Visible ="on";
                %Todas las graficas 
                app.UIAxes_4.Visible ="off";
                app.UIAxes_5.Visible ="off";
                app.UIAxes_6.Visible ="off";
                app.Panel_13.Visible ="off";
            end
            
            if strcmp(selectedButton.Text, 'All the Grahs')
               %grafiaca 1 input a
                app.UIAxes_2.Visible ="off";
                app.Panel_7.Visible ="off";
                %grafica 2 input b
                app.UIAxes.Visible ="off";
                app.Panel_11.Visible ="off";
                %grafica 3 output
                app.UIAxes_3.Visible ="off";
                app.Panel_12.Visible ="off";
                %Todas las graficas 
                app.UIAxes_4.Visible ="on";
                app.UIAxes_5.Visible ="on";
                app.UIAxes_6.Visible ="on";
                app.Panel_13.Visible ="on";
            end
        end

        % Button pushed function: StarRecordingButton
        function StarRecordingButtonPushed(app, event)

            % Crear una instancia de GetAudioInputs
            audioRecorder = GetAudioInputs();

            % Obtener el valor del Switch de la aplicación
            switchValue = app.ChannelsSwitch.Value;

            Fs = app.NumberofSamplesEditField.Value;

            % Llamar al método GrabarAudio de GetAudioInputs
             a = audioRecorder.GrabarAudio(switchValue,Fs);
             if a == 0
                 return;
             end
                
            if  app.INPUTAButton.Value == true
               
                app.timeA = audioRecorder.time;

                app.input_A= audioRecorder.audioData;

                app.freqA = audioRecorder.fs;

                app.fA = audioRecorder.f;
                
                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))
                    % Mostrar mensaje de error
                    msgbox('La señal contiene datos no válidos (NaN o Inf).', 'Error', 'error');

                    % Preguntar al usuario si desea corregir los datos
                    choice = questdlg('¿Desea interpolar los datos no válidos?', ...
                        'Datos no válidos', 'Sí', 'No', 'Sí');

                    % Si elige "Sí", interpolar los datos
                    if strcmp(choice, 'Sí')
                        % Encontrar índices con NaN o Inf
                        invalidIdx = isnan(app.input_A) | isinf(app.input_A);
                        validIdx = ~invalidIdx;

                        % Interpolación lineal para los valores no válidos
                        app.input_A(invalidIdx) = interp1(app.timeA(validIdx), app.input_A(validIdx), app.timeA(invalidIdx), 'linear');
                    else
                        return; % Si elige "No", detener la ejecución
                    end
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
                
                
            else
               
                app.timeB = audioRecorder.time;

                app.input_B = audioRecorder.audioData;
                
                app.freqB = audioRecorder.fs;

                app.fB = audioRecorder.f;

                 % Verificar si hay NaN o Inf
                 if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))
                     msgbox('La señal contiene datos no válidos (NaN o Inf).', 'Error', 'error');

                     % Preguntar al usuario si desea corregir los datos
                     choice = questdlg('¿Desea interpolar los datos no válidos?', ...
                         'Datos no válidos', 'Sí', 'No', 'Sí');

                     if strcmp(choice, 'Sí')
                         invalidIdx = isnan(app.input_B) | isinf(app.input_B);
                         validIdx = ~invalidIdx;
                         app.input_B(invalidIdx) = interp1(app.timeB(validIdx), app.input_B(validIdx), app.timeB(invalidIdx), 'linear');
                     else
                         return;
                     end
                 end
                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end

        end

        % Button pushed function: AudioFileButton
        function AudioFileButtonPushed(app, event)
            
            % Crear una instancia de GetAudioInputs
            audioRecorder = GetAudioInputs();

            % Obtener el valor del Switch de la aplicación
            switchValue = app.ChannelsSwitch.Value;

            % Llamar al método AbrirArchivo de GetAudioInputs
            a = audioRecorder.AbrirArchivo(switchValue);
            
            if a == 0
                 return;
            end
            
            if  app.INPUTAButton.Value == true
               
                app.timeA = audioRecorder.time;

                app.input_A= audioRecorder.audioData;

                app.freqA = audioRecorder.fs;
                
                app.fA = audioRecorder.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end


                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
            else
                
                app.timeB = audioRecorder.time;

                app.input_B = audioRecorder.audioData;
                
                app.freqB = audioRecorder.fs;

                app.fB = audioRecorder.f;

                 if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end

        end

        % Button pushed function: Clear1Button_6
        function Clear1Button_6Pushed(app, event)
            app.input_A = [];
            app.timeA = [];
            app.freqA = 0;
            app.fA = 0;

            cla(app.UIAxes_4);
        end

        % Button pushed function: Clear1Button_5
        function Clear1Button_5Pushed(app, event)
            app.input_B = [];
            app.timeB = [];
            app.freqB = 0;
            app.fB = 0;

            cla(app.UIAxes_5);
        end

        % Button pushed function: Clear1Button_4
        function Clear1Button_4Pushed(app, event)
            app.output = [];
            app.timeout = [];
            app.freqout = 0;
            app.fout = 0;

            cla(app.UIAxes_6);
        end

        % Button pushed function: Clear1Button_2
        function Clear1Button_2Pushed(app, event)
            app.output = [];
            app.timeout = [];
            app.freqout = 0;
            app.fout = 0;

            cla(app.UIAxes_3);
        end

        % Button pushed function: Clear1Button_3
        function Clear1Button_3Pushed(app, event)
            app.input_B = [];
            app.timeB = [];
            app.freqB = 0;
            app.fB = 0;

            cla(app.UIAxes);
        end

        % Button pushed function: Clear1Button
        function Clear1ButtonPushed(app, event)
            app.input_A = [];
            app.timeA = [];
            app.freqA = 0;
            app.fA = 0;

            cla(app.UIAxes_2);
        end

        % Button pushed function: SaveInputAButton_2
        function SaveInputAButton_2Pushed(app, event)
            dlg_title = 'Record Input A';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.input_A));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.input_A), app.freqA,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.');
            msgbox('File saved successfully.');
        end

        % Button pushed function: SaveInputBButton_3
        function SaveInputBButton_3Pushed(app, event)
            dlg_title = 'Record Input B';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.input_B));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.input_B), app.freqB,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.')
        end

        % Button pushed function: SaveOutputButton_2
        function SaveOutputButton_2Pushed(app, event)
            dlg_title = 'Record Output';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.output));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.output), app.freqout,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.');
        end

        % Button pushed function: SaveOutputButton
        function SaveOutputButtonPushed(app, event)
            dlg_title = 'Record Output';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.output));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.output), app.freqout,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.');
        end

        % Button pushed function: SaveInputBButton
        function SaveInputBButtonPushed(app, event)
             dlg_title = 'Record Input B';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.input_B));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.input_B), app.freqB,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.')
        end

        % Button pushed function: SaveInputAButton
        function SaveInputAButtonPushed(app, event)
            dlg_title = 'Record Input A';
            prompt = {'Audio Name:','Bits per sample'};
            dims = [1 25];
            defaultans = {'' '16'};
            answer = inputdlg(prompt,dlg_title,dims,defaultans);
            disp(size(app.input_A));
            audiowrite(strcat(answer{1}, '.wav'), transpose(app.input_A), app.freqA,'BitsPerSample',str2double(answer{2})); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
            msgbox('File saved sucessfully.');
        end

        % Button pushed function: Transfer1Button_6
        function Transfer1Button_6Pushed(app, event)
            

        end

        % Button pushed function: PlayButton_3
        function PlayButton_3Pushed(app, event)
            clear sound;  % Limpiar cualquier sonido en reproducción previa

% Elige la señal y frecuencia de muestreo correctas
if app.OUTPUTCheckBox.Value == true
    x = app.output;
    fs = app.freqout;
else
    if app.INPUTAButton.Value == true
        x = app.input_A;
        fs = app.freqA;
    elseif app.INPUTBButton.Value == true
        x = app.input_B;
        fs = app.freqB;
    end
end

% Validar que la tasa de muestreo es un número válido
if isempty(fs) || ~isnumeric(fs) || fs <= 0
    msgbox('Error: Invalid sample rate.', 'Error', 'Error');
    return;
end

% Definir un rango típico de sample rate (8000 Hz a 192000 Hz)
fs_min = 8000;
fs_max = 192000;

% Si la tasa de muestreo está fuera de este rango, advertir al usuario
if fs < fs_min || fs > fs_max
    msgbox('Error: Sample rate out of range (8000 Hz to 192000 Hz).', 'Error', 'Error');
    
end

% Definir una tasa de muestreo deseada (por ejemplo, 44100 Hz)
fs_deseado = 48000;

% Verificar si la tasa de muestreo actual no es la deseada
if fs ~= fs_deseado
    % Advertir al usuario y preguntar si desea remuestrear la señal
    choice = questdlg(...
        ['La señal tiene una tasa de muestreo de ', num2str(fs), ' Hz. ¿Deseas remuestrear la señal a ', num2str(fs_deseado), ' Hz para reproducirla?'], ...
        'Advertencia: Remuestrear la señal', ...
        'Sí', 'No', 'No');
    
    % Tomar acción según la elección del usuario
    switch choice
       case 'Sí'
    % Remuestrear la señal
    x = resample(x, fs_deseado, fs);  
    fs = fs_deseado;  % Actualizar la tasa de muestreo a la nueva
    disp(['Se ha remuestreado la señal a ', num2str(fs_deseado), ' Hz.']);

    % Recalcular el tiempo de la señal remuestreada
    new_time = (0:length(x)-1) / fs;  % Crear nuevo vector de tiempo para la señal remuestreada

    if app.OUTPUTCheckBox.Value == true
        app.output = x;
        app.freqout = fs;
        app.timeout = new_time;  % Actualizar el tiempo de salida con el nuevo vector de tiempo
        
        % Graficar la señal remuestreada en los ejes correspondientes
        app.plot_Graphs(app.UIAxes_3, app.output, app.timeout, app.freqout);
        app.plot_Graphs(app.UIAxes_6, app.output, app.timeout, app.freqout);

    else
        if app.INPUTAButton.Value == true
            app.input_A = x;
            app.freqA = fs;
            app.timeA = new_time;  % Actualizar el tiempo de la señal A
            
            % Graficar la señal remuestreada en los ejes correspondientes
            app.plot_Graphs(app.UIAxes_2, app.input_A, app.timeA, app.freqA);
            app.plot_Graphs(app.UIAxes_4, app.input_A, app.timeA, app.freqA);

        elseif app.INPUTBButton.Value == true
            app.input_B = x;
            app.freqB = fs;
            app.timeB = new_time;  % Actualizar el tiempo de la señal B
            
            % Graficar la señal remuestreada en los ejes correspondientes
            app.plot_Graphs(app.UIAxes, app.input_B, app.timeB, app.freqB);
            app.plot_Graphs(app.UIAxes_5, app.input_B, app.timeB, app.freqB);
        end
    end

           
           

            

        case 'No'
            % No hacer nada y continuar con la tasa de muestreo actual
            disp('No se ha remuestreado la señal.');
        otherwise
            % Si el usuario cierra la ventana, cancelar la reproducción
            return;
    end
end

% Verificar que la señal no esté vacía
if isempty(x)
    msgbox('Error: No signal to play.', 'Error', 'Error');
    return;
end

% Aplicar ajustes de volumen
signalsize = size(x);
logVolume = 10^((app.ConstantVolumeSlider.Value / 99 * 2) - 2);

% Inicializa el vector de volumen basado en la longitud de la señal
app.volumen = ones(size(x, 2), 1);

% Ajusta el vector de volumen según la opción seleccionada
if app.IncreasingVolumeButton_2.Value == 1
    app.volumen = linspace(0.01, logVolume, size(x, 2));
elseif app.DecreasingVolumeButton_2.Value == 1
    app.volumen = linspace(logVolume, 0.01, size(x, 2));
else
    app.volumen = linspace(logVolume, logVolume, size(x, 2));
end

% Asegurarse de que volumen sea una matriz de la misma forma que x para la multiplicación
if size(x, 1) == 2  % Para audio estéreo
    app.volumen = [app.volumen; app.volumen];
end

% Aplica el ajuste de volumen
x = x .* app.volumen;

% Validar si la señal debe reproducirse en reversa
if app.ReverseButton_2.Value == 1
    x = fliplr(x);
    disp('Reverse');
end

% Aplicar eco si está activado
if app.EchoButton_2.Value == 1
    dlg_title = 'Echo Amplitude';
    prompt = {'A1:', 'A2:'};
    num_lines = 1;
    defaultans = {'0.5', '0.5'};
    answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

    if isempty(answer)
        return;
    end

    A1 = str2double(answer{1});
    A2 = str2double(answer{2});

    if A1 > 1 || A2 > 1
        msgbox('Select a correct amplitude', 'Error', 'Warning');
        return;
    end

    y2 = zeros(1, signalsize(2));
    y3 = zeros(1, signalsize(2));

    for i = 1:signalsize(2) - floor(0.1 * fs)
        y2(i + floor(0.1 * fs)) = x(i);
    end

    for i = 1:signalsize(2) - floor(0.05 * fs)
        y3(i + floor(0.05 * fs)) = x(i);
    end

    x = x + A1 * y2 + A2 * y3;
end

% Reproducir el sonido ajustado
try
    sound(x, fs);  % Intentar reproducir sonido
catch ME
    msgbox(['Error playing sound: ' ME.message], 'Error', 'Error');
end

        end

        % Button pushed function: StopButton_2
        function StopButton_2Pushed(app, event)
            clear sound;
        end

        % Button pushed function: ChirpButton
function ChirpButtonPushed(app, event)

% Crear una instancia de GetSyntheticSignals
synthetic = GetSyntheticSignals();

% Llamar al método GenChirp de GetSyntheticSignals
a = synthetic.GenChirp();


if a == 0
    return;
end


if  app.INPUTAButton.Value == true

    app.timeA = synthetic.time;

    app.input_A= synthetic.audioData;

    app.freqA = synthetic.fs;

    app.fA = synthetic.f;

    if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
        msgbox(['The signal contains invalid data (NaN or Inf). Please review the data entered'], 'Error', 'error');
    end

    app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
    app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

else

    app.timeB = synthetic.time;

    app.input_B = synthetic.audioData;

    app.freqB = synthetic.fs;

    app.fB = synthetic.f;

    if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
        msgbox(['The signal contains invalid data (NaN or Inf). Please review the data entered'], 'Error', 'error');
    end

    app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
    app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

end

end

        % Button pushed function: SincButton_2
        function SincButton_2Pushed(app, event)
            % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenSinc();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                 if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end
                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
           else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: sinButton
        function sinButtonPushed(app, event)
            % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenSin();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: CosButton_2
        function CosButton_2Pushed(app, event)
             % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenCos();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;
                
                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: SawtoothButton_2
        function SawtoothButton_2Pushed(app, event)
            % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenSaw();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: RampButton_2
        function RampButton_2Pushed(app, event)
             % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenRamp();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                 if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: StepButton
        function StepButtonPushed(app, event)
            % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenStep();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end
              
                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Button pushed function: TriangularButton_2
        function TriangularButton_2Pushed(app, event)
             % Crear una instancia de GetSyntheticSignals
            synthetic = GetSyntheticSignals();

            % Llamar al método GenChirp de GetSyntheticSignals
            a = synthetic.GenTrig();

            if a == 0
                return;
            end


            if  app.INPUTAButton.Value == true

                app.timeA = synthetic.time;

                app.input_A= synthetic.audioData;

                app.freqA = synthetic.fs;

                app.fA = synthetic.f;

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else

                app.timeB = synthetic.time;

                app.input_B = synthetic.audioData;

                app.freqB = synthetic.fs;

                app.fB = synthetic.f;

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))    %Validacion NaN o Inf
                    msgbox(['La señal contiene datos no válidos (NaN o Inf). Por favor, revise los datos ingresados.'], 'Error', 'error');
                end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);

            end
        end

        % Value changed function: Button_19
function Button_19ValueChanged(app, event) %#ok<*INUSD>
    value = app.Button_19.Value; %#ok<*NASGU>
    %normalizacion 0,1
    if  app.INPUTAButton.Value == true
        minVal = min(app.input_A);
        maxVal = max(app.input_A);
        app.input_A = (app.input_A - minVal) / (maxVal - minVal);


        app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
        app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
    end
    if app.INPUTBButton.Value == true
        minVal = min(app.input_B);
        maxVal = max(app.input_B);
        app.input_B = (app.input_B - minVal) / (maxVal - minVal);

        app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
        app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
    end
end

% Value changed function: Button_20
function Button_20ValueChanged(app, event)
    value = app.Button_20.Value;
    %normalizacion 1,-1
    if  app.INPUTAButton.Value == true
        maxAbsVal = max(abs(app.input_A));
        app.input_A = app.input_A / maxAbsVal;

        app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
        app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
    end
    if app.INPUTBButton.Value == true
        maxAbsVal = max(abs(app.input_B));
        app.input_B = app.input_B / maxAbsVal;

        app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
        app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
    end
end

% Value changed function: StandardNormalizationButton_2
function StandardNormalizationButton_2ValueChanged(app, event)
    value = app.StandardNormalizationButton_2.Value;
    %standard normalization
    if  app.INPUTAButton.Value == true
        meanSignal = mean(app.input_A);
        stdSignal = std(app.input_A);
        app.input_A= (app.input_A - meanSignal) / stdSignal;

        app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
        app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
    end
    if app.INPUTBButton.Value == true
        meanSignal = mean(app.input_B);
        stdSignal = std(app.input_B);
        app.input_B= (app.input_B - meanSignal) / stdSignal;

        app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
        app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
    end
end

        % Button pushed function: DAQButton_2
        function DAQButton_2Pushed(app, event)
           
    % Crear una instancia de GetAudioInputs
    audioRecorder = GetAudioInputs();

    % Llamar al método DAQ de GetAudioInputs
    try
        a = audioRecorder.DAQ();
    catch ME
        % Mostrar mensaje de error si falla la adquisición
        msgbox(['Error in data adquisition: ', ME.message], 'Error', 'error');
        return;
    end

    if a == 0
        msgbox('Could not establish communication with the Arduino COM port.', 'Error', 'error');
        return;
    end

    % Lógica para graficar y validar los datos capturados
    if app.INPUTAButton.Value == true
        app.timeA = audioRecorder.time;
        app.input_A = audioRecorder.audioData;
        app.freqA = audioRecorder.fs;
        app.fA = audioRecorder.f;
        if app.fA == 0
        app.fA = app.calculateDominantFrequency(app.input_A, app.freqA);  %Una función para calcular la frecuencia dominante
        end



        if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:))) %Validación de NaN o Inf
            msgbox(['Signal contains invalid data (NaN o Inf). Please, check input data.'], 'Error', 'error');
            return;
        end

        app.plot_Graphs(app.UIAxes_2, app.input_A, app.timeA, app.freqA);
        app.plot_Graphs(app.UIAxes_4, app.input_A, app.timeA, app.freqA);
    else
        app.timeB = audioRecorder.time;
        app.input_B = audioRecorder.audioData;
        app.freqB = audioRecorder.fs;
        app.fB = audioRecorder.f;

        if app.fB == 0
        app.fB = app.calculateDominantFrequency(app.input_B, app.freqB);
        end


        if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:))) %Validación de NaN o Inf
            msgbox(['Signal contains invalid data (NaN o Inf). Please, check input data.'], 'Error', 'error');
            return;
        end

        app.plot_Graphs(app.UIAxes, app.input_B, app.timeB, app.freqB);
        app.plot_Graphs(app.UIAxes_5, app.input_B, app.timeB, app.freqB);
    end


        end

        % Button pushed function: MonoStereoButton_3
        function MonoStereoButton_3Pushed(app, event)


            if  app.INPUTAButton.Value == true

                if app.freqA==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end
                if size(app.input_A,1)==2
                    msgbox('The signal is already stereo','Error')
                    return;
                end
                aux = [app.input_A;app.input_A];
                app.input_A = aux;
                app.plot_Graphs(app.UIAxes_2,app.input_A,app.time1,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
            else
                if app.freqB==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end
                if size(app.input_B,1)==2
                    msgbox('The signal is already stereo','Error')
                    return;
                end
                aux = [app.input_B;app.input_B];
                app.input_B = aux;
                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end
        end

        % Button pushed function: StereoMonoButton_3
        function StereoMonoButton_3Pushed(app, event)

            if  app.INPUTAButton.Value == true

                if app.freqA==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end
                if size(app.input_A,1)==1
                    msgbox('The signal is already mono','Error')
                    return;
                end
                aux = (app.input_A(1,:)+app.input_A(2,:))/2;
                app.input_A = aux;
                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
            else
                if app.freqA==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end
                if size(app.input_B,1)==1
                    msgbox('The signal is already mono','Error')
                    return;
                end
                aux = (app.input_B(1,:)+app.input_B(2,:))/2;
                app.input_B = aux;
                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end
        end

        % Button pushed function: MergeMonoStereoButton_3
        function MergeMonoStereoButton_3Pushed(app, event)
           
            if app.freqA==0 || app.freqB==0
               msgbox('Cannot perform operation, missing data.', 'Error');
               return;
           end
            
           if size(app.input_A,1)==2
               msgbox('Input 1 signal is stereo.', 'Error');
               return;
           end
           if size(app.input_B,1)==2
               msgbox('Input 2 signal is stereo.', 'Error');
               return;
           end

           [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB]=  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
           [app.input_B,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
           app.output = [app.input_A;app.input_B];
           %app.timeout = app.timeA;
           total_time_A = length(app.input_A) / app.freqA; % Tiempo total de input_A
           total_time_B = length(app.input_B) / app.freqB; % Tiempo total de input_B
           app.timeout = 0:1/app.freqA:(total_time_A + total_time_B); % Actualiza el timeout para cubrir ambos
           app.output = [app.input_A;app.input_B];
           app.timeout = 0:1/app.freqA:(length(app.output)-1)/app.freqA;
           app.freqout = app.freqA;
           app.fout = app.fA;
           app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
           app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
          
        end

        % Button pushed function: UndersampleOversampleButton_3
        function UndersampleOversampleButton_3Pushed(app, event)

    % Crear ventana de diálogo para pedir parámetros
    dlg_title = 'Parámetros de Muestreo';
    prompt = {'Factor de Submuestreo/Sobremuestreo (k):', 'Tipo ( 1-submuestreo/ 2-sobremuestreo):'};
    dims = [1 50];
    defaultans = {'2', '1'}; % Valores por defecto
    answer = inputdlg(prompt, dlg_title, dims, defaultans);

    if isempty(answer)
        msgbox('Operación cancelada. No se ingresaron parámetros.', 'Error', 'error');
        return;
    end

    % Obtener los parámetros del diálogo
    k = str2double(answer{1});
    tipo = answer{2};

    % Verificar si el valor de k es válido
    if isnan(k) || k <= 0 || mod(k,1) ~= 0
        msgbox('El valor de k debe ser un entero positivo.', 'Error', 'error');
        return;
    end

    % Verificar si el tipo es válido
    if ~strcmp(tipo, '1') && ~strcmp(tipo, '2')
        msgbox('El tipo debe ser "1-submuestreo" o "2-sobremuestreo".', 'Error', 'error');
        return;
    end

    % Llamar a la función cambiarFrecuenciaMuestreo
             if  app.INPUTAButton.Value == true
               [app.input_A, app.timeA] = app.cambiarFrecuenciaMuestreo(app.input_A, app.timeA, k, tipo);
               
                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
                
                
            else
               
                [app.input_B, app.timeB] = app.cambiarFrecuenciaMuestreo(app.input_B, app.timeB, k, tipo);

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end


    % Mostrar un mensaje de éxito
    msgbox('Proceso completado con éxito.');


       % if  app.INPUTAButton.Value == true
            % 
            %     if app.freqA==0
            %         msgbox('Cannot perform operation, missing data.', 'Error');
            %         return;
            %     end
            % 
            %     a = menu('Oversample/Undersample','Oversample','Undersample');
            %     fs = app.freqA;
            %     x1 = app.input_A;
            %     n1 = app.timeA;
            % 
            %     dlg_title = 'n/m resample';
            %     prompt = { 'n' 'm'};
            %     dims = [1 50];
            %     defaultans = {'10' '3'};
            % 
            %     answer = inputdlg(prompt,dlg_title,dims,defaultans);
            % 
            %     n = str2double(answer{1});
            %     m = str2double(answer{2});
            % 
            %     if n/m > 1 && a ~= 1
            %         warning('No se puede realizar la operación')
            %         return ;
            %     end
            % 
            %     if n/m <1 && a == 1
            %         warning('No se puede realizar la operación')
            %         return ;
            %     end
            %     if size(app.input_A,1)==2
            %         [x2(1,:),~] = resample(x1(1,:),n,m);
            %         [x2(2,:),~] = resample(x1(2,:),n,m);
            % 
            %         lx2 = (length(x1) - length(x2));
            %         temp = x2;
            %         if mod(lx2,2)==0
            %             x2 = [zeros(2,lx2/2) temp zeros(2,lx2/2)];
            %         else
            %             lx2 = lx2/2;
            %             x2 = [zeros(2,floor(lx2)) temp zeros(2,floor(lx2)+1)];
            %         end
            %     else
            % 
            %         [x2,~] = resample(x1,n,m);
            % 
            %         lx2 = (length(x1) - length(x2));
            %         temp = x2;
            %         if mod(lx2,2)==0
            %             x2 = [zeros(1,lx2/2) temp zeros(1,lx2/2)];
            %         else
            %             lx2 = lx2/2;
            %             x2 = [zeros(1,floor(lx2)) temp zeros(1,floor(lx2)+1)];
            %         end
            %     end
            %     n2 = 0:m/fs:length(n1)-1/fs;
            % 
            %     disp(n1(1));
            %     disp(n1(end));
            % 
            %     disp(n2(1));
            %     disp(n2(end));
            % 
            %     disp(size(x2));
            %     disp(size(n2));
            %     app.input_A = x2;
            %     app.timeA = n2;
            %     app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
            %     app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);
            % 
            % else
            %     if app.freqB==0
            %         msgbox('Cannot perform operation, missing data.', 'Error');
            %         return;
            %     end
            %     a = menu('Oversample/Undersample','Oversample','Undersample');
            %         fs = app.freqB;
            %         x1 = app.input_B;
            %         n1 = app.timeB;
            % 
            %         dlg_title = 'n/m resample';
            %         prompt = { 'n' 'm'};
            %         dims = [1 50];
            %         defaultans = {'3' '5'};
            % 
            %         answer = inputdlg(prompt,dlg_title,dims,defaultans);
            % 
            %         n = str2double(answer{1});
            %         m = str2double(answer{2});
            % 
            %         if n/m > 1 && a ~= 1
            %             warning('No se puede realizar la operación')
            %             return ;
            %         end
            % 
            %         if n/m <1 && a == 1
            %             warning('No se puede realizar la operación')
            %             return ;
            %         end
            % 
            %         if size(app.input_B,1)==2
            %             [x2(1,:),~] = resample(x1(1,:),n,m);
            %             [x2(2,:),~] = resample(x1(2,:),n,m);
            % 
            %             lx2 = (length(x1) - length(x2));
            %             temp = x2;
            %             if mod(lx2,2)==0
            %                 x2 = [zeros(2,lx2/2) temp zeros(2,lx2/2)];
            %             else
            %                 lx2 = lx2/2;
            %                 x2 = [zeros(2,floor(lx2)) temp zeros(2,floor(lx2)+1)];
            %             end
            %         else
            % 
            %             [x2,~] = resample(x1,n,m);
            % 
            %             lx2 = (length(x1) - length(x2));
            %             temp = x2;
            %             if mod(lx2,2)==0
            %                 x2 = [zeros(1,lx2/2) temp zeros(1,lx2/2)];
            %             else
            %                 lx2 = lx2/2;
            %                 x2 = [zeros(1,floor(lx2)) temp zeros(1,floor(lx2)+1)];
            %             end
            %         end
            %         n2 = 0:m/fs:length(n1)-1/fs;
            %     disp(size(x2));
            %     disp(size(n2));
            %     app.input_B = x2;
            %     app.timeB = n2;
            %     app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
            %     app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            % end
        end

        % Button pushed function: TimeShiftingButton
        function TimeShiftingButtonPushed(app, event)
                
            if  app.INPUTAButton.Value == true

                if app.freqA==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end

                dlg_title = 'Delay Value';
                prompt = { 'Delay(s):'};
                num_lines = 1;
                defaultans = {'5'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

                if isempty(answer)
                    return;
                end

                n0 = str2double(answer{1});

                if n0==0
                    msgbox('Wrong delay magnitude', 'Error');
                    return;
                end

                x = app.input_A;
                fs = app.freqA;
                f = app.fA;

                txi = app.timeA(1);
                txf = app.timeA(end);
                tx = txi+n0:1/fs:txf+n0;

                app.input_A = x;
                app.timeA = tx;
                app.freqA = fs;
                app.fA = f;

                app.plot_Graphs(app.UIAxes_2,x,tx,fs);
                app.plot_Graphs(app.UIAxes_4,x,tx,fs);
            else
                if app.freqB==0
                    msgbox('Cannot perform operation, missing data.', 'Error');
                    return;
                end


                dlg_title = 'Delay Value';
                prompt = { 'Delay(s):'};
                num_lines = 1;
                defaultans = {'5'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

                if isempty(answer)
                    return;
                end

                n0 = str2double(answer{1});

                if n0==0
                    msgbox('Wrong delay magnitude', 'Error');
                    return;
                end

                x = app.input_B;
                fs = app.freqB;
                f = app.fB;

                txi = app.timeB(1);
                txf = app.timeB(end);
                tx = txi+n0:1/fs:txf+n0;
                app.input_B = x;
                app.timeB = tx;
                app.freqB = fs;
                app.fB = f;

                app.plot_Graphs(app.UIAxes,x,tx,fs);
                app.plot_Graphs(app.UIAxes_5,x,tx,fs);
            end
        end

        % Button pushed function: TemporalReflectionButton_3
        function TemporalReflectionButton_3Pushed(app, event)
             if  app.INPUTAButton.Value == true
                if app.freqA==0
                    msgbox('Cannot perform operation, no data.', 'Error');
                    return;
                end

                x = app.input_A;
                fs = app.freqA;
                f = app.fA;
                %parte importante
                txi = -app.timeA(end);
                txf = -app.timeA(1);
                tx = txi:1/fs:txf;

                app.input_A = x;
                app.timeA = tx;
                app.freqA = fs;
                app.fA = f;

                app.plot_Graphs(app.UIAxes_2,x,tx,fs);
                app.plot_Graphs(app.UIAxes_4,x,tx,fs);

            else
                if app.freqB==0
                    msgbox('Cannot perform operation, no data.', 'Error');
                    return;
                end
                x = app.input_B;
                fs = app.freqB;
                f = app.fB;

                txi = -app.timeB(end);
                txf = -app.timeB(1);
                tx = txi:1/fs:txf;

                app.input_B = x;
                app.timeB = tx;
                app.freqB = fs;
                app.fB = f;

                app.plot_Graphs(app.UIAxes,x,tx,fs);
                app.plot_Graphs(app.UIAxes_5,x,tx,fs);
            end
        end

        % Button pushed function: AdditionButton
        function AdditionButtonPushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            disp(size(app.input_A));
            disp(size(app.input_B));
            % Verificar si las longitudes coinciden
            if length(app.input_A) ~= length(app.input_B)
                msgbox('Signals have different lengths after processing.', 'Error');
                return;
            end
            app.output = app.input_A + app.input_B;
            app.timeout = app.timeA;
            app.freqout = app.freqA;
            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Button pushed function: SubstractButton_2
        function SubstractButton_2Pushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            disp(size(app.input_A));
            disp(size(app.input_B));
            app.output = app.input_A-app.input_B;
            app.timeout = app.timeA;
            app.freqout = app.freqA;
            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Button pushed function: ProductButton_2
        function ProductButton_2Pushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            disp(size(app.input_A));
            disp(size(app.input_B));
            app.output = app.input_A .* app.input_B;
            app.timeout = app.timeA;
            app.freqout = app.freqA;
            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Button pushed function: DivisionButton
        function DivisionButtonPushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            disp(size(app.input_A));
            disp(size(app.input_B));
            app.output = app.input_A ./ app.input_B;
            app.output(isnan(app.output)|isinf(app.output))=0;
            app.timeout = app.timeA;
            app.freqout = app.freqA;
            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Button pushed function: ConvolutionButton_3
        function ConvolutionButton_3Pushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            app.output = conv(app.input_A,app.input_B);
            app.output(isnan(app.output)|isinf(app.output))=0;
            % Ajustar el vector de tiempo
            output_length = length(app.output);
            % Usar el tiempo inicial de la primera señal y la duración de la salida
            start_time = app.timeA(1) + app.timeB(1);
            end_time = start_time + (output_length - 1) / app.freqA; % La duración de la señal de salida
            app.timeout = start_time:1/app.freqA:end_time;
            app.freqout = app.freqA;
            % Verifica que el tamaño del vector de tiempo coincida con la salida
            if length(app.timeout) ~= length(app.output)
                msgbox('The output length and time vector length do not match.', 'Error', 'error');
                return;
            end

            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Button pushed function: PotentiationButton_2
        function PotentiationButton_2Pushed(app, event)
            if app.freqA==0 || app.freqB==0
                msgbox('Cannot perform operation, missing data.', 'Error');
                return;
            end

            [app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB] =  app.resampleSignals(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA,app.freqB);
            [app.input_A,app.input_B,app.timeA,app.timeB] = app.pad(app.input_A,app.input_B,app.timeA,app.timeB,app.freqA);
            disp(size(app.input_A));
            disp(size(app.input_B));
            app.output = abs(app.input_A) .^ app.input_B;
            app.output(isnan(app.output)|isinf(app.output))=0;
            app.timeout = app.timeA;
            app.freqout = app.freqA;
            app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
            app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
        end

        % Value changed function: Filters_2
        function Filters_2ValueChanged(app, event)
            value = app.Filters_2.Value;
          %% parte de manzano
  % incion de filtro FIR Y IIR para la entrada a
  if app.INPUTAButton.Value == 1
                     
                     if strcmp(app.Filters_2.Value,'FIR1')
                         h = [-0.000779024944408449136062921880352405424
                             0.000999731876142587964212404827435420884
                             0.003831402877460619166294364390523696784
                             0.006331215011077922713567112111832102528
                             0.004203293161063792818255002714522561291
                            -0.006457016133813054704904121905428837636
                            -0.023107528738289289521601332921818539035
                            -0.033457559743661630102806725517439190298
                            -0.019904587044922469885044336024293443188
                             0.028977468595656773187041466144364676438
                             0.107219132728444313529259090955747524276
                             0.189581211105815133910468262001813855022
                             0.242562261249433724907831333439389709383
                             0.242562261249433724907831333439389709383
                             0.189581211105815133910468262001813855022
                             0.107219132728444313529259090955747524276
                             0.028977468595656773187041466144364676438
                            -0.019904587044922469885044336024293443188
                            -0.033457559743661630102806725517439190298
                            -0.023107528738289289521601332921818539035
                            -0.006457016133813054704904121905428837636
                             0.004203293161063792818255002714522561291
                             0.006331215011077922713567112111832102528
                             0.003831402877460619166294364390523696784
                             0.000999731876142587964212404827435420884
                            -0.000779024944408449136062921880352405424
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR1 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'Select')
                         return
                     elseif strcmp(app.Filters_2.Value,'FIR2')
                         h = [0.001884975120330046596470174868898084242
                             0.002419010748058414775329971746486990014
                             0.001590599013349243634821150727987060236
                            -0.002628390871960992539180246296837140108
                            -0.010170538297814175898725608249151264317
                            -0.015623780536385257697307160640320944367
                            -0.009593042963005444026536316926012659678
                             0.01388983701776957437046267074265415431 
                             0.048162323465281500767254385664273286238
                             0.070115607651607611550126364363677566871
                             0.044511802121696433798270220449921907857
                            -0.078704249325595526665111378861183766276
                            -0.586918083776343801893915497203124687076
                             0.586918083776343801893915497203124687076
                             0.078704249325595526665111378861183766276
                            -0.044511802121696433798270220449921907857
                            -0.070115607651607611550126364363677566871
                            -0.048162323465281500767254385664273286238
                            -0.01388983701776957437046267074265415431 
                             0.009593042963005444026536316926012659678
                             0.015623780536385257697307160640320944367
                             0.010170538297814175898725608249151264317
                             0.002628390871960992539180246296837140108
                            -0.001590599013349243634821150727987060236
                            -0.002419010748058414775329971746486990014
                            -0.001884975120330046596470174868898084242
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR2 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'FIR3')
                       h = [-0.000000000000000001922929161567581996809
                            -0.001569897304044256787849964673853264685
                             0.001799201069399839880008640236042083416
                             0.008990994433931174395491225936893897597
                             0.010792541725120091689049672822875436395
                            -0.000000000000000009354657166252341208504
                            -0.015835346743411488762420802345332049299
                            -0.019167509184833779467549419450733694248
                            -0.005459965985622059446435994090052190586
                             0.006568207170954977990040557500606155372
                            -0.000000000000000005317090080755072521592
                            -0.009469026014475398489445012728538131341
                             0.01140226611204093673523907170874736039 
                             0.058552268255241124439436362081323750317
                             0.071837724542372083380570302324485965073
                            -0.000000000000000021466635397073866283563
                            -0.116300660559505344981978680607426213101
                            -0.15922324040717714122017412137211067602 
                            -0.057452898173589478370004712814989034086
                             0.116208912714097029050819287476770114154
                             0.201752768291331929795262567495228722692
                             0.116208912714097029050819287476770114154
                            -0.057452898173589478370004712814989034086
                            -0.15922324040717714122017412137211067602 
                            -0.116300660559505344981978680607426213101
                            -0.000000000000000021466635397073866283563
                             0.071837724542372083380570302324485965073
                             0.058552268255241124439436362081323750317
                             0.01140226611204093673523907170874736039 
                            -0.009469026014475398489445012728538131341
                            -0.000000000000000005317090080755072521592
                             0.006568207170954977990040557500606155372
                            -0.005459965985622059446435994090052190586
                            -0.019167509184833779467549419450733694248
                            -0.015835346743411488762420802345332049299
                            -0.000000000000000009354657166252341208504
                             0.010792541725120091689049672822875436395
                             0.008990994433931174395491225936893897597
                             0.001799201069399839880008640236042083416
                            -0.001569897304044256787849964673853264685
                            -0.000000000000000001922929161567581996809
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR3 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'FIR4')
                       h = [-0.000000000000000007663637615541334456606
                             0.001564166308394175945381543435530602437
                            -0.001792632987859845957248627890123771067
                            -0.008958172318842353917811571761831146432
                            -0.010753142963480786728491800374740705593
                             0.000000000000000054897901653096777820552
                             0.015777538947276896486027197852308745496
                             0.019097537148141507240550751589580613654
                             0.005440034082284975550103389707601309055
                            -0.006544229572784126917894820252286081086
                            -0.000000000000000021190719014775402505014
                             0.009434458819054492365041575396844564239
                            -0.011360641518304093536273491338306484977
                            -0.058338520009537278565581175371335120872
                            -0.071575477014584637514005294178787153214
                             0.000000000000000057035387399646394219074
                             0.115876098661058898886544454853719798848
                             0.158641987292200414838561073338496498764
                             0.057243163238271957404812440017849439755
                            -0.11578468574615942421335290646311477758 
                             0.804065035269738492829105780401732772589
                            -0.11578468574615942421335290646311477758 
                             0.057243163238271957404812440017849439755
                             0.158641987292200414838561073338496498764
                             0.115876098661058898886544454853719798848
                             0.000000000000000057035387399646394219074
                            -0.071575477014584637514005294178787153214
                            -0.058338520009537278565581175371335120872
                            -0.011360641518304093536273491338306484977
                             0.009434458819054492365041575396844564239
                            -0.000000000000000021190719014775402505014
                            -0.006544229572784126917894820252286081086
                             0.005440034082284975550103389707601309055
                             0.019097537148141507240550751589580613654
                             0.015777538947276896486027197852308745496
                             0.000000000000000054897901653096777820552
                            -0.010753142963480786728491800374740705593
                            -0.008958172318842353917811571761831146432
                            -0.001792632987859845957248627890123771067
                             0.001564166308394175945381543435530602437
                            -0.000000000000000007663637615541334456606
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR4 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     % IIR   
                     elseif strcmp(app.Filters_2.Value,'IIR1')
                         dlg_title='IIR1';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.474416174978162796804781464743427932262 2.811006311911582233875606107176281511784 -1.703772240915468749733463482698425650597 0.544432694888534296495663511450402438641 -0.072315669102958557434845943134860135615', '0.003279216306360205161751775193579305778 0.016396081531801026676120613956300076097 0.032792163063602053352241227912600152194 0.032792163063602053352241227912600152194 0.016396081531801026676120613956300076097 0.003279216306360205161751775193579305778'}; 
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR2')
                         dlg_title='IIR2';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.474416174978162796804781464743427932262 2.811006311911582233875606107176281511784 -1.703772240915468305644253632635809481144 0.544432694888534185473361048934748396277 -0.072315669102958543557058135320403380319','0.268935721618647038955174366492428816855 -1.344678608093235139264720601204317063093 2.689357216186470722618651052471250295639 -2.689357216186470722618651052471250295639 1.344678608093235139264720601204317063093 -0.268935721618647038955174366492428816855'};
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR3')
                         dlg_title='IIR3';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.941867669925036121725270277238450944424 4.702317288464325173436009208671748638153 -4.634109656210405603360413806512951850891 3.077447793658273589301188621902838349342 -1.24661968102405040781377465464174747467 0.278059917634545739062446045863907784224','0.01809893300751444500384934599424013868 0 -0.054296799022543335011548037982720416039 0 0.054296799022543335011548037982720416039 0 -0.01809893300751444500384934599424013868'};
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR4')
                         dlg_title='IIR4';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.941867669925037009903689977363683283329 4.702317288464326949792848608922213315964 -4.634109656210408267895672906888648867607 3.077447793658276253836447722278535366058 -1.246619681024051962126009129860904067755 0.278059917634546294173958358442177996039','0.527624382501943212098183266789419576526 -1.956538810076257295378354683634825050831 4.001288117376634367872156872181221842766 -4.90951938700698864437299562268890440464 4.001288117376634367872156872181221842766 -1.956538810076257295378354683634825050831 0.527624382501943212098183266789419576526'};			
                        answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     end
                    

                      if app.control==1  %Convolución bloques
                       x=app.input_A;
                       app.h=transpose(str2num(answer{1}));
                       ni_h=str2num(answer{2}); % Instante inicial de h(n)
                       ni_x=app.timeA(1);
                       Lx= length(x);  %Longitud de X
                       ni_h=0;
                       Lh=length(app.h); %Longitud de h
                       Ly=Lx+Lh-1;       %tamaño del vector salida
                       
                       y=zeros(1,Ly);
                       Lxb=(str2num(answer{3}));  %Tamaño del bloque)           %Se define el tamaño del bloque 
                       Lyb=Lxb+Lh-1 ;    %Tamaño de salida del bloque
                       xb=x(1 : Lxb);    %Primer bloque
                       [y(1 : Lyb),nyi,nyf]=app.Convolucion(xb,ni_x,app.h,ni_h);
                       Nb= floor(Lx/Lxb);  %Numero de bloques totales
                       for k=2:Nb
                           xb=x(1+(k-1)* Lxb :  Lxb*k );   %Bloque siguiente
                           [yb,nyi1,nyf1]=app.Convolucion(xb,0,app.h,0);
                           y( 1+(k-1)*Lxb : (k-1)*Lxb+ Lyb)= y( 1+(k-1)*Lxb : (k-1)*Lxb+ Lyb) + yb;
                           app.salidaconv=y;
                           l1=length(app.salidaconv)/app.freqA; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                           l2=length(app.salidaconv);
                           time1=linspace(0,l1,l2);  
                           app.output=app.salidaconv;
                           app.timeout=time1;
                           app.freqout=app.freqA;
                           app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                           app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                           drawnow limitrate
                          % pause(0.00001);
                       end
                       %app.salidaconv=y;

                       Lxr= mod(Lx,Lxb); %Verificación si es bloque incompleto y adicionarlo a la salida
                       if Lxr ~= 0       %Si hay bloque incompleto
                           Lyr=Lxr+Lh-1 ;
                           xbr=x(1+ Nb*Lxb : Lx) ;
                           [ybr,nyi2,nyf2]=app.Convolucion(xbr,0,app.h,0); 
                           y( 1+Nb*Lxb: Ly)= y( 1+Nb*Lxb: Ly) + ybr;
                       end
                        app.salidaconv=y;
                        
                        l1=length(app.salidaconv)/app.freqA; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                        l2=length(app.salidaconv);
                        time1=linspace(0,l1,l2);  % transpuesta con los intervalos de tiempo
                        app.output=app.salidaconv;
                        app.timeout=time1;
                        app.freqout=app.freqA;
                        app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                        app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                        drawnow limitrate
                     
                       % plot(time1,app.salidaconv)
                        %freqz(app.salidaconv)
                        
                
                     end
                                         
                      if app.control==2 %%Ecuación en diferencia
                        dato=app.input_A;
                        
                        app.an=str2num(answer{1});
                        app.bn=str2num(answer{2});
                        
                          app.salidaedf=ecuaciondif(app,dato,app.an,app.bn);
                          l1=length(app.salidaedf)/app.freqA; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                          l2=length(app.salidaedf);
                          time1=linspace(0,l1,l2);  % transpuesta con los intervalos de tiempo
                          app.output=app.salidaedf;
                          app.timeout=time1;
                          app.freqout=app.freqA;
                          app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                          app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                          drawnow limitrate
                        
                        %freqz(app.salidaedf)
                       % plot(time1,app.salidaedf)
                      end
  end
% terminacion de la convolucion de la entrada a

%inicio entrada B
if app.INPUTBButton.Value == 1
                     
                     if strcmp(app.Filters_2.Value,'FIR1')
                         h = [-0.000779024944408449136062921880352405424
                             0.000999731876142587964212404827435420884
                             0.003831402877460619166294364390523696784
                             0.006331215011077922713567112111832102528
                             0.004203293161063792818255002714522561291
                            -0.006457016133813054704904121905428837636
                            -0.023107528738289289521601332921818539035
                            -0.033457559743661630102806725517439190298
                            -0.019904587044922469885044336024293443188
                             0.028977468595656773187041466144364676438
                             0.107219132728444313529259090955747524276
                             0.189581211105815133910468262001813855022
                             0.242562261249433724907831333439389709383
                             0.242562261249433724907831333439389709383
                             0.189581211105815133910468262001813855022
                             0.107219132728444313529259090955747524276
                             0.028977468595656773187041466144364676438
                            -0.019904587044922469885044336024293443188
                            -0.033457559743661630102806725517439190298
                            -0.023107528738289289521601332921818539035
                            -0.006457016133813054704904121905428837636
                             0.004203293161063792818255002714522561291
                             0.006331215011077922713567112111832102528
                             0.003831402877460619166294364390523696784
                             0.000999731876142587964212404827435420884
                            -0.000779024944408449136062921880352405424
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR1 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'Select')
                         return
                     elseif strcmp(app.Filters_2.Value,'FIR2')
                         h = [0.001884975120330046596470174868898084242
                             0.002419010748058414775329971746486990014
                             0.001590599013349243634821150727987060236
                            -0.002628390871960992539180246296837140108
                            -0.010170538297814175898725608249151264317
                            -0.015623780536385257697307160640320944367
                            -0.009593042963005444026536316926012659678
                             0.01388983701776957437046267074265415431 
                             0.048162323465281500767254385664273286238
                             0.070115607651607611550126364363677566871
                             0.044511802121696433798270220449921907857
                            -0.078704249325595526665111378861183766276
                            -0.586918083776343801893915497203124687076
                             0.586918083776343801893915497203124687076
                             0.078704249325595526665111378861183766276
                            -0.044511802121696433798270220449921907857
                            -0.070115607651607611550126364363677566871
                            -0.048162323465281500767254385664273286238
                            -0.01388983701776957437046267074265415431 
                             0.009593042963005444026536316926012659678
                             0.015623780536385257697307160640320944367
                             0.010170538297814175898725608249151264317
                             0.002628390871960992539180246296837140108
                            -0.001590599013349243634821150727987060236
                            -0.002419010748058414775329971746486990014
                            -0.001884975120330046596470174868898084242
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR2 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'FIR3')
                         h = [-0.000000000000000001922929161567581996809
                            -0.001569897304044256787849964673853264685
                             0.001799201069399839880008640236042083416
                             0.008990994433931174395491225936893897597
                             0.010792541725120091689049672822875436395
                            -0.000000000000000009354657166252341208504
                            -0.015835346743411488762420802345332049299
                            -0.019167509184833779467549419450733694248
                            -0.005459965985622059446435994090052190586
                             0.006568207170954977990040557500606155372
                            -0.000000000000000005317090080755072521592
                            -0.009469026014475398489445012728538131341
                             0.01140226611204093673523907170874736039 
                             0.058552268255241124439436362081323750317
                             0.071837724542372083380570302324485965073
                            -0.000000000000000021466635397073866283563
                            -0.116300660559505344981978680607426213101
                            -0.15922324040717714122017412137211067602 
                            -0.057452898173589478370004712814989034086
                             0.116208912714097029050819287476770114154
                             0.201752768291331929795262567495228722692
                             0.116208912714097029050819287476770114154
                            -0.057452898173589478370004712814989034086
                            -0.15922324040717714122017412137211067602 
                            -0.116300660559505344981978680607426213101
                            -0.000000000000000021466635397073866283563
                             0.071837724542372083380570302324485965073
                             0.058552268255241124439436362081323750317
                             0.01140226611204093673523907170874736039 
                            -0.009469026014475398489445012728538131341
                            -0.000000000000000005317090080755072521592
                             0.006568207170954977990040557500606155372
                            -0.005459965985622059446435994090052190586
                            -0.019167509184833779467549419450733694248
                            -0.015835346743411488762420802345332049299
                            -0.000000000000000009354657166252341208504
                             0.010792541725120091689049672822875436395
                             0.008990994433931174395491225936893897597
                             0.001799201069399839880008640236042083416
                            -0.001569897304044256787849964673853264685
                            -0.000000000000000001922929161567581996809
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR3 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     elseif strcmp(app.Filters_2.Value,'FIR4')
                        h = [-0.000000000000000007663637615541334456606
                             0.001564166308394175945381543435530602437
                            -0.001792632987859845957248627890123771067
                            -0.008958172318842353917811571761831146432
                            -0.010753142963480786728491800374740705593
                             0.000000000000000054897901653096777820552
                             0.015777538947276896486027197852308745496
                             0.019097537148141507240550751589580613654
                             0.005440034082284975550103389707601309055
                            -0.006544229572784126917894820252286081086
                            -0.000000000000000021190719014775402505014
                             0.009434458819054492365041575396844564239
                            -0.011360641518304093536273491338306484977
                            -0.058338520009537278565581175371335120872
                            -0.071575477014584637514005294178787153214
                             0.000000000000000057035387399646394219074
                             0.115876098661058898886544454853719798848
                             0.158641987292200414838561073338496498764
                             0.057243163238271957404812440017849439755
                            -0.11578468574615942421335290646311477758 
                             0.804065035269738492829105780401732772589
                            -0.11578468574615942421335290646311477758 
                             0.057243163238271957404812440017849439755
                             0.158641987292200414838561073338496498764
                             0.115876098661058898886544454853719798848
                             0.000000000000000057035387399646394219074
                            -0.071575477014584637514005294178787153214
                            -0.058338520009537278565581175371335120872
                            -0.011360641518304093536273491338306484977
                             0.009434458819054492365041575396844564239
                            -0.000000000000000021190719014775402505014
                            -0.006544229572784126917894820252286081086
                             0.005440034082284975550103389707601309055
                             0.019097537148141507240550751589580613654
                             0.015777538947276896486027197852308745496
                             0.000000000000000054897901653096777820552
                            -0.010753142963480786728491800374740705593
                            -0.008958172318842353917811571761831146432
                            -0.001792632987859845957248627890123771067
                             0.001564166308394175945381543435530602437
                            -0.000000000000000007663637615541334456606
                            ];
                    hstr = num2str(transpose(h));

                    dlg_title = 'FIR4 Parameters';
                    prompt = { 'h(n)' 'Initial Instant(ni)' 'Block Size(samples)'};
                    dims = [1 100];
                    defaultans = {hstr '0' '10'};
                    answer = inputdlg(prompt,dlg_title,dims,defaultans);
                    app.control=1;
                     % IIR   
                     elseif strcmp(app.Filters_2.Value,'IIR1')
                         dlg_title='IIR1';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.474416174978162796804781464743427932262 2.811006311911582233875606107176281511784 -1.703772240915468749733463482698425650597 0.544432694888534296495663511450402438641 -0.072315669102958557434845943134860135615', '0.003279216306360205161751775193579305778 0.016396081531801026676120613956300076097 0.032792163063602053352241227912600152194 0.032792163063602053352241227912600152194 0.016396081531801026676120613956300076097 0.003279216306360205161751775193579305778'}; 
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR2')
                         dlg_title='IIR2';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.474416174978162796804781464743427932262 2.811006311911582233875606107176281511784 -1.703772240915468305644253632635809481144 0.544432694888534185473361048934748396277 -0.072315669102958543557058135320403380319','0.268935721618647038955174366492428816855 -1.344678608093235139264720601204317063093 2.689357216186470722618651052471250295639 -2.689357216186470722618651052471250295639 1.344678608093235139264720601204317063093 -0.268935721618647038955174366492428816855'};
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR3')
                         dlg_title='IIR3';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.941867669925036121725270277238450944424 4.702317288464325173436009208671748638153 -4.634109656210405603360413806512951850891 3.077447793658273589301188621902838349342 -1.24661968102405040781377465464174747467 0.278059917634545739062446045863907784224','0.01809893300751444500384934599424013868 0 -0.054296799022543335011548037982720416039 0 0.054296799022543335011548037982720416039 0 -0.01809893300751444500384934599424013868'};
                         answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     elseif strcmp(app.Filters_2.Value,'IIR4')
                         dlg_title='IIR4';
                         dlg={'Input ak:', 'Input bk:'};
                         zize=[8 130];
                         defaultans = {'1 -2.941867669925037009903689977363683283329 4.702317288464326949792848608922213315964 -4.634109656210408267895672906888648867607 3.077447793658276253836447722278535366058 -1.246619681024051962126009129860904067755 0.278059917634546294173958358442177996039','0.527624382501943212098183266789419576526 -1.956538810076257295378354683634825050831 4.001288117376634367872156872181221842766 -4.90951938700698864437299562268890440464 4.001288117376634367872156872181221842766 -1.956538810076257295378354683634825050831 0.527624382501943212098183266789419576526'};			
                          answer=inputdlg(dlg,dlg_title,zize,defaultans);
                         app.control=2;
                     end
                    

                      if app.control==1  %Convolución bloques
                       x=app.input_B;
                       app.h=transpose(str2num(answer{1}));
                       ni_h=str2num(answer{2}); % Instante inicial de h(n)
                       ni_x=app.timeB(1);
                       Lx= length(x);  %Longitud de X
                       ni_h=0;
                       Lh=length(app.h); %Longitud de h
                       Ly=Lx+Lh-1;       %tamaño del vector salida
                       
                       y=zeros(1,Ly);
                       Lxb=(str2num(answer{3}));  %Tamaño del bloque)           %Se define el tamaño del bloque 
                       Lyb=Lxb+Lh-1 ;    %Tamaño de salida del bloque
                       xb=x(1 : Lxb);    %Primer bloque
                       [y(1 : Lyb),nyi,nyf]=app.Convolucion(xb,ni_x,app.h,ni_h);
                       Nb= floor(Lx/Lxb);  %Numero de bloques totales
                       for k=2:Nb
                           xb=x(1+(k-1)* Lxb :  Lxb*k );   %Bloque siguiente
                           [yb,nyi1,nyf1]=app.Convolucion(xb,0,app.h,0);
                           y( 1+(k-1)*Lxb : (k-1)*Lxb+ Lyb)= y( 1+(k-1)*Lxb : (k-1)*Lxb+ Lyb) + yb;
                           app.salidaconv=y;
                           l1=length(app.salidaconv)/app.freqB; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                           l2=length(app.salidaconv);
                           time1=linspace(0,l1,l2);  
                           app.output=app.salidaconv;
                           app.timeout=time1;
                           app.freqout=app.freqB;
                           app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                           app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                           drawnow limitrate
                          % pause(0.00001);
                       end
                       %app.salidaconv=y;

                       Lxr= mod(Lx,Lxb); %Verificación si es bloque incompleto y adicionarlo a la salida
                       if Lxr ~= 0       %Si hay bloque incompleto
                           Lyr=Lxr+Lh-1 ;
                           xbr=x(1+ Nb*Lxb : Lx) ;
                           [ybr,nyi2,nyf2]=app.Convolucion(xbr,0,app.h,0); 
                           y( 1+Nb*Lxb: Ly)= y( 1+Nb*Lxb: Ly) + ybr;
                       end
                        app.salidaconv=y;
                        
                        l1=length(app.salidaconv)/app.freqB; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                        l2=length(app.salidaconv);
                        time1=linspace(0,l1,l2);  % transpuesta con los intervalos de tiempo
                        app.output=app.salidaconv;
                        app.timeout=time1;
                        app.freqout=app.freqB;
                        app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                        app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                        drawnow limitrate
                     
                       % plot(time1,app.salidaconv)
                        %freqz(app.salidaconv)
                        
                
                     end
                                         
                      if app.control==2 %%Ecuación en diferencia
                        dato=app.input_B;
                        
                        app.an=str2num(answer{1});
                        app.bn=str2num(answer{2});
                        
                          app.salidaedf=ecuaciondif(app,dato,app.an,app.bn);
                          l1=length(app.salidaedf)/app.freqB; %dividir la longitud del vector sobre la frecuencia de muestreo se obtiene el total de tiempo
                          l2=length(app.salidaedf);
                          time1=linspace(0,l1,l2);  % transpuesta con los intervalos de tiempo
                          app.output=app.salidaedf;
                          app.timeout=time1;
                          app.freqout=app.freqB;
                          app.plot_Graphs(app.UIAxes_3,app.output,app.timeout,app.freqout);
                          app.plot_Graphs(app.UIAxes_6,app.output,app.timeout,app.freqout);
                          drawnow limitrate
                        
                        %freqz(app.salidaedf)
                       % plot(time1,app.salidaedf)
                      end
end
             
        end

        % Button pushed function: DataSerieButton
        function DataSerieButtonPushed(app, event)
            % Abrir una ventana de diálogo para ingresar los datos del usuario
            dlg_title = 'Discrete data series';
            prompt = {'Magnitude Vector (ej. [1, 2, 3]):', ...
                'Sample Frequency (Hz):', ...
                'Frequency (Hz):', ...
                'Time Vector (ej. [0, 1, 2]):'};
            num_lines = 1;
            defaultans = {'[1, 2, 3]', '1000', '1', '[0, 1, 2]'};

            % Recibir los datos a través de la ventana de diálogo
            answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

            % Validar si el usuario cerró la ventana o no ingresó datos
            if isempty(answer)
                return;  % Detener si se cierra la ventana
            end

            % Convertir los datos ingresados en vectores numéricos
            magnitudes = str2num(answer{1}); %#ok<*ST2NM>
            fs = str2double(answer{2});
            signal_frequency = str2double(answer{3});
            time = str2num(answer{4});

            % Validar los datos ingresados (NaN, Inf, y longitud de vectores)
            if isempty(magnitudes) || isempty(time) || isnan(fs) || isnan(signal_frequency)
                msgbox('The output length and time vector length do not match', 'Error', 'error');
                return;
            end

            if length(magnitudes) ~= length(time)
                msgbox('The magnitude vector and the Time vector must be of equal length.', 'Error', 'error');
                return;
            end

            % Asignar los datos ingresados a las variables correspondientes
            %input_signal = magnitudes;
            %sampling_frequency = fs;
            %time_vector = time;
            if  app.INPUTAButton.Value == true


                app.input_A = magnitudes;
                app.timeA = time;
                app.freqA = fs;
                app.fA = signal_frequency;
                if app.freqA==0
                    msgbox('Cannot perform operation, no data.', 'Error');
                    return;
                end

                if any(isnan(app.input_A(:))) || any(isinf(app.input_A(:)))
                     msgbox('The signal contains invalid data (NaN or Inf).', 'Error', 'error');

                     % Preguntar al usuario si desea corregir los datos
                     choice = questdlg('Do you want to interpolate invalid data?', ...
                         'Invalid data', 'Sí', 'No', 'Sí');

                     if strcmp(choice, 'Sí')
                         invalidIdx = isnan(app.input_A) | isinf(app.input_A);
                         validIdx = ~invalidIdx;
                         app.input_A(invalidIdx) = interp1(app.timeA(validIdx), app.input_A(validIdx), app.timeA(invalidIdx), 'linear');
                         % Rellenar con el primer o último valor válido en los extremos
                         app.input_A(1:find(validIdx, 1, 'first')-1) = app.input_A(find(validIdx, 1, 'first'));
                         app.input_A(find(validIdx, 1, 'last')+1:end) = app.input_A(find(validIdx, 1, 'last'));
                     else
                         return;
                     end
                 end

                app.plot_Graphs(app.UIAxes_2,app.input_A,app.timeA,app.freqA);
                app.plot_Graphs(app.UIAxes_4,app.input_A,app.timeA,app.freqA);

            else


                app.input_B = magnitudes;
                app.timeB = time;
                app.freqB = fs;
                app.fB = signal_frequency;

                if app.freqB==0
                    msgbox('Cannot perform operation, no data.', 'Error');
                    return;
                end

                if any(isnan(app.input_B(:))) || any(isinf(app.input_B(:)))
                     msgbox('The signal contains invalid data (NaN or Inf).', 'Error', 'error');

                     % Preguntar al usuario si desea corregir los datos
                     choice = questdlg('Do you want to interpolate invalid data?', ...
                         'Invalid data', 'Sí', 'No', 'Sí');

                     if strcmp(choice, 'Sí')
                         invalidIdx = isnan(app.input_B) | isinf(app.input_B);
                         validIdx = ~invalidIdx;
                         app.input_B(invalidIdx) = interp1(app.timeB(validIdx), app.input_B(validIdx), app.timeB(invalidIdx), 'linear');
                     else
                         return;
                     end
                 end

                app.plot_Graphs(app.UIAxes,app.input_B,app.timeB,app.freqB);
                app.plot_Graphs(app.UIAxes_5,app.input_B,app.timeB,app.freqB);
            end
            
        end

        % Button pushed function: Transfer1Button_5
        function Transfer1Button_5Pushed(app, event)
        
        end

        % Button pushed function: Transfer1Button_3
        function Transfer1Button_3Pushed(app, event)
        
        end

        % Button pushed function: Transfer1Button
        function Transfer1ButtonPushed(app, event)
         
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create Project_Dsp_2024UIFigure and hide until all components are created
            app.Project_Dsp_2024UIFigure = uifigure('Visible', 'off');
            app.Project_Dsp_2024UIFigure.AutoResizeChildren = 'off';
            app.Project_Dsp_2024UIFigure.Color = [0.0549 0.1373 0.2706];
            app.Project_Dsp_2024UIFigure.Position = [300 100 1094 559];
            app.Project_Dsp_2024UIFigure.Name = 'Project_Dsp_2024';

            % Create Panel_3
            app.Panel_3 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_3.AutoResizeChildren = 'off';
            app.Panel_3.BusyAction = 'cancel';
            app.Panel_3.Interruptible = 'off';
            app.Panel_3.Position = [144 13 942 531];

            % Create Image4
            app.Image4 = uiimage(app.Panel_3);
            app.Image4.ScaleMethod = 'stretch';
            app.Image4.Position = [1 0 939 530];
            app.Image4.ImageSource = fullfile(pathToMLAPP, 'señales_app_Disapositiva', 'AI_Generated_Image_2024-09-24_464842516034201.png');

            % Create DIGITALSIGNALSPROCESSINGPROJECT1Label
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label = uilabel(app.Panel_3);
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.HorizontalAlignment = 'center';
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.FontName = 'Times New Roman';
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.FontSize = 18;
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.FontWeight = 'bold';
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.FontAngle = 'italic';
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.FontColor = [0.6353 0.0784 0.1843];
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.Position = [206 441 518 44];
            app.DIGITALSIGNALSPROCESSINGPROJECT1Label.Text = 'DIGITAL SIGNALS PROCESSING - PROJECT 1';

            % Create JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel = uilabel(app.Panel_3);
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.FontName = 'Times New Roman';
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.FontSize = 18;
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.FontWeight = 'bold';
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.FontAngle = 'italic';
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.Position = [231 5 443 26];
            app.JrGomezAlejandroMuozJhoanSernaLuisManzanoLabel.Text = 'Jr.Gomez, Alejandro Muñoz, Jhoan Serna, Luis Manzano';

            % Create Panel_2
            app.Panel_2 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_2.AutoResizeChildren = 'off';
            app.Panel_2.Position = [12 14 121 531];

            % Create PresentProjectButton
            app.PresentProjectButton = uibutton(app.Panel_2, 'push');
            app.PresentProjectButton.ButtonPushedFcn = createCallbackFcn(app, @PresentProjectButtonPushed, true);
            app.PresentProjectButton.FontName = 'Times New Roman';
            app.PresentProjectButton.FontSize = 13;
            app.PresentProjectButton.FontWeight = 'bold';
            app.PresentProjectButton.Position = [12 494 100 24];
            app.PresentProjectButton.Text = 'Present Project';

            % Create Image2
            app.Image2 = uiimage(app.Panel_2);
            app.Image2.Position = [8 48 100 100];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'señales_app_Disapositiva', 'Imagen2.png');

            % Create UniversidadDelValleLabel
            app.UniversidadDelValleLabel = uilabel(app.Panel_2);
            app.UniversidadDelValleLabel.HorizontalAlignment = 'center';
            app.UniversidadDelValleLabel.FontName = 'Times New Roman';
            app.UniversidadDelValleLabel.FontSize = 18;
            app.UniversidadDelValleLabel.FontWeight = 'bold';
            app.UniversidadDelValleLabel.FontAngle = 'italic';
            app.UniversidadDelValleLabel.FontColor = [0.6353 0.0784 0.1843];
            app.UniversidadDelValleLabel.Position = [8 5 100 44];
            app.UniversidadDelValleLabel.Text = {'Universidad '; 'Del Valle'};

            % Create DropDown_2
            app.DropDown_2 = uidropdown(app.Panel_2);
            app.DropDown_2.Items = {'Input', 'Output'};
            app.DropDown_2.ValueChangedFcn = createCallbackFcn(app, @DropDown_2ValueChanged, true);
            app.DropDown_2.HandleVisibility = 'off';
            app.DropDown_2.Interruptible = 'off';
            app.DropDown_2.FontName = 'Times New Roman';
            app.DropDown_2.FontSize = 14;
            app.DropDown_2.FontAngle = 'italic';
            app.DropDown_2.Position = [89 424 23 28];
            app.DropDown_2.Value = 'Input';

            % Create DropDown_4
            app.DropDown_4 = uidropdown(app.Panel_2);
            app.DropDown_4.Items = {'Playback', 'Volume'};
            app.DropDown_4.Editable = 'on';
            app.DropDown_4.ValueChangedFcn = createCallbackFcn(app, @DropDown_4ValueChanged, true);
            app.DropDown_4.FontName = 'Times New Roman';
            app.DropDown_4.FontSize = 14;
            app.DropDown_4.FontAngle = 'italic';
            app.DropDown_4.BackgroundColor = [1 1 1];
            app.DropDown_4.Position = [89 389 23 28];
            app.DropDown_4.Value = 'Playback';

            % Create DropDown_5
            app.DropDown_5 = uidropdown(app.Panel_2);
            app.DropDown_5.Items = {'Signal Inputs', 'Synthetic Signals', 'Processing Functions', 'Others'};
            app.DropDown_5.Editable = 'on';
            app.DropDown_5.ValueChangedFcn = createCallbackFcn(app, @DropDown_5ValueChanged, true);
            app.DropDown_5.FontName = 'Times New Roman';
            app.DropDown_5.FontSize = 14;
            app.DropDown_5.FontAngle = 'italic';
            app.DropDown_5.BackgroundColor = [1 1 1];
            app.DropDown_5.Position = [88 461 23 23];
            app.DropDown_5.Value = 'Signal Inputs';

            % Create DataButton
            app.DataButton = uibutton(app.Panel_2, 'push');
            app.DataButton.ButtonPushedFcn = createCallbackFcn(app, @DataButtonPushed, true);
            app.DataButton.FontName = 'Times New Roman';
            app.DataButton.FontSize = 13;
            app.DataButton.FontWeight = 'bold';
            app.DataButton.Position = [13 424 80 28];
            app.DataButton.Text = 'Data';

            % Create AudioButton
            app.AudioButton = uibutton(app.Panel_2, 'push');
            app.AudioButton.ButtonPushedFcn = createCallbackFcn(app, @AudioButtonPushed, true);
            app.AudioButton.FontName = 'Times New Roman';
            app.AudioButton.FontSize = 13;
            app.AudioButton.FontWeight = 'bold';
            app.AudioButton.Position = [13 389 80 28];
            app.AudioButton.Text = 'Audio ';

            % Create DropDown_6
            app.DropDown_6 = uidropdown(app.Panel_2);
            app.DropDown_6.Items = {'Time', 'Frequency', 'Spectogram'};
            app.DropDown_6.Editable = 'on';
            app.DropDown_6.ValueChangedFcn = createCallbackFcn(app, @DropDown_6ValueChanged, true);
            app.DropDown_6.FontName = 'Times New Roman';
            app.DropDown_6.FontSize = 14;
            app.DropDown_6.FontAngle = 'italic';
            app.DropDown_6.BackgroundColor = [1 1 1];
            app.DropDown_6.Position = [90 356 23 24];
            app.DropDown_6.Value = 'Time';

            % Create GraphsButton
            app.GraphsButton = uibutton(app.Panel_2, 'push');
            app.GraphsButton.ButtonPushedFcn = createCallbackFcn(app, @GraphsButtonPushed, true);
            app.GraphsButton.FontName = 'Times New Roman';
            app.GraphsButton.FontSize = 13;
            app.GraphsButton.FontWeight = 'bold';
            app.GraphsButton.Position = [13 356 80 24];
            app.GraphsButton.Text = 'Graphs';

            % Create GraphsPanel
            app.GraphsPanel = uipanel(app.Panel_2);
            app.GraphsPanel.AutoResizeChildren = 'off';
            app.GraphsPanel.TitlePosition = 'centertop';
            app.GraphsPanel.Title = 'Graphs';
            app.GraphsPanel.FontName = 'Times New Roman';
            app.GraphsPanel.FontAngle = 'italic';
            app.GraphsPanel.FontWeight = 'bold';
            app.GraphsPanel.Position = [6 147 107 188];

            % Create ButtonGroup_7
            app.ButtonGroup_7 = uibuttongroup(app.GraphsPanel);
            app.ButtonGroup_7.AutoResizeChildren = 'off';
            app.ButtonGroup_7.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_7SelectionChanged, true);
            app.ButtonGroup_7.BorderType = 'none';
            app.ButtonGroup_7.TitlePosition = 'centertop';
            app.ButtonGroup_7.FontName = 'Times New Roman';
            app.ButtonGroup_7.Position = [6 20 100 135];

            % Create GraphInputAButton
            app.GraphInputAButton = uiradiobutton(app.ButtonGroup_7);
            app.GraphInputAButton.Text = 'Graph Input A';
            app.GraphInputAButton.FontName = 'Times New Roman';
            app.GraphInputAButton.Position = [4 112 91 22];
            app.GraphInputAButton.Value = true;

            % Create GraphInputBButton
            app.GraphInputBButton = uiradiobutton(app.ButtonGroup_7);
            app.GraphInputBButton.Text = 'Graph Input B';
            app.GraphInputBButton.FontName = 'Times New Roman';
            app.GraphInputBButton.Position = [4 79 91 22];

            % Create GraphOutputButton
            app.GraphOutputButton = uiradiobutton(app.ButtonGroup_7);
            app.GraphOutputButton.Text = 'Graph Output ';
            app.GraphOutputButton.FontName = 'Times New Roman';
            app.GraphOutputButton.Position = [4 49 91 22];

            % Create AlltheGrahsButton
            app.AlltheGrahsButton = uiradiobutton(app.ButtonGroup_7);
            app.AlltheGrahsButton.Text = 'All the Grahs';
            app.AlltheGrahsButton.FontName = 'Times New Roman';
            app.AlltheGrahsButton.Position = [4 21 87 22];

            % Create MenuButton
            app.MenuButton = uibutton(app.Panel_2, 'push');
            app.MenuButton.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton.FontName = 'Times New Roman';
            app.MenuButton.FontSize = 13;
            app.MenuButton.FontWeight = 'bold';
            app.MenuButton.Position = [12 458 80 28];
            app.MenuButton.Text = 'Menu';

            % Create VolumeControlPanel
            app.VolumeControlPanel = uipanel(app.Panel_2);
            app.VolumeControlPanel.AutoResizeChildren = 'off';
            app.VolumeControlPanel.TitlePosition = 'centertop';
            app.VolumeControlPanel.Title = 'Volume Control';
            app.VolumeControlPanel.FontName = 'Times New Roman';
            app.VolumeControlPanel.FontAngle = 'italic';
            app.VolumeControlPanel.FontWeight = 'bold';
            app.VolumeControlPanel.Position = [8 156 107 179];

            % Create ConstantVolumeSlider
            app.ConstantVolumeSlider = uislider(app.VolumeControlPanel);
            app.ConstantVolumeSlider.Limits = [1 100];
            app.ConstantVolumeSlider.MajorTicks = [];
            app.ConstantVolumeSlider.Position = [12 18 80 3];
            app.ConstantVolumeSlider.Value = 50;

            % Create ButtonGroup_2
            app.ButtonGroup_2 = uibuttongroup(app.VolumeControlPanel);
            app.ButtonGroup_2.AutoResizeChildren = 'off';
            app.ButtonGroup_2.BorderType = 'none';
            app.ButtonGroup_2.TitlePosition = 'centertop';
            app.ButtonGroup_2.FontName = 'Times New Roman';
            app.ButtonGroup_2.Position = [6 31 100 115];

            % Create IncreasingVolumeButton_2
            app.IncreasingVolumeButton_2 = uiradiobutton(app.ButtonGroup_2);
            app.IncreasingVolumeButton_2.Text = {'Increasing'; ' Volume'};
            app.IncreasingVolumeButton_2.FontName = 'Times New Roman';
            app.IncreasingVolumeButton_2.Position = [7 82 72 30];
            app.IncreasingVolumeButton_2.Value = true;

            % Create DecreasingVolumeButton_2
            app.DecreasingVolumeButton_2 = uiradiobutton(app.ButtonGroup_2);
            app.DecreasingVolumeButton_2.Text = {'Decreasing'; ' Volume'};
            app.DecreasingVolumeButton_2.FontName = 'Times New Roman';
            app.DecreasingVolumeButton_2.Position = [7 43 76 30];

            % Create ConstantVolumeButton
            app.ConstantVolumeButton = uiradiobutton(app.ButtonGroup_2);
            app.ConstantVolumeButton.Text = {'Constant'; ' Volume'};
            app.ConstantVolumeButton.FontName = 'Times New Roman';
            app.ConstantVolumeButton.Position = [7 4 65 30];

            % Create Image5
            app.Image5 = uiimage(app.ButtonGroup_2);
            app.Image5.Position = [82 81 15 29];
            app.Image5.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'volumen.png');

            % Create Image5_2
            app.Image5_2 = uiimage(app.ButtonGroup_2);
            app.Image5_2.Position = [82 42 15 29];
            app.Image5_2.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'bajar-volumen.png');

            % Create AudioPlaybackPanel_2
            app.AudioPlaybackPanel_2 = uipanel(app.Panel_2);
            app.AudioPlaybackPanel_2.AutoResizeChildren = 'off';
            app.AudioPlaybackPanel_2.TitlePosition = 'centertop';
            app.AudioPlaybackPanel_2.Title = 'Audio Playback';
            app.AudioPlaybackPanel_2.FontName = 'Times New Roman';
            app.AudioPlaybackPanel_2.FontAngle = 'italic';
            app.AudioPlaybackPanel_2.FontWeight = 'bold';
            app.AudioPlaybackPanel_2.Position = [8 154 107 179];

            % Create PlayButton_3
            app.PlayButton_3 = uibutton(app.AudioPlaybackPanel_2, 'push');
            app.PlayButton_3.ButtonPushedFcn = createCallbackFcn(app, @PlayButton_3Pushed, true);
            app.PlayButton_3.Position = [10 121 63 23];
            app.PlayButton_3.Text = 'Play';

            % Create StopButton_2
            app.StopButton_2 = uibutton(app.AudioPlaybackPanel_2, 'push');
            app.StopButton_2.ButtonPushedFcn = createCallbackFcn(app, @StopButton_2Pushed, true);
            app.StopButton_2.Position = [10 83 63 24];
            app.StopButton_2.Text = 'Stop';

            % Create ReverseButton_2
            app.ReverseButton_2 = uibutton(app.AudioPlaybackPanel_2, 'state');
            app.ReverseButton_2.Text = 'Reverse';
            app.ReverseButton_2.Position = [10 44 63 23];

            % Create EchoButton_2
            app.EchoButton_2 = uibutton(app.AudioPlaybackPanel_2, 'state');
            app.EchoButton_2.Text = 'Echo';
            app.EchoButton_2.Position = [10 10 63 21];

            % Create Image5_3
            app.Image5_3 = uiimage(app.AudioPlaybackPanel_2);
            app.Image5_3.Position = [78 118 15 29];
            app.Image5_3.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'jugar.png');

            % Create Image5_4
            app.Image5_4 = uiimage(app.AudioPlaybackPanel_2);
            app.Image5_4.Position = [78 81 15 29];
            app.Image5_4.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'prohibicion.png');

            % Create Image5_5
            app.Image5_5 = uiimage(app.AudioPlaybackPanel_2);
            app.Image5_5.Position = [78 41 15 29];
            app.Image5_5.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'girar-hacia-atras.png');

            % Create Image5_6
            app.Image5_6 = uiimage(app.AudioPlaybackPanel_2);
            app.Image5_6.Position = [78 6 15 29];
            app.Image5_6.ImageSource = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'png-transparent-computer-icons-avatar-amazon-echo-sound-audio-heroes-text-logo.png');

            % Create DataInputPanel
            app.DataInputPanel = uipanel(app.Panel_2);
            app.DataInputPanel.AutoResizeChildren = 'off';
            app.DataInputPanel.TitlePosition = 'centertop';
            app.DataInputPanel.Title = 'DataInput';
            app.DataInputPanel.FontName = 'Times New Roman';
            app.DataInputPanel.FontAngle = 'italic';
            app.DataInputPanel.FontWeight = 'bold';
            app.DataInputPanel.Position = [8 224 107 111];

            % Create ButtonGroup_6
            app.ButtonGroup_6 = uibuttongroup(app.DataInputPanel);
            app.ButtonGroup_6.AutoResizeChildren = 'off';
            app.ButtonGroup_6.BorderType = 'none';
            app.ButtonGroup_6.TitlePosition = 'centertop';
            app.ButtonGroup_6.FontName = 'Times New Roman';
            app.ButtonGroup_6.Position = [6 21 100 57];

            % Create INPUTAButton
            app.INPUTAButton = uiradiobutton(app.ButtonGroup_6);
            app.INPUTAButton.Text = 'INPUT A';
            app.INPUTAButton.FontName = 'Times New Roman';
            app.INPUTAButton.Position = [14 35 68 22];
            app.INPUTAButton.Value = true;

            % Create INPUTBButton
            app.INPUTBButton = uiradiobutton(app.ButtonGroup_6);
            app.INPUTBButton.Text = 'INPUT B';
            app.INPUTBButton.FontName = 'Times New Roman';
            app.INPUTBButton.Position = [14 4 68 22];

            % Create ShowOutputPanel
            app.ShowOutputPanel = uipanel(app.Panel_2);
            app.ShowOutputPanel.AutoResizeChildren = 'off';
            app.ShowOutputPanel.TitlePosition = 'centertop';
            app.ShowOutputPanel.Title = 'Show Output';
            app.ShowOutputPanel.FontName = 'Times New Roman';
            app.ShowOutputPanel.FontAngle = 'italic';
            app.ShowOutputPanel.FontWeight = 'bold';
            app.ShowOutputPanel.Position = [8 156 107 69];

            % Create ButtonGroup_5
            app.ButtonGroup_5 = uibuttongroup(app.ShowOutputPanel);
            app.ButtonGroup_5.AutoResizeChildren = 'off';
            app.ButtonGroup_5.BorderType = 'none';
            app.ButtonGroup_5.TitlePosition = 'centertop';
            app.ButtonGroup_5.FontName = 'Times New Roman';
            app.ButtonGroup_5.Position = [6 6 100 30];

            % Create OUTPUTCheckBox
            app.OUTPUTCheckBox = uicheckbox(app.ButtonGroup_5);
            app.OUTPUTCheckBox.Text = 'OUTPUT';
            app.OUTPUTCheckBox.FontName = 'Times New Roman';
            app.OUTPUTCheckBox.Position = [14 13 69 22];

            % Create Panel_7
            app.Panel_7 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_7.AutoResizeChildren = 'off';
            app.Panel_7.BusyAction = 'cancel';
            app.Panel_7.Interruptible = 'off';
            app.Panel_7.Position = [144 13 944 421];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.Panel_7);
            title(app.UIAxes_2, 'INPUT A')
            xlabel(app.UIAxes_2, 'X')
            ylabel(app.UIAxes_2, 'Y')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.FontName = 'Times New Roman';
            app.UIAxes_2.FontWeight = 'bold';
            app.UIAxes_2.XGrid = 'on';
            app.UIAxes_2.YGrid = 'on';
            app.UIAxes_2.Position = [10 34 911 375];

            % Create SaveInputAButton
            app.SaveInputAButton = uibutton(app.Panel_7, 'push');
            app.SaveInputAButton.ButtonPushedFcn = createCallbackFcn(app, @SaveInputAButtonPushed, true);
            app.SaveInputAButton.FontName = 'Times New Roman';
            app.SaveInputAButton.FontWeight = 'bold';
            app.SaveInputAButton.Position = [771 7 82 23];
            app.SaveInputAButton.Text = 'Save Input A';

            % Create Clear1Button
            app.Clear1Button = uibutton(app.Panel_7, 'push');
            app.Clear1Button.ButtonPushedFcn = createCallbackFcn(app, @Clear1ButtonPushed, true);
            app.Clear1Button.FontName = 'Times New Roman';
            app.Clear1Button.FontWeight = 'bold';
            app.Clear1Button.Position = [728 6 40 24];
            app.Clear1Button.Text = 'Clear';

            % Create CopyButton
            app.CopyButton = uibutton(app.Panel_7, 'push');
            app.CopyButton.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'copy.png');
            app.CopyButton.BackgroundColor = [1 1 1];
            app.CopyButton.Position = [894 6 27 24];
            app.CopyButton.Text = '';

            % Create Transfer1Button
            app.Transfer1Button = uibutton(app.Panel_7, 'push');
            app.Transfer1Button.ButtonPushedFcn = createCallbackFcn(app, @Transfer1ButtonPushed, true);
            app.Transfer1Button.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'girar-hacia-atras.png');
            app.Transfer1Button.BackgroundColor = [1 1 1];
            app.Transfer1Button.Position = [860 6 28 23];
            app.Transfer1Button.Text = '';

            % Create Panel_4
            app.Panel_4 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_4.AutoResizeChildren = 'off';
            app.Panel_4.Position = [146 448 945 95];

            % Create ChirpButton
            app.ChirpButton = uibutton(app.Panel_4, 'push');
            app.ChirpButton.ButtonPushedFcn = createCallbackFcn(app, @ChirpButtonPushed, true);
            app.ChirpButton.FontName = 'Times New Roman';
            app.ChirpButton.FontWeight = 'bold';
            app.ChirpButton.FontAngle = 'italic';
            app.ChirpButton.Position = [10 10 82 23];
            app.ChirpButton.Text = 'Chirp';

            % Create SincButton_2
            app.SincButton_2 = uibutton(app.Panel_4, 'push');
            app.SincButton_2.ButtonPushedFcn = createCallbackFcn(app, @SincButton_2Pushed, true);
            app.SincButton_2.FontName = 'Times New Roman';
            app.SincButton_2.FontWeight = 'bold';
            app.SincButton_2.FontAngle = 'italic';
            app.SincButton_2.Position = [101 10 82 23];
            app.SincButton_2.Text = 'Sinc';

            % Create Button
            app.Button = uibutton(app.Panel_4, 'push');
            app.Button.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Chirp_Waveform.png');
            app.Button.Position = [10 40 82 42];
            app.Button.Text = '';

            % Create Button_2
            app.Button_2 = uibutton(app.Panel_4, 'push');
            app.Button_2.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Sinc_Waveform.png');
            app.Button_2.Position = [100 40 82 42];
            app.Button_2.Text = '';

            % Create sinButton
            app.sinButton = uibutton(app.Panel_4, 'push');
            app.sinButton.ButtonPushedFcn = createCallbackFcn(app, @sinButtonPushed, true);
            app.sinButton.FontName = 'Times New Roman';
            app.sinButton.FontWeight = 'bold';
            app.sinButton.FontAngle = 'italic';
            app.sinButton.Position = [192 10 82 23];
            app.sinButton.Text = 'sin';

            % Create Button_3
            app.Button_3 = uibutton(app.Panel_4, 'push');
            app.Button_3.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Sin_Waveform.png');
            app.Button_3.Position = [191 40 82 42];
            app.Button_3.Text = '';

            % Create CosButton_2
            app.CosButton_2 = uibutton(app.Panel_4, 'push');
            app.CosButton_2.ButtonPushedFcn = createCallbackFcn(app, @CosButton_2Pushed, true);
            app.CosButton_2.FontName = 'Times New Roman';
            app.CosButton_2.FontWeight = 'bold';
            app.CosButton_2.FontAngle = 'italic';
            app.CosButton_2.Position = [281 10 82 23];
            app.CosButton_2.Text = 'Cos';

            % Create Button_4
            app.Button_4 = uibutton(app.Panel_4, 'push');
            app.Button_4.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Cos_Waveform.png');
            app.Button_4.Position = [280 40 82 42];
            app.Button_4.Text = '';

            % Create SawtoothButton_2
            app.SawtoothButton_2 = uibutton(app.Panel_4, 'push');
            app.SawtoothButton_2.ButtonPushedFcn = createCallbackFcn(app, @SawtoothButton_2Pushed, true);
            app.SawtoothButton_2.FontName = 'Times New Roman';
            app.SawtoothButton_2.FontWeight = 'bold';
            app.SawtoothButton_2.FontAngle = 'italic';
            app.SawtoothButton_2.Position = [375 10 82 23];
            app.SawtoothButton_2.Text = 'Sawtooth';

            % Create Button_5
            app.Button_5 = uibutton(app.Panel_4, 'push');
            app.Button_5.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'ST_Waveform.png');
            app.Button_5.Position = [374 40 82 42];
            app.Button_5.Text = '';

            % Create RampButton_2
            app.RampButton_2 = uibutton(app.Panel_4, 'push');
            app.RampButton_2.ButtonPushedFcn = createCallbackFcn(app, @RampButton_2Pushed, true);
            app.RampButton_2.FontName = 'Times New Roman';
            app.RampButton_2.FontWeight = 'bold';
            app.RampButton_2.FontAngle = 'italic';
            app.RampButton_2.Position = [465 10 82 23];
            app.RampButton_2.Text = 'Ramp';

            % Create Button_6
            app.Button_6 = uibutton(app.Panel_4, 'push');
            app.Button_6.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Ramp_Waveform.png');
            app.Button_6.Position = [464 40 82 42];
            app.Button_6.Text = '';

            % Create StepButton
            app.StepButton = uibutton(app.Panel_4, 'push');
            app.StepButton.ButtonPushedFcn = createCallbackFcn(app, @StepButtonPushed, true);
            app.StepButton.FontName = 'Times New Roman';
            app.StepButton.FontWeight = 'bold';
            app.StepButton.FontAngle = 'italic';
            app.StepButton.Position = [560 10 82 23];
            app.StepButton.Text = 'Step';

            % Create Button_7
            app.Button_7 = uibutton(app.Panel_4, 'push');
            app.Button_7.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'US_Waveform.png');
            app.Button_7.Position = [559 40 82 42];
            app.Button_7.Text = '';

            % Create TriangularButton_2
            app.TriangularButton_2 = uibutton(app.Panel_4, 'push');
            app.TriangularButton_2.ButtonPushedFcn = createCallbackFcn(app, @TriangularButton_2Pushed, true);
            app.TriangularButton_2.FontName = 'Times New Roman';
            app.TriangularButton_2.FontWeight = 'bold';
            app.TriangularButton_2.FontAngle = 'italic';
            app.TriangularButton_2.Position = [649 10 82 23];
            app.TriangularButton_2.Text = 'Triangular';

            % Create Button_8
            app.Button_8 = uibutton(app.Panel_4, 'push');
            app.Button_8.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Trig_Waveform.png');
            app.Button_8.Position = [648 40 82 42];
            app.Button_8.Text = '';

            % Create DataSerieButton
            app.DataSerieButton = uibutton(app.Panel_4, 'push');
            app.DataSerieButton.ButtonPushedFcn = createCallbackFcn(app, @DataSerieButtonPushed, true);
            app.DataSerieButton.FontName = 'Times New Roman';
            app.DataSerieButton.FontWeight = 'bold';
            app.DataSerieButton.FontAngle = 'italic';
            app.DataSerieButton.Position = [784 10 82 23];
            app.DataSerieButton.Text = 'Data Serie';

            % Create Button_29
            app.Button_29 = uibutton(app.Panel_4, 'push');
            app.Button_29.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'image.png');
            app.Button_29.Position = [784 40 82 42];
            app.Button_29.Text = '';

            % Create Panel_9
            app.Panel_9 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_9.AutoResizeChildren = 'off';
            app.Panel_9.Position = [147 447 944 97];

            % Create AdditionButton
            app.AdditionButton = uibutton(app.Panel_9, 'push');
            app.AdditionButton.ButtonPushedFcn = createCallbackFcn(app, @AdditionButtonPushed, true);
            app.AdditionButton.FontName = 'Times New Roman';
            app.AdditionButton.FontWeight = 'bold';
            app.AdditionButton.FontAngle = 'italic';
            app.AdditionButton.Position = [88 11 64 23];
            app.AdditionButton.Text = 'Addition';

            % Create Button_21
            app.Button_21 = uibutton(app.Panel_9, 'push');
            app.Button_21.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'mas.png');
            app.Button_21.Position = [97 41 45 42];
            app.Button_21.Text = '';

            % Create OperationsOnSignalsPanel_2
            app.OperationsOnSignalsPanel_2 = uipanel(app.Panel_9);
            app.OperationsOnSignalsPanel_2.AutoResizeChildren = 'off';
            app.OperationsOnSignalsPanel_2.TitlePosition = 'centertop';
            app.OperationsOnSignalsPanel_2.Title = 'Operations On Signals';
            app.OperationsOnSignalsPanel_2.FontName = 'Times New Roman';
            app.OperationsOnSignalsPanel_2.FontWeight = 'bold';
            app.OperationsOnSignalsPanel_2.Position = [759 1 184 95];

            % Create TemporalReflectionButton_3
            app.TemporalReflectionButton_3 = uibutton(app.OperationsOnSignalsPanel_2, 'push');
            app.TemporalReflectionButton_3.ButtonPushedFcn = createCallbackFcn(app, @TemporalReflectionButton_3Pushed, true);
            app.TemporalReflectionButton_3.WordWrap = 'on';
            app.TemporalReflectionButton_3.FontName = 'Times New Roman';
            app.TemporalReflectionButton_3.FontWeight = 'bold';
            app.TemporalReflectionButton_3.Position = [15 15 76 51];
            app.TemporalReflectionButton_3.Text = 'Temporal Reflection';

            % Create TimeShiftingButton
            app.TimeShiftingButton = uibutton(app.OperationsOnSignalsPanel_2, 'push');
            app.TimeShiftingButton.ButtonPushedFcn = createCallbackFcn(app, @TimeShiftingButtonPushed, true);
            app.TimeShiftingButton.WordWrap = 'on';
            app.TimeShiftingButton.FontName = 'Times New Roman';
            app.TimeShiftingButton.FontWeight = 'bold';
            app.TimeShiftingButton.Position = [106 15 76 51];
            app.TimeShiftingButton.Text = {'Time'; 'Shifting'};

            % Create OperationsBetweenSignalsButton
            app.OperationsBetweenSignalsButton = uibutton(app.Panel_9, 'push');
            app.OperationsBetweenSignalsButton.FontName = 'Times New Roman';
            app.OperationsBetweenSignalsButton.FontWeight = 'bold';
            app.OperationsBetweenSignalsButton.FontAngle = 'italic';
            app.OperationsBetweenSignalsButton.Position = [1 1 73 95];
            app.OperationsBetweenSignalsButton.Text = {'Operations '; 'Between '; 'Signals'};

            % Create SubstractButton_2
            app.SubstractButton_2 = uibutton(app.Panel_9, 'push');
            app.SubstractButton_2.ButtonPushedFcn = createCallbackFcn(app, @SubstractButton_2Pushed, true);
            app.SubstractButton_2.FontName = 'Times New Roman';
            app.SubstractButton_2.FontWeight = 'bold';
            app.SubstractButton_2.FontAngle = 'italic';
            app.SubstractButton_2.Position = [163 11 64 23];
            app.SubstractButton_2.Text = 'Substract';

            % Create Button_22
            app.Button_22 = uibutton(app.Panel_9, 'push');
            app.Button_22.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'resta.png');
            app.Button_22.Position = [172 41 45 42];
            app.Button_22.Text = '';

            % Create ProductButton_2
            app.ProductButton_2 = uibutton(app.Panel_9, 'push');
            app.ProductButton_2.ButtonPushedFcn = createCallbackFcn(app, @ProductButton_2Pushed, true);
            app.ProductButton_2.FontName = 'Times New Roman';
            app.ProductButton_2.FontWeight = 'bold';
            app.ProductButton_2.FontAngle = 'italic';
            app.ProductButton_2.Position = [239 11 64 23];
            app.ProductButton_2.Text = 'Product';

            % Create Button_23
            app.Button_23 = uibutton(app.Panel_9, 'push');
            app.Button_23.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'simbolo-de-multiplicacion.png');
            app.Button_23.Position = [248 41 45 42];
            app.Button_23.Text = '';

            % Create DivisionButton
            app.DivisionButton = uibutton(app.Panel_9, 'push');
            app.DivisionButton.ButtonPushedFcn = createCallbackFcn(app, @DivisionButtonPushed, true);
            app.DivisionButton.FontName = 'Times New Roman';
            app.DivisionButton.FontWeight = 'bold';
            app.DivisionButton.FontAngle = 'italic';
            app.DivisionButton.Position = [312 11 64 23];
            app.DivisionButton.Text = 'Division';

            % Create Button_24
            app.Button_24 = uibutton(app.Panel_9, 'push');
            app.Button_24.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'division (1).png');
            app.Button_24.Position = [321 41 45 42];
            app.Button_24.Text = '';

            % Create PotentiationButton_2
            app.PotentiationButton_2 = uibutton(app.Panel_9, 'push');
            app.PotentiationButton_2.ButtonPushedFcn = createCallbackFcn(app, @PotentiationButton_2Pushed, true);
            app.PotentiationButton_2.FontName = 'Times New Roman';
            app.PotentiationButton_2.FontWeight = 'bold';
            app.PotentiationButton_2.FontAngle = 'italic';
            app.PotentiationButton_2.Position = [389 11 76 23];
            app.PotentiationButton_2.Text = 'Potentiation';

            % Create Button_25
            app.Button_25 = uibutton(app.Panel_9, 'push');
            app.Button_25.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Potentiation.png');
            app.Button_25.Position = [389 41 76 42];
            app.Button_25.Text = '';

            % Create ConvolutionButton_3
            app.ConvolutionButton_3 = uibutton(app.Panel_9, 'push');
            app.ConvolutionButton_3.ButtonPushedFcn = createCallbackFcn(app, @ConvolutionButton_3Pushed, true);
            app.ConvolutionButton_3.FontName = 'Times New Roman';
            app.ConvolutionButton_3.FontWeight = 'bold';
            app.ConvolutionButton_3.FontAngle = 'italic';
            app.ConvolutionButton_3.Position = [473 11 76 23];
            app.ConvolutionButton_3.Text = 'Convolution';

            % Create Button_27
            app.Button_27 = uibutton(app.Panel_9, 'push');
            app.Button_27.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'Convolution.png');
            app.Button_27.Position = [473 41 76 42];
            app.Button_27.Text = '';

            % Create FIRIIRFiltersPanel_2
            app.FIRIIRFiltersPanel_2 = uipanel(app.Panel_9);
            app.FIRIIRFiltersPanel_2.AutoResizeChildren = 'off';
            app.FIRIIRFiltersPanel_2.TitlePosition = 'centertop';
            app.FIRIIRFiltersPanel_2.Title = 'FIR / IIR Filters';
            app.FIRIIRFiltersPanel_2.FontName = 'Times New Roman';
            app.FIRIIRFiltersPanel_2.FontWeight = 'bold';
            app.FIRIIRFiltersPanel_2.Position = [577 1 183 95];

            % Create Filters_2
            app.Filters_2 = uidropdown(app.FIRIIRFiltersPanel_2);
            app.Filters_2.Items = {'Select', 'FIR1', 'FIR2', 'FIR3', 'FIR4', 'IIR1', 'IIR2', 'IIR3', 'IIR4'};
            app.Filters_2.ValueChangedFcn = createCallbackFcn(app, @Filters_2ValueChanged, true);
            app.Filters_2.Position = [25 29 100 22];
            app.Filters_2.Value = 'Select';

            % Create Panel_8
            app.Panel_8 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_8.AutoResizeChildren = 'off';
            app.Panel_8.Position = [147 448 943 95];

            % Create StarRecordingButton
            app.StarRecordingButton = uibutton(app.Panel_8, 'push');
            app.StarRecordingButton.ButtonPushedFcn = createCallbackFcn(app, @StarRecordingButtonPushed, true);
            app.StarRecordingButton.FontName = 'Times New Roman';
            app.StarRecordingButton.FontWeight = 'bold';
            app.StarRecordingButton.FontAngle = 'italic';
            app.StarRecordingButton.Position = [211 12 78 22];
            app.StarRecordingButton.Text = 'Star Recording';

            % Create Button_12
            app.Button_12 = uibutton(app.Panel_8, 'push');
            app.Button_12.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'onda-sonora.png');
            app.Button_12.Position = [211 40 78 43];
            app.Button_12.Text = '';

            % Create AudioFileButton
            app.AudioFileButton = uibutton(app.Panel_8, 'push');
            app.AudioFileButton.ButtonPushedFcn = createCallbackFcn(app, @AudioFileButtonPushed, true);
            app.AudioFileButton.FontName = 'Times New Roman';
            app.AudioFileButton.FontWeight = 'bold';
            app.AudioFileButton.FontAngle = 'italic';
            app.AudioFileButton.Position = [298 12 82 23];
            app.AudioFileButton.Text = 'Audio File';

            % Create Button_13
            app.Button_13 = uibutton(app.Panel_8, 'push');
            app.Button_13.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'wav.png');
            app.Button_13.Position = [298 40 82 42];
            app.Button_13.Text = '';

            % Create DAQButton_2
            app.DAQButton_2 = uibutton(app.Panel_8, 'push');
            app.DAQButton_2.ButtonPushedFcn = createCallbackFcn(app, @DAQButton_2Pushed, true);
            app.DAQButton_2.FontName = 'Times New Roman';
            app.DAQButton_2.FontWeight = 'bold';
            app.DAQButton_2.FontAngle = 'italic';
            app.DAQButton_2.Position = [463 13 81 23];
            app.DAQButton_2.Text = 'DAQ';

            % Create Button_15
            app.Button_15 = uibutton(app.Panel_8, 'push');
            app.Button_15.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'NI.png');
            app.Button_15.Position = [463 40 82 42];
            app.Button_15.Text = '';

            % Create PreprocessingNormalizationPanel
            app.PreprocessingNormalizationPanel = uipanel(app.Panel_8);
            app.PreprocessingNormalizationPanel.AutoResizeChildren = 'off';
            app.PreprocessingNormalizationPanel.TitlePosition = 'centertop';
            app.PreprocessingNormalizationPanel.Title = 'Preprocessing Normalization';
            app.PreprocessingNormalizationPanel.FontName = 'Times New Roman';
            app.PreprocessingNormalizationPanel.FontWeight = 'bold';
            app.PreprocessingNormalizationPanel.Position = [783 1 161 93];

            % Create Button_19
            app.Button_19 = uibutton(app.PreprocessingNormalizationPanel, 'state');
            app.Button_19.ValueChangedFcn = createCallbackFcn(app, @Button_19ValueChanged, true);
            app.Button_19.Text = '0 , 1';
            app.Button_19.WordWrap = 'on';
            app.Button_19.FontName = 'Times New Roman';
            app.Button_19.FontWeight = 'bold';
            app.Button_19.Position = [32 38 39 23];

            % Create Button_20
            app.Button_20 = uibutton(app.PreprocessingNormalizationPanel, 'state');
            app.Button_20.ValueChangedFcn = createCallbackFcn(app, @Button_20ValueChanged, true);
            app.Button_20.Text = '-1, 1';
            app.Button_20.WordWrap = 'on';
            app.Button_20.FontName = 'Times New Roman';
            app.Button_20.FontWeight = 'bold';
            app.Button_20.Position = [94 38 39 23];

            % Create StandardNormalizationButton_2
            app.StandardNormalizationButton_2 = uibutton(app.PreprocessingNormalizationPanel, 'state');
            app.StandardNormalizationButton_2.ValueChangedFcn = createCallbackFcn(app, @StandardNormalizationButton_2ValueChanged, true);
            app.StandardNormalizationButton_2.Text = 'Standard Normalization';
            app.StandardNormalizationButton_2.WordWrap = 'on';
            app.StandardNormalizationButton_2.FontName = 'Times New Roman';
            app.StandardNormalizationButton_2.FontWeight = 'bold';
            app.StandardNormalizationButton_2.Position = [16 5 136 24];

            % Create RecordingParamatersPanel
            app.RecordingParamatersPanel = uipanel(app.Panel_8);
            app.RecordingParamatersPanel.AutoResizeChildren = 'off';
            app.RecordingParamatersPanel.TitlePosition = 'centertop';
            app.RecordingParamatersPanel.Title = 'Recording Paramaters';
            app.RecordingParamatersPanel.FontName = 'Times New Roman';
            app.RecordingParamatersPanel.FontWeight = 'bold';
            app.RecordingParamatersPanel.Position = [1 1 199 95];

            % Create ChannelsSwitchLabel
            app.ChannelsSwitchLabel = uilabel(app.RecordingParamatersPanel);
            app.ChannelsSwitchLabel.HorizontalAlignment = 'center';
            app.ChannelsSwitchLabel.FontName = 'Times New Roman';
            app.ChannelsSwitchLabel.FontWeight = 'bold';
            app.ChannelsSwitchLabel.Position = [18 41 53 22];
            app.ChannelsSwitchLabel.Text = 'Channels';

            % Create ChannelsSwitch
            app.ChannelsSwitch = uiswitch(app.RecordingParamatersPanel, 'slider');
            app.ChannelsSwitch.Items = {'1', '2'};
            app.ChannelsSwitch.FontName = 'Times New Roman';
            app.ChannelsSwitch.FontWeight = 'bold';
            app.ChannelsSwitch.Position = [23 14 45 20];
            app.ChannelsSwitch.Value = '1';

            % Create NumberofSamplesEditFieldLabel
            app.NumberofSamplesEditFieldLabel = uilabel(app.Panel_8);
            app.NumberofSamplesEditFieldLabel.HorizontalAlignment = 'center';
            app.NumberofSamplesEditFieldLabel.FontName = 'Times New Roman';
            app.NumberofSamplesEditFieldLabel.FontWeight = 'bold';
            app.NumberofSamplesEditFieldLabel.Position = [115 37 63 30];
            app.NumberofSamplesEditFieldLabel.Text = {'Number of '; 'Samples'};

            % Create NumberofSamplesEditField
            app.NumberofSamplesEditField = uieditfield(app.Panel_8, 'numeric');
            app.NumberofSamplesEditField.ValueDisplayFormat = '%d ';
            app.NumberofSamplesEditField.HorizontalAlignment = 'center';
            app.NumberofSamplesEditField.FontName = 'Times New Roman';
            app.NumberofSamplesEditField.FontWeight = 'bold';
            app.NumberofSamplesEditField.Position = [105 13 81 22];
            app.NumberofSamplesEditField.Value = 44100;

            % Create MonoStereoSignalsPanel_3
            app.MonoStereoSignalsPanel_3 = uipanel(app.Panel_8);
            app.MonoStereoSignalsPanel_3.AutoResizeChildren = 'off';
            app.MonoStereoSignalsPanel_3.TitlePosition = 'centertop';
            app.MonoStereoSignalsPanel_3.Title = 'Mono/Stereo Signals';
            app.MonoStereoSignalsPanel_3.FontName = 'Times New Roman';
            app.MonoStereoSignalsPanel_3.FontWeight = 'bold';
            app.MonoStereoSignalsPanel_3.Position = [551 1 233 94];

            % Create MonoStereoButton_3
            app.MonoStereoButton_3 = uibutton(app.MonoStereoSignalsPanel_3, 'push');
            app.MonoStereoButton_3.ButtonPushedFcn = createCallbackFcn(app, @MonoStereoButton_3Pushed, true);
            app.MonoStereoButton_3.WordWrap = 'on';
            app.MonoStereoButton_3.FontName = 'Times New Roman';
            app.MonoStereoButton_3.FontWeight = 'bold';
            app.MonoStereoButton_3.Position = [7 43 76 22];
            app.MonoStereoButton_3.Text = 'Mono/Stereo';

            % Create StereoMonoButton_3
            app.StereoMonoButton_3 = uibutton(app.MonoStereoSignalsPanel_3, 'push');
            app.StereoMonoButton_3.ButtonPushedFcn = createCallbackFcn(app, @StereoMonoButton_3Pushed, true);
            app.StereoMonoButton_3.WordWrap = 'on';
            app.StereoMonoButton_3.FontName = 'Times New Roman';
            app.StereoMonoButton_3.FontWeight = 'bold';
            app.StereoMonoButton_3.Position = [7 14 76 22];
            app.StereoMonoButton_3.Text = 'Stereo/Mono';

            % Create MergeMonoStereoButton_3
            app.MergeMonoStereoButton_3 = uibutton(app.MonoStereoSignalsPanel_3, 'push');
            app.MergeMonoStereoButton_3.ButtonPushedFcn = createCallbackFcn(app, @MergeMonoStereoButton_3Pushed, true);
            app.MergeMonoStereoButton_3.WordWrap = 'on';
            app.MergeMonoStereoButton_3.FontName = 'Times New Roman';
            app.MergeMonoStereoButton_3.FontWeight = 'bold';
            app.MergeMonoStereoButton_3.Position = [99 42 124 23];
            app.MergeMonoStereoButton_3.Text = 'Merge Mono/Stereo';

            % Create UndersampleOversampleButton_3
            app.UndersampleOversampleButton_3 = uibutton(app.MonoStereoSignalsPanel_3, 'push');
            app.UndersampleOversampleButton_3.ButtonPushedFcn = createCallbackFcn(app, @UndersampleOversampleButton_3Pushed, true);
            app.UndersampleOversampleButton_3.WordWrap = 'on';
            app.UndersampleOversampleButton_3.FontName = 'Times New Roman';
            app.UndersampleOversampleButton_3.FontWeight = 'bold';
            app.UndersampleOversampleButton_3.Position = [89 13 145 23];
            app.UndersampleOversampleButton_3.Text = 'Undersample/Oversample';

            % Create Button_28
            app.Button_28 = uibutton(app.Panel_8, 'push');
            app.Button_28.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'archivo-excel.png');
            app.Button_28.Position = [386 41 73 41];
            app.Button_28.Text = '';

            % Create ExcelButton
            app.ExcelButton = uibutton(app.Panel_8, 'push');
            app.ExcelButton.FontName = 'Times New Roman';
            app.ExcelButton.FontWeight = 'bold';
            app.ExcelButton.FontAngle = 'italic';
            app.ExcelButton.Position = [386 12 73 23];
            app.ExcelButton.Text = 'Excel';

            % Create Panel_11
            app.Panel_11 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_11.AutoResizeChildren = 'off';
            app.Panel_11.BusyAction = 'cancel';
            app.Panel_11.Interruptible = 'off';
            app.Panel_11.Position = [146 9 944 421];

            % Create UIAxes
            app.UIAxes = uiaxes(app.Panel_11);
            title(app.UIAxes, 'INPUT B')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Times New Roman';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [17 32 911 378];

            % Create SaveInputBButton
            app.SaveInputBButton = uibutton(app.Panel_11, 'push');
            app.SaveInputBButton.ButtonPushedFcn = createCallbackFcn(app, @SaveInputBButtonPushed, true);
            app.SaveInputBButton.FontName = 'Times New Roman';
            app.SaveInputBButton.FontWeight = 'bold';
            app.SaveInputBButton.Position = [774 10 82 23];
            app.SaveInputBButton.Text = 'Save Input B';

            % Create Clear1Button_3
            app.Clear1Button_3 = uibutton(app.Panel_11, 'push');
            app.Clear1Button_3.ButtonPushedFcn = createCallbackFcn(app, @Clear1Button_3Pushed, true);
            app.Clear1Button_3.FontName = 'Times New Roman';
            app.Clear1Button_3.FontWeight = 'bold';
            app.Clear1Button_3.Position = [731 9 40 24];
            app.Clear1Button_3.Text = 'Clear';

            % Create CopyButton_3
            app.CopyButton_3 = uibutton(app.Panel_11, 'push');
            app.CopyButton_3.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'copy.png');
            app.CopyButton_3.BackgroundColor = [1 1 1];
            app.CopyButton_3.Position = [896 9 27 24];
            app.CopyButton_3.Text = '';

            % Create Transfer1Button_3
            app.Transfer1Button_3 = uibutton(app.Panel_11, 'push');
            app.Transfer1Button_3.ButtonPushedFcn = createCallbackFcn(app, @Transfer1Button_3Pushed, true);
            app.Transfer1Button_3.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'girar-hacia-atras.png');
            app.Transfer1Button_3.BackgroundColor = [1 1 1];
            app.Transfer1Button_3.Position = [863 9 28 23];
            app.Transfer1Button_3.Text = '';

            % Create Panel_12
            app.Panel_12 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_12.AutoResizeChildren = 'off';
            app.Panel_12.BusyAction = 'cancel';
            app.Panel_12.Interruptible = 'off';
            app.Panel_12.Position = [146 9 944 421];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.Panel_12);
            title(app.UIAxes_3, 'OUTPUT')
            xlabel(app.UIAxes_3, 'X')
            ylabel(app.UIAxes_3, 'Y')
            zlabel(app.UIAxes_3, 'Z')
            app.UIAxes_3.FontName = 'Times New Roman';
            app.UIAxes_3.XGrid = 'on';
            app.UIAxes_3.YGrid = 'on';
            app.UIAxes_3.Position = [17 30 911 380];

            % Create SaveOutputButton
            app.SaveOutputButton = uibutton(app.Panel_12, 'push');
            app.SaveOutputButton.ButtonPushedFcn = createCallbackFcn(app, @SaveOutputButtonPushed, true);
            app.SaveOutputButton.FontName = 'Times New Roman';
            app.SaveOutputButton.FontWeight = 'bold';
            app.SaveOutputButton.Position = [840 8 82 23];
            app.SaveOutputButton.Text = 'Save Output';

            % Create Clear1Button_2
            app.Clear1Button_2 = uibutton(app.Panel_12, 'push');
            app.Clear1Button_2.ButtonPushedFcn = createCallbackFcn(app, @Clear1Button_2Pushed, true);
            app.Clear1Button_2.FontName = 'Times New Roman';
            app.Clear1Button_2.FontWeight = 'bold';
            app.Clear1Button_2.Position = [796 7 40 24];
            app.Clear1Button_2.Text = 'Clear';

            % Create Panel_13
            app.Panel_13 = uipanel(app.Project_Dsp_2024UIFigure);
            app.Panel_13.AutoResizeChildren = 'off';
            app.Panel_13.BusyAction = 'cancel';
            app.Panel_13.Interruptible = 'off';
            app.Panel_13.Position = [144 13 944 421];

            % Create UIAxes_4
            app.UIAxes_4 = uiaxes(app.Panel_13);
            title(app.UIAxes_4, 'INPUT A')
            xlabel(app.UIAxes_4, 'X')
            ylabel(app.UIAxes_4, 'Y')
            zlabel(app.UIAxes_4, 'Z')
            app.UIAxes_4.FontName = 'Times New Roman';
            app.UIAxes_4.XGrid = 'on';
            app.UIAxes_4.YGrid = 'on';
            app.UIAxes_4.Position = [31 230 416 180];

            % Create UIAxes_5
            app.UIAxes_5 = uiaxes(app.Panel_13);
            title(app.UIAxes_5, 'INPUT B')
            xlabel(app.UIAxes_5, 'X')
            ylabel(app.UIAxes_5, 'Y')
            zlabel(app.UIAxes_5, 'Z')
            app.UIAxes_5.FontName = 'Times New Roman';
            app.UIAxes_5.XGrid = 'on';
            app.UIAxes_5.YGrid = 'on';
            app.UIAxes_5.Position = [497 229 416 181];

            % Create UIAxes_6
            app.UIAxes_6 = uiaxes(app.Panel_13);
            title(app.UIAxes_6, 'OUTPUT')
            xlabel(app.UIAxes_6, 'X')
            ylabel(app.UIAxes_6, 'Y')
            zlabel(app.UIAxes_6, 'Z')
            app.UIAxes_6.FontName = 'Times New Roman';
            app.UIAxes_6.XGrid = 'on';
            app.UIAxes_6.YGrid = 'on';
            app.UIAxes_6.Position = [33 17 882 187];

            % Create SaveOutputButton_2
            app.SaveOutputButton_2 = uibutton(app.Panel_13, 'push');
            app.SaveOutputButton_2.ButtonPushedFcn = createCallbackFcn(app, @SaveOutputButton_2Pushed, true);
            app.SaveOutputButton_2.FontName = 'Times New Roman';
            app.SaveOutputButton_2.FontSize = 9;
            app.SaveOutputButton_2.FontWeight = 'bold';
            app.SaveOutputButton_2.Position = [833 5 82 18];
            app.SaveOutputButton_2.Text = 'Save Output';

            % Create Clear1Button_4
            app.Clear1Button_4 = uibutton(app.Panel_13, 'push');
            app.Clear1Button_4.ButtonPushedFcn = createCallbackFcn(app, @Clear1Button_4Pushed, true);
            app.Clear1Button_4.FontName = 'Times New Roman';
            app.Clear1Button_4.FontSize = 9;
            app.Clear1Button_4.FontWeight = 'bold';
            app.Clear1Button_4.Position = [782 4 40 17];
            app.Clear1Button_4.Text = 'Clear';

            % Create SaveInputBButton_3
            app.SaveInputBButton_3 = uibutton(app.Panel_13, 'push');
            app.SaveInputBButton_3.ButtonPushedFcn = createCallbackFcn(app, @SaveInputBButton_3Pushed, true);
            app.SaveInputBButton_3.FontName = 'Times New Roman';
            app.SaveInputBButton_3.FontSize = 9;
            app.SaveInputBButton_3.FontWeight = 'bold';
            app.SaveInputBButton_3.Position = [751 212 82 17];
            app.SaveInputBButton_3.Text = 'Save Input B';

            % Create Clear1Button_5
            app.Clear1Button_5 = uibutton(app.Panel_13, 'push');
            app.Clear1Button_5.ButtonPushedFcn = createCallbackFcn(app, @Clear1Button_5Pushed, true);
            app.Clear1Button_5.FontName = 'Times New Roman';
            app.Clear1Button_5.FontSize = 9;
            app.Clear1Button_5.FontWeight = 'bold';
            app.Clear1Button_5.Position = [708 212 40 17];
            app.Clear1Button_5.Text = 'Clear';

            % Create CopyButton_5
            app.CopyButton_5 = uibutton(app.Panel_13, 'push');
            app.CopyButton_5.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'copy.png');
            app.CopyButton_5.BackgroundColor = [1 1 1];
            app.CopyButton_5.FontSize = 9;
            app.CopyButton_5.Position = [879 211 36 17];
            app.CopyButton_5.Text = '';

            % Create Transfer1Button_5
            app.Transfer1Button_5 = uibutton(app.Panel_13, 'push');
            app.Transfer1Button_5.ButtonPushedFcn = createCallbackFcn(app, @Transfer1Button_5Pushed, true);
            app.Transfer1Button_5.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'girar-hacia-atras.png');
            app.Transfer1Button_5.BackgroundColor = [1 1 1];
            app.Transfer1Button_5.FontSize = 9;
            app.Transfer1Button_5.Position = [840 211 35 17];
            app.Transfer1Button_5.Text = '';

            % Create SaveInputAButton_2
            app.SaveInputAButton_2 = uibutton(app.Panel_13, 'push');
            app.SaveInputAButton_2.ButtonPushedFcn = createCallbackFcn(app, @SaveInputAButton_2Pushed, true);
            app.SaveInputAButton_2.FontName = 'Times New Roman';
            app.SaveInputAButton_2.FontSize = 9;
            app.SaveInputAButton_2.FontWeight = 'bold';
            app.SaveInputAButton_2.Position = [283 211 82 17];
            app.SaveInputAButton_2.Text = 'Save Input A';

            % Create Clear1Button_6
            app.Clear1Button_6 = uibutton(app.Panel_13, 'push');
            app.Clear1Button_6.ButtonPushedFcn = createCallbackFcn(app, @Clear1Button_6Pushed, true);
            app.Clear1Button_6.FontName = 'Times New Roman';
            app.Clear1Button_6.FontSize = 9;
            app.Clear1Button_6.FontWeight = 'bold';
            app.Clear1Button_6.Position = [240 210 40 18];
            app.Clear1Button_6.Text = 'Clear';

            % Create CopyButton_6
            app.CopyButton_6 = uibutton(app.Panel_13, 'push');
            app.CopyButton_6.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'copy.png');
            app.CopyButton_6.BackgroundColor = [1 1 1];
            app.CopyButton_6.FontSize = 9;
            app.CopyButton_6.Position = [411 211 36 17];
            app.CopyButton_6.Text = '';

            % Create Transfer1Button_6
            app.Transfer1Button_6 = uibutton(app.Panel_13, 'push');
            app.Transfer1Button_6.ButtonPushedFcn = createCallbackFcn(app, @Transfer1Button_6Pushed, true);
            app.Transfer1Button_6.Icon = fullfile(pathToMLAPP, 'inconos_app_diapositiva', 'girar-hacia-atras.png');
            app.Transfer1Button_6.BackgroundColor = [1 1 1];
            app.Transfer1Button_6.FontSize = 9;
            app.Transfer1Button_6.Position = [375 211 35 17];
            app.Transfer1Button_6.Text = '';

            % Show the figure after all components are created
            app.Project_Dsp_2024UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Project_of_Dsp_App_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.Project_Dsp_2024UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.Project_Dsp_2024UIFigure)
        end
    end
end