classdef GetAudioInputs < handle
    properties
        audioData
        time
        fs %Frecuencia de muestreo
        f
        serialObj % Objeto de la conexión serial
    end
    
    methods
        function audio = GrabarAudio(obj, switchValue,Fs)
            % Button pushed function: MicButton
            % Este método es un callback que se ejecuta cuando se presiona un botón llamado MicButton

            answer = inputdlg({'Audio Name:','Duration:','Bits per sample:'},'Recording Properties',1,{'Recording_1','3','16'});
            
            if isempty(answer) %En caso de que se cierre la pestaña
                audio = 0;
                return;
            end
            
            duration = str2double(answer{2}); %La duración será la segunda respuesta
            bps = str2double(answer{3});
            rec = audiorecorder(Fs,bps, str2double(regexp(switchValue,'\d*','match'))); %La frecuencia de muestreo será la que el usuario elija, 16 bits, Mono
            recordblocking(rec, duration); %Ingresa los parametros y la duración para grabar

            y = (getaudiodata(rec))'; %Almacena los valores de la grabación en y recuperando los datos de audio de rec
            audiowrite(strcat(answer{1}, '.wav'), transpose(y), Fs,'BitsPerSample',bps); %Nombra el audio, lo guarda en un .wav, le ingresa los datos de grabación y la frec muestreo
                 
            msgbox(strcat('Sucessfully Saved', " ", answer{1}, '.wav'), 'Success', 'info'); %Imprime un aviso que dice que grabó, el nombre del archivo

            ts = (0:length(y)-1) / Fs; %Recupera los datos del tiempo a partir de la frec y la longitud de y
            
            obj.audioData = y;
            obj.time = ts;
            obj.fs = Fs;
            obj.f = 0;
            audio = 1;
        end
        function arch = AbrirArchivo(obj, ~)
           [filename, pathname] = uigetfile('*.mp3;*.m4a;*.wav', 'Select a audio file'); % Seleccionar el archivo de audio
           
           if filename == 0
               arch = 0;
               return;
           end

           [x,Fs] = audioread(fullfile(pathname,filename)); % Leer el archivo de audio
           
           obj.audioData = transpose(x);
           obj.time = transpose((0:length(x)-1) / Fs);
           obj.fs = Fs;
           obj.f = 0;
           arch = 1;
        end
        function dq = DAQ(obj)
            % Configurar la conexión serial con el Arduino
            try
                obj.serialObj = serialport('COM3', 115200); % Ajusta el puerto COM
                configureTerminator(obj.serialObj, 'LF');
            catch ME
                msgbox(['Cant acces to port ', ME.message], 'Error', 'error');
                dq = 0;
                return;
            end

            % Parámetros de adquisición
            answer = inputdlg({'Time (s):', 'Sampling Frequency (Hz):'}, 'Arduino Parameters', 1, {'5', '1000'});
            
            if isempty(answer) % En caso de que se cierre la ventana de diálogo
                dq = 0;
                return;
            end

            duration = str2double(answer{1});
            sampleRate = str2double(answer{2});
            numSamples = duration * sampleRate;

            % Verificar la validez de los parámetros
            if isnan(duration) || isnan(sampleRate) || duration <= 0 || sampleRate <= 0
                msgbox('Invalid duration of signal or sampling frequency.', 'Error', 'error');
                dq = 0;
                return;
            end

            % Inicializar vectores
            obj.time = zeros(1, numSamples);
            obj.audioData = zeros(1, numSamples);
            obj.fs = sampleRate;
            obj.f = 0;

            % Captura de datos desde el Arduino
            try
                for i = 1:numSamples
                    newValue = readline(obj.serialObj);
                    newValue = str2double(newValue);

                    if isnan(newValue) || isinf(newValue)
                        msgbox('Invalid data from Arduino.', 'Error', 'error');
                        dq = 0;
                        return;
                    end

                    obj.audioData(i) = newValue;
                    obj.time(i) = (i - 1) / sampleRate;
                end
                dq = 1;
            catch ME
                msgbox(['Error in input data reading: ', ME.message], 'Error', 'error');
                dq = 0;
            end
        end


        



    end


end
