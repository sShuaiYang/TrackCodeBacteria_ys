classdef LEDpwm_control_system_96
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        sr= serial('COM3','baudrate',115200,'parity','none','databits',8,'stopbits',1);%serial connect
        num_LED =7;
        luminosity_intensity_matrix_example = zeros(96,3);%0-255
        luminosity_PWM_intensity_example = zeros(96,1);%0-100
%         status_matrix_example = zeros(96,1);% 0 or 1
        delay_in_interrupt_recover;
        seq;
    end
    
    methods
        function obj = turn_port_on(obj)
%             obj.sr.OutputBufferSize = 10240 ;
            fopen(obj.sr);
        end
        function obj = send_LED_matrix(obj,luminosity_intensity_matrix)
            obj.seq=char();
            obj.seq=[obj.seq,dec2hex(160,2)];
            %first byte can be any non-zero number
            for i=1:length(luminosity_intensity_matrix)
                    %status
                    for j=1:3
                        obj.seq=[obj.seq,dec2hex(luminosity_intensity_matrix(i,j),2)];
                    end
            end
            obj.seq = [obj.seq,dec2hex(13,2),dec2hex(10,2)];
            D = sscanf(obj.seq, '%2x'); %convert char into byte
            fwrite(obj.sr, D, 'uint8') %send byte to target UART port
        end
        function obj = send_PWM_intensity(obj,luminosity_PWM_intensity)
            obj.seq=char();
            obj.seq=[obj.seq,dec2hex(161,2)];
            %first byte can be any non-zero number
            for i=1:length(luminosity_PWM_intensity)
                        obj.seq=[obj.seq,dec2hex(luminosity_PWM_intensity(i,:),2)];
            end
            obj.seq = [obj.seq,dec2hex(13,2),dec2hex(10,2)];
            D = sscanf(obj.seq, '%2x'); %convert char into byte
            fwrite(obj.sr, D, 'uint8') %send byte to target UART port
        end
        function obj = send_LED_PWM_matrix(obj,luminosity_intensity_matrix)
            obj.seq=char();
            obj.seq=[obj.seq,dec2hex(162,2)];
            %first byte can be any non-zero number
            for i=1:length(luminosity_intensity_matrix)
                    %status
                    for j=1:3
                        obj.seq=[obj.seq,dec2hex(luminosity_intensity_matrix(i,j),2)];
                    end
            end
            obj.seq = [obj.seq,dec2hex(13,2),dec2hex(10,2)];
            D = sscanf(obj.seq, '%2x'); %convert char into byte
            fwrite(obj.sr, D, 'uint8') %send byte to target UART port
        end

        function obj = turn_port_off(obj)
            fclose(obj.sr);
        end
    end
end