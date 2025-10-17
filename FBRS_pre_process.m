%% 

function Data_E = FBRS_pre_process(data,frame_length,frame_stride)

len = length(data);
% t = linspace(0,(len-1)/fs,len);

%% Optional. Not used here.

% alpha = 0.97;  
% for i = 2:len
%     y(i) = data(i) - alpha * data(i-1);
% end
% data = y';

%% 

frame_num = ceil((len-frame_length)/frame_stride); 
if mod((len-frame_length), frame_stride) == 0 
    frame_num = frame_num+1;
end
S = enframe(data,frame_length,frame_stride); 


%% 
a = 0.46; 
n=1:frame_length;
W = (1-a) - a*cos((2*pi.*n)/frame_length);

C=zeros(frame_num,frame_length);
for i=1:frame_num
    C(i,:)=W;
end

SC=S.*C; 

%% 
nfft = 1024;
for i=1:frame_num
    len = length(SC(i,:));
    fftamp1 = abs(fft(SC(i,:), nfft*2)) / len;
    fftamp2 = fftamp1(1:nfft*2/2);
    Data_E(i,:) = fftamp2.^2;
end

end
