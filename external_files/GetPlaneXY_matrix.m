function [ tdom_dataX, tdom_dataY ] = GetPlaneXY_matrix(fft_data,General,ix,iy,KXi,KYj)

% loop over fft_data{k}

idx_intern = 1:length(KXi);

u_planes = zeros(General.ifftNx, General.ifftNy,length(length(fft_data)));
for i = 1:size(fft_data,1)
    
    u_vol = zeros(General.ifftNx, General.ifftNy);
    
    % Generate the top and bottom left corners
    u_vol(General.idx(KXi),General.idy(KYj)) = fft_data(i,idx_intern,KYj);
    % complex conjugate the zero'th mode and map to the top right side
    u_vol(1,General.idy_reverse(KYj(1:end-1))) = conj(fft_data(i,1,General.idy_intern_map(KYj(1:end-1))));
    % complex conjugate the remaining wave numbers and copy to the top and
    % bottom right corner
    u_vol(General.idx_reverse(fliplr(length(General.Kx) + 1 - KXi(2:end))) ,General.idy_reverse(KYj(1:end-1))) = conj( fft_data(i,idx_intern(2:end),General.idy_intern_map(KYj(1:end-1))));
    
    %  General.idx_reverse    General.idx_intern_map 66 - KXi(2:end)
    
    u_planes(:,:,i) = ifft2(u_vol(:,:));
    
end

tdom_dataX = squeeze(u_planes(ix,:,:));
tdom_dataY = squeeze(u_planes(1:(end/8*3),iy,:));

end