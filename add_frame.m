function framed_phi = add_frame(phi,k, value)

% Helper method to add a frame of with k into a gray-scale image
% Input: 
% ------- phi       :   a gray-scale image of size M * N
% ------- k         :   the width of the frame
% ------- value     :   the value of the boundary
% 
% Output:
% ------- framed_phi:   a grapy-scale image of size (M+k) * (N+k)



framed_phi                      = value.*ones(size(phi)+[k*2,k*2]);
framed_phi(k+1:end-k,k+1:end-k) = phi;

end

