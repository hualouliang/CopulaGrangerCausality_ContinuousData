function[out]=bootrsp(in,B)
%   out=bootrsp(in,B)
%    
%   Bootstrap  resampling  procedure. 
%
%     Inputs:
%        in - input data 
%         B - number of bootstrap resamples (default B=1)        
%     Outputs:
%       out - B bootstrap resamples of the input data  
%
%   For a vector input data of size [N,1], the  resampling 
%   procedure produces a matrix of size [N,B] with columns 
%   being resamples of the input vector.
%
%   For a matrix input data of size  [N,M], the resampling
%   procedure produces a 3D matrix of  size  [N,M,B]  with 
%   out(:,:,i), i = 1,...,B, being a resample of the input 
%   matrix.
%
%   Example:
%
%   out=bootrsp(randn(10,1),10);

%  Created by A. M. Zoubir and D. R. Iskander
%  May 1998
%
%  References:
% 
%  Efron, B.and Tibshirani, R.  An Introduction to the Bootstrap.
%               Chapman and Hall, 1993.
%
%  Zoubir, A.M. Bootstrap: Theory and Applications. Proceedings 
%               of the SPIE 1993 Conference on Advanced  Signal 
%               Processing Algorithms, Architectures and Imple-
%               mentations. pp. 216-235, San Diego, July  1993.
%
%  Zoubir, A.M. and Boashash, B. The Bootstrap and Its Application
%               in Signal Processing. IEEE Signal Processing Magazine, 
%               Vol. 15, No. 1, pp. 55-76, 1998.

if (exist('B')~=1), B=1;  end;
if (exist('in')~=1), error('Provide input data'); end;

s=size(in);     
if length(s)>2, 
  error('Input data can be a vector or a 2D matrix only'); 
end;
if min(s)==1,  
  out=in(ceil(max(s)*rand(max(s),B)));    
else         
  out=in(ceil(s(1)*s(2)*rand(s(1),s(2),B))); 
end;





