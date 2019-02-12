function [p,ellipse]=lei_Simulate_Space_myphantom(varargin)
%PHANTOM Create head phantom image.
%   P = PHANTOM(DEF,N) generates an image of a head phantom that can   
%   be used to test the numerical accuracy of RADON and IRADON or other  
%   2-D reconstruction algorithms.  P is a grayscale intensity image that
%   consists of one large ellipse (representing the brain) containing
%   several smaller ellipses (representing features in the brain).
%
%   DEF is a string that specifies the type of head phantom to generate.
%   Valid values are: 
%         
%      'Shepp-Logan'            A test image used widely by researchers in
%                               tomography
%      'Modified Shepp-Logan'   (default) A variant of the Shepp-Logan phantom
%                               in which the contrast is improved for better  
%                               visual perception.
%
%   N is a scalar that specifies the number of rows and columns in P.
%   If you omit the argument, N defaults to 256.
% 
%   P = PHANTOM(E,N) generates a user-defined phantom, where each row
%   of the matrix E specifies an ellipse in the image.  E has six columns,
%   with each column containing a different parameter for the ellipses:
%   
%     Column 1:  A    the additive intensity value of the ellipse
%     Column 2:  a    the length of the horizontal semi-axis of the ellipse 
%     Column 3:  b    the length of the vertical semi-axis of the ellipse
%     Column 4:  x0   the x-coordinate of the center of the ellipse
%     Column 5:  y0   the y-coordinate of the center of the ellipse
%     Column 6:  phi  the angle (in degrees) between the horizontal semi-axis 
%                     of the ellipse and the x-axis of the image        
%
%   For purposes of generating the phantom, the domains for the x- and 
%   y-axes span [-1,1].  Columns 2 through 5 must be specified in terms
%   of this range.
%
%   [P,E] = PHANTOM(...) returns the matrix E used to generate the phantom.
%
%   Class Support
%   -------------
%   All inputs must be of class double.  All outputs are of class double.
%
%   Remarks
%   -------
%   For any given pixel in the output image, the pixel's value is equal to the
%   sum of the additive intensity values of all ellipses that the pixel is a 
%   part of.  If a pixel is not part of any ellipse, its value is 0.  
%
%   The additive intensity value A for an ellipse can be positive or negative;
%   if it is negative, the ellipse will be darker than the surrounding pixels.
%   Note that, depending on the values of A, some pixels may have values outside
%   the range [0,1].
%    
%   Example
%   -------
%        ph = phantom(256);
%        figure, imshow(ph)
%
%   See also RADON, IRADON, FANBEAM, IFANBEAM.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.13.4.5 $  $Date: 2006/06/15 20:09:16 $

%   References: 
%      A. K. Jain, "Fundamentals of Digital Image Processing", p. 439.
%      P. A. Toft, "The Radon Transform, Theory and Implementation" (unpublished
%      dissertation), p. 199.

n = varargin{1};            % a scalar is the image size
flag=varargin{2};
ellipse = modified_shepp_logan(flag);

p = zeros(n);

xax =  ( (0:n-1)-(n-1)/2 ) / ((n-1)/2); 
xg = repmat(xax, n, 1);   % x coordinates, the y coordinates are rot90(xg)

for k = 1:size(ellipse,1)    
   asq = ellipse(k,2)^2;       % a^2
   bsq = ellipse(k,3)^2;       % b^2
   phi = ellipse(k,6)*pi/180;  % rotation angle in radians
   x0 = ellipse(k,4);          % x offset
   y0 = ellipse(k,5);          % y offset
   A = ellipse(k,1);           % Amplitude change for this ellipse
   x=xg-x0;                    % Center the ellipse
   y=rot90(xg)-y0;  
   cosp = cos(phi); 
   sinp = sin(phi);
   idx=find(((x.*cosp + y.*sinp).^2)./asq + ((y.*cosp - x.*sinp).^2)./bsq <= 1); 
   p(idx) = p(idx) + A;
end

      

function toft=modified_shepp_logan(flag)
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%     Column 1:  A    the additive intensity value of the ellipse
%     Column 2:  a    the length of the horizontal semi-axis of the ellipse 
%     Column 3:  b    the length of the vertical semi-axis of the ellipse
%     Column 4:  x0   the x-coordinate of the center of the ellipse
%     Column 5:  y0   the y-coordinate of the center of the ellipse
%     Column 6:  phi  the angle (in degrees) between the horizontal semi-axis 
%                     of the ellipse and the x-axis of the image    
%         A    a     b    x0    y0    phi
%        ---------------------------------

switch flag
    case 0
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*6/7)  .08 .15  .72    0    -18%听觉
        -round(256*6/7)  .08 .15  -.72    0     18%听觉
        -round(256*5/7)  .11 .08  -.3    0.6     0%左前叶
        -round(256*4/7)  .14 .12   0    .2      0%运动
        -round(256*3/7)  .07 .1   0    .6    0%default
        -round(256*3/7)  .1 .15     0   -.45      0%default
        -round(256*3/7)  .04 .07  .52    -.55    -38%default
        -round(256*3/7)  .04 .07  -.52    -.55    38%default
        -round(256*2/7)  .09 .09  .5    0.5    0%右前叶
        -round(256*1/7)  .25 .05    0    -.75   0%视觉
       ];     
    case 1
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*1/7)  .25 .05    0    -.75   0%视觉
       ];           
    case 2
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*3/7)  .07 .1   0    .6    0%default
        -round(256*3/7)  .1 .15     0   -.45      0%default
        -round(256*3/7)  .04 .07  .52    -.55    -38%default
        -round(256*3/7)  .04 .07  -.52    -.55    38%default
       ];  
       case 3
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*6/7)  .08 .15  .72    0    -18%听觉
        -round(256*6/7)  .08 .15  -.72    0     18%听觉
       ];   
       case 4
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*4/7)  .14 .12   0    .2      0%运动
       ];   
       case 5
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*5/7)  .11 .08  -.3    0.6     0%左前叶
       ];  
       case 6
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*2/7)  .09 .09  .5    0.5    0%右前叶
       ];  
       case 7%fMRI的空间模式
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
        -round(256*6/7)  .08 .15  .72    0    -18%听觉
        -round(256*6/7)  .08 .15  -.72    0     18%听觉
        -round(256*5/7)  .11 .08  -.3    0.6     0%左前叶
        -round(256*4/7)  .14 .12   0    .2      0%运动
        -round(256*3/7)  .07 .1   0    .6    0%default
        -round(256*3/7)  .1 .15     0   -.45      0%default
        -round(256*3/7)  .04 .07  .52    -.55    -38%default
        -round(256*3/7)  .04 .07  -.52    -.55    38%default
%         -round(256*2/7)  .09 .09  .5    0.5    0%右前叶
        -round(256*1/7)  .25 .05    0    -.75   0%视觉
       ]; 
      case 8%只有白质和灰质，显示结果时使用
        toft = [255   .88    .89    0     0     0 
            -1   .85    .85    0     0     0 
        1  .1 .35000  .32    0    0%白质
        1  .1 .3500  -.32    0     0%白质
       ]; 
end     
