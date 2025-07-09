% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G2
% Author: Salvi Niccolò (10773726)

% read the provided image
image = im2double(imread('images\Look-outCat.jpg'));

l_infty = [-0.0001; -0.0011; 1.0000]; % line at infinity obtained from G1

H_aff = [eye(2), zeros(2, 1); l_infty(:)']; % affine matrix

% visualize the affine rectified image
tform_aff = projective2d(H_aff');

% load conic params obtained from F2
fileID = fopen('data\conic_C_parameters.txt', 'r');
data = textscan(fileID, '%c: %f');
fclose(fileID);
params = data{2};
a = params(1);
b = params(2);
c = params(3);
d = params(4);
e = params(5);
f = params(6);

C = [a b/2 d/2; b/2 c e/2; d/2 e/2 f];
C = C ./ C(3, 3);

syms 'x';
syms 'y';

% every conic provide a 2n ddegree equation
eq1 = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f;
eq2 = l_infty(1)*x + l_infty(2)*y + l_infty(3);
% solve a system with a line and the conic꞉ this gives a 2th degree equation
eqns = [eq1 == 0, eq2 == 0];
S12 = solve(eqns, [x,y], 'IgnoreAnalyticConstraints', true, 'Maxdegree', 4);
% hence you get 2 pairs of complex conjugate solution
s1 = [double(S12.x(1)); double(S12.y(1)); 1];
s2 = [double(S12.x(2)); double(S12.y(2)); 1];

alpha = -l_infty(2) / l_infty(1);
beta = -l_infty(3) / l_infty(1);
% substitute in eq2 to get a second degree equation in y
poli = [a*alpha*alpha + b*alpha + c, 2*a*alpha*beta + b*beta + d*alpha + e, a*beta*beta + d*beta + f];
sol = roots(poli);

y1 = sol(1);
y2 = sol(2);
x1 = alpha*y1 + beta;
x2 = alpha*y2 + beta;
q1 = [x1; y1; 1];
q2 = [x2; y2; 1];

II = s1;
JJ = s2;

imDCCP = II*JJ' + JJ*II'; % image of the circular points
imDCCP = imDCCP ./ norm(imDCCP);

[V, D] = eigs(imDCCP);
[d, ind] = sort(abs(diag(D)), 'descend');
Ds = D(ind, ind);
Vs = V(:, ind);

% Compute Q matrix and normalize
Q = inv(H_aff)' * C * inv(H_aff);
Q = Q / Q(3, 3); % Normalize Q

par_geo = AtoG([Q(1, 1), 2*Q(1, 2), Q(2, 2), 2*Q(1, 3), 2*Q(2, 3), Q(3, 3)]);
center = par_geo(1:2);
axes = par_geo(3:4);
angle = par_geo(5);

% define rotation angle
alpha = angle;

% define axes lengths
a = axes(1);
b = axes(2);

% rotation matrix
U = [
    cos(alpha), -sin(alpha); 
    sin(alpha),  cos(alpha)
];

% rescaling the axis to be unitary
S = diag([1, a/b]); % scale matrix with adjusted axes ratio
% compute the combined transformation matrix
K = U * S * U'; % rotate, scale, then reverse rotation

% extend the transformation to 3D homogeneous coordinates
A = [K, zeros(2, 1); zeros(1, 2), 1];

% Apply the transformation to H
H_met = A * H_aff;

disp("H_met");
disp(H_met);

tform_met = projective2d(H_met');
metric_rectified_image = imwarp(image, tform_met);
figure(1);
imshow(metric_rectified_image);
imwrite(metric_rectified_image, 'metric_rectified_image_ora.jpg');

% AtoG
% Convert conic section equation from ABCDEF form to {center, axes, angle}.
% Minor modifications made to the original AtoG posted by Hui Ma
%  Note: Ma claimed copyright, but that is automatically released to the 
% public upon his posting to FileExchange
%
%  ParA = [A,B,C,D,E,F]'-  parameter vector of the conic:  Ax^2 + Bxy + Cy^2 +Dx + Ey + F = 0
%  to geometric parameters  ParG = [Center(1:2), Axes(1:2), Angle]'
%
% The Angle value is in radians
%  code is:  1 - ellipse  
%            2 - hyperbola 
%            3 - parabola
%           -1 - degenerate cases  
%            0 - imaginary ellipse  
%            4 - imaginary parelle lines
function [ParG,code] = AtoG(ParA)
    tolerance1 = 1.e-10;
    tolerance2 = 1.e-20;
    % ParA = ParA/norm(ParA);
    if (abs(ParA(1)-ParA(3)) > tolerance1)
        Angle = atan(ParA(2)/(ParA(1)-ParA(3)))/2;
    else
        Angle = pi/4;
    end
    c = cos(Angle);  s = sin(Angle);
    Q = [c s; -s c];
    M = [ParA(1)  ParA(2)/2;  ParA(2)/2  ParA(3)];
    D = Q*M*Q';
    N = Q*[ParA(4); ParA(5)];
    O = ParA(6);
    if (D(1,1) < 0) && (D(2,2) < 0)
        D = -D;
        N = -N;
        O = -O;
    end
    UVcenter = [-N(1,1)/2/D(1,1); -N(2,1)/2/D(2,2)];
    free = O - UVcenter(1,1)*UVcenter(1,1)*D(1,1) - UVcenter(2,1)*UVcenter(2,1)*D(2,2);
    % if the determinant of [A B/2 D/2;B/2 C E/2;D/2 E/2 F]is zero 
    % and if K>0,then it is an empty set;
    % otherwise the conic is degenerate
    Deg =[ParA(1),ParA(2)/2,ParA(4)/2;...
         ParA(2)/2,ParA(3),ParA(5)/2;...
         ParA(4)/2,ParA(5)/2,ParA(6)];
    K1=[ParA(1),ParA(4)/2;ParA(4)/2 ParA(6)];
    K2=[ParA(3),ParA(5)/2;ParA(5)/2 ParA(6)];
    K = det(K1)+ det(K2);
    if (abs(det(Deg)) < tolerance2)
        if (abs(det(M))<tolerance2) &&(K > tolerance2)
            code = 4;  % empty set(imaginary parellel lines)
        else
            code = -1; % degenerate cases
        end
    else
        if (D(1,1)*D(2,2) > tolerance1)
            if (free < 0)
                code = 1; % ellipse
            else
                code = 0; % empty set(imaginary ellipse)
            end
        elseif (D(1,1)*D(2,2) < - tolerance1)
            code = 2;  % hyperbola
        else
            code = 3;  % parabola
        end
    end
    XYcenter = Q'*UVcenter;
    Axes = [sqrt(abs(free/D(1,1))); sqrt(abs(free/D(2,2)))];
    if code == 1 && Axes(1)<Axes(2)
        AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
        Angle = Angle + pi/2;
    end
    if code == 2 && free*D(1,1)>0
        AA = Axes(1); Axes(1) = Axes(2); Axes(2) = AA;
        Angle = Angle + pi/2;
    end
    % some people never learn...
    Angle = mod(Angle,pi);
    % while Angle > pi
    %     Angle = Angle - pi;
    % end
    % while Angle < 0
    %     Angle = Angle + pi;
    % end
    ParG = [XYcenter; Axes; Angle];
end