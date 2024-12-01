close all; clear; clc;

%{
Author: Anand K., Gonzalez C., Taha H.E.
Affiliation: UC Irivine
Last Updated: December 2024

This is the main configuration file used to interface with the conformal
mapping package. Note that this code is designed to handle Jordan
Star Shaped Domains with exactly one point which is singular in the first
derivative of the conformal map. Star shaped domains are domains in the
interior of which there exist atleast one point such that each point on the
boundary of the curve can be connected to this point without touching or
intersecting the boundary at any other point. 

Code Usage instructions: Refer to description of each section below. 

%}

%% Configure the results folder
%{ Creates a Results folder in the working directory to save the images /
%plots published by the code in the solving process. If a Results folder
%was already present, it will be wiped clean each time this code this run.
%}
makeResultsFolder();

%% STEP 1: CONFIGURE THE GEOMETRY
%{
NACA airfoils have been used as the standard example. The code used to
generate NACA airfoils is attributed to the NACA 4 digit airfoil generator
published by Divahar J. at https://www.mathworks.com/matlabcentral/fileexchange/19915-naca-4-digit-airfoil-generator

In general, is required to configure their custom geometry and feed it to
variable z below. The shape must be configured such that the longest axis 
of the body must be aligned with the x axis. The shape need not be centered
about the origin. However, the author suggests placing the left end of the
longest axis on the origin for easy debugging of the Karman Trefftz map. 

The user can quickly search for "USER INPUTS" to find variables where user inputs
are required. 
%}

% Refer the MATLAB article referred to in the document section above to
% obtain further instructions on using naca4gen.m
iaf.designation='2412';
iaf.n=200;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0; % 1 for Non-joined TE
af = naca4gen(iaf);
% Ensuring TE and LE points are not repeated
z = [af.xU(2:end)'+1i*af.zU(2:end)' af.xL(2:end)'+1i*af.zL(2:end)']; % End Points on the Upper and Lower surfaces are repeated
[~,n] = size(z);

% USER INPUTS: z - array of points on the contour starting at the positive
% x-axis and moving CCW on the contour. 
% z = ; % Comment out the NACA example and enter your geometry here

% Get the coordinates
x = real(z);
y = imag(z);
a = max(x)/2;

% Visualize the geometry - plots the user provided geometry and saves it in
% the Results folder as 01_UserGeometry.png
visualizeGeometry(z);

%% STEP 2-6: COMPUTE THE CONFORMAL MAP
%{
The user provides z2 (the second point in the z domain), zetaPrime1 and
zetaPrim2 (both points in the zeta Prime domain). The code identifies the
corner and turning angle at the corner. The turning angle is shown in the
terminal as output to verify if it is the intended turning angle. If it is
not, the user needs to refine the contour by increasing the density of
points near the corner point. 

List of figures generated -
-> Variation of turning angle along the contour points 
-> Intermediate Domain Contour
-> Radial vs Angular Distribution (Intermediate Domain)
-> Boundary Correspondence Function 
-> Angular Distortion Function 
-> Psi vs Phi
-> All Three Domains
%}

% Empirical Guess for z2: Let z3 be point where the angle bisector of the
% corner point intersects with the contour. z2 should be located about 0.01*D
% away from z3 along the angle bisector and be inside the contour. D is the
% total length from the corner point to z3 along the angle bisector. 

% USER INPUTS
% Example for NACA 2412
z2 = 0.01;
zetaPrime1 = 1;
zetaPrime2 = 0;
iterMax = 1e04; % Nax iterations for Conjugation Method (Theodorsen)
accuracy = 1e-08; % Accuracy threshold for Conjugation Method
q = 100; % Max order of Fourier Coefficients

interDomData = computeConformalKTMap(z2,zetaPrime1,zetaPrime2,n,a,z,iterMax,accuracy,q);

% Unzip the outputs
minAngle = interDomData.minAngle;
z1 = interDomData.z1;
beta = interDomData.beta;
A = interDomData.A;
B = interDomData.B;
psi_0 = interDomData.psi_0;
zetaPrimeOffset = interDomData.zetaPrimeOffset;
n = interDomData.n;
zRetrace = interDomData.zRetrace;
zeta = interDomData.zeta;

% Print the results
fprintf("The corner point in the Z Domain is (%4.2f,%4.2f). \n ",real(z1),imag(z1));
fprintf("The turning angle at the corner point is %4.2f degrees. \n",beta*180/pi)

%% STEPS 7-12: FLOW FIELD COMPUTATION USING PMPG
%{
After, the conformal map has been successfully computed. The flow fields in
all the three domains are computed and the Appellian in the Z domain is
minimized. Images 09-19 are generated in the Results folder with
self-explanatory names. The plots in images 18 and 19 show the final
results of minimizing the appellian and the flow field for the
corresponding circulation.

The relevant results - normalized circulation and appellians for the
Kutta's and PMPG's solution are shown in the terminal. They can can be
unzipped from the flowFieldData as required by the user by uncommenting the
code at the end.
%}

% USER INPUTS
% Note: In shapes with corners, it is essential to capture the diminishing
% of velocity at the corner. Hence, the grid finess is crucial to the
% correct calculation of the Appellian. The autor suggests using atleast 10
% times more density in the redial direction and 5 times more density in
% the angular direction as compared to the smooth Jordan curve counterpart.
radialPts = 101;
angularPts = 501;
terminalR = 10;
U = 1; % Far Field Conditions
alpha = 5; % In degrees
kuttaNormGDev = 0.05; % Bandwidth around equivalent Kutta circulation to look for minimizing solution

flowFieldData = computeKTFlowField(radialPts,angularPts,terminalR,a,psi_0, ...
    zetaPrimeOffset,zetaPrime1,zetaPrime2,z1,z2,U,alpha,beta,kuttaNormGDev,A,B,n,zRetrace,zeta);

% Unzipping flowFieldData
% minNormG = flowFieldData.minNormG;
% minApp = flowFieldData.minApp;
% kuttaApp = flowFieldData.kuttaApp;
% kuttaNormG = flowFieldData.kuttaNormG;


