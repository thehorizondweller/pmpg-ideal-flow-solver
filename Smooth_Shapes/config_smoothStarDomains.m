close all; clear; clc;

%{
Author: Anand K., Gonzalez C., Taha H.E.
Affiliation: UC Irivine
Last Updated: December 2024

This is the main configuration file used to interface with the conformal
mapping package. Note that this code is designed to handle smooth Jordan
Star Shaped Domains i.e. there should be no corners/cusps and there must
exist atleast one point in the interior of the geometry such that each
point on the boundary of the curve can be connected to this point without
touching or intersecting the boundary at any other point. 

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
The smooth rounded TE Kutta-Zhoukhovsky airfoil developed by Gonzalez C. &
Taha H.E. is used as an example. The user can read the descriptions of the
variables in this section and adjust them according to their geometry. Note
that the geometry must be preconditioned such that the longest axis of the
geometry is aligned with the horizontal axis (referred as the x-axis) and
the origin must lie midway on the longest axis.

The user can quickly search for "USER INPUTS" to find variables where user inputs
are required. 
%}

% Code for generating Smooth TE Kutta-Zhoukhovsky airfoil
n = 1000; % number of points on the surface
C = 1; % scaling factor
e = 0.1; % eccentricity for the x offset
beta = 0; % Input in degrees
beta = beta * pi / 180; % Convert to radians
D = 0.25; % TE Smoothing Parameter
% Grid controls
terminalR = 10; % Extent of the domain 
radialPts = 41; % Radial spacing
angularPts = 101; % Angular spacing
% To get the surface points only
geoData = kjSmoothTEAirfoilGrid(n, C, e, beta, D, radialPts, angularPts, terminalR);

% USER INPUTS: z - array of points on the contour starting at the positive
% x-axis and moving CCW on the contour.
% beta = 0; % User must provide beta = 0 for a generic shape. 
z = geoData.z;

% Get the coordinates and the scaling factor
x = real(z);
y = imag(z);
a = max(x)/2;

% Visualize the geometry - plots the user provided geometry and saves it in
% the Results folder as 01_UserGeometry.png
visualizeGeometry(z);

%% STEPS 2-4: COMPUTATION OF THE INTERMEDIATE DOMAIN
%{
The radial and angular distributions in the Theodorsen's method for the
geometry is computed here. The intermediate domain zetaPrime is
constructed. The number of points on the curve are increased by
interpolation to ensure higher accuracy and faster convergence in the
conjugation iteration. 

Two plots are saved in the ./Results - Radial vs Angular distribution and
Intermediate Domain contour. 
%}
interDomData = computeIntermediateDomain(x,y,z,a,n);

% Unzipping the data
theta = interDomData.theta; % Angular distribution
psi = interDomData.psi; % Radial distribution
n = interDomData.n; % Number of points on the contour after interpolation
zetaPrime = interDomData.zetaPrime; % contour of the intermediate domain

thetaExt = interDomData.thetaExt;
psiExt = interDomData.psiExt;

%% STEPS 5-6: COMPUTATION OF THE CONFORMAL MAP
%{
The conformal map from the intermediate domain to the circular domain is
computed here.

List of figures generated in the Results folder -
-> Boundary Correspondence Function
-> Angular Distortion Function
-> Psi (Radial Dist. of Intermediate Domain) vs Phi (Angular Dist. of Circular
    Domain)
-> All three domains - Z, Intermediate (Zeta Prime) and Circular (Zeta)
%}

% USER INPUTS
iterMax = 2e03; % Maximum number of iterations
accuracy = 1e-08; % Accuracy for convergence
q = n/200; % Max order of the Laurent Series Coefficients

conformalMapData = computeConformalMap(a,psi,theta,psiExt,thetaExt,n,z,zetaPrime,q);

A = conformalMapData.A;
B = conformalMapData.B;
psi_0 = conformalMapData.psi_0;

%% STEPS 7-12: FLOW FIELD COMPUTATION USING PMPG
%{
After, the conformal map has been successfully computed. The flow fields in
all the three domains are computed and the Appellian in the Z domain is
minimized. Images 08-18 are generated in the Results folder with
self-explanatory names. The plots in images 17 and 18 show the final
results of minimizing the appellian and the flow field for the
corresponding circulation.

The relevant results - normalized circulation and appellians for the
Kutta's and PMPG's solution are shown in the terminal. They can can be
unzipped from the flowFieldData as required by the user by uncommenting the
code at the end.
%}


% USER INPUTS
radialPts = 41;
angularPts = 101;
terminalR = 10;
U = 1; % Far Field Conditions
alpha = 5; % In degrees
kuttaNormGDev = 0.1; % Bandwidth around equivalent Kutta circulation to look for minimizing solution

flowFieldData = computeFlowField(radialPts,angularPts,terminalR,a,psi_0,U,alpha,beta,kuttaNormGDev,A,B);

% Unzipping flowFieldData
% minNormG = flowFieldData.minNormG;
% minApp = flowFieldData.minApp;
% kuttaApp = flowFieldData.kuttaApp;
% kuttaNormG = flowFieldData.kuttaNormG;