%% MXB103 Project Group 69:  BUNGEE!

%% 1  Introduction
% As part of Brisbane's "New World City" transformation, the Brisbane City
% Council is investigating a proposal to allow bungee jumping off the
% Story Bridge. This report addresses several key questions about the
% proposal.
%
% In Section 2 of this report, the proposal and the key questions relating 
% to it are discussed in more detail.  In Section 3, the model ...
% In Section 4, etc.

%% 2  The proposal
% What is the proposal exactly?
% What questions have you been asked to answer?
% What's your approach to doing this?
% A few paragraphs should be enough.
% Finish with a sentence linking to the next section on the model.

%% 3 The model
%
% The equation of motion for bungee jumping is
%%
% $$\frac{dv}{dt} = g - C|v|v - \max(0,K(y-L))$$
%
% Write a paragraph or two about the equation, including where it comes
% from and what the symbols mean.

%% 3.1 Assumptions and limitations
%
% Write a paragraph on the assumptions and limitations on the model,
% and a paragraph on the limitations on the study you have conducted.

%% 3.2 Parameters
H = 74;             % Height of jump point (m)
D = 31;             % Deck height (m)
c = 0.9;            % Drag coefficient (kg/m)
m = 80;             % Mass of the jumper (kg)
L = 25;             % Length of bungee cord (m)
k = 90;             % Spring constant of bungee cord (N/m)
g = 9.8;            % Gravitational acceleration (m/s^2)
C = c/m;            % Scaled drag coefficient
K = k/m;            % Scaled spring constant

%% 4 The numerical method
%
% How do you formulate the model to solve numerically? (hint: you write it
% as two equations)

%% 4.1 Parameters

T = 60;             % Final time in simulation (s)
n = 10000;          % Number of subintervals (you decide how many you need)

%% 4.2 Solution
%
% The ordinary differential equations are solved using the Classical Fourth Order Runge-Kutta method.
[t, y, v, h] = rk4_bungee(T, n, g, C, K, L);

figure
plot(t, y);
xlabel('time (s)');
ylabel('distance fallen (m)');
title('Figure 1: Distance between Jumper and Platform.');

%% 5 Analysis
%
% In this section, the model predictions are analysed with respect to the
% key questions being asked about the proposal.
%

%% 5.1 Timing and bounces
%
% The bungee jump company suggests that the standard jump will consist of
% 10 "bounces" which should take approximately 60 seconds.
% Our model, Figure 1, does agree with the timing.
% 

%% 5.2 Maximum speed experienced by the jumper
%
% The "thrill factor" of bungee jumping is partly determined by the maximum
% speed experienced by the jumper.
%
% We can find the maximum velocity by:
[t, y, v, h] = rk4_bungee(T, n, g, C, K, L);

figure
plot(t,v);
xlabel('time (s)');
ylabel('velocity (m/s)');
title('Figure 2: Velocity of Jumper Vs Time');

% According to the graph, it seems that the maximum velocity achieved by
% the jumper is approximately 20 m/s.

%% 5.3 Maximum acceleration experienced by the jumper
%
% Describe the question, and then answer it.  In this case, you will call
% a function to numerically differentiate $v$ to find $a$.  Then plot $a$ 
% versus $t$ and refer to it to answer the question.

%% 5.4 Distance travelled by the jumper
%
% Describe the question, and then answer it.  In this case, you will call
% a function to numerically integrate $|v|$ to find the distance.

%% 5.5 Automated camera system
%
% Describe the question, and then answer it.  In this case, you will
% fit an interpolating polynomial through the four points in your solution
% $y$ that lie either side of the camera location.  Then use that
% polynomial to solve for when the jumper passes the camera.

%% 5.6 Water touch option
%
% Describe the question, and then answer it.  In this case, you will
% re-solve the equations with different parameters for $L$ and $k$.  
% Experiment to find which values work best for the water touch option, but 
% include only the best combination that you found in the submitted code.

%% 6 Conclusion
%
% Conclude your report by summing up your findings and making any
% recommendations.