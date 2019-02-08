%% PART A
%{
* Unobserved variable: rain
* Observable variable: umbrella
* P(Xt|Xt-1) = 0.7, P(E1|Xt) = 0.9

* Dynamic model: 
===================
    |Xt-1=t|Xt-1=f|
===================
Xt=t| 0.7  | 0.3  |
Xt=f| 0.3  | 0.7  |
===================

* Observational model: 
===================
    |Xt=t  |Xt=f  |
===================
Et=t| 0.9  | 0.1  |
Et=f| 0.2  | 0.8  |
===================

* Assumptions:
==============
-Stationary process: Probability of raining on a certain day is only 
dependant ONLY on whether or it rained the day before. 
-Sensor Markov assumption: P(Et|X1:t,E1:t-1) = P(Et|Xt) and the probability 
is constant and not change over time.
-1st order Markov process: P(Xt|X0:t-1) = P(Xt|Xt-1)

Are assumptions reasonable?
===========================
The stationary assumption is not accurate. Weather is highly dependant on 
the time of the year. 1st order Markov process modelling is too simple. 
The information about whether it rained the day before is not enough to 
accurately predict weather conditions the day after. A higher order Markov 
process could be more accurate. Sensor Markov assumption is not very good 
assumption as well. Since our sensor a person, and his habit of bringing 
an umbrella could be highly weather dependant, but it could also be
influenced by other factors. 
%}

%% PART B
T = [0.7 0.3;0.3 0.7];   % transistion model P(xt|xt-1)
Ot = [0.9 0;0 0.2];      % sensor model for Umbrella = true
Of = eye(length(Ot))-Ot; % sensor model for Umbrella = false

f0 = [0.5;0.5];          % initial rain probability

f1 = Forward(f0,T,Ot);
f2 = Forward(f1,T,Ot);
%{ 
Probability of rain at day 2:
f0 = [0.5000,0.5000]
f1 = [0.8182,0.1818]
f2 = [0.8834,0.1166]  

P(x2|e1,e2) = 0.8834
%}

f1 = Forward(f0,T,Ot);   % Umbrella
f2 = Forward(f1,T,Ot);   % Umbrella
f3 = Forward(f2,T,Of);   % No umbrella
f4 = Forward(f3,T,Ot);   % Umbrella
f5 = Forward(f4,T,Ot);   % Umbrella
%{ 
Probability of rain at day 5:
f0 = [0.5000,0.5000]
f1 = [0.8182,0.1818]
f2 = [0.8834,0.1166]   
f3 = [0.1907,0.8093]
f4 = [0.7308,0.2692]
f5 = [0.8673,0.1327]
%}

%% PART C
ev = [Ot,Ot];
Forward_Backward(ev,x0,T);
%{
sv = 
    0.8834    0.8834
    0.1166    0.1166

P(x1|e1:2) = 0.883
%}

ev = [Ot,Ot,Of,Ot,Ot];
Forward_Backward(ev,x0,T);
%{
b = 
1   0.6900  0.4593  0.0906  0.0661  0.0444    
1   0.4100  0.2437  0.1503  0.0455  0.0242

sv = 
    0.8673    0.8204    0.3075    0.8204    0.8673
    0.1327    0.1796    0.6925    0.1796    0.1327

P(x1|e1:5) = 0.866
%}

%% Functions

% PART B
% Normalize vector elements so the sum of all elements add up to 1
function normalized_vector = normalize(vector)
    normalized_vector = vector/sum(vector);
end

% Based on equation 15.12
function out = Forward(f,T,O)
    out = normalize(O*T'*f);
end

% PART C
% Based on equation 15.13
function out = Backward(b,T,O)
    out = T*O*b;
end

% Based on figure 15.4
function out = Forward_Backward(ev,prior,T)
    % ev: a vector of evidence values for n-steps
    % prior: prior distribution on the initial state P(x0)
    
    % n: number of steps    
    n = length(ev)/length(prior);
    
    % RUN FORWARD GIVEN PRIOR AND EVIDENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fv: a vector of forward messages
    fv = zeros(length(prior),n+1);
    fv(:,1) = prior;  
    
    for i = 1:n
        fv(:,i+1) = Forward(fv(:,i),T,ev(:,i*2-1:i*2));
    end
    
    % RUN BACKWARD ON FORWARD DATA 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % b: a representation of backward messages
    b = [1;1];  % initialize to 1s

    % sv: a vector of smoothed estimates
    sv = zeros(length(prior),n);

    for i=n:-1:1
        sv(:,i) = normalize(fv(:,i+1).*b);
        b = Backward(b,T,ev(:,i*2-1:i*2));
    end
    out = sv;
end







