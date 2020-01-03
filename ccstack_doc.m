%% Stacking cross-correlation traces
% Stack cross-correlation functions with variety of methods
% using MATLAB(R).


%% Arguments
% 
%   method::String      : choose "linear", "selective", "robust".
%   A::Array            : 2D Array of traces [Ntimelag, Ntraces].     %
% 
% *Options:*
% 
%   ref::Vector         : reference trace used for selective stack [Ntimelag].
%   ccthreshold::Float  : threshold of selective stack.
%   eps::Float          : Threshold for convergence of robust stack.
%   maxiter::Int        : Maximum number of iterations to converge to robust stack.
%   v::Int              : if v = 1, output process status.
% 
%% Return
%   X::Vector           : stacked trace [Ntimelag]
%   Stats::Struct       : status associated with stacking method.
%% Usage
A = rand(101,3);
reference = rand(101,1);
%%%
% * linear stack
[X_linear, Stats_linear] = ccstack("linear", A);
%%%
% * selective stack
[X_selective, Stats_selective] = ccstack("selective", A, "ref", reference);
%%%
% * robust stack
[X_robust, Stats_robust] = ccstack("robust", A, "eps", 1e-4, "maxiter", 100);
%%

% figure()
% hold on;
% plot(X_linear);
% plot(X_selective);
% plot(X_robust);
% legend({'linear','selective', 'robust'}) 