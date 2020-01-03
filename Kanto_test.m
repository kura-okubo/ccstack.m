clear all;
%load data
D = load("E.ABHM_E.AYHM.mat");

data = D.data'; %size = (Ntimelag, Ntrace)
timestamp = datetime(D.timestamp(:),'ConvertFrom','epochtime','Format','dd-MM-yyyy HH:mm:ss');

%%%
% linear stack
[X_linear, Stats_linear] = ccstack("linear", data);

%%%
% selective stack
reference = X_linear;
[X_selective, Stats_selective] = ccstack("selective", data, "ref", reference);

% 2nd iterated selective stack
reference = X_selective;
[X_selective_2nd, Stats_selective_2nd] = ccstack("selective", data, "ref", reference);

%%%
% robust stack
[X_robust, Stats_robust] = ccstack("robust", data, "eps", 1e-4, "maxiter", 100);

% plot result
fs = 1.0; %[Hz]
npts = size(data, 1);
halfnpts = (npts-1)/2;
tvec = -halfnpts/fs:1/fs:halfnpts/fs;

RMS_selective = rms(X_selective - X_linear, 2);
RMS_selective_2nd = rms(X_selective_2nd - X_linear, 2);
RMS_robust = rms(X_robust - X_linear, 2);


set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1)
set(0,'defaulttextinterpreter','none')


figure(1);
h = pcolor(tvec, datenum(timestamp), data');
datetick('y', 'yyyy-mm-dd',  'keepticks');
set(h, 'EdgeColor', 'none');
colormap(jet);
colorbar;
xlabel('Time lag');
saveas(gcf, "allccfunction_Kanto", "jpg");

fig = figure(2);
clf;
subplot(2,1,1)
hold on;
plot(tvec, X_linear, "k");
plot(tvec, X_selective, "g--");
plot(tvec, X_selective_2nd, "b--");
plot(tvec, X_robust, "r-");
legend({'Linear', '1st selective', '2nd selective', 'Robust'});
xlabel('Time lag');
ylabel('Coherence');
box on;
subplot(2,1,2)
plot(Stats_robust.weight, "k");
xlabel('Trace id');
ylabel('Robust stack weights');
box on;
% subplot(3,1,3)
% plot(tvec, RMS_selective);
% plot(tvec, RMS_selective_2nd);
% plot(tvec, RMS_robust);

% 
set(gcf, 'Units', 'Normalized', 'Position',  [0.2, 0.8, 0.4 0.5])
saveas(gcf, "stackedtrace_Kanto", "jpg");


