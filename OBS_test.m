clear all;
%load data
D = load("datafortest.mat");

%%
data = D.egfraw; %size = (Ntimelag, Ntrace)

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
tvec = D.timeflag(:,1);

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


% figure(1);
% date = 1:1:155;
% h = pcolor(tvec, date, data');
% %datetick('y', 'yyyy-mm-dd',  'keepticks');
% set(h, 'EdgeColor', 'none');
% colormap(jet);
% colorbar;
% xlabel('Time lag');
% saveas(gcf, "allccfunction_OBS", "jpg");

fig = figure(2);
clf;
subplot(2,1,1)
hold on;
plot(tvec, X_linear, "k");
plot(tvec, X_selective, "m-");
plot(tvec, X_selective_2nd, "b--");
%plot(tvec, X_robust, "r-");
legend({'Linear', '1st selective', '2nd selective', 'Robust'});
xlabel('Time lag');
ylabel('Coherence');
xlim([-500, 500]);
box on;
subplot(2,1,2)
plot(Stats_robust.weight, "k");
xlabel('Trace id');
ylabel('Robust stack weights');
box on;
% subplot(3,1,3)
% %plot(tvec, RMS_selective);
% plot(tvec, RMS_selective_2nd);
% %plot(tvec, RMS_robust);


set(gcf, 'Units', 'Normalized', 'Position',  [0.2, 0.8, 0.4 0.5])
saveas(gcf, "stackedtrace_OBS", "jpg");


