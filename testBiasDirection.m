dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
    'projects/0/neurodec/Data/MEG-PL'));

close;
subplot(441); scatter(dat.bias, ...
    dat.z__stimcodingdcprevresp, 10, dat.session);
xlabel('P(resp 1)'); ylabel('z'); title('Stimcoding');
axis square;
subplot(442); scatter(dat.bias, ...
    (dat.dc__stimcodingzprevresp), 10, dat.session);
xlabel('P(resp 1)'); ylabel('dc'); title('Stimcoding');
axis square;

subplot(443); scatter(dat.criterion, ...
    dat.z__regressdcprevresp, 10, dat.session);
xlabel('P(resp 1)'); ylabel('z'); title('Regression');
axis square;

subplot(444); scatter(dat.bias, ...
    (dat.dc__stimcodingzprevresp), 10, dat.session);
xlabel('P(resp 1)'); ylabel('dc'); title('Regression');
axis square;

% against each other
subplot(445); scatter(dat.z__stimcodingdczprevresp, ...
    dat.dc__stimcodingdczprevresp, 10, dat.session);
xlabel('z'); ylabel('dc'); title('Stimcoding');
axis square;

% what is the correct link function?
z_link_func = @(x) log(1 ./ (1-x));
z_link_func = @(x) 1./(1+exp(-x))';

subplot(446); 
scatter(dat.v_prevresp__regressdczprevresp, ...
   dat.dc_1__stimcodingdczprevresp - dat.dc_2__stimcodingdczprevresp, 10, dat.session);
xlabel('z'); ylabel('dc'); title('Regression');
axis square;


