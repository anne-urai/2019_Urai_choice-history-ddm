
close;
subplot(221); scatter(dat.criterion, ...
    z_link_func(dat.z_Intercept__regressdczprevresp), 10, dat.session);
subplot(221); scatter(dat.criterion, ...
    (dat.z__stimcodingdcprevresp), 10, dat.session);

lsline;
subplot(222); scatter(dat.criterion(dat.session == 0), ...
    dat.v_Intercept__regressdczprevresp(dat.session == 0));
lsline;