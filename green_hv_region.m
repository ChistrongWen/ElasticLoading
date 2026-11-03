function Udis = green_hv_region(LLZ,boundfile)
    % Get Green's functions of PREM model
    all_green = load('cf_PREM.txt');
    Gfun.theta = all_green(:,1);
    Gfun.ur2   = all_green(:,4);
    Gfun.v2    = all_green(:,5);
    % Calculate the loading deformation of TP's points
    N = length(boundfile(:,1));
    % p = parpool(10);
    for ii = 1:N
        disp(['finish_',num2str(ii),'_point!']);
        Udis(ii,:) = mass_load_hv(Gfun,LLZ,[boundfile(ii,1),boundfile(ii,2)],1e3);
    end
    % delete(p);
end


