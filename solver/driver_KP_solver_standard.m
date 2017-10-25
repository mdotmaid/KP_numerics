% Driver script for running the KP equation solver
% KP_solver_periodic.m
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 1; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver

%% Numerical Parameters
tmax     = 50;    % Solver will run from t=0 to t=tmax
xmax     = 250;   % Solver will solve on domain z=0 to z=zmax
ymax     = 250;   % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax);           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dxinit =  1/10;     % Spatial Discretization 
dyinit =  1/10;     % KP Solver not valdiated yet                  
Nx       = round(xmax/dxinit);
Ny       = round(ymax/dyinit);
if periodic
    dx       = xmax/Nx;    % Spatial  discretization
	dy       = ymax/Ny;    % Spatial  discretization
else
    dx       = xmax/(Nx+1);    % Spatial  discretization
    dy       = ymax/(Ny+1);    % Spatial  discretization
end
%% KP PARAMETERS
epsilon = 1; % u_xxx scaling factor
lambda  = 1; % Sign on the mixed term (+1 is KP2, -1 is KP1)

%% PDE Initial and Boundary Conditions
u0 = @(X,Y) zeros(size(X));
    ic_type = '';
maindir = '/Users/dhl/Documents/MATLAB';
%% Create directory run will be saved to
data_dir = [maindir,'/Numerics/KP/',...
            '_tmax_',   num2str(round(tmax)),...
            '_xmax_',   num2str(round(xmax)),...
            '_Nx_',     num2str(Nx),...
            '_ymax_',   num2str(round(ymax)),...
            '_Ny_',     num2str(Ny),...
            '_epsilon_',num2str(epsilon),...
            '_lambda_', num2str(lambda),...
            '_bndry_condns_','periodic',...
            '_init_condns_',ic_type,...
            '/'];
% Create the data directory if necessary
if ~exist(data_dir,'dir')
    mkdir(data_dir);
else
    disp(['Warning, directory ',data_dir]);
    disp('already exists, possibly overwriting data');
end

savefile = sprintf('%sparameters.mat',data_dir);

%% If chosen, run the solver using the parameters and conditions above
if save_on
    % Load initial data
      xplot  = dx*[1:Nx];
      yplot  = dy*[1:Ny];
      [XPLOT,YPLOT] = meshgrid(xplot,yplot);
      tplot  = linspace(0,tmax,floor(tmax*10));
      u_init = u0(XPLOT,YPLOT);
    if plot_on
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        if ~periodic
            subplot(3,1,1);
        end
            contourf(XPLOT,YPLOT,u_init,'edgecolor','none'); xlabel('x'); ylabel('y'); 
            title('Initial Conditions');
%         if ~periodic
%         subplot(3,1,2);
%             plot(tplot,g0(tplot),tplot,g1(tplot));
%             xlabel('t'); ylabel('A(t)'); title('Boundary Conditions');
%             legend('At z=0','At z=zmax');
%         subplot(3,1,3);
%             plot(tplot,dg0(tplot),tplot,dg1(tplot));
%             xlabel('t'); ylabel('A(t)'); title('Derivs of Boundary Conditions');
%             legend('At z=0','At z=zmax');
%         end
        set(gca,'fontsize',fontsize,'fontname','times');
        pause(0.25);
        if check_IC
            legend(ic_type);
            drawnow; 
            return;
        end
    end
    
    % Save parameters
        save(savefile,'t','Nx','dx','xmax',...
                          'Ny','dy','ymax',...
                          'u0','periodic','epsilon','lambda');
    % Run timestepper
        KP_solver_periodic( t, xmax, Nx,...
                               ymax, Ny,...
                               epsilon, lambda,...
                               u0, data_dir );     
else
    load(savefile);
end

% %% If chosen, plot data associated with the parameters and conditions above
% if plot_on
%     disp('Calculating maximum time increment in saved data files...');
%     for ii=1:length(t)+1
%         [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
%         if fid == -1 % File does not exist
%             tm = ii-1;
%             disp(['Maximum time = ',num2str(t(tm))]);
%             break;
%         end
%         fclose(fid);
%     end
%     if ii == length(t)
%         disp(['Maximum time = ',num2str(t(tm))]);
%     end
%     % Get rid of larger t values
%     t = t(1:tm);
%     if periodic
%         xplot  = dx:dx:xmax;
%     else
%         xplot  = dx:dx:xmax-dx;
%     end
%     A_full = zeros(length(t)-1,length(xplot));
%     % Load first time step
%     load(strcat(data_dir,num2str(0,'%05d')),'A_init');
%         fontsize = 12;
%         fig=figure(2); clf;
%         plot(xplot,A_init);
%         qaxis = axis;
%         title(['Time: ',num2str(t(1))]);
%         set(gca,'fontsize',fontsize,'fontname','times');
%         drawnow
% %         input('Return');
%     %Plot subsequent time steps
%     for tind=2:tm
%         load(strcat(data_dir,num2str(tind,'%05d')),'A','tnow');
%         fontsize = 12;
%         fig=figure(2); clf;
%         plot(xplot,A);
%         A_full(tind-1,:) = A;
%         axis(qaxis)
%         hold off;
%         title(['Time: ',num2str(tnow)]);
%         set(gca,'fontsize',fontsize,'fontname','times');
%         drawnow;
%     end
% 
% end
