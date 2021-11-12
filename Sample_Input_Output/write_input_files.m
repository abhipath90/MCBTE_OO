% This is to write input files for the test case 1
% These are the test cases that have a proper simulation with their
% description

%{ 
Problem decription goes here
100nm cube with x and y boundaries periodic z boundaries isothermal
Steady state simulation with wall temp 305 and 295 K 
equilibrium temperature 300K

%}


%% Define parameters
ds = 100e-9;
t = 100e-9;
N = 10000000;
max_scat = 10000;
Teq = 300;
% Interface conditions
Trans_type = 5; % transmission probability: 1 for low pass, 2 for high pass, 3 for normal incidence and 4 for grazing incidence 5 for constant
Spec_type = 1; % specularity type: 1 for fixed and 2 for ziman specular
% Rest of the interface conditions are to be defined directly in interface_implement.m

FF_list = 100; % filling fraction in percentage
L_list = 100;
Ly = 100e-9; % through the thickness of paper

for mm=1:length(FF_list)

    for ll = 1:length(L_list)
        L = 1e-9*L_list(ll);
        P = L*100/FF_list(mm);
        
%         dir_name = [num2str(FF_list(mm)) '_' num2str(1e9*L)];
%         % write something to make directory and then copy standard files into it
%         mkdir(dir_name);
%         
%         copyfile('*.txt', dir_name);
%         copyfile('*.m', dir_name);
%         
%         cd(dir_name);
        
        %% Boundary.txt
        % Defining points to construct Boundary.txt
        Pa = [0 0 0]; Pb = [0 0 ds]; Pd = [(P-L)/2 0 ds+t];
        Pe = [(P+L)/2 0 ds+t];  Pg = [P 0 ds]; Ph = [P 0 0];
        Pi = [P Ly 0]; Pj = [P Ly ds];  Pl = [(P+L)/2 Ly ds+t];
        Pm = [(P-L)/2 Ly ds+t]; Po = [0 Ly ds]; Pp = [0 Ly 0];
        
        mat_bnd = [Pa Pb Pg Ph;
            Ph Pg Pj Pi;
            Pi Pj Po Pp;
            Pb Pa Pp Po;
            Pa Ph Pi Pp;
            Pd Pe Pg Pb;
            Po Pj Pl Pm;
            Pd Pm Pl Pe;
            Pb Pg Pj Po;
            Po Pm Pd Pb;
            Pe Pl Pj Pg];
        
        writematrix(mat_bnd, 'Boundary.txt');
        
        %% Boundary_prop.txt
        % BndID, type, matID, otherProperties 
        % matID 1 for substrate, 2 for heater
        % for interface the other material is defined as first property. 
        % in the reverse direction
        mat_bnd_prop = [1 3 1 0 Ly 0 ;
            2 3 1 -P 0 0 ;
            3 3 1 0 -Ly 0 ;
            4 3 1 P 0 0 ;
            5 1 1 295 0 0 ;
            6 3 2 0 Ly 0 ;
            7 3 2 0 -Ly 0 ;
            8 1 2 305 0 0 ;
            9 5 2 1 0 0 ;
            10 3 2 P 0 0 ;
            11 3 2 -P 0 0 ];
        
        writematrix(mat_bnd_prop,'Boundary_prop.txt');
        
        
        %% Periodic boundary pairs
        % list periodic bounday pairs. For each periodic boundary list its
        % conjugate ID. Leave blank if there is not periodic boundary
        mat_peri_pair = [1 3;
                         2 4;
                         3 1;
                         4 2;
                         6 7;
                         7 6;
                         10 11;
                         11 10];
       writematrix(mat_peri_pair, 'PeriBnd_pairs.txt');
        %% Measure_regions.txt
        % Detectors also contain information about the material and no
        % refinement information
        % Region1 Si: matID 2
        N_heater = 25;
        xmin = (P-L)/2; xmax = (P+L)/2; % xmin = P/2 - 2.5e-9; xmax = P/2 + 2.5e-9;
        ymin = 0; ymax = Ly;
        zmin = ds; zmax = ds+t;
        ZZ = linspace(zmin,zmax,N_heater+1);
        Heater_measure = zeros(N_heater,7);
        for heater_index=1:N_heater
            
            Heater_measure(heater_index,:) = [xmin xmax ymin ymax...
                ZZ(heater_index) ZZ(heater_index+1) 2];
        end
        
        %Region2 - substrate Ge center: matID1
        N_subs_top = 25;
        xmin = 0; xmax = P; %xmin = P-5e-9; xmax = P;
        ymin = 0; ymax = Ly;
        zmin = 0; zmax = ds;
        
        ZZ = linspace(zmin,zmax,N_subs_top+1);
        Subs_cent_measure = zeros(N_subs_top,7);
        for subs_index=1:N_subs_top
            
            Subs_cent_measure(subs_index,:) = [xmin xmax ymin ymax...
                ZZ(subs_index) ZZ(subs_index+1) 1];
        end
        %Region2 Ge
        
        mat_measure_reg = [Heater_measure; Subs_cent_measure];
        
        writematrix(mat_measure_reg,'Measure_regions.txt');
        
       
        %%
        % Sim_param.txt
        volume = P*Ly*ds + L*Ly*t;
        mat_sim_param = [N; max_scat; volume; Teq; Trans_type; Spec_type];
        
        writematrix(mat_sim_param,'Sim_param.txt');
        
        
        % Writing the initial conditions to the file in text format
        %         Generating initial conditions data on a grid based on whatever method we choose.
        %         finally the cells are to be reported in the same way as detectors because we
        %         intend to imploy same sampling mechanism for particle emmission as other cases
        %
        % Domain extent for mat1 (substrate)
        xmin=0; xmax=P; ymin=0; ymax=Ly; zmin=0; zmax=ds;
        
        Sample_points = [];%[xmin xmax ymin ymax zmin zmax -0.5];
        
        writematrix(Sample_points,'Initial_conditions1.txt');
        
        
        % Domain extent for mat2 (heater)
        xmin=(P-L)/2; xmax=(P+L)/2; ymin=0; ymax=Ly; zmin=ds; zmax=ds+t;
        
        Sample_points = [];%[xmin xmax ymin ymax zmin zmax 0.5];
        
        writematrix(Sample_points,'Initial_conditions2.txt');
        
        cd ..
        
    end
end


