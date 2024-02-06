source('mystartdefaults.m');
pos = get(groot,'DefaultFigurePosition');  #To get default matplotlib graphics
tic

# GRAPHICAL OUTPUT COMMANDS (To be solved one at a time)
# ***(Note: Please set only one graphical outcome to 'yes' and the rest to 'no'. 
#     Incorrect settings may lead to unintended executions or errors.
#     Ensure that only the desired graphical output command is set to 'yes'.****


    
Part1_Question2_plot = 'yes';      % Takes nearly 3 min for execution
Part1_Question4_plot = 'no';
Part1_Question5_plot = 'no';
Part2_Question1_plot = 'no';
Part2_Question3_plot = 'no';      % Takes nearly 7 min for execution
Part2_Question4_plot = 'no';
Part2_Question5_plot = 'no';      % Takes nearly 10 min for execution


#Declaring units
recipunit = 1.0E+10;
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;
#-----------------------------------------------
# PART-A
#-----------------------------------------------
#-----------------------------------------------
#Question 1
#-----------------------------------------------
#Energy E0
Ea = 0;
# To calculate Eb
%E=hbar^2*k^2/2
lambda_min = 1;
k=2*pi/lambda_min;
Eb=ceil((hbar*recipunit)^2*k^2/(2*elm)/qel);

#-----------------------------------------------
#Question 2
#-----------------------------------------------
#space Discretization inside potentials

#Inside U1
x1_min = 0;
x1_max = 80;
n = 160;

#Inside U2
step = (x1_max-x1_min)/n;

U = zeros(n,1);    #Potential barrier;
x1 = zeros(n,1);

#Energy E0
E_min = 0; % in eV
E_max = 0.3; % in eV
deltaE = 0.0005;
divide_E = round((E_max-E_min)/deltaE);
for i = 1:divide_E
    E0(i) = deltaE/2 + E_min+(i-1) * deltaE;
end
Ref = zeros(length(E0),1);
Tra = zeros(length(E0),1);
Ab = zeros(length(E0),1);


if strcmp(Part1_Question2_plot, 'yes') == 1
  tic
    for i = 1:n
        x1(i) = x1_min + step/2 + (i-1)*step;
        if x1(i) <= 15
            U(i) = 0.2;
        elseif x1(i) >= 65 && x1(i) <= 80
            U(i) = 0.2;
        else
            U(i) = 0;
        end
        V(i) = U(i)/ekinscale;
    end
   
    for i = 1:length(E0)
        k_0(i) = sqrt(E0(i)/ekinscale);
        lmbda(i) = 2*pi/k_0(i);
    end

    m = 2;
    x = zeros(m, 1);
    Phis = zeros(m, 1);
    Phi = zeros(m, 1);
    x(1) = x1_min - 1;
    x(2) = x1_max + 1;

    fprintf('    n           k             E0                R             T              R+T\n');
    for k = 1:length(E0)
        if E0(k) == 0
            Ref(k) = 1;
            Tra(k) = 0;
            RpT(k) = 1;
            Ab(k) = 0;
        else
            Phi0p = exp(1i*k_0(k)*x1);
            G0 = zeros(n, n);
            for j = 1:n
                for i = 1:n
                    G0(i, j) = step*exp(1i*k_0(k)*abs(x1(i)-x1(j)))/(2i*k_0(k));
                endfor
            endfor

            T = eye(n, n) - G0*diag(V);
            Phip = T\Phi0p;

            for i = 1:m
                Phis(i) = 0;
                for j = 1:n
                    Phis(i) += step*(exp(1i*k_0(k)*abs(x(i)-x1(j))))/(2*1i*k_0(k)) * V(j) * Phip(j);
                endfor
                Phi(i) = exp(1i*k_0(k)*x(i)) + Phis(i);
            endfor
            Ref(k) = abs(Phis(1)/exp(1i*k_0(k)*x(1)))^2;
            Tra(k) = abs(Phi(2))^2;
            Ab(k) = 1 - (Ref(k) + Tra(k));
        endif
        fprintf('%4d %15.3E %15.3E %15.3E %15.3E %15.3E\n', k, k_0(k), E0(k), Ref(k), Tra(k), Ab(k));
    endfor

    % Plotting
    spec_graph = plot(E0, Ref, 'b', E0, Tra, 'r', E0, Ab, 'g');
    
    % Customizing the plot
    xlabel('Energy (eV)');
    ylabel('Coefficient');
    title('Reflection, Transmission, and Absorption vs. Energy');
    legend('Reflection', 'Transmission', 'Absorption', 'Location', 'west');
    grid on;
    
    xlim([0, 0.3]);
    ylim([0,1]);
  toc
  return
end

#-----------------------------------------------------------------
# Question 3 (Resonance Theory Included in pdf)
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# Question 4
#-----------------------------------------------------------------

if strcmp(Part1_Question4_plot, 'yes') == 1
 tic
  mystartdefaults
  step_size = 0.5;
  x_min = -20;
  x_max = 100;
  x1_min = 0;
  x1_max = 80;
  n = (x_max - x_min)/step_size;
  n1 = (x1_max - x1_min)/step_size;
  tau = 1e-9;
  gam = (hbar * 2 * pi/tau)/qel;
  recipunit=1.0E+10;
  ekinscale=(hbar * recipunit)^2/(2 * elm)/qel;
  
  function [Gf]= Green(step,x,xp,E,gam,ekinscale) 
    ko=sqrt((E+1i*gam)/ekinscale);
    k1=sqrt((E+1i*gam-step)/ekinscale);
    if(x>=0 && xp>=0)
        Gf=exp(1i*k1*abs(x-xp))/(2i*k1) +  exp(1i*k1*(x+xp))*((k1-ko)/(k1+ko))/(2i*k1); % Green's function [A.76] for both positions greater than 0
    end

    if (x<0 && xp>=0)
        Gf=exp(-1i*ko*x+1i*k1*xp)/(1i*(ko+k1)); % Green's function [A.77] for x and xp  on different sides
    end
    if (xp<0 && x>=0)
        Gf=exp(-1i*ko*xp+1i*k1*x)/(1i*(ko+k1)); % Green's function [A.77] for xp and x  on different sides
    end
  end
  %% The Phi_K part added to the solution due to the potential 
  function[Phi_K]= sol_out(x,pos,V,dx,phisol,step,E,gam,ekinscale)
    Phi_K=0;
    [~,nco]=size(pos);
    for ii=1:nco
        Phi_K=Phi_K+Green(step,x,pos(ii),E,gam,ekinscale)*dx*(V(ii)/ekinscale)*phisol(ii); %equation [4.36]
    end
  end
  
  U_total = zeros(1,n); % Potential Barriers
  x = zeros(1,n);
  for i = 1:n
    x(i) = x_min + step_size/2 + (i - 1) * step_size;
    if x(i) >= 0 && x(i) <= 15
        U_total(i) = 0.2;
    elseif x(i) >= 65 && x(i) <= 80
        U_total(i) = 0.2;
    else
        U_total(i) = 0;
    end
  end
  U = zeros(1,n1); % Potential barrier
  x1 = zeros(1,n1);
  for i = 1:n1
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
    if x1(i) <= 15
        U(i) = 0.2;
    elseif x1(i) >= 65 && x1(i) <= 80
        U(i) = 0.2;
    else
        U(i) = 0;
    end
  end
  resoener=[0.01075,0.04325,0.09525,0.16325,0.07]; %Resonance energies found by looking at the R-T-A graph AND a non resonat enegry(the last one)
  [~,n]=size(resoener);
  wavefunctions=zeros(n,length(U_total));
  W = zeros(length(U),length(U));
  
  wavefunctions = zeros(n-1, length(U_total));
  for i = 1:n
    k1 = sqrt((resoener(i) + 1i * gam)/ekinscale);
    for j = 1:length(U)
        phi0k(j) = exp(1i * k1 * x1(j));
        W(j,j) = step_size * U(j)/ekinscale;
        for k = 1:length(U)
            G0(j,k) = Green(0,x1(j),x1(k),resoener(i),gam,ekinscale);
        end
    end
    T = eye(length(U)) - G0 * W;
    phisol = T\(phi0k.');
    for j = 1:(length(U_total))
        if i < n
            wavefunctions(i,j) = exp(x(j) * 1i * k1) + sol_out(x(j), x1, U, step_size, phisol, 0, resoener(i), gam, ekinscale);
        else
            nonresophi(j) = exp(x(j) * 1i * k1) + sol_out(x(j),x1,U,step_size,phisol,0,resoener(i),gam,ekinscale);
        end
        end
    end
    resoprobas = abs(wavefunctions).^2;
    nonresoproba = abs(nonresophi).^2;
    for jj=1:length(U_total)    
      resoprobas(1,jj)= 0.0065*resoprobas(1,jj)+resoener(1);
      resoprobas(2,jj)= 0.002*resoprobas(2,jj)+resoener(2);
      resoprobas(3,jj)= 0.0005*resoprobas(3,jj)+resoener(3);
      resoprobas(4,jj)= 0.0027*resoprobas(4,jj)+resoener(4); 
      nonresoproba(jj)= 0.01*nonresoproba(jj)+resoener(n);
    endfor
    
    plot(x, nonresoproba, 'r', x, resoprobas, 'g', x, U_total, 'k', 'LineWidth', 2);
    title('Probability densities against Barriers');
    ylim([0, 0.25]);
    xlabel('x (Angstrom)');
    ylabel('E (eV)'); 
    grid on;
 toc
 return;
end


#-----------------------------------------------------------------
# Question 5
#-----------------------------------------------------------------
#Part1_Question5_plot ='yes';

if strcmp(Part1_Question5_plot, 'yes') == 1
  tic
    for i = 1:n
        x1(i) = x1_min + step/2 + (i-1)*step;
        if x1(i) <= 15
            U(i) = 0.2;
        elseif x1(i) >= 65 && x1(i) <= 80
            U(i) = 0.1;
        else
            U(i) = 0;
        end
        V(i) = U(i)/ekinscale;
    end

    k_0 = zeros(length(E0), 1);
    lmbda = zeros(length(E0), 1);
    for i = 1:length(E0)
        k_0(i) = sqrt(E0(i)/ekinscale);
        lmbda(i) = 2*pi/k_0(i);
    end

    m = 2;
    x = zeros(m, 1);
    Phis = zeros(m, 1);
    Phi = zeros(m, 1);
    x(1) = x1_min - 1;
    x(2) = x1_max + 1;

    fprintf('    n           k             E0                R             T              R+T\n');
    for k = 1:length(E0)
        if E0(k) == 0
            Ref(k) = 1;
            Tra(k) = 0;
            RpT(k) = 1;
        else
            Phi0p = exp(1i*k_0(k)*x1);
            G0 = zeros(n, n);
            for j = 1:n
                for i = 1:n
                    G0(i, j) = step*exp(1i*k_0(k)*abs(x1(i)-x1(j)))/(2i*k_0(k));
                end
            end

            T = eye(n, n) - G0*diag(V);
            Phip = T\Phi0p;

            for i = 1:m
                Phis(i) = 0;
                for j = 1:n
                    Phis(i) = Phis(i) + step*(exp(1i*k_0(k)*abs(x(i)-x1(j))))/(2*1i*k_0(k)) * V(j) * Phip(j);
                end
                Phi(i) = exp(1i*k_0(k)*x(i)) + Phis(i);
            end
            Ref(k) = abs(Phis(1)/exp(1i*k_0(k)*x(1)))^2;
            Tra(k) = abs(Phi(2))^2;
            Ab(k) = 1 - (Ref(k) + Tra(k));
        end
        fprintf('%4d %15.3E %15.3E %15.3E %15.3E %15.3E\n', k, k_0(k), E0(k), Ref(k), Tra(k), Ab(k));
    end

    % Plotting
    spec_graph = plot(E0, Ref, 'b', E0, Tra, 'r', E0, Ab, 'g');
    
    % Customizing the plot
    xlabel('Energy (eV)');
    ylabel('Coefficient');
    title('Reflection, Transmission, and Absorption vs. Energy');
    legend('Reflection', 'Transmission', 'Absorption', 'Location', 'west');
    grid on;
    
    xlim([0, 0.3]);
    ylim([0,1]);
    

    waitfor(spec_graph);
  toc  
  return;
end
#----------------------------------------------------------------------------
# PART 2
#----------------------------------------------------------------------------
#--------------------------------
# Question 1 and 2
#--------------------------------
# Declaring units
recipunit = 1.0E+10;
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;


Ea = 0;
lambda_min = 1;
k = 2 * pi / lambda_min;
Eb = ceil((hbar * recipunit)^2 * k^2 / (2 * elm) / qel);

x1_min = 0;
x1_max = 80;
step = 0.5;
n = (x1_max - x1_min) / step;
x1=zeros(n,1);
for i = 1:n
  x1(i) = x1_min + step/2 + (i-1)*step;
end
x_min = -20;
x_max = 100;
m = (x_max - x_min) / step;
x = zeros(m, 1);
U = zeros(m, 1);
W = zeros(m, 1);
U_0 = 0;
U_1 = -0.10;


E = -U_1 / (x1_max - x1_min);

for i = 1:m
    x(i) = x_min + step/2 + (i-1)*step;
    if x(i) >= 0 && x(i) <= 15
        U(i) = 0.2;
    elseif x(i) >= 65 && x(i) <= 80
        U(i) = 0.2;
    else
        U(i) = 0;
    end
    V(i) = U(i) / ekinscale;
end

% Modified potential W

for i = 1:m
    if x(i) >= 0 && x(i) <= 80
        W(i) = U(i) - E * x(i);
    else
        W(i) = 0;
    end
    V(i) = W(i)/ekinscale;
end

% Plot the potential profile

if strcmp(Part2_Question1_plot, 'yes') == 1
 tic
  figure;
  plot(x, V, 'b','LineWidth', 2);
  title('Potential Profile');
  xlabel('x');
  ylabel('Potential (U)');
  grid on;
 toc
 return;
end

#--------------------------------
# Question 3
#--------------------------------

function [Gf]= Green(step,x,xp,E,gam,ekinscale) 
    ko=sqrt((E+1i*gam)/ekinscale);
    k1=sqrt((E+1i*gam-step)/ekinscale);
    if(x>=0 && xp>=0)
        Gf=exp(1i*k1*abs(x-xp))/(2i*k1) +  exp(1i*k1*(x+xp))*((k1-ko)/(k1+ko))/(2i*k1); % Green's function [A.76] for both positions greater than 0
    end

    if (x<0 && xp>=0)
        Gf=exp(-1i*ko*x+1i*k1*xp)/(1i*(ko+k1)); % Green's function [A.77] for x and xp  on different sides
    end
    if (xp<0 && x>=0)
        Gf=exp(-1i*ko*xp+1i*k1*x)/(1i*(ko+k1)); % Green's function [A.77] for xp and x  on different sides
    end
end


   %% Barrier potential
   function Bar_V =Barriers(pos,bar1,bar2,u1,u2)
     [~,nco]=size(pos);
     Vbar=zeros(1,nco);
     for ii=1:nco
        if(pos(ii)<=bar1)
            Vbar(ii)=u1;  
        elseif(bar2<=pos(ii))
            Vbar(ii)=u2;   
     end
    end    
    Bar_V=Vbar;
   end
   
   function[Phi_K]= sol_out(x,pos,V,dx,phisol,step,E,gam,ekinscale)
    Phi_K=0;
    [~,nco]=size(pos);
    for ii=1:nco
        Phi_K=Phi_K+Green(step,x,pos(ii),E,gam,ekinscale)*dx*(V(ii)/ekinscale)*phisol(ii); %equation [4.36]
    end
   end
   %% Reflection, transmission and absorption coefficients
   function[refl,tra,absor] = spectra(step,E,gam,pos,V,dx,ekinscale)
     [~,nco]=size(pos);
     ko=sqrt((E+1i*gam)/ekinscale);
     k1=sqrt((E+1i*gam-step)/ekinscale); 
     tb=2*ko/(ko+k1);
     rb=(ko-k1)/(k1+ko);
     Wmat=zeros(nco,nco);
     for ii=1:nco
       phiok(ii)=tb*exp(1i*k1*pos(ii));% Function to the right of the step [A.65]
       Wmat(ii,ii)=dx*V(ii)/ekinscale; % Equation[4.39]
       for jj=1:nco
        Go(ii,jj)= Green(step,pos(ii),pos(jj),E,gam,ekinscale); %Equation [4.39]
       end
     end
     scatt=eye(nco)-Go*Wmat;
     phisol=scatt\(phiok.'); % Solution inside the perturbated region [4.51]
   
     phirefl=rb*exp(-1i*ko*(pos(1)-1))+sol_out(pos(1)-1,pos,V,dx,phisol,step,E,gam,ekinscale); %The first part comes from a pre-existing reflection in equation [A.65]
     phitrans=tb*exp(1i*k1*(pos(nco)+1))+sol_out(pos(nco)+1,pos,V,dx,phisol,step,E,gam,ekinscale); %The first part comes from an already exinsting transmission in equation [A.65]
     
     tra=(real(k1)/real(ko))*abs(phitrans)^2;
     refl=(abs(phirefl/(exp(1i*ko*(pos(1)-1)))))^2; %comparing the reflected part with just an incoming wave
     absor=1-(refl+tra);
   end
   
   
   dx= 0.5; % discretization of horizontal position
   xmin=-20; %lower bound
   xmax=100; %upper bound
   xpmin=0; %lower bound  of the perturbated region
   xpmax=80; %upper bound of the perturbated region
   bar1=15; %position of the end of the first barrier
   bar2=65; %position of the beggining of the second barrier
   u1=0.2; %heights of Barriers
   u2=0.1;
   step=-0.1; % Potential bias (for the second part)
   dE=0.0005; %discretization of energy steps
   Emax=0.3; %Upper bound for incoming energy
   Emin=0.00; %lower bound for incoming energy
   liftim=1.0E-9; % Lifetime
   gam=(hbar*2*pi/liftim)/qel; % damping factor
   stepmin=-0.2; %Lower bound for energy bias
   stepmax=0.2; %Upper bound for energy bias
   %% Setting up the x grid and the two equal height Barriers
   piecesin=round((xpmax-xpmin)/dx);
   for ii=1:piecesin
    xin(ii) = dx/2+xpmin+(ii-1)*dx; %positions inside the perturbated zone
   end
   piecesL=round((xpmin-xmin)/dx);
   for ii=1:piecesL
    xL(ii) = dx/2+xmin+(ii-1)*dx; %positions to the left
   end
   piecesR=round((xmax-xpmax)/dx);
   for ii=1:piecesR
    xR(ii) = dx/2+xpmax+(ii-1)*dx; %positions to the right
   end
   x=[xL xin xR]; %all positions
   Vbar=Barriers(xin,bar1,bar2,u1,u1); %potential with two equal height Barriers inside
   Vbartot=[zeros(1,piecesL) Vbar zeros(1,piecesR)];

   
   %% Creating the potential for the biased junction
   for ii=1:piecesin
    efield(ii)=step*xin(ii)/(xpmax-xpmin);%-step; %Taking the reference
   end
   steppot=Vbar+efield;
   for ii=1:piecesR
    Vright(ii)=step;
   end
   Vbartotstep=[zeros(1,piecesL) steppot Vright];

   %% Setting up the x grid and the two equal height Barriers
   piecesin=round((xpmax-xpmin)/dx);
   for ii=1:piecesin
    xin(ii) = dx/2+xpmin+(ii-1)*dx; %positions inside the perturbated zone
   end
   piecesL=round((xpmin-xmin)/dx);
   for ii=1:piecesL
    xL(ii) = dx/2+xmin+(ii-1)*dx; %positions to the left
   end
   piecesR=round((xmax-xpmax)/dx);
   for ii=1:piecesR
    xR(ii) = dx/2+xpmax+(ii-1)*dx; %positions to the right
   end
   x=[xL xin xR]; %all positions
   Vbar=Barriers(xin,bar1,bar2,u1,u1); %potential with two equal height Barriers inside
   Vbartot=[zeros(1,piecesL) Vbar zeros(1,piecesR)];
   
 
if strcmp(Part2_Question3_plot, 'yes') == 1
  tic
   %% Taking the step as a reference(IMPORTANT!)
   for ii=1:piecesin
    steppot(ii)=steppot(ii)-step; 
   end
   %% Calculate minimum and maximum wavelengths on both sides
   wavelengths = true;
   if(wavelengths)
    k0m=real(sqrt((Emin+1i*gam)*2*elm*qel/hbar^2));       % Dispersion relation for the left media [A.58]
    k0M=real(sqrt((Emax+1i*gam)*2*elm*qel/hbar^2));
    k1m=real(sqrt((Emin-step+1i*gam)*2*elm*qel/hbar^2));  % Dispersion relation for the right media [A.58]
    k1M=real(sqrt((Emax-step+1i*gam)*2*elm*qel/hbar^2));

    L0m= 2*pi*(1E+10)/k0M;
    L0M= 2*pi*(1E+10)/k0m;
    L1m= 2*pi*(1E+10)/k1M;
    L1M= 2*pi*(1E+10)/k1m;

    fprintf('Longest left wavelenght= %#5.5G [angstrom]\n', L0M);
    fprintf('Shortest left wavelenght= %#5.5G [angstrom]\n', L0m);
    fprintf('Longest right wavelenght= %#5.5G [angstrom]\n', L1M);
    fprintf('Shortest right wavelenght= %#5.5G [angstrom]\n', L1m);
   end
   %% R-T-A plot for step of -0.1eV
   reftrastep = true;
   if(reftrastep) 
   edivs=round((Emax-Emin)/dE);
   for ii=1:edivs
     Eners(ii) =dE/2 +Emin+(ii-1) * dE; %Energy range we´ll look at
   end
   for ii=1:edivs
    [Reff(ii),Traa(ii),absor(ii)]=spectra(step,Eners(ii),gam,xin,steppot,dx,ekinscale) ; 
   end
   compalabel{1}=sprintf('Transmission');
   compalabel{2}=sprintf('Reflection');
   compalabel{3}=sprintf('Absorption');

   figure
   plot(Eners,Traa,'r', Eners, Reff,'b', Eners,absor,'g');

   xlabel('Energy(eV)',"fontsize", 15);
   ylabel('Coefficient',"fontsize", 15);
   title('Biased junction(-0.1eV)');
   legend(compalabel,"fontsize",8,"location",'west');
   xlim([0, 0.3]);
   ylim([0,1]);
   grid on;
   end
  toc 
  return
end
#----------------------------------------------------------------------------
# Question 4
#----------------------------------------------------------------------------
%% Profiles  for the limits of the varying biased junction

if strcmp(Part2_Question4_plot, 'yes') == 1
 tic
  for ii=1:piecesin
    efield1(ii)=stepmin*xin(ii)/(xpmax-xpmin);
    efield2(ii)=stepmax*xin(ii)/(xpmax-xpmin);
  end
  steppot1=Vbar+efield1;
  steppot2=Vbar+efield2;
  for ii=1:piecesR
    Vright1(ii)=stepmin;
    Vright2(ii)=stepmax;
  end
  Vbartotstep1=[zeros(1,piecesL) steppot1 Vright1];
  Vbartotstep2=[zeros(1,piecesL) steppot2 Vright2];
  figure
  subplot(1,2,1)
  plot(x,Vbartotstep1,'b');
  title('U1=-0.2eV');
  ylabel('E(eV)');
  xlabel('x(Angstrom)');
  grid on;

  subplot(1,2,2)
  plot(x,Vbartotstep2,'b');
  title('U1=0.2eV');
  ylabel('E(eV)');
  xlabel('x(Angstrom)');
  grid on;
 toc
 return;
end

#----------------------------------------------------------------------------
# Question 4
#----------------------------------------------------------------------------
%% Profiles  for the limits of the varying biased junction

if strcmp(Part2_Question5_plot, 'yes') == 1
  tic
   %% R-T-A plot for varying step
   fixedE=0.01; % We fix this energy for the incoming  electron

   edivs=round((stepmax-stepmin)/dE);
   for ii=1:edivs
     Enerstep(ii) = dE/2 + stepmin+(ii-1)*dE; 
   end

   for ii=1:edivs
     for jj=1:piecesin
        efield(jj)=Enerstep(ii)*xin(jj)/(xpmax-xpmin)-Enerstep(ii); % The previous is to take the step as a reference
     end
     steppot=Vbar+efield;
     [Reff(ii),Traa(ii),absor(ii)]=spectra(Enerstep(ii),fixedE,gam,xin,steppot,dx,ekinscale) ; 
   end

   compalabel{1}=sprintf('Transmission');
   compalabel{2}=sprintf('Reflection');
   compalabel{3}=sprintf('Absorption');
   figure
   plot(-Enerstep,Traa,'r', -Enerstep, Reff,'b',-Enerstep,absor,'g'); % The minus sign indicates the fact that electron direction is opposite to conventional direction of current

   xlabel('Bias(eV)');
   ylabel('Coefficient');
   title('Varying bias');
   legend(compalabel,"location",'west');
   grid on;
   hold on;
  toc
  return;
end

  
  
  
  
  
  