%%% Dichiarazione Variabili %%%
% Variabili
a = 2.8;
b = 10 ;
Len = 3;
espr = 2.3;
%Costanti
epsylon0 = 8.85e-12;
mu0 = 4 * pi * 1e-7;

C = (epsylon0 * espr * 2 * pi)/log(b/a); % [Fm]
L = mu0 / (2 * pi) * log(b/a); % [H/m] 
v = 1/sqrt(L*C);
fmax =1.3e9; 
dz = v/ (fmax*50);
N = ceil(Len/dz);   %% Arrotonda all'intero piu' prossimoo il rapporto Len dz

dz = Len /(N-1); %% ???

% Time
 T = 50e-9;
 dt = (0.5 *dz/v)/500;
 t0 = T/8; 
 deltaT = 300e-12; 
 
 %% Funzioni di Ingresso
 Vin = @(t) normpdf(t,t0,deltaT); % funzione gaussiana ( Normal probability density function )
 %Vin = @(t) sin(t*2*pi*1.4e9); % sinusoide
 
 %% Carico
 %R_load = 5e1;
 R_load =sqrt(L/C);
 C_load = 1;
 I_load = 1;
 %% Grafici
 figure; %% Crea la finestra
 
 Z = linspace(0,Len,N);  % Segmenti di percorso
 V = zeros(1,N); %% Inizializzo la tensione
 I = zeros(1,N); %% Inizializzo la corrente
 
%subplot(3,1,2);
fplot(@(x) Vin(x),[0,T]); %% Crea il grafico
title('$V_{in}\;\;in\;\;T$','Interpreter','latex'); %% Titolo del grafico
xlabel('time[s]'); %% Titolo Asse X
ylabel('tensione[V]'); %% Titolo Asse Y
%grid(); %% Crea una griglia nel grafico
 
%subplot(3,1,2);
render(Z,V,I);
 
 pause();
 
 tensioni = [];
 correnti = [];
 

%% Loop
tmp = 0;
i = 0;
for t = 0:dt:T
   i=i+1;
   
   V(1) = Vin(t);
   %%calcolo condizioni di contorno in base al carico
   [I,V,tmp] = contorni(I,V,tmp);   
   [V,I] = steps(V,I,L,C,dt,dz);
   
   
   if mod(i,10000) == 0
       disp(t/T);
       %subplot(3,1,1);
       render(Z,V,I);
       
       drawnow();
   end
   stride = 100;
    if (mod(i, stride) == 0)
        tensioni(1, i / stride) = V(1);
        correnti(1, i / stride) = I(1);
        tensioni(2, i / stride) = V(floor(end/2));
        correnti(2, i / stride) = I(floor(end/2));
        tensioni(3, i / stride) = V(end);
        correnti(3, i / stride) = I(end);
    end
   
end
 
 %% Funzioni
 function [V, I] = steps(V, I, L, C, dt, dz)
    N = numel(V);
    
    dv_dz = diff(V) / dz;
    di_dz = diff(I) / dz;
    
    I(1:end-1) = I(1:end-1) - (dv_dz / L) * dt;
    V(2:end) = V(2:end) - (di_dz / C) * dt;
end
 function [] = render(xx, V, I)
    subplot(2,1,1);
    xlabel("Lunghezza [m]")
    %yyaxis left
    plot(xx, V,'b');
    ylabel('Tensione [ V ]')
    legend('V (t)');
    subplot(2,1,2);
    %ylabel("tensione [V]")
    %yyaxis right
    plot(xx(1:end-1) + (xx(2) - xx(1))*0.5, I(1:end-1),'r'); % derivata
    legend('I (T)');
    ylabel('Corrente [ A ]');
    %ylabel("corrente [A]")
    %yyaxis left
 end
 function [I,V,tmp] = contorni(I,V,tmp,scelta)
    switch scelta       
        case 0 %% Resistenza a carico
             I(end) = V(end)/R_load; %% Condizione di contorno della resistenza
        case 1 %% Condensatore a carico
            I(end) =I(end-1);
            contorno = tmp + I(end) /C_load * delta_t;   % Condizione di contorno del condensatore
            V(end) = contorno;
            tmp = contorno;
        case 2 % Induttore a carico
            I(end) = I(end)+ V(end)/L_load*delta_t;
        otherwise
            error('Scelta non valida!! Deve essere un numero compreso tra 0 e 2');
    end

end
 