n = 100;
h = 0.5;
z = (0:(n-1))'*h;
[XX,YY,ZZ] = ndgrid(z,z,z);
e = ones(n,1);
Df = spdiags([-e e], 0:1, n, n)/h;
Df(n,1) = Df(1,2);
Db = spdiags([-e e], -1:0, n, n)/h;
Db(1,n) = Db(2,1);
I= speye(n);
Z = 0*I;

Dfx = kron(kron(Df,I),I);
Dfy = kron(kron(I,Df),I);
Dfz = kron(kron(I,I),Df);

Dbx = kron(kron(Db,I),I);
Dby = kron(kron(I,Db),I);
Dbz = kron(kron(I,I),Db);

Z3 = kron(kron(Z,Z),Z);

Ex = zeros(n,n,n);
Ex(:) =exp(-(XX-5).^2-(YY-5).^2-(ZZ-5).^2);
%Ex(1) =1;
Ey = zeros(n,n,n);
Ez = zeros(n,n,n);
E = [Ex(:);Ey(:);Ez(:)];
size(E);

Bx = zeros(n,n,n);
By = zeros(n,n,n);
Bz = zeros(n,n,n);
B = [Bx(:);By(:);Bz(:)];

Yx = zeros(n,n,n);
Yy = zeros(n,n,n);
Yz = zeros(n,n,n);
Y = [Yx(:);Yy(:);Yz(:)];

kappa = zeros(n,n,n);
kappa(n/2:n,:,:) = 2;
kappa = [kappa(:);kappa(:);kappa(:)];

Px = zeros(n,n,n);
Py = zeros(n,n,n);
Pz = zeros(n,n,n);
P = [Px(:);Py(:);Pz(:)];

nabla_cross_f= [Z3 -Dfz Dfy  ; ...%[ Ex ]
                Dfz Z3 -Dfx  ; ...% [Ey]
                -Dfy Dfx Z3 ]; ...%[Ez]
                
nabla_cross_b= [Z3 -Dbz Dby  ; ...%[ Ex ]
                Dbz Z3 -Dbx  ; ...% [Ey]
                -Dby Dbx Z3 ]; ...%[Ez]
                            
epsilon0 = 1;
mu0 = 1;
c=1/sqrt(epsilon0*mu0);
%Ex(1) = 1;

Nsteps = 1000;
dt = 0.03;
data = zeros(n,Nsteps);
omega = 0;
rho = 1;
%kappa = 1;

for i=1:Nsteps
    Er = reshape(E,n,n,n,3);
    temp = fftn(Er(:,:,:,1));
    data(:,i) = squeeze(temp(1,:,1));
    %data(:,i)=Ex;
    
    dBdt = -nabla_cross_f * E;
    B = B + dBdt*dt;
    Y = Y - omega^2 * P*dt - 1/rho *(kappa.* E)*dt;
        
    dEdt = c^2*nabla_cross_b*B;
    E = E + dEdt*dt + 1/epsilon0*Y*dt;
    P = P + Y * dt;
    
   if mod(i,50)==0
        Pr = reshape(P,n,n,n,3);
        subplot(2,2,1);
        imagesc(squeeze(Pr(:,1,:,1)));
        colorbar
        title('P');
        Br = reshape(B,n,n,n,3);
        subplot(2,2,2);
        imagesc(Br(:,:,5,2));
        colorbar
        pause(0.01);
        hold off
        
        subplot(2,2,3);
        dataG = fft(data,[],2);
        imagesc(flipud(abs(dataG(:,1:20))'));
        xlabel('k');
        ylabel('\omega');
        colorbar
   end
end


