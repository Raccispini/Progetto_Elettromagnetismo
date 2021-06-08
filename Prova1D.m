n = 100;
h = 0.1;
z = (0:(n-1))'*h;
Ex = zeros(n,1);
Px = zeros(n,1);
Yx = zeros(n,1);
By = zeros(n,1);
epsilon0 = 1;
mu0 = 1;

c=1/sqrt(epsilon0*mu0);

e = ones(n,1);
Df = spdiags([-e e], 0:1, n, n)/h;
Df(n,1) = Df(1,2);
Db = spdiags([-e e], -1:0, n, n)/h;
Db(1,n) = Db(2,1);
%Ex(1) = 1;
Ex(:) = exp(-(z-5).^2);
Nsteps = 1000;
dt = 0.03;
data = zeros(n,Nsteps);
omega = 1;
rho = 1;

for i=1:Nsteps
    
    
    
    data(:,i)=Ex;
    dBydt = -Df*Ex;
    By = By + dBydt*dt;
    Yx = Yx - omega^2 * Px*dt + 1/rho * Ex*dt;
    
    
    dExdt = -c^2*Db*By;
    Ex = Ex + dExdt*dt - 1/epsilon0*Yx*dt;
    Px = Px + Yx * dt;
    
    if mod(i,1)==0
        hold off
        plot(z-h/2,Ex,'r');
        
        
        hold on
        plot(z,By,'b--');   
        plot(z-h/2,Px,'g');
        plot(z-h/2,Yx,'k');
        ylim([-1 1]);
        legend({'E_x','B_y','P_x','Y_x'});
        pause(0.001);
    end
end

dataG = fftn(data);
clf
imagesc(flipud(abs(dataG(:,1:100))'))
xlabel('k');
ylabel('\omega');
