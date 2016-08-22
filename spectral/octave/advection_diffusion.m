clear all; close all;
graphics_toolkit("gnuplot");

nu=0.001;
Lx=20; Ly=20; nx=64; ny=64; N=nx*ny;
x2=linspace(-Lx/2,Lx/2,nx+1); x=x2(1:nx);
y2=linspace(-Ly/2,Ly/2,ny+1); y=y2(1:ny);


% INITIAL CONDITIONS
[X,Y]=meshgrid(x,y);
w=1*exp(-0.25*X.^2-Y.^2);
figure(1), pcolor(abs(w)); 


% SPECTRAL K VALUES
kx=(2*pi/Lx)*[0:(nx/2-1) (-nx/2):-1]'; kx(1)=10^(-6);
ky=(2*pi/Ly)*[0:(ny/2-1) (-ny/2):-1]'; ky(1)=10^(-6);

% RHS of the system of differential equations which results from
% Fourier transforming is contained within the function spc_rhsrhs
function rhs=spc_rhs(tspan,wt2,dummy,psix,psiy,kx,ky,nu,nx,ny,N);
    wt=reshape(wt2,nx,ny);
    % w_x
    for j=1:nx
        wtx(j,:)=i*(wt(j,:).*kx');
    end
    wx=ifft2(wtx*N);
    
    % w_y
   for j=1:ny
       wty(:,j)=i*(wt(:,j).*ky);
    end
    wy=ifft2(wty*N);
    
    % transform w_x*psi_y and w_y*psi_x and reshape
    wxpyt=fft2(wx.*psiy)/N;
    wypxt=fft2(wy.*psix)/N;
    wxpyt2=reshape(wxpyt,N,1);
    wypxt2=reshape(wypxt,N,1);
    
    % Laplacian Terms
    for j=1:nx
        wtxx(j,:)=-wt(j,:).*(kx.^2)';
    end
    for j=1:ny
        wtyy(:,j)=-wt(:,j).*(ky.^2);
    end
    wtxx2=reshape(wtxx,N,1);
    wtyy2=reshape(wtyy,N,1);
    rhs=(nu*(wtxx2+wtyy2)-wypxt2+wxpyt2);

endfunction    

% Streamfunction
for stream_loop=1:30
    wt=fft2(w)/N;
    for j=1:ny
        psit(:,j)=-wt(:,j)./(kx.^2+ky(j)^2);
    end
    for j=1:nx
        psitx(j,:)=i*(psit(j,:).*kx');
    end
    psix=real(ifft2(psitx*N));
    for j=1:ny
        psity(:,j)=i*(psit(:,j).*ky);
    end
    psiy=real(ifft2(psity*N));

    wt2=reshape(wt,N,1);
    [t,wsol]=ode23(@spc_rhs,[0 0.5],wt2,[],psix,psiy,kx,ky,nu,nx,ny,N);

    wt2=wsol(end,:)';
    wt=reshape(wt2,nx,ny);
    w=ifft2(wt*N);

    figure(2)
    pcolor(abs(w)); 
    pause(1)
end



