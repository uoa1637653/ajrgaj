function plot_eigs(vec,lam)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Plots the given eigenvectors against the discretised, periodic 
% domain pre-initialised by init_domain().
% Groups the eigenvalues according to wavenumber, and plots the
% corresponding eigenvectors for each wavenumber group separately.
% The eigenvalues are shown in the legend.
%----------------------------------------------------------------
global W L x;
vec=[vec(L,:);vec];
wavenums=round(real(sqrt(-real(lam)))/(2*pi/W));
for k=0:floor(L/2)
  idx=find(wavenums==k);
  if ~isempty(idx)
     figure;
     plot(x,vec(:,idx));
     legend(num2str(lam(idx),3));
     xlabel('x'); ylabel('v_k(x)');
	 title(strcat('Eigenvectors/values for wavenumber k=',num2str(k)));
	 vmin=min(min(vec(:,idx)));
	 vmax=max(max(vec(:,idx)));
	 axis([0 W vmin-0.1 vmax+0.1]);
  end
end
