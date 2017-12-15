function [uEF, duEF, xEF] = postProcFibre(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,u, phir, dphir, r, p, fig)
	%plot the solution u and its derivative
	if fig
		img = figure;
	end
	for ne=1:nEls
	    for i=1:elDof(ne)
	        c=dFreedom(ne,i);
	        uL(ne,i)=u(c);
	    end
	end
	duu=zeros(nEls-1,1);
	for ne=1:nEls
	    nx=10;
	    u=zeros(nx,1);
	    du=zeros(nx,1);
	    x1=nodes(connect(ne,1));
	    x2=nodes(connect(ne,2));
	    h=x2-x1;
	    jac=2/h;
	    x=(x1:(x2-x1)/(nx-1):x2);
	    xi=(-1:2/(nx-1):1);
	    for j=1:nx
	        [N,dN]=shape(xi(j),ne,pDeg,pType);
	        dN=dN*jac;
	        for i=1:elDof(ne)
	            u(j)=u(j)+N(i)*uL(ne,i);
	            du(j)=du(j)+dN(i)*uL(ne,i);
	            %             if i==1 && j==1 && ne==1
	            %                 duu(1)=du(1);
	            %                 test=du
	            %             end
	            %             if ne>1
	            %                 duu(ne)=du(end);
	            %             end
	        end
	    end

		% if fig
		%     subplot(2,2,1)
		% 	title('EF U')
		%     hold on
		%     plot(x,u)
		%     subplot(2,2,2)
		% 	title('EF dU')
		%     hold on
		%     plot(x,du)
		% end

		uEF(ne) = u(1);
	    duEF(ne) = du(1);
		xEF(ne) = x(1);
	end

	if fig
		[PHI, u] = normalise(phir, uEF, false);
		[dPHI, du] = normalise(dphir, duEF, false);
		for i = 2:length(uEF)
		    subplot(2,2,1)
		    hold on
		    plot(xEF(i-1:i),u(i-1:i))
		    subplot(2,2,2)
		    hold on
		    plot(xEF(i-1:i),du(i-1:i))
		end

		subplot(2,2,1)
		title('EF U')
		subplot(2,2,2)
		title('EF dU')
		subplot(2,2,3)
		title('Exacte U')
		plot(r,PHI)
		subplot(2,2,4)
		title('Exacte dU')
		plot(r, dPHI)

		saveImg(img, strcat('ef_vs_exact_p', num2str(p)));
	end
end
