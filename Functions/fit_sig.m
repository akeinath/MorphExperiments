function fit_sig(vals)

    x = 1:length(vals(1,:));
    for k = 1:length(vals(:,1))

        [params fval] = fmincon(@(p)sig_error(vals(k,:),p),[(x(end)./2)+0.5 1 0 2] ,[],[],[],[],...
                [-inf 1e-10 -1 1e-10],[inf 10 1 inf],[],optimoptions('fmincon','Display','none'));
    end

end

function error = sig_error(data,p)
    pred = doSig(1:length(data(1,:)),p);
    error = nanmean((data-pred).^2);
end

function pred = doSig(x,p)
    pred = ([1./(exp((-x+p(1)).*p(2))+1)].*p(4))+p(3);
end