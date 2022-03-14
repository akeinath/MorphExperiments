function [params error] = fitContextPattern(v,nsims)
    if nargin < 2 || isempty(nsims)
        nsims = 1;
    end
    
    params = nan(length(v(1,1,:)),5);
    error = nan(length(v(1,1,:)),1);
    parfor i = 1:length(v(1,1,:))
        
        x = [1:length(v(:,1,1))];
        
        pparams = nan(nsims,5);
        perror = nan(nsims,1);
        for si = 1:nsims
            seed = [(0.5+(length(v(:,:,1)))./4)+rand.*(0.5+(length(v(:,:,1)))./2) ...
                rand.*10 1.*(rand-0.5) 2.*(rand) 2.*(rand-0.5)];
            [pparams(si,:) perror(si)] = fmincon(@(p)sig_error(v(:,:,i),p),seed ,[],[],[],[],...
                    [0 0 -0.5 0 -5],[x(end)+1 10 0.5 10 5],[],optimoptions('fmincon','Display','none'));
        end
        [a b] = nanmin(perror);
        params(i,:) = pparams(b(1),:);
        error(i) = perror(b(1));
    end

end

function error = sig_error(data,p)
    pred = DOS(1:length(data(1,:)),p);
%     error = nanmean(sqrt((data(:)-pred(:)).^2));
    
%     error = nanmean(((atanh(data(:))-pred(:)).^2));
    error = nanmean((((data(:))-pred(:)).^2));
end

function error = sig_error2(data,p)
    pred = DO2S(1:length(data(1,:)),p);
    error = nanmean(sqrt((data(:)-pred(:)).^2));
end