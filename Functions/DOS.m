function pred = DOS(x,p)
    pred = ([1./(exp((-x+p(1)).*p(2))+1)]-0.5)+p(3);
    pred = p(4).*bsxfun(@times,pred,pred')+p(5);
end