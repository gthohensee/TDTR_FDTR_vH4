function Z = knife_edge(X,xpos,knife_data)

A = max(knife_data);
mid = X(1);
w0 = X(2);

knife_model = (A/2) * (1-erf( (sqrt(2)*(xpos-mid))/w0));

figure(10)
plot(xpos,knife_data,'ob',xpos,knife_model,'r-'); 
xlabel('Position (microns)','FontSize',16);
pause(0.1)

res = ((knife_model-knife_data)).^2;
Z = sum(abs(res));
end