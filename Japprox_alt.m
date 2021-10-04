
function J=Japprox_alt(alp,bet,the)
Ma=Mapprox(alp,alp+bet+1,the); 
Mb=Mapprox(alp,bet+1,the);
J = ((alp+1)/(2*bet))*Ma + ((the-1)/(2*bet))*Mb;
end