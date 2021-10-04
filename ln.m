function l=ln(x)
% ln(x) outputs the natural log of x for all x~=0, and when x==0 ln(x) 
% outputs the limit of x*log(x) as x->0, ie 0. 
if x==0
    l=0;
else
    l=log(x);
end
end