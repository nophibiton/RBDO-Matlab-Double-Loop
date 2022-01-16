function [xIP,weight] = GaussLaguerrePoints(nIP)
%GAUSS locations and weights of Gauss Laguerre integration scheme     
% [xIP,weight] = GaussLaguerre(nIP)
%
% xIP    : Gauss-Laguerre integration point locations 
% weight : Gauss-Laguerre integration weights
% nIP    : number of integration points
%
xIP    = zeros(nIP,1);
weight = zeros(nIP,1);

switch nIP
case 1
   xIP (1)   = 1.0;
   weight(1) = 1.0;
   
case 2
   xIP (1)   = 0.585786437626905;
   weight(1) = 0.853553390593274;
   xIP (2)   = 3.41421356237309;
   weight(2) = 0.146446609406726;
   
case 3
   xIP (1)   = 0.415774556659002;
   weight(1) = 0.711093009954157;
   xIP (2)   = 2.29428036046204;
   weight(2) = 0.278517733544999;
   xIP (3)   = 6.28994508287896;
   weight(3) = 0.010389256500845;
   
case 4
   xIP (1)   = 0.322547688627706;
   weight(1) = 0.603154025190419;
   xIP (2)   = 1.7457609373338;
   weight(2) = 0.357418778179241;
   xIP (3)   = 4.53662055670391;
   weight(3) = 0.038887901973767;
   xIP (4)   = 9.39507081733458;
   weight(4) = 0.000539294656572;
   
case 5
   xIP(1)    = 0.263560317199723;
   weight(1) = 0.521755613197929;
   xIP(2)    = 1.41340306420443;
   weight(2) = 0.398666808306289;
   xIP(3)    = 3.59642576749151;
   weight(3) = 0.075942449840528;
   xIP(4)    = 7.08581000694528;
   weight(4) = 0.003611758682921;
   xIP(5)    = 12.6408008441591;
   weight(5) = 0.000023369972334;
   
case 6
   xIP(1)    =  0.222846599957007;
   weight(1) =  0.458964681100265;
   xIP(2)    =  1.18893211096868;
   weight(2) =  0.417000822772701;
   xIP(3)    =  2.99273631836904;
   weight(3) =  0.113373382929006;
   xIP(4)    =  5.77514357233998;
   weight(4) =  0.010399197447799;
   xIP(5)    =  9.83746741771355;
   weight(5) =  0.000261017202320;
   xIP(6)    = 15.9828739806518;
   weight(6) =  0.000000898547910;
   
   
case 8
   xIP(1)    = -0.960289856497536;
   weight(1) =  0.101228536290376;
   xIP(2)    = -0.796666477413627;
   weight(2) =  0.222381034453374;
   xIP(3)    = -0.525532409916329;
   weight(3) =  0.313706645877887;
   xIP(4)    = -0.183434642495650;
   weight(4) =  0.362683783378362;
   xIP(5)    = -xIP(4);
   weight(5) =  weight(4);  
   xIP(6)    = -xIP(3);
   weight(6) =  weight(3);
   xIP(7)    = -xIP(2);
   weight(7) =  weight(2);
   xIP(8)    = -xIP(1);
   weight(8) =  weight(1);
   
case 9
   xIP (1)   = -0.9681602395;
   weight(1) =  0.0812743884;
   xIP (2)   = -0.8360311073;
   weight(2) =  0.1806481607;
   xIP (3)   = -0.6133714327;
   weight(3) =  0.2606106964;
   xIP (4)   = -0.3242534234;
   weight(4) =  0.3123470770;
   xIP (5)   =  0.;
   weight(5) =  0.3302393550;
   xIP (6)   = -xIP (4);
   weight(6) =  weight(4);  
   xIP (7)   = -xIP (3);
   weight(7) =  weight(3);
   xIP (8)   = -xIP (2);
   weight(8) =  weight(2);
   xIP (9)   = -xIP (1);
   weight(9) =  weight(1);
   
case 10
   xIP (1)   = -0.9739065285;
   weight(1) =  0.0666713443;
   xIP (2)   = -0.8650633667;
   weight(2) =  0.1494513492;
   xIP (3)   = -0.6794095683;
   weight(3) =  0.2190863625;
   xIP (4)   = -0.4333953941;
   weight(4) =  0.2692667193;
   xIP (5)   = -0.1488743390;
   weight(5) =  0.2955242247;
   xIP (6)   = -xIP (5);
   weight(6) =  weight(5);  
   xIP (7)   = -xIP (4);
   weight(7) =  weight(4);
   xIP (8)   = -xIP (3);
   weight(8) =  weight(3);
   xIP (9)   = -xIP (2);
   weight(9) =  weight(2);
   xIP (10)  = -xIP (1);
   weight(10)=  weight(1);
   
case 16
   xIP(1)    =  0.9894009349916499;
   weight(1) =  0.0271524594117541;
   xIP(2)    =  0.9445750230732326;
   weight(2) =  0.0622535239386479;
   xIP(3)    =  0.8656312023878317;
   weight(3) =  0.0951585116824928;
   xIP(4)    =  0.7554044083550030;
   weight(4) =  0.1246289712555339;
   xIP(5)    =  0.6178762444026437;
   weight(5) =  0.1495959888165767;
   xIP(6)    =  0.4580167776572274;
   weight(6) =  0.1691565193950025;
   xIP(7)    =  0.2816035507792589;
   weight(7) =  0.1826034150449236;
   xIP(8)    =  0.0950125098376374;
   weight(8) =  0.1894506104550685;
   
   xIP(1:8)  = -xIP(1:8);
   xIP(9:16) = -xIP(8:-1:1);
   weight(9:16) =  weight(8:-1:1);
   
otherwise
   error ('integration order not supported')
   stop
end
return
