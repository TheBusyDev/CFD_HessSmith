function [xd, yd, xv, yv, xe, ye] = NACA_4d_a(M, P, SS, n, varargin)
%
% NACA_4d.m
%
% Generazione dei punti sul dorso e sul ventre di un profilo alare NACA a 4
% cifre (four digits) a partire dalla linea media  linea_media_4d(x, M, P) 
% e dallo spessore  spessore(x, SS)  del profilo simmetrico:
%
% M  : ordinata massima della linea media             (una cifra)
%
% P  : posizione sulla corda dell'ordinata massima M  (una cifra)
%
% SS : spessore massimo del profilo simmetrico        (due cifre)
%
% n  : numero di punti lungo la corda distribuiti uniformemente,
%      inclusi i punti del bordo di attacco e del bordo di uscita,  
%      in corrispondenza dei quali si calcolano i punti sul dorso
%      e sul ventre
% 
%
% M, P e SS possono essere forniti come NUMERI INTERI, compresi, 
%      rispettivamente,
%
%      M  fra 1 e 9,  in  %  della corda      
%      P  fra 1 e 9   in 1/10  della corda   
%      SS fra 1 e 99  in  %  della corda
%
%      oppure, se si preferisce, come VALORI REALI  < 1  
%      frazione della corda
%
% REF:
% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections, pp 113
% Dover, New York, 1949, 1959
%

%

% Input
if M >= 1
   m = M/100; 
else  
   m = M;
end

if P >= 1
   p = P/10;  % in 1/10 
else  
   p = P;
end

if SS >= 1
   s = SS/100; 
else  
   s = SS;
end

if (~isempty (varargin))
    spacing = varargin{1};
else
    spacing = 'constant';
end

xd = zeros(n+1,1);
yd = zeros(n+1,1);

xv = zeros(n+1,1);
yv = zeros(n+1,1);

k = n+1;

switch spacing
    
   case 'constant'
        
      for i = 1 : k
         
         x = (i-1)/(k-1);
         
         [y_lm, Dy_lm] = linea_media_4d(x, m, p);
         
         y_sp = spessore(x, s);
         
         theta = atan(Dy_lm);
         
         xd(i) = x  -  y_sp * sin(theta); 
         xv(i) = x  +  y_sp * sin(theta);
         
         yd(i) = y_lm  +  y_sp * cos(theta); 
         yv(i) = y_lm  -  y_sp * cos(theta); 

      end
        
   case 'halfcos'
        
      n = n*2+1;
      k = n+1;
      i = 1:k;
      x = 1-0.5*(1+cos(((i-1)*pi)/n));
      x = x(1:ceil(length(x)/2))/0.5;
      
      for i = 1 : length(x)
         
         [y_lm, Dy_lm] = linea_media_4d(x(i), m, p);
         
         y_sp = spessore(x(i), s);
         
         theta = atan(Dy_lm);
         
         xd(i) = x(i)  -  y_sp * sin(theta); 
         xv(i) = x(i)  +  y_sp * sin(theta);
         
         yd(i) = y_lm  +  y_sp * cos(theta); 
         yv(i) = y_lm  -  y_sp * cos(theta);
         
      end
    
   case 'cos'
        
      for i = 1 : k
         
         x = 1-0.5*(1+cos(((i-1)*pi)/n));
         
         [y_lm, Dy_lm] = linea_media_4d(x, m, p);
         
         y_sp = spessore(x, s);
         
         theta = atan(Dy_lm);
         
         xd(i) = x  -  y_sp * sin(theta); 
         xv(i) = x  +  y_sp * sin(theta);
         
         yd(i) = y_lm  +  y_sp * cos(theta); 
         yv(i) = y_lm  -  y_sp * cos(theta); 

      end
        
end
    
   if (xd(end,1) ~= 1)
      xd(end+1) = 1;
      yd(end+1) = 0;
   end

   if (xv(end,1) ~= 1)
      xv(end+1) = 1;
      yv(end+1) = 0;
   end

   
   xe=[flipud(xv); xd(2:end)];
   ye=[flipud(yv); yd(2:end)];

return

%==========================================================================
%==========================================================================
% linea_media_4d del profilo NACA a 4 cifre (four digits)

function [y, Dy] = linea_media_4d(x, M, P)

% Calcolo dell'ordinata della linea media  y(x)  e della
% sua derivata  dy(x)/dx  nel punto  x  lungo la corda
% del profilo NACA a 4 cifre con curvatura:  NACA-MP__
%
% Le coordinate x e y sono adimensionali rispetto alla corda.
%
% M :  ordinata massima della linea media             (una cifra)
%
% P :  posizione sulla corda dell'ordinata massima M  (una cifra)
%
%
% M e P possono essere forniti come NUMERI INTERI compresi 
%            
%      M fra 1 e 9,  in  %  della corda      
%
%      P fra 1 e 9   in  1/10  della corda   
%
%      oppure, se si preferisce, come VALORI REALI  < 1  
%      frazione della corda
%           

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 114


if M >= 1
   m = M/100; 
else  
   m = M;
end


if P >= 1
   p = P/10;  % in 1/10
else  
   p = P;
end


if x < p

    y = (m/p^2) * (2*p - x) * x;

   Dy = (m/p^2) * 2*(p - x); 

else   

    y = (m/(1-p)^2) * (1 - 2*p + (2*p - x) * x);
  
   Dy = (m/(1-p)^2) * 2*(p - x); 
   
end

return

%==========================================================================
%==========================================================================
% spessore di un profilo NACA (con 4 o 5 cifre)

function [y] = spessore(x, SS)

% Valore dello spessore  y  nel punto  x  lungo la corda
% per un determinato spessore massimo SS (due cifre)
%
% Le coordinate x e y sono adimensionali rispetto alla corda.
%
% Lo spessore massimo SS del profilo puo' essere fornito
% come NUMERO INTERO compreso fra 1 e 99  in  %  della corda
% oppure come VALORE REALE  < 1  come frazione della corda

% Ad esempio, entrambe le specificazioni  SS = 12  e  SS = 0.12 
% producono lo spessore del profilo simmetrico  NACA0012   

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 113


if SS >= 1
   s = SS/100; 
else  
   s = SS;
end


y = 5 * s * (0.29690 * sqrt(x)  -  0.12600 * x  -  0.35160 * x^2 ...
                                +  0.28430 * x^3  -  0.10360 * x^4); %0.10150 * x^4);


return

%==========================================================================
%==========================================================================
% Derivata dello spessore di un profilo NACA (con 4 o 5 cifre)

function [Dy] = d_spessore_dx(x, SS)

% Valore della derivata dello spessore  y  nel punto  x  lungo 
% la corda per un determinato spessore massimo SS (due cifre)
%
% Le coordinate x e y sono adimensionali rispetto alla corda.
%
% Lo spessore massimo SS del profilo puo' essere fornito
% come NUMERO INTERO compreso fra 1 e 99  in  %  della corda
% oppure come VALORE REALE  < 1  come frazione della corda

% Ad esempio, entrambe le specificazioni  SS = 12  e  SS = 0.12 
% producono la derivata dello spessore del profilo simmetrico  
% NACA0012   

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 113


if SS >= 1
   s = SS/100; 
else  
   s = SS;
end


Dy = 5 * s * (0.14845 / sqrt(x)  -  0.12600  -  0.7032 * x  ...
                                 +  0.8529 * x^2  - 0.4144*x^3);  %0.406 * x^3);

return