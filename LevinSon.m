
function [a]=LevinSon(message_code,p)

gamma=autocorrelation(message_code,0:p); %fonction d'autocorr�lation du message
%avec "p" la taille du message p=2060
a=gamma(2)/gamma(1);%premier coefficient de r�gression du filtre 

for i=2:p   % on d�marre l'it�ration � 2 jusqu'a 2060 
    beta=gamma(i+1)-a'*flipud(gamma(2:i)');    
    alpha=gamma(1)-a'*(gamma(2:i)');
    Kn=beta/alpha;     
    a_fin=a-Kn*flipud(a); %Coefficiant de regression 
   
    a=[a_fin;Kn]; %concat�nation des coefficiants de regregression qui donne le  vecteur de regression 
end

