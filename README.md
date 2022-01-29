# Traitement-du-signal
# Rapport TP1: Analyse spectrale d’un signal et Transformée de Fourier discrète

## réalisé par Wafae Hamdaoui et Oumayma EL Bakkali

## Objectifs du TP:

      • Représentation de signaux et applications de la transformée de Fourier discrète
        (TFD) sous Matlab. 

       • Evaluation de l’intérêt du passage du domaine temporel au domaine fréquentiel
          dans l’analyse et l’interprétation des signaux physiques réels.

## Représentation temporelle et fréquentielle 

 ***Q1,Q2,Q3:***

 ```Matlab
     clear all
     close all
     clc

    % on déclare l'intervalle de temps
t=0:1/50:10-1/50;

% Puis on déclare l'intervalle de fréquence
f=0:50:500*50-50;

% on déclare notre signal périodique qui est constitué d'une somme de deux sinusoïdes
x=sin(2*pi*15*t)+sin(2*pi*20*t);

% maintenant on veux tracer le signal x(t) Avec un Pas de temps : Te = 1/50s, et un  Intervalle : 0, 10-Te.
% puis on veux  tracer une approximation de la TF continue du signal x(t), et la représenté suivant un pas Te,
% on utilise les deux commandes : fft et fftshif.

% Calculons la TFD du signal x(t) en utilisant la commande fft
y = fft (x);

%  Pour mieux visualiser le contenu fréquentiel du signal, on peux utiliser la fonction fftshift, 
% qui effectue un décalage circulaire centré sur zéro du spectre en amplitude obtenu par 
% la commande fft.
z= fftshift(y);

%tracer le signal x(t):
 subplot(3,3,1),plot(t,x);
% tracer le spectre en amplitude de la TFD du signal:
 subplot(3,3,2), plot(f,abs(y));
 subplot(3,3,3), plot(f,abs(z));

```
> Output:
![image](https://user-images.githubusercontent.com/75392302/151665752-1d2d1160-f746-4cbe-84df-95f41ba0b4cd.png)

Pour bien voir le signal on clique sur zoom in:

![image](https://user-images.githubusercontent.com/75392302/151665766-060aa038-af14-40ef-8e5f-9cc9c9c0c90a.png)




▪  On remarquera que la TF est une fonction complexe et que la fonction ainsi
obtenue décrit la TF de x(t) entre –1/(2Te) et 1/(2Te) par pas de 1/(nTe) où n
est le nombre de points constituant le signal x(t).

▪ La commande fft codant les fréquences positives sur les n/2 premières valeurs
du signal et les valeurs négatives entre n/2+1 et n, la commande fftshift permet
de les inverser.

▪ On remarquera que la fonction fftshift, qui effectue un décalage circulaire centré 
sur zéro du spectre en amplitude obtenu par la commande fft.


## L'ajout d'un bruit


*Un bruit correspond à tout phénomène perturbateur gênant la transmission ou
l'interprétation d'un signal. Dans les applications scientifiques, les signaux sont
souvent corrompus par du bruit aléatoire, modifiant ainsi leurs composantes
fréquentielles. La TFD peut traiter le bruit aléatoire et révéler les fréquences qui y
correspond.*

```Matlab

close all 
clc
% Le Pas de temps
Te=1/50;

% Le Pas de fréquence
fe= 50; 

%l'intervalle de temps
t= 0:Te:2-Te;
% le nombre de points constituant le signal
N=length(t);
% le signal qui contient l'information
x=sin(2*pi*15*t)+sin(2*pi*20*t);
% La TFD du signal
y = fft (x);
% L'intervalle de la fréquence
f=[0:N-1]*(fe/N);
% plot(f,abs(y));
fshift=[-N/2:N/2-1]*(fe/N);

% Créons un nouveau signal xnoise, en introduisant un bruit blanc gaussien dans 
% le signal d’origine x
xnoise=x+randn(size(x));
% La TFD de xnoise
ynoise=fft(xnoise);

% Traçons tous les signaux résultants:
subplot(2,2,1)
plot(t,x);
subplot(2,2,2)
plot(f,abs(y));
subplot(2,2,3)
plot(t,xnoise);

% Calculons puis traçons le spectre de puissance du signal bruité centré à la fréquence zéro
subplot(2,2,4)
plot(fshift,abs((ynoise).^2)/500);

```

la commande randn pour générer le bruit. Il est à noter qu’un bruit blanc est une réalisation d'un processus aléatoire dans lequel la densité spectrale de puissance est la même pour toutes les fréquences de la bande passante. Ce bruit suit une loi normale de moyenne et variance données.

> Output:

***Signal original:***

![image](https://user-images.githubusercontent.com/75392302/151665985-1c2ffc14-2927-4b70-a3ed-064bbed210fb.png)

***Son Transformée :***

![image](https://user-images.githubusercontent.com/75392302/151666006-96eeb4e5-0495-4eda-bcdb-f04da37fe207.png)

	Les deux  pics de puissance perment de distinguer les fréquenses du signal.

***Signal bruité:***

![image](https://user-images.githubusercontent.com/75392302/151666037-3cd80924-fd22-4999-8ba7-7dcdd2961dcf.png)

***Son Transformée :***

![image](https://user-images.githubusercontent.com/75392302/151666053-ee68006b-bd73-49a6-9462-1af1ad934eec.png)

	Il est possible de distinguer les fréquences du signal en raison des pics de puissance.

**Le spectre de puissance du signal bruité:***

![image](https://user-images.githubusercontent.com/75392302/151666095-fe9dea30-79ea-49d5-9f64-e8ae7bad6ede.png)

	Malgré le bruit on voit qu’il est toujours possible de distinguer les fréquences du signal en raison des pics de puissance.
En plus, on voit que, entre 0Hz et 15Hz, la densité spectrale de puissance est relativement constante. Ceci est dû au bruit blanc gaussien.

## Analyse fréquentielle du chant du rorqual bleu 

-*Il existe plusieurs signaux dont l’information est encodée dans des sinusoïdes. Les
ondes sonores est un bon exemple. Considérons maintenant des données audios
collectées à partir de microphones sous - marins au large de la Californie. On cherche
à détecter à travers une analyse de Fourier le contenu fréquentiel d’une onde sonore
émise pas un rorqual bleu.*

Depuis le fichier ‘bluewhale.au’, le sous-ensemble de données qui
correspond au chant du rorqual bleu du Pacifique. En effet, les appels de rorqual bleu
sont des sons à basse fréquence, ils sont à peine audibles pour les humains. Utiliser
la commande audioread pour lire le fichier. Le son à récupérer correspond aux indices
allant de 2.45e4 à 3.10e4.



```Matlab

whaleFile = fullfile(matlabroot,'examples','matlab','data','bluewhale.au');
[w,ft] = audioread(whaleFile);
chant=w(2.45e4:3.10e4);
subplot(2,1,1);
soundsc(w,fs);

Nchant=length(chant);
t=[0:Nchant-1]*1/ft;
plot(t,chant);
title('courbe du signal bluewhale');


```

La TFD peut être utilisée pour identifier les composantes fréquentielles de ce signal
audio. Dans certaines applications qui traitent de grandes quantités de données avec
fft, il est courant de redimensionner l'entrée de sorte que le nombre d'échantillons soit
une puissance de 2. fft remplit automatiquement les données avec des zéros pour
augmenter la taille de l'échantillon. Cela peut accélérer considérablement le calcul de
la transformation.

```Matlab

%la densité spectrale de puissance du signal
fshift=[-Nchant/2:Nchant/2-1]*(ft/Nchant)/10;
DSP=abs(fft(chant).^2/Nchant);
subplot(2,1,2);
plot(fshift,fftshift(DSP));
title('la densité spectrale de puissance du signal')
```


* On peut Déterminer à partir du tracé, la fréquence fondamentale du gémissement de rorqual
bleu.*

![image](https://user-images.githubusercontent.com/85129301/149664604-28481b50-bfbc-49a2-8cbd-1ade3023551f.png)
