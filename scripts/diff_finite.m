                       %Anno Accademico 2012/2013
                   %Corso di Elaborazione di Immagini
            
        %Gruppo 4: Chiara Capurro, Giulia Marconi, Mara Scussolini

%Esercitazione 4:
%Integrazione di informazione da bioimmagini funzionali ed anatomiche. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               ROUTINE DIFFERENZE FINITE DI UN'IMMAGINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: - A matrice di cui si vogliono calcolare le derivate parziali con 
%          le differenze finite e il rispettivo modulo del gradiente
%        - type stringa per scegliere se operare con:
%          differenze finite centrali --> 'centrali' 
%          differenze finite avanti --> 'avanti'
%          differenze finite indietro --> 'indietro'

% OUTPUT: - Dx matrice della derivata a riga fissata rispetto alle colonne 
%         - Dy matrice della derivata a colonna fissata rispetto alle righe
%         - G matrice del rispettivo modulo del gradiente

%DIFFERENZE FINITE <--> Filtro spaziale con maschere Wx Wy di
%                       dimensione 3x3 con entrate specifiche.
%                       (--> vedi commento DIFFERENZE FINITE -line 145-)

function [Dx,Dy,G]=diff_finite(A,type)

         Dx=zeros(size(A)); %inizializzazione matrice derivata parziale Dx (rispetto alle colonne)
         Dy=zeros(size(A)); %inizializzazione matrice derivata parziale Dy (rispetto alle righe)
         G=zeros(size(A)); %inizializzazione matrice modulo gradiente G
         
         L=3; %dimensione maschere fissata LxL=3x3

         %-----------------------------------------------------------------
         %OPERAZIONE DI PERIODICITA':
         %si replica la matrice A dell'immagine originale in tutte le
         %direzioni e si considera quella centrale a cui si aggiungono
         %(L-1)/2 righe e (L-1)/2 colonne per parte; in questo modo la
         %finestra LxL agisce su tutti i pixel dell'immagine di partenza
         %==> si aggiungono righe/colonne della stessa matrice e pertanto
         %si mantengono i livelli di grigio originali, non si altera il
         %contrasto dell'immagine e la trasformazione rimane più fedele 
         %rispetto all'operazione di ZERO-PADDING (--> vedi NB -line 111-).
         %-----------------------------------------------------------------
         
         n_plus=(L-1)/2; %righe/colonne da aggiungere alla matrice originale    
         A_rep=repmat(A,3,3); %si replica la matrice A in tutte le direzioni 
                              %<--> 3 volte in riga, 3 volte in colonna
         
         %Per isolare nella replicazione la matrice centrale con (L-1)/2 
         %righe e (L-1)/2 colonne aggiunte per parte, si calcola il numero 
         %di righe/colonne da scartare della matrice replicata:                   
         scarto_righe=size(A,1)-n_plus+1; 
         scarto_colonne=size(A,2)-n_plus+1;
         
         %Si seleziona la regione della matrice replicata di interesse:
         A_grande=A_rep(scarto_righe:size(A_rep,1)-scarto_righe+1, scarto_colonne:size(A_rep,2)-scarto_colonne+1);
        
         %Scelta del metodo alle differenze finite da utilizzare (passato in input): 
         if strcmp(type,'centrali')==true
             fai=true;
             %Maschere Wx e Wy delle differenze finite centrali:
             Wx=[0 0 0 ; -1/2 0 1/2 ; 0 0 0]; 
             Wy=[0 -1/2 0 ; 0 0 0 ; 0 1/2 0];
         elseif strcmp(type,'avanti')==true
             fai=true;
             %Maschere Wx e Wy delle differenze finite avanti:
             Wx=[0 0 0 ; 0 -1 1 ; 0 0 0];
             Wy=[0 0 0 ; 0 -1 0 ; 0 1 0];
         elseif strcmp(type,'indietro')==true
             fai=true;
             %Maschere Wx e Wy delle differenze finite indietro:
             Wx=[0 0 0 ; -1 1 0 ; 0 0 0];
             Wy=[0 -1 0 ; 0 1 0 ; 0 0 0]; 
         else 
             fai=false;
             Wx=[]; Wy=[]; %nel caso di controllo fallito, le maschere Wx e Wy sono matrici vuote
             Dx=[]; Dy=[]; G=[]; %nel caso di controllo fallito, le matrici di output Dx, Dy e G sono vuote
             disp('WARNING function ''diff_finite''!')
             disp('input errato: [Dx,Dy,G]=diff_finite(A,type)')
             disp('              ''type'' stringa per la definizione del tipo di differenze finite:')
             disp('              input accettati: - centrali')
             disp('                               - avanti')
             disp('                               - indietro') 
         end
         
         if fai==true %se il controllo del tipo di differenza finita passato in input è andato a buon fine:
             
             for i=1:size(A_grande,1)-L+1     %|==> si sposta la maschera su tutta la matrice A_grande
                 for j=1:size(A_grande,2)-L+1 %|    nelle due direzioni colonne/righe
                     
                     submatrix=A_grande(i:i+L-1,j:j+L-1); %sottomatrice di A_grande di dimensione LxL su cui agisce la maschera
                     Dx(i,j)=sum(sum(submatrix.*Wx)); %| ==> l'entrata (i,j) della derivata parziale è la somma delle entrate
                     Dy(i,j)=sum(sum(submatrix.*Wy)); %|     della matrice corrispondente al prodotto puntale (entrata per entrata,
                                                      %|     non l'usuale prodotto matriciale righe/colonne) tra la sottomatrice 
                                                      %|     di A_grande e la maschera stessa
                 
                 end
             end

             %Modulo del gradiente:
             G=sqrt(Dx.^2+Dy.^2);
             
         end

end

%--------------------------------------------------------------------------
%NB - Modo alternativo per espandere la matrice - :
%OPERAZIONE DI ZERO-PADDING:
%LxL dimensione maschera ==> si aggiungono alla matrice A dell'immagine 
%                            originale (L-1)/2 righe e (L-1)/2 colonne 
%                            per parte di ZERI in modo tale che la
%                            finestra agisca su tutti i pixel 
%                            dell'immagine di partenza. 
%==> aggiungendo righe/colonne di zeri non si mantengono i livelli di 
%    grigio originali, si altera il contrato dell'immagine e la
%    trasformazione è meno fedele rispetto all'operazione di PERIODICITA'.

%A_grande=[zeros(n_plus, size(A,2)+2*n_plus) ; zeros(size(A,1), n_plus) A zeros(size(A,1), n_plus); zeros(n_plus, size(A,2)+2*n_plus)]; 
%--------------------------------------------------------------------------

%__________________________________________________________________________
%OPERATORI LOCALI:
%Gli operatori locali sono trasformazioni che cambiano il livello di grigio
%dell'immagine ma il risultato dell'operazione su ciascun pixel dipende dal 
%pixel stesso e da quelli a lui adiacenti.
%Per determinare quanto locale è l'operazione bisogna definire una finestra 
%di dimensione LxL con L dispari che agisce attorno allo specifico pixel.
%L'operatore locale muove questa finestra (detta "maschera") su tutta 
%l'immagine e modifica solo il pixel centrale della finestra basandosi sui
%restanti valori dei pixel in essa presenti. 
%Allo scorrere della finestra sull'immagine si incontrano pixel sia 
%modificati sia non modificati; noi lavoriamo con operatori locali in 
%parallelo, ovvero si considerano i valori originali dei pixel.
%Partendo da un'immagine MxM, dopo la trasformazione l'immagine risulta più 
%piccola (M-L+1) x (M-L+1) in quanto non tutti i pixel sono raggiunti dalla 
%trasformazione; infatti il centro della finestra non può raggiungere il 
%bordo dell'immagine. Per risolvere questo problema si applica zero-padding 
%o periodicità. 

%__________________________________________________________________________
%DIFFERENZE FINITE: 
%Per il calcolo delle derivate parziali di un immagine f si ricorre ai
%metodi alle differenze finite.
%Posto f(i,j) il valore dell'immagine nel pixel (i,j), il metodo delle:

%DIFERENZE FINITE CENTRALI è definito come:
% - derivata parziale f_x a riga fissata i:
%   f_x(i,j)= (1/2)*f(i,j+1) - (1/2)*f(i,j-1) ;
% - derivata parziale f_y a colonna fissata j:
%   f_y(i,j)= (1/2)*f(i+1,j) - (1/2)*f(i-1,j) .

%Implementativamente operare le differenze finite centrali su un'immagine 
%equivale ad applicare il filtro spaziale con maschere (=matrici dei pesi):
%Wx=[0 0 0 ; -1/2 0 1/2 ; 0 0 0] e Wy=[0 -1/2 0 ; 0 0 0 ; 0 1/2 0].

%DIFERENZE FINITE AVANTI è definito come:
% - derivata parziale f_x a riga fissata i:
%   f_x(i,j)= f(i,j+1) - f(i,j) ;
% - derivata parziale f_y a colonna fissata j:
%   f_y(i,j)= f(i+1,j) - f(i,j) .

%Implementativamente operare le differenze finite avanti su un'immagine 
%equivale ad applicare il filtro spaziale con maschere (=matrici dei pesi):
%Wx=[0 0 0 ; 0 -1 1 ; 0 0 0] e Wy=[0 0 0 ; 0 -1 0 ; 0 1 0].

%DIFERENZE FINITE INDIETRO è definito come:
% - derivata parziale f_x a riga fissata i:
%   f_x(i,j)= f(i,j) - f(i,j-1) ;
% - derivata parziale f_y a colonna fissata j:
%   f_y(i,j)= f(i,j) - f(i-1,j) .

%Implementativamente operare le differenze finite indietro su un'immagine 
%equivale ad applicare il filtro spaziale con maschere (=matrici dei pesi):
%Wx=[0 0 0 ; -1 1 0 ; 0 0 0] e Wy=[0 -1 0 ; 0 1 0 ; 0 0 0].
%__________________________________________________________________________

