//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

//a = read("/home/naudin/Documents/FichierScilab/Fourre tout/Fig1A_AIYCurrentClampTrace2.txt",-1,12);
a = read("/scilab-scripts/Fig 1A_AIY Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:10)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:25];

//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];
InfSD=[2.88 2.55 1.47 1.31 1.04 0.809 0.7 0.658 0.638 0.889 1.94 3.5 5.36 7.63 10.6 13.3 16 17.9]

//Mono-objective solution
bM=[0.123860430031553
   0.157354280722996
   0.214036806923923
   0.139693180935831
   128.6006770757747
  -98.14825263991114
  -41.07483204558734
  -19.09574559613729
  -21.92789650381577
  -89.84199258265514
  -17.71413240455412
   4.565037462232217
  -15.11575568039892
  -30.
   7.395161906785207
   0.0001
   10.59767882678912
   0.001632023755844
   0.279407407379098
   0.804888377484526
   0.660772380318157
   0.042294111350645]

//////////////////////////////////////////////////
///////////////    Cost function    //////////////
//////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Ca,t+Kir+K,p+L-model
function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*x(3)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(10),pa(14))*(x(1)-pa(6)) - pa(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(9),pa(13))-x(3))/pa(17)
    Hdot(4)=(xinf(x(1),pa(11),pa(15))-x(4))/pa(18)
endfunction

//Function that computes the standard deviation
function y=sigma(v)
    s=0;
    moy=mean(v);
    for i=1:length(v)
        s=s+(v(i)-moy)^2
    end
    y=sqrt(s/(length(v)-1));
endfunction

//Noise level (standard deviation) for each I
dev=[]
dev1=sigma(A(500:$,1));
dev=[dev1]
for i=2:length(stim)
    dev=[dev sigma(A(5000:$,i))]
end

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Cost function voltage
function y=fct11(pa)
    tmp=0;
    condini = [-53; pa(19); pa(20); pa(21)]
    for i=1:length(stim)
        c=0;
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i));
        end
        c=sqrt(c/length(t))/dev(i);
        tmp=tmp+c;
    end
    y=tmp/length(stim);
endfunction

///////////////////////////////////////////////////////////////////////
///////////////    Steady-state current cost function    //////////////
///////////////////////////////////////////////////////////////////////

function y=WSS(pa)
    e=0;
    for i=1:length(vecV)
        tmp=0;
        tmp=(Inf(i)-(pa(1)*xinf(vecV(i),pa(8),pa(12))*xinf(vecV(i),pa(9),pa(13))*(vecV(i)-pa(5)) + pa(2)*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(6)) + pa(3)*xinf(vecV(i),pa(11),pa(15))*(vecV(i)-pa(6)) + pa(4)*(vecV(i)-pa(7))))^2
        tmp=tmp/InfSD(i)
        e=e+tmp;
    end
    y=e/length(vecV) 
endfunction

///////////////////////////////////////////////////////////////
/////////    Crowding Sorting and Domination Front    /////////
///////////////////////////////////////////////////////////////

function [Front]=NDS(A)
    dominationCount = zeros(size(A,'r'),1)
    S=list(); // S(1) is composed of all solutions' indexes dominated by i
    Front=list(); // F(1) is composed of solutions' indexes of the front 1, F(2), is composed of solutions' indexes of the front 2, etc.
    for i=1:size(A,'r')
        Stmp=[]; // set of solutions dominated by the solution i
        for j=1:size(A,'r')
            if i~=j then 
                // number of solution which dominates the solution i 
                if A(i,1)>A(j,1) & A(i,2)>A(j,2) then
                    dominationCount(i) = dominationCount(i) + 1;
                end
                // set of solution dominated by the solution i
                if A(i,1)<A(j,1) & A(i,2)<A(j,2) then
                    Stmp=[Stmp j]
                end
            end
        end
        S(i)=Stmp
    end
    Front(1)=find(0==dominationCount); // index of solutions belonging to F(1)
    
    // determining all next fronts
    m=1
    while Front(m)~=[] // while the front m is non-empty 
        Q=[];
        for i=Front(m)
            for j=S(i)
                dominationCount(j) = dominationCount(j) - 1;
                if dominationCount(j)==0 then
                    Q=[Q j];
                end
            end
        end
        m=m+1;
        Front(m)=Q;
    end
endfunction

function [d]=crowdingSorting(A)
    l=size(A, 1); // l=number of individual in A (=set of objective functions of the last acceptable front)
    M=size(A,'c'); // M=number of objective function
    d = zeros(l, 1);
    for m=1:M
        [tmp, Index] = gsort(A(:, m)); // Step C2 : sort the set in ascendant order of magnitude
//        pause;
        d(Index(1)) = %inf;
        d(Index(l)) = %inf;
        fmax = max(A(:, m));
        fmin = min(A(:, m));
//        pause;
        for j=2:l-1
            d(Index(j)) = d(Index(j)) + abs(tmp(j+1) - tmp(j-1)) / (fmax - fmin);
        end
//        pause;
    end
endfunction

//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

function [popInit, valInit, pop2500, val2500, pop5000, val5000, popFinal, valFinal]=simulation(NP,itermax,F,CR)
  
    D=22; 
    pop=zeros(D,NP);

    ///////////////////////////
    //// Bound constraints ////
    ///////////////////////////
    
    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 -30 1  0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  -1  30 20     20     20     0.999 0.999 0.999 10];
    
    ////////////////////////////////////
    ////  Population initialization ////
    ////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    // Save popInit
    popInit=pop;
    //Integration of mono-objective solution in random position 
    pop(:,floor(1+NP*rand()))=bM
    disp(pop);
    
    ///////////////////////////////////////
    //// Initial population evaluation ////
    ///////////////////////////////////////
    
    val=zeros(NP,2); // tableau avec le co??t de chacun des individus. 1??re colonne = cout voltage. 2??me colonne = cout SS.
    
    for j=1:NP
        val(j,1)=fct11(pop(:,j))
        val(j,2)=WSS(pop(:,j))
    end
    
    disp(val);
    
    // Save valInit
    valInit=val;

    ///////////////////
    //// Next step ////
    ///////////////////
     
    iter=1; // number of iteration
    U=zeros(D,NP); // matrix resulting of mutation + crossover
    tempvalVol=0;
    tempvalSS=0;
    while iter<itermax
        for j=1:NP
            // Building of the matrix U

            // Random selection of two different integers r1 and r2, different from j as well
            r1=j; r2=j;//////////////////////////////////////
            while (r1==r2 | r1==j | r2==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
            end
            // best member on the voltage cost function
            [tmp, Index] = gsort(val(:,1),'g','i');
                        
            // Differential variation
            V=pop(:,Index(1)) + F*(pop(:,r1)-pop(:,r2));
            
            // Constraints
            for i=1:length(Xmin)
                if V(i)<=Xmin(i) then V(i)=Xmin(i);
                elseif V(i)>Xmax(i) then V(i)=Xmax(i);
                end
            end
            // Crossover
            for i=1:D
                if rand()<CR then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // end for j=1:NP
        // Adding of children U in the pop if they dominate parents or if they are not dominated. 
        tempPop=pop;        
        tempval=val;
        
        for j=1:NP
            tempvalVol = fct11(U(:,j));
            tempvalSS = WSS(U(:,j));
            if tempvalVol<tempval(j,1) & tempvalSS<=tempval(j,2) then
                tempPop(:,j) = U(:,j);
                tempval(j,1) = tempvalVol;
                tempval(j,2) = tempvalSS;
            end
            if (tempvalVol>tempval(j,1) & tempvalSS<=tempval(j,2)) | (tempvalVol<tempval(j,1) & tempvalSS>=tempval(j,2)) then
                tempPop=[tempPop U(:,j)]
                tempval=[tempval; [tempvalVol tempvalSS]]
            end
        end
        
        // Front ranking of tempPop > NP
        [Front]=NDS(tempval);
        
        // Integration of fronts in the pop until its size is exceeded
        pop=[];
        val=[];
        k=1;
        while (size(pop,2)+length(Front(k)))<NP
            for i=1:length(Front(k))
                pop=[pop tempPop(:,Front(k)(i))];
                val=[val; tempval(Front(k)(i),:)];
            end
            k=k+1;
        end
       
        // Compute crowding distance of the last considered front
        lastFront=[];
        for i=1:length(Front(k))
            lastFront=[lastFront; tempval(Front(k)(i),:)];
        end
        
        cs=crowdingSorting(lastFront);//Asignation d'une distance de crowding
        
        // Integration of individuals according to their crowding distance
        [osef, indice]=gsort(cs);
        n=1;
        while size(pop,2)<NP
            pop=[pop tempPop(:,Front(k)(indice(n)))];
            val=[val; tempval(Front(k)(indice(n)),:)];
            n=n+1;
        end
        
        if iter==700 then
            disp(pop);
            disp(val);
            pop2500=pop;
            val2500=val;
        end
        if iter==1500 then
            disp(pop);
            disp(val);
            pop5000=pop;
            val5000=val;
        end

        if (iter==3 | iter==200 | iter==400 | iter==1000 | iter==1200 | iter==1800) then
            disp(pop);
            disp(val);
        end

        disp(iter);
        iter = iter + 1;
    end  //end while
    
    popFinal=pop;
    valFinal=val;
    disp(pop);
    disp(val);
endfunction




