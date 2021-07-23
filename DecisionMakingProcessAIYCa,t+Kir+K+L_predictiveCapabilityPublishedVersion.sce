//////////////////////////////////////////////////////////
///////////////    Voltage cost function    //////////////
//////////////////////////////////////////////////////////

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
for i=2:11
    dev=[dev sigma(A(5000:$,i))]
end

//Cost function voltage for validation
stim=[30];
B=A(:,10);
devV=[sigma(A(5000:$,10))]

function y=fct11V(pa)
    tmp=0;
    condini = [-53; pa(19); pa(20); pa(21)]
    for i=1:length(stim)
        c=0;
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-B(k,i))*(V(k)-B(k,i));
        end
        c=sqrt(c/length(t))/devV(i);
        tmp=tmp+c;
    end
    y=tmp/length(stim);
endfunction

////////////////////////////////////////////////////
///////////////   Results reading    ///////////////
////////////////////////////////////////////////////

result=[]; //concatenation of all final results
pop=[]; // concatenation of all final population

//Results for DEMO/rand/best/biased
for i=14:24;
    tmpResult=csvRead('/home/naudin/Documents/article-2/RIMAIYAFDmonoObjRandBestNP=600F=1.5CR=0.3/Results-20210513BiObjectifSSVoltageAIYCa,t+Kir+K+LparalleleUntil25pAfirstIterationRandBestNP=600F=1.5CR=0.3/Results/run_EstimationBiObjectifSSVoltageAIYCa,t+Kir+K+LparalleleUntil25pAfirstIterationRandBest.sce_' + string(i) + '/valFinal.csv');
    tmpPop=csvRead('/home/naudin/Documents/article-2/RIMAIYAFDmonoObjRandBestNP=600F=1.5CR=0.3/Results-20210513BiObjectifSSVoltageAIYCa,t+Kir+K+LparalleleUntil25pAfirstIterationRandBestNP=600F=1.5CR=0.3/Results/run_EstimationBiObjectifSSVoltageAIYCa,t+Kir+K+LparalleleUntil25pAfirstIterationRandBest.sce_' + string(i) + '/popFinal.csv');
    result=[result; tmpResult];
    pop=[pop tmpPop];
//    disp(result);
//    disp(pop);
end

///////////////////////////////////////////////////////////////////////
///////////////   Front pareto recuperation (Step 1)    ///////////////
///////////////////////////////////////////////////////////////////////

[f_pareto,X_pareto]=pareto_filter(result,pop')
pop = X_pareto';

//////////////////////////////////////////////////////////////////////////
//////   Selecting solution with a correct bifurcation structure    //////
//////////////////////////////////////////////////////////////////////////

function y=xinf(VH,V12,k)
    y = 1./(1+exp((V12-VH) ./k))
endfunction

function y=alpha(V)
    y = gCa.*xinf(V,V12mCa,kmCa).*xinf(V,V12hCa,khCa)
endfunction

function y=betta(V)
    y = gKir.*xinf(V,V12hKir,kKir)
endfunction

function y=delta(V)
    y = gK.*xinf(V,V12mK,kmK)
endfunction

function y=fmCa(V)
    y = (1/kmCa)*exp((V12mCa-V) ./kmCa).*xinf(V,V12mCa,kmCa)
endfunction

function y=fhCa(V)
    y = (1/khK)*exp((V12hCa-V) ./khCa).*xinf(V,V12hCa,khCa)
endfunction

function y=fhKir(V)
    y = (1/kKir)*exp((V12hKir-V) ./kKir).*xinf(V,V12hKir,kKir)
endfunction

function y=fmK(V)
    y = (1/kmK)*exp((V12mK-V) ./kmK).*xinf(V,V12mK,kmK)
endfunction

function y=IinfDer(V)
    y = alpha(V).*[1+(fmCa(V)+fhCa(V)).*(V-ECa)] + betta(V).*[1+fhKir(V).*(V-EK)] + delta(V).*[1+fmK(V).*(V-EK)] + gL
endfunction

vecV=[-100:0.01:50];
popGoodBifStructure=[];

for i=1:size(pop,2)
    pa=pop(:,i);
    
    //Model parameters
    gCa=pa(1); gKir=pa(2); gK=pa(3); gL=pa(4);
    ECa=pa(5); EK=pa(6); EL=pa(7);
    V12mCa=pa(8); V12hCa=pa(9); V12hKir=pa(10); V12mK=pa(11); 
    kmCa=pa(12); khCa=pa(13); kKir=pa(14); kmK=pa(15); 
    tmCa=pa(16); thCa=pa(17); tmK=pa(18); 
    mCa0=pa(19); mK0=pa(20); hK0=pa(21);
    C=pa(22);
        
    if min(IinfDer(vecV))>0 then
        popGoodBifStructure=[popGoodBifStructure pop(:,i)]
    end
end


///////////////////////////////////////////////////////////////////////////
///////   Selection of solution that best fit validation trace     ////////
/////////////////////////b//////////////////////////////////////////////////

valBest=fct11V(popGoodBifStructure(:,1))
bM=popGoodBifStructure(:,1)
for i=2:size(popGoodBifStructure,2)
    if fct11V(popGoodBifStructure(:,i))<valBest then
        valBest=fct11V(popGoodBifStructure(:,i))
        bM=popGoodBifStructure(:,i);
    end
    disp(i)
end




