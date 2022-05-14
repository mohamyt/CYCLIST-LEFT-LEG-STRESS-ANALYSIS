% Knowns & Assumptions
%Angles and Forces
aAnkle=readmatrix('Ankle Angle.xlsx');
aCrank=readmatrix('Crank Angle.xlsx');
aHip=readmatrix('Hip Angle.xlsx');
fHorizontal=readmatrix('Horizontal Pedal Force.xlsx');
fVertical=readmatrix('Vertical Pedal Force.xlsx');
aKnee=readmatrix('Knee Angle.xlsx');
aPedal=readmatrix('Pedal Angle.xlsx');
%Constant assumptions
gHip=9.81*9.411;
gFoot=9.81*1.2817;
gShank=9.81*5.54;
lFoot=0.273;
lShank=0.51051;
lHip=0.52961;
bFoot=0.04914;
hFoot=0.07098;
rShank=0.039585;
rHip=0.031395;
%Second moments of area
ixFoot=bFoot*power(hFoot,3)/12;
iyFoot=hFoot*power(bFoot,3)/12;
ixShank=pi*power(rShank,4)/4;
iyShank=ixShank;
ixHip=pi*power(rHip,4)/4;
iyHip=ixHip;


%Coefficients of F to get maximum shear stress (Q/It) or (4/3*(Area of circle)) or (4/3*(Area of ellipse))
fcFoot=0.5*hFoot*bFoot/(ixFoot*bFoot);
fcShank=4/(3*pi*rShank*rShank);
fcHip=4/(3*pi*rHip*rHip);

%calculating stresses

n=length(aCrank);
maxNSfoot=zeros(n,1);
maxTfoot=zeros(n,1);
maxNSshank=zeros(n,1);
maxTshank=zeros(n,1);
maxNShip=zeros(n,1);
maxThip=zeros(n,1);
maxSfoot=zeros(n,1);
maxSshank=zeros(n,1);
maxShip=zeros(n,1);
USfoot=0;
indFoot=0;
indShank=0;
indHip=0;
USshank=0;
UShip=0;
tfoot=zeros(n,1);
tshank=zeros(n,1);
thip=zeros(n,1);
Ffoot=zeros(n,1);
Fshank=zeros(n,1);
Fhip=zeros(n,1);
mzfoot=zeros(n,1);
mzshank=zeros(n,1);
mzhip=zeros(n,1);



for i=1:n
    %Maximum stresses in foot
    footH=-fHorizontal(i)*cos(aPedal(i)) - fVertical(i)*sin(aPedal(i));
    footV=-fHorizontal(i)*sin(aPedal(i)) - fVertical(i)*cos(aPedal(i));
    wHfoot=gFoot*sin(aPedal(i));
    wVfoot=gFoot*cos(aPedal(i));
    mzfoot(i)=(0.5*lFoot*footV)+(0.5*lFoot*(footV-wVfoot));           %Maximum moment at section a-a of the foot
    RVfoot=wVfoot-footV;
    RHfoot=-wHfoot-footH;
    Ffoot(i)=sqrt(RHfoot*RHfoot+RVfoot*RVfoot);
    tfoot(i)=0;                                                     %Shear stress at section a-a
    sxfoot=mzfoot(i)*(hFoot/2)/ixFoot -RHfoot/(bFoot*hFoot);        %Normal stress foot in x direction
    syfoot=-RVfoot/(bFoot*hFoot);                                   %Normal stress foot in y direction
    maxNSfoot(i)=((sxfoot+syfoot)/2)+sqrt(power((sxfoot-syfoot)/2,2)+power(tfoot(i),2)); %Maximum normal stress
    maxTfoot(i)=sqrt(power(syfoot-sxfoot,2)+power(2*tfoot(i),2))/2;   %Maximum shear stress
    stress=((sxfoot+syfoot)/2)-sqrt(power((sxfoot-syfoot)/2,2)+power(tfoot(i),2));
    if abs(stress)>abs(maxNSfoot(i))
        maxNSfoot(i)=stress;
    end
    if abs(maxNSfoot(i))<abs(maxTfoot(i))
        maxSfoot(i)=abs(maxTfoot(i));
    else
        maxSfoot(i)=abs(maxNSfoot(i));
    end
    if USfoot < maxNSfoot(i)
        USfoot=maxNSfoot(i);
        indFoot = i;
        sigmaxfoot=sxfoot;
        sigmayfoot=syfoot;
    end
    %Maximum stresses in shank
    shankH=RVfoot*cos(90-aAnkle(i)) + RHfoot*sin(90-aAnkle(i));
    shankV=RVfoot*sin(90-aAnkle(i)) + RHfoot*cos(90-aAnkle(i));
    wHshank=gShank*sin(90-aAnkle(i)-aPedal(i));
    wVshank=gShank*cos(90-aAnkle(i)-aPedal(i));
    mzshank(i)=-mzfoot(i)+(0.5*lShank*shankV)+(0.5*lShank*(shankV-wVshank));           %Maximum moment at section b-b of the foot
    RVshank=wVshank-shankV;
    RHshank=-wHshank-shankH;
    Fshank(i)=sqrt(RHshank*RHshank+RVshank*RVshank);
    tshank(i)=0;                                                            %Shear stress at section b-b
    sxshank=mzshank(i)*(rShank)/ixShank -RHshank/(pi*rShank*rShank);        %Normal stress foot in x direction
    syshank=-RVshank/(pi*rShank*rShank);                                %Normal stress foot in y direction
    maxNSshank(i)=((sxshank+syshank)/2)+sqrt(power((sxshank-syshank)/2,2)+power(tshank(i),2)); %Maximum normal stress
    maxTshank(i)=sqrt(power(syshank-sxshank,2)+power(2*tshank(i),2))/2;   %Maximum shear stress
    stress=((sxshank+syshank)/2)-sqrt(power((sxshank-syshank)/2,2)+power(tshank(i),2));
    if abs(stress)>abs(maxNSshank(i))
        maxNSshank(i)=stress;
    end
    if abs(maxNSshank(i))<abs(maxTshank(i))
        maxSshank(i)=abs(maxTshank(i));
    else
        maxSshank(i)=abs(maxNSshank(i));
    end
    if USshank < maxNSshank(i)
        USshank=maxNSshank(i);
        indShank = i;
        sigmaxshank=sxshank;
        sigmayshank=syshank;
    end
    %Maximum stresses in hip
    hipH=RVshank*sin(aKnee(i)) - RHshank*cos(aKnee(i));
    hipV=-RVshank*cos(aKnee(i)-90) - RHshank*cos(aKnee(i)-90);
    wHhip=gHip*cos(aHip(i));
    wVhip=gHip*sin(aHip(i));
    mzhip(i)=-mzshank(i)+(0.5*lHip*hipV)+(0.5*lHip*(hipV-wVhip));           %Maximum moment at section c-c of the foot
    RVhip=wVhip-hipV;
    RHhip=-wHhip-hipH;
    Fhip(i)=sqrt(RHhip*RHhip+RVhip*RVhip);
    thip(i)=0;                                         %Shear stress at section c-c
    sxhip=mzhip(i)*(rHip)/ixHip -RHhip/(pi*rHip*rHip);        %Normal stress foot in x direction
    syhip=-RVhip/(pi*rHip*rHip);                                %Normal stress foot in y direction
    maxNShip(i)=((sxhip+syhip)/2)+sqrt(power((sxhip-syhip)/2,2)+power(thip(i),2)); %Maximum normal stress
    maxThip(i)=sqrt(power(syhip-sxhip,2)+power(2*thip(i),2))/2;   %Maximum shear stress
    stress=((sxhip+syhip)/2)-sqrt(power((sxhip-syhip)/2,2)+power(thip(i),2)); %Maximum normal stress
    if abs(stress)>abs(maxNShip(i))
        maxNShip(i)=stress;
    end
    if abs(maxNShip(i))<abs(maxThip(i))
        maxShip(i)=abs(maxThip(i));
    else
        maxShip(i)=abs(maxNShip(i));
    end
    if UShip < maxNShip(i)
        UShip=maxNShip(i);
        indHip = i;
        sigmaxhip=sxhip;
        sigmayhip=syhip;
    end
end
%plotting data
figure('Name','Maximum Normal σ Stress with Crank Angle');
title('Maximum Normal σ Stress with Crank Angle');
hold on
plot (aCrank,maxNSfoot/power(10,6),'r');
plot (aCrank,maxNSshank/power(10,6),'g');
plot (aCrank,maxNShip/power(10,6),'b');
ylabel('Maximum Normal Stress σ (MPa)');
xlabel('Crank Angle (degrees)');
legend ('Foot','Shank','Hip');
figure('Name','Maximum Shear Stress τ with Crank Angle');
title('Maximum Shear Stress τ with Crank Angle');
hold on
plot (aCrank,maxTfoot/power(10,6),'r');
plot (aCrank,maxTshank/power(10,6),'g');
plot (aCrank,maxThip/power(10,6),'b');
legend ('Foot','Shank','Hip');
ylabel('Maximum Shear Stress τ (MPa)');
xlabel('Crank Angle (degrees)');
figure('Name','Reaction at Joints with respect to Crank Angle');
title('Reaction at Joints with respect to Crank Angle');
hold on
plot (aCrank,Ffoot,'r');
plot (aCrank,Fshank,'g');
plot (aCrank,Fhip,'b');
legend ('Ankle','Knee','Hip');
ylabel('Reaction Force at joint (N)');
xlabel('Crank Angle (degrees)');
figure('Name','Bending Moment at Joints with respect to Crank Angle');
title('Bending Moment at Joints with respect to Crank Angle');
hold on
plot (aCrank,mzfoot,'r');
plot (aCrank,mzshank,'g');
plot (aCrank,mzhip,'b');
legend ('Foot','Shank','Hip');
ylabel('Moment reaction at joint (N.m)');
xlabel('Crank Angle (degrees)');
figure('Name','Mohr Circle for Maximum Stress in Foot');
title('Mohr Circle for Maximum Stress in Foot');
hold on
pCircle(((sigmayfoot+sigmaxfoot)/2)/power(10,6),0,maxTfoot(indFoot)/power(10,6))
plot([sigmayfoot/power(10,6),sigmaxfoot/power(10,6)],[tfoot(indFoot)/power(10,6),-tfoot(indFoot)/power(10,6)],'r')
ylabel('τ (MPa)');
xlabel('σ (MPa)');
figure('Name','Mohr Circle for Maximum Stress in Shank');
title('Mohr Circle for Maximum Stress in Shank');
hold on
pCircle(((sigmayshank+sigmaxshank)/2)/power(10,6),0,maxTshank(indShank)/power(10,6))
plot([sigmayshank/power(10,6),sigmaxshank/power(10,6)],[tshank(indShank)/power(10,6),-tshank(indShank)/power(10,6)],'r')
ylabel('τ (MPa)');
xlabel('σ (MPa)');
figure('Name','Mohr Circle for Maximum Stress in Hip');
title('Mohr Circle for Maximum Stress in Hip');
hold on
pCircle(((sigmayhip+sigmaxhip)/2)/power(10,6),0,maxThip(indHip)/power(10,6))
plot([sigmayhip/power(10,6),sigmaxhip/power(10,6)],[thip(indHip)/power(10,6),-thip(indHip)/power(10,6)],'r')
ylabel('τ (MPa)');
xlabel('σ (MPa)');