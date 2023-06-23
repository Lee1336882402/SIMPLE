clear
clc
M=129;
Re=1000;
dx=1/M;
dy=1/M;
dt=1e-2;
P=ones(M);
U=zeros(M,M+1);%U(i,j)对应U(i-1/2,j)
V=zeros(M+1,M);%V(i,j)对应V(i,j-1/2)
step=0;
Uin=3;
eps=1e-7;
while true
    step=step+1;
    U(:,1)=zeros(M,1);
    U(:,end)=zeros(M,1);
    V(1,:)=zeros(1,M);
    V(end,:)=zeros(1,M);
    alphap=max(zeros(M),(U(:,1:end-1)+U(:,2:end))/2);
    alpham=min(zeros(M),(U(:,1:end-1)+U(:,2:end))/2);
    betap=max(zeros(M),(V(1:end-1,:)+V(2:end,:))/2);
    betam=min(zeros(M),(V(1:end-1,:)+V(2:end,:))/2);
    alphaph=max(zeros(M),(V(2:end,:)+[V(2:end,2:end),-V(2:end,end)])/2);
    alphamh=min(zeros(M),(V(2:end,:)+[V(2:end,2:end),-V(2:end,end)])/2);
    betaph=max(zeros(M),(U(:,2:end)+[U(2:end,2:end);ones(1,M)*2*Uin-U(end,2:end)])/2);
    betamh=min(zeros(M),(U(:,2:end)+[U(2:end,2:end);ones(1,M)*2*Uin-U(end,2:end)])/2);
    a00=dx*dy/dt*ones(M,M)+dy*([alphap(:,2:M),-alpham(:,M)]-alpham+2/Re/dx*ones(M,M))+dx*(alphaph-[zeros(1,M);alphamh(1:M-1,:)]+2/Re/dy*ones(M,M));
    am0=dy*(alphap+ones(M,M)/Re/dx);
    ap0=dy*(-[alpham(:,2:M),-alpham(:,M)]+ones(M,M)/Re/dx);
    a0m=dx*([zeros(M,1),alphaph(:,1:end-1)]+ones(M,M)/Re/dy);
    a0p=dx*(-alphamh+ones(M,M)/Re/dy);
    b00=dx*dy/dt*ones(M,M)+dx*([betap(2:M,:);-betam(M,:)]-betam+2/Re/dy*ones(M,M))+dy*(betaph-[zeros(M,1),betamh(:,1:M-1)]+2/Re/dx*ones(M,M));
    b0m=dx*(betap+ones(M,M)/Re/dy);
    b0p=dx*(-[betam(2:M,:);-betam(M,:)]+ones(M,M)/Re/dy);
    bm0=dy*([zeros(1,M);betaph(1:end-1,:)]+ones(M,M)/Re/dx);
    bp0=dy*(-betamh+ones(M,M)/Re/dx);
    Ppr=P;
    Upr=U(:,2:end);
    Vpr=V(2:end,:);
    instep=0;
    ineps=1e-1;
    if step > 2000
        ineps = 1e-2;
    end
    if step > 10000
        ineps = 1e-3;
    end
    %ineps=1e-1;
    %速度估算
    while true
        instep=instep+1;
        if instep >10000
            error('内迭代1步数过大！');
        end
        Upr0=Upr;
        Vpr0=Vpr;
        for i = 1:M
            for j = 1:M-1
                if i == 1
                    U3=-Upr(1,j);
                else
                    U3=Upr(i-1,j);
                end
                if i == M
                    U4=2*Uin-Upr(M,j);
                else
                    U4=Upr(i+1,j);
                end
                if j == 1
                    U1=0;
                else
                    U1=Upr(i,j-1);
                end
%                 if j==M
%                     Upr(i,j)=0;
%                     continue
%                 else
%                     U2=Upr(i,j+1);
%                 end
                U2=Upr(i,j+1);
                %Upr(i,j)=am0(i,j)*Upr(i,j-1)+ap0(i,j)*Upr(i,j+1)+a0m(i,j)*Upr(i-1,j)+a0p(i,j)*Upr(i+1,j)-deltay*(Ppr(i,j+1)-Ppr(i,j))+deltax*deltay/dt*Upr(i,j);
                Upr(i,j)=(am0(i,j)*U1+ap0(i,j)*U2+a0m(i,j)*U3+a0p(i,j)*U4-dy*(Ppr(i,j+1)-Ppr(i,j))+dx*dy/dt*U(i,j))/a00(i,j);

            end
        end
        for i =1:M-1
            for j =1:M
                if i == 1
                    V1=0;
                else
                    V1=Vpr(i-1,j);
                end
%                 if i == M
%                     Vpr(i,j)=0;
%                     continue
%                 else
%                     V2=Vpr(i+1,j);
%                 end
                V2=Vpr(i+1,j);
                if j == 1
                    V3=-Vpr(i,1);
                else
                    V3=Vpr(i,j-1);
                end
                if j ==M
                    V4 = -Vpr(i,M);
                else
                    V4 = Vpr(i,j+1);
                end
                Vpr(i,j)=(b0m(i,j)*V1+b0p(i,j)*V2+bm0(i,j)*V3+bp0(i,j)*V4-dx*(Ppr(i+1,j)-Ppr(i,j))+dx*dy/dt*V(i,j))/b00(i,j);
                
            end
        end
        if max(max(max(abs(Upr-Upr0))),max(max(abs(Vpr-Vpr0))))<ineps
            break
        end
    end 
    Upr = [zeros(M,1),Upr];
    Vpr = [zeros(1,M);Vpr];
    %a00=deltax*deltay/dt*ones(M,M)+deltay*([alphap(:,2:M),zeros(M,1)]-alpham+2/Re/deltax*ones(M,M))+deltax*(alphaph-[zeros(1,M);alphamh(1:M-1,:)]+2/Re/deltay*ones(M,M))
    cp0=dy^2*ones(M,M)./a00;
    cm0=dy^2*ones(M,M)./[dx*dy/dt*ones(M,1)+dy*(2*alphap(:,1)+2/Re/dx*ones(M,1))+dx*(2/Re/dy*ones(M,1)),a00(:,1:M-1)];
    c0p=dx^2*ones(M,M)./b00;
    c0m=dx^2*ones(M,M)./[dx*dy/dt*ones(1,M)+dx*(2*betap(1,:)+2/Re/dy*ones(1,M))+dy*(2/Re/dx*ones(1,M));b00(1:M-1,:)];
%     cp0=dy^2*ones(M,M)./(dx*dy/dt);cm0=dy^2*ones(M,M)./(dx*dy/dt);
%     c0p=dx^2*ones(M,M)./(dx*dy/dt);c0m=dx^2*ones(M,M)./(dx*dy/dt);
    c00=cp0+cm0+c0p+c0m;
    c00hat=-dy*(Upr(:,2:end)-Upr(:,1:end-1))-dx*(Vpr(2:end,:)-Vpr(1:end-1,:));
    instep=0;
    P_rev=zeros(M,M);
    while true
        instep=instep+1;
        if instep >100000
            error('内迭代2步数过大！');
        end
        P_rev0=P_rev;
        for i = 1:M
            for j = 1:M
                if i ==1
                    if j == 1
                        P_rev(i,j)=(cp0(i,j)*P_rev(i,j+1)+c0p(i,j)*P_rev(i+1,j)+c00hat(i,j))/(c00(i,j)-c0m(i,j)-cm0(i,j));
                        continue
                    end
                    if j == M
                        P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+c0p(i,j)*P_rev(i+1,j)+c00hat(i,j))/(c00(i,j)-c0m(i,j)-cp0(i,j));
                        continue
                    end
                    P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+cp0(i,j)*P_rev(i,j+1)+c0p(i,j)*P_rev(i+1,j)+c00hat(i,j))/(c00(i,j)-c0m(i,j));
                    continue
                end
                if i == M
                    if j == 1
                        P_rev(i,j)=(cp0(i,j)*P_rev(i,j+1)+c0m(i,j)*P_rev(i-1,j)+c00hat(i,j))/(c00(i,j)-c0p(i,j)-cm0(i,j));
                        continue
                    end
                    if j == M
                        P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+c0m(i,j)*P_rev(i-1,j)+c00hat(i,j))/(c00(i,j)-c0p(i,j)-cp0(i,j));
                        continue
                    end
                    P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+cp0(i,j)*P_rev(i,j+1)+c0m(i,j)*P_rev(i-1,j)+c00hat(i,j))/(c00(i,j)-c0p(i,j));
                    continue
                end
                if j == 1
                    P_rev(i,j)=(c0m(i,j)*P_rev(i-1,j)+cp0(i,j)*P_rev(i,j+1)+c0p(i,j)*P_rev(i+1,j)+c00hat(i,j))/(c00(i,j)-cm0(i,j));
                    continue
                end
                if j == M
                    P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+c0m(i,j)*P_rev(i-1,j)+c0p(i,j)*P_rev(i+1,j)+c00hat(i,j))/(c00(i,j)-cp0(i,j));
                    continue
                end
                P_rev(i,j)=(cm0(i,j)*P_rev(i,j-1)+cp0(i,j)*P_rev(i,j+1)+c0p(i,j)*P_rev(i+1,j)+c0m(i,j)*P_rev(i-1,j)+c00hat(i,j))/c00(i,j);
            end
        end
        if max(max(abs(P_rev-P_rev0)))<ineps
            break
        end
    end
    U_new=Upr-dy*[zeros(M,1),(P_rev(:,2:end)-P_rev(:,1:end-1))./a00(:,1:end-1),zeros(M,1)];
    V_new=Vpr-dx*[zeros(1,M);(P_rev(2:end,:)-P_rev(1:end-1,:))./b00(1:end-1,:);zeros(1,M)];
    P_new=Ppr+P_rev;
    if max(max(max(abs(U_new-U))),max(max(abs(V_new-V))))<eps
        break
    end
    U=U_new;
    V=V_new;
    P=P_new;
    %测试
%     if step==20
%         break
%     end
    if step > 500000
        error('整体迭代步数过大！');
    end
end
%%
%展示结果
U_show=(U(:,1:end-1)+U(:,2:end))/2;
V_show=(V(1:end-1,:)+V(2:end,:))/2;
[X ,Y]=meshgrid(linspace(0,1,M),linspace(0,1,M));
figure(1)
streamslice(U_show,V_show)
figure(2)
contour(X,Y,P)
figure(3)
quiver(X,Y,U_show,V_show);