
###CC0293-An?lise Multivariada- Prof. Mauricio#


##X1=comprimento do osso aos 8 anos de idade
##X2=comprimento do osso aos 8 anos e meio de idade
##X2=comprimento do osso aos 9 anos   de idade
##X4=comprimento do osso aos 9 anos e meio de idade


###X=(X1,X2,X3,X4)'~N_4(mu,Sigma), mu=(mu1,mu2,mu3,mu4)


####Dados amostrais


X_1=c(478,488,490,497)/10
X_2=c(464,473,477,484)/10
X_3=c(463,468,478,485)/10
X_4=c(451,453,461,472)/10
X_5=c(476,485,489,493)/10
X_6=c(525,532,533,537)/10
X_7=c(512,530,543,545)/10
X_8=c(498,500,503,527)/10
X_9=c(481,508,523,544)/10
X_10=c(450,470,473,483)/10
X_11=c(512,514,516,519)/10
X_12=c(485,492,530,555)/10
X_13=c(521,528,537,550)/10
X_14=c(482,489,493,498)/10
X_15=c(496,504,512,518)/10
X_16=c(507,517,527,533)/10
X_17=c(472,477,484,495)/10
X_18=c(533,546,551,553)/10
X_19=c(462,475,481,484)/10
X_20=c(463,476,513,518)/10
dad=rbind(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,X_15,X_16,
X_17,X_18,X_19,X_20)
dad
colnames(dad)=c("t1","t2","t3","t4");dad

p=ncol(dad)###n?mero de vari?veis
p
n=nrow(dad)### tamanho da amostra
n


C=matrix(c(1,1,1,1,1,-1,1,-2,0,0,0,-3),ncol=4);C

postoC=qr(C)$rank;postoC




plot(dad)
##item a. Estudo longitudinal


##item b.

Xb=matrix(colMeans(dad),ncol=1);Xb

S_Y=C%*%S%*%t(C);S_Y


n=nrow(dad);n
p=ncol(dad);p


T2_cal=n*t((C%*%Xb))%*%solve(C%*%S%*%t(C))%*%(C%*%Xb);T2_cal
aux= (n-p+1)/((p-1)*(n-1));aux

F_cal= aux*T2_cal;F_cal

alfa=0.05

F_tab=qf(1-alfa,p-1,n-p+1);F_tab
T2_tab=(1/aux)*F_tab;T2_tab

nd=1-pf(F_cal,p-1,n-p+1);nd



S=cov(dad);round(S,4)

R=cor(dad);round(R,4)

##item c.
VTE=sum(diag(S));VTE
VGE=det(S);VGE




Yb=C%*%Xb;Yb


##item d. Z=CY; EZ=C*E(Y),CovZ=C*Sigma*t(C)


C=matrix(c(1,-2,2,3,1,-1,-3,2),ncol=4);C

Zb=C%*%Yb;Zb
SZ=C%*%S%*%t(C);SZ

RZ=sqrt(solve(diag(diag(SZ))))%*%SZ%*%sqrt(solve(diag(diag(SZ))));RZ

##item e.    W=C1*Y

C1=matrix(c(2,-2,3,3,-1,-2,-1,4,-1,4,-2,3),ncol=4);C1

Wb=C1%*%Yb;Wb
SW=C1%*%S%*%t(C1);SW
RW=sqrt(solve(diag(diag(SW))))%*%SW%*%sqrt(solve(diag(diag(SW))));RW






###Cria??o da fun??o transforma


Trans=function(dad,C,d){
dadT=dad%*%t(C) + kronecker(matrix(1,nrow(dad),1),t(d))
return(dadT)
}D_12=t((Y_1 -Yb))%*%solve(S)%*%(Y_1 -Yb);D_12


	
dadZ=Trans(dad,C,d);dadZ


matrix(colMeans(dadZ),ncol=1);Zb

cov(dadZ);SZ


cor(dadZ);RZ




d1=matrix(0,3,1);d1


dadW=Trans(dad,C1,d1);dadW



matrix(colMeans(dadW),ncol=1);Wb

cov(dadW);SW


cor(dadW);RW


S
Xb=matrix(colMeans(dad),ncol=1);Xb
####

D_12=t((X_1 -Xb))%*%solve(S)%*%(X_1 -Xb);D_12
D_22=t((Y_2 -Yb))%*%solve(S)%*%(Y_2 -Yb);D_22
D_32=t((Y_3 -Yb))%*%solve(S)%*%(Y_3 -Yb);D_32
D_42=t((Y_4 -Yb))%*%solve(S)%*%(Y_4 -Yb);D_42
D_52=t((Y_5 -Yb))%*%solve(S)%*%(Y_5 -Yb);D_52
D_62=t((Y_6 -Yb))%*%solve(S)%*%(Y_6 -Yb);D_62
D_72=t((Y_7 -Yb))%*%solve(S)%*%(Y_7 -Yb);D_72
D_82=t((Y_8 -Yb))%*%solve(S)%*%(Y_8 -Yb);D_82
D_92=t((Y_9 -Yb))%*%solve(S)%*%(Y_9 -Yb);D_92
D_102=t((Y_10 -Yb))%*%solve(S)%*%(Y_10 -Yb);D_102
D_112=t((Y_11 -Yb))%*%solve(S)%*%(Y_11 -Yb);D_112
D_122=t((Y_12 -Yb))%*%solve(S)%*%(Y_12 -Yb);D_122
D_132=t((Y_13 -Yb))%*%solve(S)%*%(Y_13 -Yb);D_132
D_142=t((Y_14 -Yb))%*%solve(S)%*%(Y_14 -Yb);D_142
D_152=t((Y_15 -Yb))%*%solve(S)%*%(Y_15 -Yb);D_152
D_162=t((Y_16 -Yb))%*%solve(S)%*%(Y_16 -Yb);D_162
D_172=t((Y_17 -Yb))%*%solve(S)%*%(Y_17 -Yb);D_172
D_182=t((Y_18 -Yb))%*%solve(S)%*%(Y_18 -Yb);D_182
D_192=t((Y_19 -Yb))%*%solve(S)%*%(Y_19 -Yb);D_192
D_202=t((Y_20 -Yb))%*%solve(S)%*%(Y_20 -Yb);D_202


D2=matrix(c(D_12,D_22,D_32,D_42,D_52,D_62,D_72,D_82,D_92,D_102,
D_112,D_122,D_132,D_142,D_152,D_162,D_172,D_182,D_192,D_202),ncol=1);D2

D2_n=max(D2);D2_n


U=n*D2/(n-1)^2;U
alfa=(p-2)/(2*p)
beta=(n-p-2)/(2*(n-p-1))
i=1:n

V=(i-alfa)/(n-alfa-beta +1)
u=sort(U);u

plot(u,V)


##Fazer Direto

aux=kronecker(matrix(1,20,1),t(Yb));aux

D2=(dad -aux)%*%solve(S)%*%t((dad -aux));D







#####Vamos transformar os dados originais usando a matriz C




C=matrix(c(1,1,1,1,1,-1,1,-2,0,-3,0,0),ncol=4);C


Yb=C%*%Xb;Yb
S_Y=C%*%S%*%t(C);S_Y



T2_cal=n*t((C%*%Xb))%*%solve(C%*%S%*%t(C))%*%(C%*%Xb);T2_cal

aux= (n-p+1)/((p-1)*(n-1));aux

F_cal= aux*T2_cal;F_cal

alfa=0.05

F_tab=qf(1-alfa,p-1,n-p+1);F_tab
T2_tab=(1/aux)*F_tab;T2_tab

nd=1-pf(F_cal,p-1,n-p+1);nd





dad_Y=dad%*%t(C)
rownames(dad_Y)=c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9","Y10",
"Y11","Y12","Y13","Y14","Y15","Y16","Y17","Y18","Y19","Y20")


dad_Y


matrix(colMeans(dad_Y),ncol=1);Yb

cov(dad_Y);S_Y


n=nrow(dad_Y);n
p=ncol(dad_Y);p ###### (  s? temos 3 vari?veis agora)


T2_cal=n*t(Yb)%*%solve(S_Y)%*%(Yb);T2_cal

aux= (n-p)/(p*(n-1));aux

F_cal= aux*T2_cal;F_cal

alfa=0.05

F_tab=qf(1-alfa,p,n-p);F_tab
T2_tab=(1/aux)*F_tab;T2_tab

nd=1-pf(F_cal,p,n-p);nd






























C_1=matrix(c(1,1,1,-1,0,0,0,-1,0,0,0,-1),ncol=4);C_1

postoC_1=qr(C_1)$rank;postoC_1




Xb=matrix(colMeans(dad),ncol=1);Xb

mu_est_Y=C_1%*%Xb;mu_est_Y

SC1_Y=C_1%*%S%*%t(C_1);SC1_Y


n=nrow(dad);n
p=ncol(dad);p


T2_cal=n*t((C_1%*%Xb))%*%solve(C_1%*%S%*%t(C_1))%*%(C_1%*%Xb);T2_cal

aux= (n-p+1)/((p-1)*(n-1));aux

F_cal= aux*T2_cal;F_cal

alfa=0.05

F_tab=qf(1-alfa,p-1,n-p+1);F_tab
T2_tab=(1/aux)*F_tab;T2_tab

nd=1-pf(F_cal,p-1,n-p+1);nd



####################C_2##########################################



C_2=matrix(c(1,0,0,-1,1,0,0,-1,1,0,0,-1),ncol=4);C_2

postoC_2=qr(C_2)$rank;postoC_2




Xb=matrix(colMeans(dad),ncol=1);Xb

C2mu_est_Y=C_2%*%Xb;C2mu_est_Y

SC2_Y=C_2%*%S%*%t(C_2);SC2_Y


n=nrow(dad);n
p=ncol(dad);p


T2_cal=n*t((C_2%*%Xb))%*%solve(C_2%*%S%*%t(C_2))%*%(C_2%*%Xb);T2_cal

aux= (n-p+1)/((p-1)*(n-1));aux

F_cal= aux*T2_cal;F_cal

alfa=0.05

F_tab=qf(1-alfa,p-1,n-p+1);F_tab
T2_tab=(1/aux)*F_tab;T2_tab

nd=1-pf(F_cal,p-1,n-p+1);nd






