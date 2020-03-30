
function [ uout,xout,xextout,yout,yextout] = sim_NM_network(Ttot,deltat,Nout,NM_network)
% Simulates a model of coupled (Wendling) neural masses using the settings
% in the struct NM_network. External input to a neural mass is added to all
% populations, where alpha, beta and gamma scale the input to the local
% excitatory, slow inhibitory and fast inhibitory population w.r.t. the
% pyramidal cells, respectively. Differential equations are numerically
% solved using the Euler-Maruyama method.
%
%Input:
%   Ttot: (positive double) Total simulation time (seconds)
%   deltat: (positive double) Time step size (seconds)
%   Nout: (positive integer) output every Nout-th timestep
%   NM_network: struct containing networkstructure and parameters of all
%   NMs.
%
%Output:
%   uout: (4*Nnodes) x NTSTout matrix containing the simulated potential of
%   the neural masses. Here Nnodes is the number of NMs and NTSTout the
%   number of output timesteps (=ceil(NTST/Nout)+1);
%   Rows contain in order the potential of the following
%   cells: [u_py^(1),u_ex^(1),u_is^(1),u_if^(1),u_py^(2),u_ex^(2),etc].
%   So rows 1:4:end contain the potential of the pyramidal cells.
%   xout: (4*Nnodes) x NTSTout matrix containing the PSPs originating from
%   py, ex, is, if (x_1,..,x_4 in the manuscript). Rows contain in order the potential of the following
%   cells: [x_1^(1),x_2^(1),x_3^(1),x_4^(1),x_1^(2),x_2^(2),etc].
%   xextout: Nnodes x NTSTout matrix containing the PSP originating from
%   timevarying input and noise. Note that the input from pyramidal cells
%   of other NMs influencing a NMs is not contained in this variable. These
%   are incorporated by adding the PSP of the pyramidal cells of the other
%   NMs directly when calculating the potential of the populations.
%   yout: (4*Nnodes) x NTSTout matrix containing derivatives of PSPs. The
%   order of the elements is equal to xout.
%   yextout: Nnodes x NTSTout matrix containing the derivatives of xext.
%
% This code is part of the simulation code of the manuscript 'Pathological responses to single pulse electrical stimuli in epilepsy: the role of feedforward inhibition'
% (c) 2019 Jurgen Hebbink (University of Twente, University Medical Centre Utrecht)


NTST=ceil(Ttot/deltat);             %number of timestepes
Nnodes=size(NM_network.nodes,1);    %number of nodes


%% Pre-allocate state variables
x=zeros(4*Nnodes,1);
y=zeros(4*Nnodes,1);
xext=zeros(Nnodes,1);
yext=zeros(Nnodes,1);

%% Write DVs of the derivatives of the PSPs (y's) in quasi-linear representation dy=Q1x+Q2y+Q3*sigm(u)
% Due to the specific form of the ODE's Q1, Q2 and Q3 are diagonal
% matrices. We represent these as vectors and do elementwise multiplication
% rather than matrix multiplication as matrix multiplication is
% computationally more expensive.
Q1=-deltat*reshape([NM_network.nodes.a;NM_network.nodes.a;NM_network.nodes.b;NM_network.nodes.g].^2,4*Nnodes,1);
Q2=1-deltat*reshape(2*[NM_network.nodes.a;NM_network.nodes.a;NM_network.nodes.b;NM_network.nodes.g],4*Nnodes,1);
tvA=[NM_network.nodes.A].*[NM_network.nodes.a];
tvB=[NM_network.nodes.B].*[NM_network.nodes.b];
tvG=[NM_network.nodes.G].*[NM_network.nodes.g];
Q3=deltat*reshape([tvA;tvA;tvB;tvG],4*Nnodes,1);

%% External input
%same matrices as for the state variables
Q1ext=-deltat*([NM_network.nodes.a].^2)';
Q2ext=1-deltat*2*[NM_network.nodes.a]';

% Vectors containing scaling factors of the input
alpha=[NM_network.nodes.alpha]';
beta=[NM_network.nodes.beta]';
gamma=[NM_network.nodes.gamma]';

% Vector containing time varying input
fvar=deltat*repmat(tvA',1,NTST).*cell2mat(cellfun(@(c) c(deltat:deltat:Ttot),{NM_network.nodes.Ivar}','UniformOutput',false));

%% Stochastic input
% normal situation independent noise
fs=sqrt(deltat)*repmat((tvA.*[NM_network.nodes.sd])',1,NTST).*randn(Nnodes,NTST);

% comment or uncomment following lines for same noise on neural masses 2:end
% disp('Warning: channels 2:end receive same noise and noise is init from seed');
% rng(0);
% tempnoise=nan(Nnodes,NTST);
% tempnoise(1,:)=randn(1,NTST);
% tempnoise(2:end,:)=repmat(randn(1,NTST),Nnodes-1,1);
% fs=sqrt(deltat)*repmat((tvA.*[NM_network.nodes.sd])',1,NTST).*tempnoise; %stochastic input
%%

%% Initialize matrix to compute potentials from PSPs.
temp=zeros([4,4,Nnodes]);
temp(1,2,:)=cat(3,NM_network.nodes.C2);
temp(1,3,:)=-cat(3,NM_network.nodes.C4);
temp(1,4,:)=-cat(3,NM_network.nodes.C7);
temp(2,1,:)=cat(3,NM_network.nodes.C1);
temp(3,1,:)=cat(3,NM_network.nodes.C3);
temp(4,1,:)=cat(3,NM_network.nodes.C5);
temp(4,3,:)=-cat(3,NM_network.nodes.C6);
temp(4,4,:)=-cat(3,NM_network.nodes.C8);

if(Nnodes>1)
    temp=mat2cell(temp,4,4,ones(Nnodes,1));    
    P=sparse(blkdiag(temp{:}));
else
    P=temp;
end

u=P*x;      % vector with potentials of the populations
su=Sigm(u); % vector with firing rates of populations

%% Pre-allocate output variables
uout=zeros(4*Nnodes,ceil(NTST/Nout)+1);
uout(:,1)=u;
xout=zeros(4*Nnodes,ceil(NTST/Nout)+1);
xextout=zeros(Nnodes,ceil(NTST/Nout)+1);
yout=zeros(4*Nnodes,ceil(NTST/Nout)+1);
yextout=zeros(Nnodes,ceil(NTST/Nout)+1);

%% Network structure
conmat=NM_network.constrength*NM_network.conmat;

indupy=1+(0:Nnodes-1)*4;
induex=2+(0:Nnodes-1)*4;
induins=3+(0:Nnodes-1)*4;
induinf=4+(0:Nnodes-1)*4;
indpyEu=1+(0:Nnodes-1)*4; %indices of the excitatory potential evoked by the pyramidal cells
indEpy=2+(0:Nnodes-1)*4;

%%
tel=0;          % counter for saving output   
indexplot=2;    % index of output

%loop performing the Euler-Maruyama iterations
for ii=1:NTST
    xn=x+deltat*y;                                      %New x
    yn=Q1.*x+Q2.*y+Q3.*su;                              %New y
    xextn=xext+deltat*yext;                             %New xextn
    yextn=Q1ext.*xext+Q2ext.*yext+fvar(:,ii)+fs(:,ii);  %New yextn
    
    %update vars
    x=xn;
    y=yn;
    xext=xextn;
    yext=yextn;

    %compute new potentials
    u=P*x;  % sum internal PSPs
    
    temp=xext+conmat*x(indpyEu); %Potentials contribution from external PSPs, consisting of PSPs due to timevarying input and PSPs from pyramidal cells of other populations
    %add external potential contribution to populations
    u(indupy)=u(indupy)+temp;           %py
    u(induex)=u(induex)+alpha.*temp;    %ex
    u(induins)=u(induins)+beta.*temp;   %is
    u(induinf)=u(induinf)+gamma.*temp;  %if
    
    su=Sigm(u);     %new firing rate
    
    tel=tel+1;      %update output counterl
    %Save output every Nout-th step
    if tel==Nout
        tel=0;
        uout(:,indexplot)=u;
        xout(:,indexplot)=x;
        xextout(:,indexplot)=xext;
        yout(:,indexplot)=y;
        yextout(:,indexplot)=yext;
        indexplot=indexplot+1;
    end
end


end

%% Sigmoid function
% The sigmoid function is taken the same for all populations across all
% NMs. The parameters of this function are hard coded to speed-up
% computation. It is very ugly but sigmoid function can be changed by
% outcommenting the desired form...

% % sigmoid shifted s.t. S(0)=0 baselevel v0=4.5 (used in manuscript)
function outp=Sigm(v)
    outp=5./(1+exp(0.56*(4.5-v)))-0.372339725830140;
end

% % % sigmoid shifted s.t. S(0)=0, v0=6
% function outp=Sigm(v)
%     outp=5./(1+exp(0.56*(6-v)))-0.1678461164074125;
% end

% % % sigmoid shifted s.t. S(0)=0 baselevel v0=5
% function outp=Sigm(v)
%     outp=5./(1+exp(0.56*(5-v)))-0.286620879494344;
% end

% % % sigmoid shifted s.t. S(0)=0 baselevel v0=4
% function outp=Sigm(v)
%     outp=5./(1+exp(0.56*(4-v)))-0.481077708553464;
% end

% % % sigmoid shifted s.t. S(0)=0 baselevel v0=0
% function outp=Sigm(v)
%     outp=5./(1+exp(0.56*(v)))-2.5;
% end

% % % Unshited sigmoid
% function outp=Sigm(v)
%     outp=5./(1+exp(0.56*(6-v)));
% end







