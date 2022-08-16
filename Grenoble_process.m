function [W, dims, parties] = Grenoble_process(V1,V2,V3,V4,V5,V6,trace_F,input_state)
%Grenoble_process Generates the tripartite Grenoble QC-QC process 
%   [W, dims, parties] = Grenoble_process()
%   [W, dims, parties] = Grenoble_process(trace_F,input_state)
%   [W, dims, parties] = Grenoble_process(V1,V2,V3,V4,V5,V6,trace_F,input_state) 
%
%   The general form is given by Eq. (D7) of arXiv:2101.08796, and is specified by 6 isometries.
%   The canonical version (default), has V1=V2=V3=V_COPY and V4=V5=V6=CNOT
%   The systems are ordered as Pc,Pt,AI,AO,BI,BO,CI,CO,Ft,Falpha,Fc
%   By default we leave P and F open, but options allow F to be traced, and a joint control-target 
%   input state provided
%
% Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 16 August 2022
 
    if nargin == 0 || nargin == 6
        trace_F = false;
    end

    if nargin == 1 || nargin == 2
        trace_F = V1;
    end
    if nargin == 2
        input_state = V2;
    end

    COPY = Tensor([1;0],[1;0],[1,0]) + Tensor([0;1],[0;1],[0,1]);
    CNOT = Tensor([1;0],[1;0],[1,0],[1,0]) + Tensor([1;0],[0;1],[1,0],[0,1]) + Tensor([0;1],[0;1],[0,1],[1,0]) + Tensor([0;1],[1;0],[0,1],[0,1]);

    % default isometries
    if nargin < 6
        V1 = COPY;
        V2 = COPY;
        V3 = COPY;
        V4 = CNOT;
        V5 = CNOT;
        V6 = CNOT;
    end

    id2 = eye(2); 
    id3 = eye(3);

	d_t = 2; % target dimension
	d_c = 3; % control dimension
    
	% Pc Pt A1I A1O A2I A2O A3I A3O Ft Falpha Fc
	dims = [d_c d_t d_t*ones(1,6) d_t d_t d_c];

    Pc = 1;
    Pt = 2;
    P = [Pc,Pt];
    AI = 3;
    AO = 4;
    BI = 5;
    BO = 6;
    CI = 7;
    CO = 8;
    Ft = 9;
    Falpha = 10;
    Fc = 11;
    F = [Ft,Falpha,Fc]; 
    
    parties = {{[],P}, {AI,AO}, {BI,BO}, {CI,CO}, {F,[]}};

	V1_CJ = pure_CJ(V1);
	V2_CJ = pure_CJ(V2);
	V3_CJ = pure_CJ(V3);
	V4_CJ = pure_CJ(V4);
	V5_CJ = pure_CJ(V5);
	V6_CJ = pure_CJ(V6);
    
	id_CJ = pure_CJ(eye(d_t));

	ketA = id3(:,1);
	ketB = id3(:,2);
	ketC = id3(:,3);
	
	ket0 = id2(:,1);
	ket1 = id2(:,2);

    % Defined following Eq. (D7)
    % These vectors help perform the link products in (D7) following Eq. (6)
    link1 = Tensor(eye(d_t^2),id_CJ');
    link2 = (Tensor(eye(d_t),ket0',eye(d_t^2),ket0')+Tensor(eye(d_t),ket1',eye(d_t^2),ket1'));
    
	wABC = PermuteSystems(Tensor(...
             ketA,id_CJ,link1*Tensor(V1_CJ,ket0),...
             link2*Tensor(V6_CJ,ket0),...
             id_CJ,ketC...
           ),[Pc, Pt, AI, AO, BI, BO, CI, Falpha, CO, Ft, Fc],dims,0,1); % Put spaces back in right order
	wACB = PermuteSystems(Tensor(...
             ketA,id_CJ,link1*Tensor(V1_CJ,ket1),...
             link2*Tensor(V5_CJ,ket1),...
             id_CJ,ketB...
           ),[Pc, Pt, AI, AO, CI, CO, BI, Falpha, BO, Ft, Fc],dims,0,1);
	wBCA = PermuteSystems(Tensor(...
             ketB,id_CJ,link1*Tensor(V2_CJ,ket0),...
             link2*Tensor(V4_CJ,ket0),...
             id_CJ,ketA...
           ),[Pc, Pt, BI, BO, CI, CO, AI, Falpha, AO, Ft, Fc],dims,0,1);
	wBAC = PermuteSystems(Tensor(...
             ketB,id_CJ,link1*Tensor(V2_CJ,ket1),...
             link2*Tensor(V6_CJ,ket1),...
             id_CJ,ketC...
           ),[Pc, Pt, BI, BO, AI, AO, CI, Falpha, CO, Ft, Fc],dims,0,1);
	wCAB = PermuteSystems(Tensor(...
             ketC,id_CJ,link1*Tensor(V3_CJ,ket0),...
             link2*Tensor(V5_CJ,ket0),...
             id_CJ,ketB...
           ),[Pc, Pt, CI, CO, AI, AO, BI, Falpha, BO, Ft, Fc],dims,0,1);
	wCBA = PermuteSystems(Tensor(...
             ketC,id_CJ,link1*Tensor(V3_CJ,ket1),...
            link2*Tensor(V4_CJ,ket1),...
             id_CJ,ketA...
           ),[Pc, Pt, CI, CO, BI, BO, AI, Falpha, AO, Ft, Fc],dims,0,1);

	w = wABC + wACB + wBCA + wBAC + wCAB + wCBA;
	W = w*w';	

    if trace_F == true
        % trace all of F
        [W, dims, parties] = trace_superop_output(W,dims,parties,[1,2,3]);
    end

    if nargin == 2 || nargin == 8
        [W, dims, parties] = insert_superop_input(W,dims,parties,input_state,[1,2]);
    end
    
end
