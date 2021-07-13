%PSO >> function for the PSO ALGORITHM
%
% USAGES:   1.) [fxmin, xmin, Swarm, history] = PSO(psoOptions);
%           2.) [fxmin, xmin, Swarm, history] = PSO;
%           3.) fxmin = PSO(psoOptions);
%           3.) PSO 
%           etc.
%
% Arguments     : psoOptions--> A Matlab stucture containing all PSO related options. (see also: get_psoOptions)
% Return Values : [fxmin, xmin, Swarm, history]
%                     |     |       |       |_The history of the algorithm. Depicts how the function value of GBest changes over the run.
%                     |     |       |_The final Swarm. (A matrix containing co-ordinates of all particles)
%                     |     |_The co-ordinates of Best (ever) particle found during the PSO's run.
%                     |__The objective value of Best (^xmin) particle.
%
%  History        :   Author      :   JAG (Jagatpreet Singh)
%                     Created on  :   05022003 (Friday. 2nd May, 2003)
%                     Comments    :   The basic PSO algorithm.
%                     Modified on :   0710003 (Thursday. 10th July, 2003)
%                     Comments    :   It uses psoOptions structure now. More organized.
%  
%                     see also: get_psoOptions
function [fxmin, xmin, Swarm,accept_fes,success, history] = CMSCLPSO_MU(psoOptions,T,initSWARM,initStep)
%     need init_exemplers ,updating_exmplersCL ,updating_exmplersCLformerged,init_DLPE,update_DLPE
			%Initializations
			if nargin == 0
				psoOptions = get_psoOptions;
			end

			close all;
			%For Displaying 
			if psoOptions.Flags.ShowViz
				global vizAxes; %Use the specified axes if using GUI or create a new global if called from command window
				vizAxes = plot(0,0, '.');
				axis([-1000 1000 -1000 1000 -1000 1000]);   %Initially set to a cube of this size
				axis square;
				grid off;
				set(vizAxes,'EraseMode','xor','MarkerSize',15); %Set it to show particles.
				pause(1);
			end
			%End Display initialization
			%% 免疫克隆初始化种群
			Swarm = init_ica(psoOptions);
			%%普通初始化种群
			% Swarm = rand(psoOptions.Vars.SwarmSize, psoOptions.Vars.Dim)*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb; 
			me = psoOptions.Vars.max_gen;
			[num,Dim]=size(Swarm);
			SwarmSize = num;
		
			% Inertia Weight
			Weight = 0.9-(1:me)*0.7/me; 

			c1 = 1.49445.*ones(me,1);
			c2 = 1.49445.*ones(me,1);
			DLNum = 1;

			[M1,M2,shifto,lambda10,lambda100]= functionparameter(num, Dim);
			success = 0; % Success Flag
			f2eval = psoOptions.Obj.f2eval;
			accuracy = psoOptions.Vars.Threshold;
			FEs = 0;

			% init global
			fSwarm = feval(f2eval,Swarm,M1,M2,shifto,lambda10,lambda100);
			% fSwarm = feval(f2eval,Swarm);
			FEs = FEs + num;

			pbest_pos = Swarm;
			pbest_val= fSwarm;
			[gbest_val,gbest_id] = min(pbest_val);
			gbest_pos = pbest_pos(gbest_id,:);
			%初始化合并群信息
			pos_m = [];
			vel_m = [];
			pbest_val_m = [];
			pbest_pos_m =[];
			obj_func_slope_m  = [];
			gbest_pos_m = [];
			fri_best_pos_m = [];
			obj_func_slope_m = [];
			gbest_val_m_r= max( fSwarm); %初始化一个最大解做为合并群最优
			history = [0, gbest_val];
			%% 划分子群
			subSwarm = divideSwarm(Swarm,psoOptions,M1,M2,shifto,lambda10,lambda100);
			% subSwarm = divideSwarm(Swarm,psoOptions);
			subSwarmSize = size(subSwarm,2);
			%% Range for initial swarm's elements
			% Initialization VRmin and VRmax
			VRmin = psoOptions.Obj.lb;
			VRmax = psoOptions.Obj.ub;
			mv=0.2*(VRmax-VRmin);
			v_min=repmat(-mv,num,Dim);
			v_max=-v_min;
			beta=1;
			interval = psoOptions.Obj.ub - psoOptions.Obj.lb;
			vmax = 0.5*interval;
			vmin = -vmax;



			for i=1:subSwarmSize
			 
				 pos_g{i} = subSwarm{i};
				 num_g(i) = size(pos_g{i},1);
				 % vel_g{i} = v_min(1:num_g(i),:) + (v_max(1:num_g(i),:) - v_min(1:num_g(i),:)).*rand(num_g(i),Dim);
				 vel_g{i}= v_min(1:num_g(i),:)+2.*v_max(1:num_g(i),:).*rand(num_g(i),Dim);    
				 % 初始化每个子群的最优值
				 pbest_pos_g{i} = pos_g{i};
				 pbest_val_g{i}= feval(f2eval,pos_g{i},M1,M2,shifto,lambda10,lambda100);
				 % pbest_val_g{i}= feval(f2eval,pos_g{i});
				 FEs = FEs + num_g(i);
				 
				 
				 [gbest_val_g(i),gbest_id(i)] = min(pbest_val_g{i});
				 gbest_pos_g{i} = pbest_pos_g{i}(gbest_id(i),:);
				 gbest_val_g_r(i) = gbest_val_g(i);	   
				 % 检查每个粒子的维度是否更新
				 obj_func_slope_g{i}= zeros(num_g(i),1); 
				 num1=max(num_g);
				 t=0:1/(num1 - 1):1;t=5.*t;
				 Pc1=0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(num1))-exp(t(1)));	
				 % 初始化每个粒子的最好邻居
				 fri_best_pos_g{i} = init_exemplers(pbest_pos_g{i},pbest_val_g{i}',Pc1,num_g(i),Dim);
				 % 检查每个种群陷入局部最优的次数
				 stop_g(i) = 0;
				 % 种群开始合并的标志
				 merge_flag(i) = 0; 
				 
				 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		 
				 noupdatetimes{i}=zeros(num_g(i),1);%记录每一个子群中的粒子的running更新状态
				 currentrunstate{i} =zeros(num_g(i),1);%记录每一个子群中的粒子的running更新状态
				 %0 表示 CL，1表示 DL
				 PreviousPE{i} = pbest_pos_g{i} ;%记录每一个粒子的上一次运行DL后的最后一个最新的PE
				 needupdatePEagain{i} = zeros(num_g(i),1);%记录每一个粒子需要更新PE，由于CL运行时导致的Pbest的改变
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				 DLparticles{i} = [];
				 DLVec{i} = [];
				 
				 DLPBestval{i} = [];
				 DLPBest{i} =[];
				 DLPe{i} = [];
				 
				 DLFlag{i} = [];
				 DLparticleID{i} = [];
				 DLAge{i} = [];
				 DLnum(i) = 0;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
			end

			
%%%%%%%%%%%%%%%%%%%%改动%%%%%%%%%%%%%%%保留global子群的最优解初始化用，不是全局群的最优解%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
        gbest_val_m = inf;
		stop_m = 0;
		
	    MaxDLNumforeachsubSwarm = min(5, min(num_g(1:end)));
	    MaxCLruntimes = 10;
		MaxDLruntimes = 10;
%%%%%%%%%%%%%%%%%%%%改动%%%%%%%%%%%%%%%保留global子群的最优解初始化用，不是全局群的最优解%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
			
			
		fitcount = 0;
		k=0;
		
		Initmutestd = (psoOptions.Obj.ub-psoOptions.Obj.lb);%初始化方差

		%初始化变异参数
		T_d = repmat(0.1,1,psoOptions.Vars.Dim); %形成一个行向量[0.1 0.1 0.1 0.1]表示为V的阈值
		F_d = repmat(0,1,psoOptions.Vars.Dim);   %%
		Initmutestd = (psoOptions.Obj.ub-psoOptions.Obj.lb);%初始化方差
		Maxm = 15;
		M = 5;
		K_1 = 5;
		K_2 = 10;
		%%%scale
		% P = num/M;
		%初始化标准方差
		InitmutestdQ = repmat(Initmutestd,[1,M]);
		stastagancym=0;
		dec = 0;
		
		
		if psoOptions.Disp.Interval & (rem(k, psoOptions.Disp.Interval) == 0)
			fprintf('Iterations\t\tfGBest\t\t\tfevals\n');
		end
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  THE  PSO  LOOP                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while fitcount < psoOptions.Vars.Iterations
        fitcount = fitcount + 1;  
        pbest_pos_g_all=[];
		pbest_val_g_all=[];
        pbest_pos = [];
		
       k = k + 1;	 
		if k > me
			 k = k -1;
		end
		
	 
%%%%%%%%%%%%%%%%%%%%%%%%这个k，感觉合并后的群和单独子群用一个k，而且每次改变个数的时候就会k=1，会影响CL子群的k    
	
		for i=1:subSwarmSize

				   if merge_flag(i) == 0
						% 更新粒子邻居，只更新CL状态的粒子 
						 [fri_best_pos_g{i},obj_func_slope_g{i}] = updating_exmplersCL(fri_best_pos_g{i},pbest_pos_g{i},pbest_val_g{i}',Pc1,num_g(i),Dim,obj_func_slope_g{i},currentrunstate{i});  
                   
				   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%下面才是合并global的进化步骤%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				   else
				        pbest_pos_whole = [];pbest_val_whole=[];
                         pbest_pos_unmerge = [];pbest_val_unmerge=[]; 
						 tempii = 1:subSwarmSize;
						 tempii(i)=[];

                         for kk = 1: size(tempii,2)
							 pbest_pos_unmerge = [pbest_pos_unmerge;pbest_pos_g{tempii(kk)}];  
							 pbest_val_unmerge = [pbest_val_unmerge;pbest_val_g{tempii(kk)}];
                         end                
                         pbest_pos_whole = [pbest_pos_whole;pbest_pos_g{i};pbest_pos_unmerge];  
						 pbest_val_whole = [pbest_val_whole;pbest_val_g{i};pbest_val_unmerge];
						 %%%%当前的i的信息必须放在前面
				         [fri_best_pos_g{i},obj_func_slope_g{i}] = updating_exmplersCLformerged(fri_best_pos_g{i},pbest_pos_whole,pbest_val_whole',SwarmSize,Pc1,num_g(i),Dim,obj_func_slope_g{i},currentrunstate{i}); 
				  end
				         
						 %只对CL粒子进化
						 CLparaticlesID = [];
						 CLparaticlesID = find(currentrunstate{i}(:) == 0);
						 numgCL = length(CLparaticlesID);

					if merge_flag(i) == 0		 
						 
						 delta= 1.49445 .*(rand(numgCL,Dim).*(fri_best_pos_g{i}(CLparaticlesID,:) - pos_g{i}(CLparaticlesID,:)));  %邻域型

				    else
				         gbest_pos_temp_g{i} = [];
				 		 gbest_pos_temp_g{i} = repmat ( gbest_pos_g{i}, numgCL ,1 );		
						 delta= (c1(k).*rand(numgCL,Dim).*(fri_best_pos_g{i}(CLparaticlesID,:)- pos_g{i}(CLparaticlesID,:)))+(c2(k).*rand(numgCL,Dim).*(gbest_pos_temp_g{i} - pos_g{i}(CLparaticlesID,:)));  %带全局型         
                    end
						
						vel_g{i}(CLparaticlesID,:)= Weight(k) .* vel_g{i}(CLparaticlesID,:) + delta;
						  %检查边界
						 vel_g{i}(CLparaticlesID,:)=(vel_g{i}(CLparaticlesID,:) < v_min(1:numgCL,:)).*v_min(1:numgCL,:)+(vel_g{i}(CLparaticlesID,:) >= v_min(1:numgCL,:)).*vel_g{i}(CLparaticlesID,:);
						 vel_g{i}(CLparaticlesID,:)= (vel_g{i}(CLparaticlesID,:) > v_max(1:numgCL,:)).*v_max(1:numgCL,:)+(vel_g{i}(CLparaticlesID,:) <= v_max(1:numgCL,:)).*vel_g{i}(CLparaticlesID,:);
						 %更新位置
						 pos_g{i}(CLparaticlesID,:) = pos_g{i}(CLparaticlesID,:) + vel_g{i}(CLparaticlesID,:);				           								
						 
						 %更新局部与全局最优,fSwarm是列向量
						 fSwarm_g{i}(CLparaticlesID,:) = feval(f2eval,pos_g{i}(CLparaticlesID,:),M1,M2,shifto,lambda10,lambda100);
						 % fSwarm_g{i} = feval(f2eval,pos_g{i});
						 FEs = FEs + numgCL;
						 
				

				         changeRowsf1 = [];
						 changeRowsf1 = find((fSwarm_g{i}(CLparaticlesID,:)-pbest_val_g{i}(CLparaticlesID,:))<0);
						 curCLparticleIDwithupdatePbest = CLparaticlesID(changeRowsf1);
						 
						 pbest_val_g{i}(curCLparticleIDwithupdatePbest, :) = fSwarm_g{i}(curCLparticleIDwithupdatePbest, :);  %更新局部最优
						 pbest_pos_g{i}(curCLparticleIDwithupdatePbest, :) = pos_g{i}(curCLparticleIDwithupdatePbest, :);    %更新局部位置            
						 needupdatePEagain{i}(curCLparticleIDwithupdatePbest,:) = 1; %%%%表示需要更新PE如果运行DL
						
						%检查没有更新的粒子，标记加1
						tem1 = find((fSwarm_g{i}(CLparaticlesID,:)-pbest_val_g{i}(CLparaticlesID,:))>=0);
						curCLparticleIDwithnoupdatePbest = CLparaticlesID(tem1);
						
						obj_func_slope_g{i}(curCLparticleIDwithnoupdatePbest)=obj_func_slope_g{i}(curCLparticleIDwithnoupdatePbest)+1;	
						%[gbest_val_g(i),id]= min(pbest_val_g{i});                                    %寻找本次跟新后的全局最优
						
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%改动%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						
						 [tempgbest_val_gi,tempiid] = min(pbest_val_g{i});
						%pos_g{i}(tempiid,:);
						 %
						 if tempgbest_val_gi< gbest_val_g(i)
								gbest_id(i) = tempiid;
								gbest_pos_g{i}=pbest_pos_g{i}(gbest_id(i),:);;
								gbest_val_g(i)=tempgbest_val_gi;
								stop_g(i) = 0; 
						 else
						
							if fitcount >= psoOptions.Vars.Iterations/5 
								stop_g(i) = stop_g(i) + 1; 
							end
								
						 end
						 
						 if  gbest_val_g(i) <  gbest_val
								gbest_pos = gbest_pos_g{i};
								gbest_val = gbest_val_g(i);
						 end    			

						 
						 
						 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%判断每一个子群中每一粒子的更新状态%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%					 
						 
						 noupdatetimes{i}(curCLparticleIDwithupdatePbest) = 0;%
						 noupdatetimes{i}(curCLparticleIDwithnoupdatePbest) = noupdatetimes{i}(curCLparticleIDwithnoupdatePbest)+1;%
						 
						
						    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%选取方法，可以考虑按照适应值排序或者随机%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
						    tempID = find(noupdatetimes{i}(CLparaticlesID,1) > MaxCLruntimes);
							stagnantCLparaticlesID = CLparaticlesID(tempID);
							tempID = stagnantCLparaticlesID;
							
							if ~isempty(tempID)
							   leftDLnumi = MaxDLNumforeachsubSwarm-DLnum(i);
							   	%最好不要选择到gbest_pos！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
                               if leftDLnumi~=0							
 							      if length(tempID) <= leftDLnumi
								   
									DLparticleID{i}= [DLparticleID{i};tempID];	
									%列向量
									DLparticles{i} = [DLparticles{i};pos_g{i}(tempID,:)];
									DLVec{i} = [DLVec{i};vel_g{i}(tempID,:)];
									DLPBestval{i} = [DLPBestval{i};pbest_val_g{i}(tempID,:)];
									DLPBest{i} = [DLPBest{i};pbest_pos_g{i}(tempID,:)];
									DLFlag{i} = [DLFlag{i};zeros(length(tempID),1)];
									
									DLAge{i} =[DLAge{i};zeros(length(tempID),1)];
								
									for seqindex =1:length(tempID)
											if needupdatePEagain{i}(tempID(seqindex),:) == 1
											   [DLPe{i}(DLnum(i)+seqindex,:),gbest_pos,gbest_val,gbest_pos_g{i},gbest_val_g(i),FEs1] = init_DLPE(gbest_val,gbest_pos,gbest_val_g(i),gbest_pos_g{i},DLPBestval{i}(DLnum(i)+seqindex),DLPBest{i}(DLnum(i)+seqindex,:),DLparticles{i}(DLnum(i)+seqindex,:),psoOptions,M1,M2,shifto,lambda10,lambda100);
										        FEs = FEs+FEs1;
												PreviousPE{i}(tempID(seqindex),:) = DLPe{i}(DLnum(i)+seqindex,:);
												needupdatePEagain{i}(tempID(seqindex),:) = 0;
											else
											    DLPe{i}(DLnum(i)+seqindex,:) = PreviousPE{i}(tempID(seqindex),:);
                                            end											
									end
									
									DLnum(i) = DLnum(i)+length(tempID);
									
									currentrunstate{i}(tempID) = 1;
								  
								  else
								%%%%%%%%%%%%%%%%选取方法，可以考虑按照适应值排序或者随机或者按照noupdatetimes大小排序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		  
									tempID1  = randperm(length(tempID));
									tempID = tempID(tempID1(1:leftDLnumi));
									
									%最好不要选择到gbest_pos！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
									
									DLparticleID{i}= [DLparticleID{i};tempID];	
									DLparticles{i} = [DLparticles{i};pos_g{i}(tempID,:)];
									DLVec{i} = [DLVec{i};vel_g{i}(tempID,:)];
									DLPBestval{i} = [DLPBestval{i};pbest_val_g{i}(tempID)];
									DLPBest{i} = [DLPBest{i};pbest_pos_g{i}(tempID,:)];
									DLFlag{i} = [DLFlag{i};zeros(length(tempID),1)];
									
									DLAge{i} =[DLAge{i};zeros(length(tempID),1)];
									
									for seqindex =1:length(tempID)
											if needupdatePEagain{i}(tempID(seqindex),:) == 1
											   [DLPe{i}(DLnum(i)+seqindex,:),gbest_pos,gbest_val,gbest_pos_g{i},gbest_val_g(i),FEs1] = init_DLPE(gbest_val,gbest_pos,gbest_val_g(i),gbest_pos_g{i},DLPBestval{i}(DLnum(i)+seqindex),DLPBest{i}(DLnum(i)+seqindex,:),DLparticles{i}(DLnum(i)+seqindex,:),psoOptions,M1,M2,shifto,lambda10,lambda100);
										        FEs = FEs+FEs1;
												PreviousPE{i}(tempID(seqindex),:) = DLPe{i}(DLnum(i)+seqindex,:);
												needupdatePEagain{i}(tempID(seqindex),:) = 0;
											else
											    DLPe{i}(DLnum(i)+seqindex,:) = PreviousPE{i}(tempID(seqindex),:);
                                            end											
									end
									
									
									DLnum(i) = DLnum(i)+length(tempID);
									
									currentrunstate{i}(tempID) = 1;
							     end %end of if length(tempID) <= leftDLnumi )
							  end	%end of if leftDLnumi~=0		
						  end	%end of ~isempy(tempID)
							   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end 选取方法，可以考虑按照适应值排序或者随机%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%							 
						if stop_g(i) >= psoOptions.Vars.stop_num 
						merge_flag(i) = 1;
					  end 
						 
			 
						 %定期的检查子群，看是否陷入局部最优
			
					
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DL 执行阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
        if  DLnum(i)~=0
		
			c3=0.5+2*(fitcount/psoOptions.Vars.Iterations);
            w1=0.9 - 0.5*(fitcount/psoOptions.Vars.Iterations); 	
		 
		 

			%%%  update vel and pos bu DL
		    for ii=1:DLnum(i)
				
								  
					  for dem = 1:Dim					                          
						DLVec{i}(ii,dem)=w1* DLVec{i}(ii,dem)+1.49445*rand*(DLPe{i}(ii,dem)-DLparticles{i}(ii,dem))+c3*rand*(gbest_pos(dem)-DLparticles{i}(ii,dem));
						DLVec{i}(ii,dem)=min(vmax,max(-vmax, DLVec{i}(ii,dem)));
						DLparticles{i}(ii,dem)=DLparticles{i}(ii,dem)+ DLVec{i}(ii,dem);
					  end 
			          
					  %%%更新我对应的真正的粒子位置和粒子vec！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
					  
						 vel_g{i}(DLparticleID{i}(ii),:) = DLVec{i}(ii,:);
						 pos_g{i}(DLparticleID{i}(ii),:)= DLparticles{i}(ii,:);
					  
					  %%%
				
						fitnessValue1 = feval(f2eval,pos_g{i}(DLparticleID{i}(ii),:),M1,M2,shifto,lambda10,lambda100);
						FEs=FEs+1;
					
					   if fitnessValue1<pbest_val_g{i}(DLparticleID{i}(ii))
					       % DLPBestval{i}(DLparticleID{i}(ii)) = fitnessValue1;
						    % DLPBest{i}(DLparticleID{i}(ii),:) = DLparticles{i}(ii,:);
						   pbest_val_g{i}(DLparticleID{i}(ii)) = fitnessValue1;
						   pbest_pos_g{i}(DLparticleID{i}(ii),:)=pos_g{i}(DLparticleID{i}(ii),:);
						   DLFlag{i}(ii) = 1;
						   DLAge{i}(ii) = 0;
						   
						   noupdatetimes{i}(DLparticleID{i}(ii)) = 0;
						   %已经改变了当前的粒子的pbest因此置零
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%改动%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%					  
					   else
					       DLAge{i}(ii) =DLAge{i}(ii)+1;
						   DLFlag{i}(ii) = 0;
					   end
					   
				       if fitnessValue1 < gbest_val_g(i)
							gbest_pos_g{i}=pos_g{i}(DLparticleID{i}(ii),:);
							gbest_val_g(i)=fitnessValue1;
							stop_g(i) = stop_g(i)-1; 
					 
					   end
					
					   if  gbest_val_g(i) <  gbest_val
							gbest_pos = gbest_pos_g{i};
							gbest_val = gbest_val_g(i);
					   end
				      
            end %end of for each DL particle in current subswarm i 
			
			stagnantDLid = find (DLAge{i}>MaxDLruntimes);
			%%%保留要退出的这些DL粒子的PE信息
			PreviousPE{i}(DLparticleID{i}(stagnantDLid,:),:) = DLPe{i}(stagnantDLid,:);
			%当前粒子回去执行CL进化去吧！！！！！！！！！！！！！！！！
			currentrunstate{i}(DLparticleID{i}(stagnantDLid,:)) = 0;
						            
			%随后，我这边的信息全部清空			
            DLparticleID{i}(stagnantDLid,:)= [];	
			DLparticles{i}(stagnantDLid,:) = [];
			DLVec{i}(stagnantDLid,:) = [];
			DLPBestval{i}(stagnantDLid) = [];
			DLPBest{i}(stagnantDLid,:) = [];
  		    DLFlag{i}(stagnantDLid) = [];
					
			DLAge{i}(stagnantDLid) =[];
			DLPe{i}(stagnantDLid,:)= [];
			DLnum(i) = DLnum(i)-length(stagnantDLid);
				
		end %end of if~isempty(DLparticlesSubSnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end DL 执行阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%					
		
		
	end	%end of all subswarms 	

	
        
	
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%改动%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%放到这里等到the global swarm也进化完了，再更新，这样可以照顾永远不合并的情况，而不是放在子群进化完成后%%%%%%%%%%%%%%%
    
   for curi = 1:subSwarmSize	
	%%%	if merge_flag(curi) == 0    %%%
			for jj = 1:DLnum(curi)
			% update learning exemplars for DL-subpopulation ,将1写死了，这样不好，考虑知道那个是DL子群就好了，label啥的
				%[pe,gbest_pos_g{curi},gbest_val_g(curi),FEs1] = update_DL(gbest_val,gbest_pos,pbest_val_g{DLNum},pbest_pos_g{DLNum},pos_g{DLNum},psoOptions,flag1,pe);
                [DLPe{curi}(jj,:),gbest_pos,gbest_val,gbest_pos_g{curi},gbest_val_g(curi),FEs1] = update_DLPE(gbest_val,gbest_pos,gbest_val_g(curi),gbest_pos_g{curi},DLPBestval{curi}(jj),DLPBest{curi}(jj,:),DLparticles{curi}(jj,:),psoOptions,DLFlag{curi}(jj),DLPe{curi}(jj,:),M1,M2,shifto,lambda10,lambda100);
				DLFlag{curi}(jj) = 0;
				FEs = FEs+FEs1;
			end	
	%%%	end		  
    end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%放到这里等到the global swarm也进化完了，再更新，这样可以照顾永远不合并的情况，而不是放在子群进化完成后%%%%%%%%%%%%%%%	
	   
    if mod(fitcount,100) == 0
      merger_num = size(find(merge_flag >= 1),2);
      fprintf('已有 %d 个种群开始合并,还有 %d 个子群未开始合并\n',merger_num ,subSwarmSize - merger_num);
    end 
    
   
     if success == 0
            if gbest_val <= accuracy
                accept_fes = FEs;
                success = 1;
            else
                accept_fes = 0;
            end
     end
        
    if psoOptions.Save.Interval & (rem(fitcount, psoOptions.Save.Interval) == 0)
        history((size(history,1)+1), :) = [fitcount, gbest_val];
    end
	
    if psoOptions.Disp.Interval & (rem(k, psoOptions.Disp.Interval) == 0)
        fprintf('%4d\t\t\t%.5g\t\t\t%5d\t\t\t%5d\n', fitcount, gbest_val, FEs,T);
    end
	
end	
	
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  END  OF PSO  LOOP                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fxmin = gbest_val;
xmin = gbest_pos;
Swarm = pbest_pos;
history = history(:,2);

function [stdq]=StdQstandard(q,W);
while(find(abs(q)>W/2)>0)
    q(find(abs(q)>W/2))=abs(W/2-q(find(abs(q)>W/2)));
end
stdq=q;
