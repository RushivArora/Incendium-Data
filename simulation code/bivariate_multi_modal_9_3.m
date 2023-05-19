

t_back_force=zeros(4,1);
t_follow_maximum=0;
t_visit_follow_maximum=0;
numberofiteration=20;
t_gradient_ascent=zeros(numberofiteration,1);
t_step_ahead=zeros(numberofiteration,10);
t_partition=zeros(10,4);
t_visit_normal=0;
t_visit_normal_order=0;
t_visit_normal_order2=0;
notcount=0;
for a=50:50:200
digitsOld = digits(5000);
for indexloop=1:100
	%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
    %Muti normal
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
    numrand=randi(3)+1;
    %numrand=1;
	rrr = rand(1, numrand); % Start with numrand random numbers that don't sum to 1.
	p = rrr / sum(rrr);  % Normalize so the sum is 1.
    mu=zeros(numrand,2);
	mu(:,1) = randi([1 2*a+1],numrand,1);
    mu(:,2) = randi([1 2*a+1],numrand,1);
	sig=randi(30,numrand,2);
	%sig=5*ones(numrand,2);
	if(numrand==1)        
		sigma=cat(3,sig(1,:));
	elseif(numrand==2)        
		sigma=cat(3,sig(1,:),sig(2,:));
	elseif(numrand==3)       
		sigma=cat(3,sig(1,:),sig(2,:),sig(3,:));
	elseif(numrand==4)       
		sigma=cat(3,sig(1,:),sig(2,:),sig(3,:),sig(4,:));
	elseif(numrand==5)
		sigma=cat(3,sig(1,:),sig(2,:),sig(3,:),sig(4,:),sig(5,:));
	end
	gm = gmdistribution(mu,sigma);
	x1 = 1:1:2*a+1;
	x2 = 1:1:2*a+1;
	[X1,X2] = meshgrid(x1,x2);
	X = [X1(:) X2(:)];
	y = pdf(gm,X);
	y = reshape(y,length(x2),length(x1));
    %figure;imshow(y*(10^10));
    %%
    %imagesc(y)
    %%colorbar
	surf(x1,x2,y);
	R = random(gm,1);
	R=ceil(R);
	R(1,1)=R(1,1);
	R(1,2)=R(1,2);
    
	if(R(1,1)<1)
		R(1,1)=1;
	elseif(R(1,1)>2*a+1)
		R(1,1)=2*a+1;
	end
	if(R(1,2)<1)
		R(1,2)=1;
	elseif(R(1,2)>2*a+1)
		R(1,2)=2*a+1;
    end
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
    %Uniform
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
%     v = rand(2*a+1,2*a+1);
%     vSum = sum(v);
%     v2v = v ./ vSum;
%     y = v2v;
%     v = v2v(:);
%     b=1:(2*a+1)*(2*a+1);
%     N = 1;
% 
%     [~,R] = histc(rand(1,N),cumsum([0;v(:)./sum(v)]));
%     [R(1,1),R(1,2)] = ind2sub(size(v2v),R);
%     
%     if(R(1,1)<1)
% 		R(1,1)=1;
% 	elseif(R(1,1)>2*a+1)
% 		R(1,1)=2*a+1;
% 	end
% 	if(R(1,2)<1)
% 		R(1,2)=1;
% 	elseif(R(1,2)>2*a+1)
% 		R(1,2)=2*a+1;
%     end
    %R = v(R)
    %%%%%%%%%%%%%%%%
    %y(R(1,2),R(1,1))
    hold on;
    % Plot cross at row 100, column 50
    plot(R(1,1),R(1,2), 'r+', 'MarkerSize', 30, 'LineWidth', 2);
    
	%%%%back-and-force
	t=0;
	loc_back_force=[1,1];
    y_temp=y;
	dir=0;
	while(~((loc_back_force(1)==R(1,1)) && (loc_back_force(2)==R(1,2))))
		t=t+1;
		if(dir==0)
			loc_back_force(1)=loc_back_force(1)+1;
			if(loc_back_force(1)==2*a+2)
				loc_back_force(1)=2*a+1;
				loc_back_force(2)=loc_back_force(2)+1;
				dir=1-dir;
				t=t+1;
			end
		else
			loc_back_force(1)=loc_back_force(1)-1;
			if(loc_back_force(1)==0)
				loc_back_force(1)=1;
				loc_back_force(2)=loc_back_force(2)+1;
				dir=1-dir;
				t=t+1;
			end		
        end
        y_temp(loc_back_force(1), loc_back_force(2))=0;
    end
    
    %
    plot(loc_back_force(1),loc_back_force(2), 'b+', 'MarkerSize', 30, 'LineWidth', 2);
    t_last_back_force=t;
	t_back_force(a/50)=t_back_force(a/50)+t
    %%%%%%%%%%%%%%%
    for i=10:-1:1
        j=1;
        k=1;
        notfound=0;
        y_temp=y;
        ttemp=0;
        while((notfound==0)&&(ttemp<t))
            for iii=1:(2*i)
                % extract region of interest
                BRegion = y_temp((j-1)*2*i+1:j*2*i+1, (k-1)*2*i+1:k*2*i+1);

                % find max value and get its index
                [value, kin] = max(BRegion(:));
                [iin, jin] = ind2sub(size(BRegion), kin);

                % move indexes to correct spot in matrix
                iin = iin + (j-1)*2*i;
                jin = jin + (k-1)*2*i;
                if(((jin==R(1,1)) && (iin==R(1,2)))||(iin==R(1,1)) && (jin==R(1,2)))
                    notfound=1;
                    break;
                end
                y_temp(iin, jin)=0;
            end
            ttemp=ttemp+2*i;
            yyy=zeros(floor(size(y,1)/(2*i)),floor(size(y,2)/(2*i)));
            maxoverall=0;
            
            for jl=1:size(y,1)/(2*i)
                for kl=1:size(y,2)/(2*i)
                    for i1=1:2*i
                        for i2=1:2*i
                            yyy(jl,kl)=yyy(jl,kl)+y_temp((jl-1)*2*i+i1,(kl-1)*2*i+i2);
                        end
                    end
                
                    if(yyy(jl,kl)>maxoverall)
                        maxoverall=yyy(jl,kl);
                        maxjl=jl;
                        maxkl=kl;
                    end
                end
            end
            if(maxjl>j)
                j=j+1;
            else
                if(maxjl<j)
                    j=j-1;
                end
            end
            if(maxkl>k)
                k=k+1;
            else
                if(maxkl<k)
                    k=k-1;
                end
            end
            
        end
        t_partition(i,a/50)=t_partition(i,a/50)+ttemp;
    end
    
    %end
	%%%%Just_folow-maximum
% 	y_temp=y;
% 	t=0;
% 	loc_follow_maximum_prev=[1,1];
% 	[max_num, max_idx]=max(y_temp(:));
% 	[maxx,maxy]=ind2sub(size(y_temp),max_idx);
% 	while(~((maxy==R(1,1)) && (maxx==R(1,2))))	
% 		y_temp(maxx, maxy)=0;
% 		y_temp=y_temp*(1-max_num);
% 		maxdis=sqrt((loc_follow_maximum_prev(1)-maxx)^2+(loc_follow_maximum_prev(2)-maxy)^2);
% 		t=t+maxdis;
% 		loc_follow_maximum_prev(1)=maxx;
% 		loc_follow_maximum_prev(2)=maxy;
% 		[max_num, max_idx]=max(y_temp(:));
% 		[maxx,maxy]=ind2sub(size(y_temp),max_idx);
% 	end
% 	t_follow_maximum=t_follow_maximum+t
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%%%%Visit-folow-maximum
% 	y_temp=y;
% 	t=0;
% 	loc_follow_maximum_prev=[1,1];
% 	loc_follow_maximum_next=[1,1];
% 	[max_num, max_idx]=max(y_temp(:));
% 	[maxx,maxy]=ind2sub(size(y_temp),max_idx);	
% 	while(~((loc_follow_maximum_prev(2)==R(1,1)) && (loc_follow_maximum_prev(1)==R(1,2))))	
% 		if((maxx==loc_follow_maximum_prev(1))&&(maxy==loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 		elseif((maxx==loc_follow_maximum_prev(1))&&(maxy>loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 		elseif((maxx>loc_follow_maximum_prev(1))&&(maxy>loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 		elseif((maxx>loc_follow_maximum_prev(1))&&(maxy==loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 		elseif((maxx>loc_follow_maximum_prev(1))&&(maxy<loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 		elseif((maxx==loc_follow_maximum_prev(1))&&(maxy<loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 		elseif((maxx<loc_follow_maximum_prev(1))&&(maxy<loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 		elseif((maxx<loc_follow_maximum_prev(1))&&(maxy==loc_follow_maximum_prev(2)))
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 		else
% 			loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 			loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 		end
% 		t=t+1;
% 		if(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0)
% 			ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
% 			y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
% 			y_temp=y_temp*(1-ptemp);
% 		end	
% 		loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
% 		loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
% 		[max_num, max_idx]=max(y_temp(:));
% 		[maxx,maxy]=ind2sub(size(y_temp),max_idx);
% 	end
% 	t_visit_follow_maximum=t_visit_follow_maximum+t
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%Visit-gradient_ascent 
	% y_temp=y;
	% t=0;
	% loc_follow_maximum_prev=[1,1];
	% loc_follow_maximum_next=returnmax(y_temp,loc_follow_maximum_prev);
	% while(~((loc_follow_maximum_next(2)==R(1,1)) && (loc_follow_maximum_next(1)==R(1,2))))
		% loc_follow_maximum_prev=loc_follow_maximum_next;
		% t=t+1;
		% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));		
		% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
		% y_temp=y_temp*(1-ptemp);
		% loc_follow_maximum_next=returnmax(y_temp,loc_follow_maximum_prev);
		% if(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))==0)
			% window=3;
			% while(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))==0)
				% loc_follow_maximum_next=returnmax2(y_temp,loc_follow_maximum_prev,window);
				% window=window+1;
			% end
		% end
	% end
	% t_gradient_ascent=t_gradient_ascent+t
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%Visit-gradient_ascent
% 	for indexi=1:numberofiteration
% 		y_temp=y;
% 		t=0;
% 		loc_follow_maximum_prev=[1,1];
% 		loc_follow_maximum=returnmax3(y_temp,loc_follow_maximum_prev,indexi);
% 		while(~((loc_follow_maximum_prev(2)==R(1,1)) && (loc_follow_maximum_prev(1)==R(1,2))))
% 			if((loc_follow_maximum(1)==loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)==loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 			elseif((loc_follow_maximum(1)==loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)>loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 			elseif((loc_follow_maximum(1)>loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)>loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 			elseif((loc_follow_maximum(1)>loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)==loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 			elseif((loc_follow_maximum(1)>loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)<loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)+1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 			elseif((loc_follow_maximum(1)==loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)<loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 			elseif((loc_follow_maximum(1)<loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)<loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)-1;
% 			elseif((loc_follow_maximum(1)<loc_follow_maximum_prev(1))&&(loc_follow_maximum(2)==loc_follow_maximum_prev(2)))
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
% 			else
% 				loc_follow_maximum_next(1)=loc_follow_maximum_prev(1)-1;
% 				loc_follow_maximum_next(2)=loc_follow_maximum_prev(2)+1;
% 			end
% 			t=t+1;
% 			if(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0)
% 				ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
% 				y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
% 				y_temp=y_temp*(1-ptemp);
% 			end	
% 			loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
% 			loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);
% 			
% 			loc_follow_maximum=returnmax3(y_temp,loc_follow_maximum_prev,indexi);
% 			if(y_temp(loc_follow_maximum(1), loc_follow_maximum(2))==0)
% 				window=indexi+1;
% 				while(y_temp(loc_follow_maximum(1), loc_follow_maximum(2))==0 && window<202)
% 					loc_follow_maximum=returnmax3(y_temp,loc_follow_maximum_prev,window);
% 					window=window+1;
% 				end
% 			end
% 		end
% 		t_gradient_ascent(indexi)=t_gradient_ascent(indexi)+t
% 	end
	%%%%Visit-gradient_ascent
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
% 	for indexi=1:numberofiteration
%         for indexj=1:10
%             y_temp=y;
%             t=0;
%             loc_follow_maximum_prev=[1,1];
%             %[loc_follow_maximum_next, maxsum]=returnmax6(y_temp,loc_follow_maximum_prev,2*indexi,(t_last_back_force-t)/500);
%             loc_follow_maximum_next(1)=loc_follow_maximum_prev(1);
%             loc_follow_maximum_next(2)=loc_follow_maximum_prev(2);
%             if(y_temp(loc_follow_maximum_prev(1), loc_follow_maximum_prev(2))~=0)
%                 ptemp=y_temp(loc_follow_maximum_prev(1), loc_follow_maximum_prev(2));
%                 y_temp(loc_follow_maximum_prev(1), loc_follow_maximum_prev(2))=0;
%                 y_temp=y_temp*(1-ptemp);
%             end
%             while(~((loc_follow_maximum_prev(2)==R(1,1)) && (loc_follow_maximum_prev(1)==R(1,2))))
%                 %t=t+1;
%                 if(t>t_last_back_force)
%                     loc_follow_maximum_prev(2)=R(1,1);
%                     loc_follow_maximum_prev(1)=R(1,2);
%                     
%                 else
%                     [loc_follow_maximum_next, maxsum, poslist]=returnmax7_1(y_temp,loc_follow_maximum_prev,2*indexi,(t_last_back_force-t)*indexj/500);
%                     if(maxsum<0)
%                         window=2*indexi+1;
%                         while(maxsum<=0 && window<202)
%                             [loc_follow_maximum_next, maxsum, poslist]=returnmax7_1(y_temp,loc_follow_maximum_prev,window,(t_last_back_force-t)/500);
%                             window=window+1;
%                         end
%                         
%                     end
%                     if(maxsum>0)
%                         t=t+sqrt((loc_follow_maximum_prev(1)-poslist(size(poslist,1),1))^2+(loc_follow_maximum_prev(2)-poslist(size(poslist,1),2))^2);
%                         checkloc=0;
%                         for indexii=1:size(poslist,1)
%                             if((poslist(indexii,2)==R(1,1)) && (poslist(indexii,1)==R(1,2)))
%                                 checkloc=1;
%                                 loc_follow_maximum_prev(1)=poslist(indexii,1);
%                                 loc_follow_maximum_prev(2)=poslist(indexii,2);
%                             end
%                             %loc_follow_maximum_next
%                             if(y_temp(poslist(indexii,1), poslist(indexii,2))~=0)
%                                 ptemp=y_temp(poslist(indexii,1), poslist(indexii,2));
%                                 y_temp(poslist(indexii,1), poslist(indexii,2))=0;
%                                 y_temp=y_temp*(1-ptemp);
%                             end
%                         end
%                         if(checkloc==0)
%                             loc_follow_maximum_prev(1)=poslist(size(poslist,1),1);
%                             loc_follow_maximum_prev(2)=poslist(size(poslist,1),2);
%                         end
%                     end
%                 end
%                 
%                 
%                 
%                 
%             end
%             %figure;imshow(y_temp*(10^30))
%             t_step_ahead(indexi,indexj)=t_step_ahead(indexi,indexj)+t;
%         end
%     end
%     
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%Visit-maximum-normal
	% p_temp=p;
	% order_p=ones(numrand,1);
	% for iloop=1:numrand
		% [max_num, max_idx]=max(p_temp(:));
		% order_p(iloop)=max_idx;
		% p_temp(max_idx)=0;
	% end
	% p_temp=p;
	% order_p=ones(numrand,1);
	% for iloop=1:numrand
		% [max_num, max_idx]=max(p_temp(:));
		% order_p(iloop)=max_idx;
		% p_temp(max_idx)=0;
	% end
	% y_temp=y;
	% t=0;
	% tfinal=-1;
	% loc_follow_maximum_prev=[1,1];
	% loc_follow_maximum_next=[1,1];
	% for iloop1=1:numrand		
		% %diff=[0 0;1 0;1 1;0 1; -1 1; -1 0;-1 -1;0 -1;1 -1];
		% %diff=zeros((2*sig(order_p(iloop1),1)+1)*(2*sig(order_p(iloop1),2)+1),2);	
		
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=-1*sig(order_p(iloop1),2):sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		
	% end
	% %2\delta
	% for iloop1=1:numrand		
		% %diff=[2 0;2 1; 2 2;1 2;0 2; -1 2;-2 2;-2 1;-2 0; -2 1; -2 -2;-1 -2; 0 -2;1 -2;2 -2;2 -1];
		% for iloop3=-2*sig(order_p(iloop1),1):-1*sig(order_p(iloop1),1)-1
			% for iloop4=-2*sig(order_p(iloop1),2):2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=sig(order_p(iloop1),1)+1:2*sig(order_p(iloop1),1)
			% for iloop4=-2*sig(order_p(iloop1),2):2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=-2*sig(order_p(iloop1),2):-1*sig(order_p(iloop1),2)-1
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=sig(order_p(iloop1),2)+1:2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end	
		
	% end
	% %3\delta
	% for iloop1=1:numrand		
		% %diff=[3 0;3 1; 3 2; 3 3; 2 3; 1 3;0 3; -1 3; -2 3; -3 3; -3 2; -3 1; -3 0; -3 -1; -3 -2;-3 -3; -2 -3; -1 -3; 0 -3; 1 -3; 2 -3; 3 -3; 3 -2; 3 -1];
		% for iloop3=-3*sig(order_p(iloop1),1):-2*sig(order_p(iloop1),1)-1
			% for iloop4=-3*sig(order_p(iloop1),2):3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=2*sig(order_p(iloop1),1)+1:3*sig(order_p(iloop1),1)
			% for iloop4=-3*sig(order_p(iloop1),2):3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-2*sig(order_p(iloop1),1):2*sig(order_p(iloop1),1)
			% for iloop4=-3*sig(order_p(iloop1),2):-2*sig(order_p(iloop1),2)-1
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-2*sig(order_p(iloop1),1):2*sig(order_p(iloop1),1)
			% for iloop4=2*sig(order_p(iloop1),2)+1:3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
	% end
	% if(tfinal==-1)
		% notcount=notcount+1;
	% end
	% t_visit_normal=t_visit_normal+tfinal
	% %%%%Visit-maximum-normal-order from p/sigma	
	% p_temp=p;
	% for iloop=1:numrand
		% p_temp(iloop)=p_temp(iloop)/((2*sig(iloop,1)+1)*(2*sig(iloop,2)+1));
	% end
	% order_p=ones(numrand,1);
	% for iloop=1:numrand
		% [max_num, max_idx]=max(p_temp(:));
		% order_p(iloop)=max_idx;
		% p_temp(max_idx)=0;
	% end
	% y_temp=y;
	% t=0;
	% tfinal=-1;
	% loc_follow_maximum_prev=[1,1];
	% loc_follow_maximum_next=[1,1];
	% for iloop1=1:numrand		
		% %diff=[0 0;1 0;1 1;0 1; -1 1; -1 0;-1 -1;0 -1;1 -1];
		% %diff=zeros((2*sig(order_p(iloop1),1)+1)*(2*sig(order_p(iloop1),2)+1),2);	
		
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=-1*sig(order_p(iloop1),2):sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		
	% end
	% %2\delta
	% for iloop1=1:numrand		
		% %diff=[2 0;2 1; 2 2;1 2;0 2; -1 2;-2 2;-2 1;-2 0; -2 1; -2 -2;-1 -2; 0 -2;1 -2;2 -2;2 -1];
		% for iloop3=-2*sig(order_p(iloop1),1):-1*sig(order_p(iloop1),1)-1
			% for iloop4=-2*sig(order_p(iloop1),2):2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=sig(order_p(iloop1),1)+1:2*sig(order_p(iloop1),1)
			% for iloop4=-2*sig(order_p(iloop1),2):2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=-2*sig(order_p(iloop1),2):-1*sig(order_p(iloop1),2)-1
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-1*sig(order_p(iloop1),1):sig(order_p(iloop1),1)
			% for iloop4=sig(order_p(iloop1),2)+1:2*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end	
		
	% end
	% %3\delta
	% for iloop1=1:numrand		
		% %diff=[3 0;3 1; 3 2; 3 3; 2 3; 1 3;0 3; -1 3; -2 3; -3 3; -3 2; -3 1; -3 0; -3 -1; -3 -2;-3 -3; -2 -3; -1 -3; 0 -3; 1 -3; 2 -3; 3 -3; 3 -2; 3 -1];
		% for iloop3=-3*sig(order_p(iloop1),1):-2*sig(order_p(iloop1),1)-1
			% for iloop4=-3*sig(order_p(iloop1),2):3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=2*sig(order_p(iloop1),1)+1:3*sig(order_p(iloop1),1)
			% for iloop4=-3*sig(order_p(iloop1),2):3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-2*sig(order_p(iloop1),1):2*sig(order_p(iloop1),1)
			% for iloop4=-3*sig(order_p(iloop1),2):-2*sig(order_p(iloop1),2)-1
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
		% for iloop3=-2*sig(order_p(iloop1),1):2*sig(order_p(iloop1),1)
			% for iloop4=2*sig(order_p(iloop1),2)+1:3*sig(order_p(iloop1),2)
				% loc_follow_maximum_next(1)=mu(order_p(iloop1),1)+a+1+iloop3;
				% loc_follow_maximum_next(2)=mu(order_p(iloop1),2)+a+1+iloop4;
				% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
				% t=t+maxdis;
				% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
					% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
					% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
					% y_temp=y_temp*(1-ptemp);
				% end	
				% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
				% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
				% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
					% tfinal=t;
				% end
			% end
		% end
	% end
	% if(tfinal==-1)
		% notcount=notcount+1;
	% end
	% t_visit_normal_order=t_visit_normal_order+tfinal
	% %%%%Visit-maximum-normal-order from p/sigma	
	% p_temp=p;
	% tstore=zeros(numrand,3);
	% ppp_temp=zeros(numrand,3);	
	% for iloop=1:numrand
		% ppp_temp(iloop,1)=0.68*p_temp(iloop)/((2*sig(iloop,1)+1)*(2*sig(iloop,2)+1));
		% ppp_temp(iloop,2)=0.27*p_temp(iloop)/((4*sig(iloop,1)+1)*(4*sig(iloop,2)+1)-(2*sig(iloop,1)+1)*(2*sig(iloop,2)+1));
		% ppp_temp(iloop,3)=0.047*p_temp(iloop)/((6*sig(iloop,1)+1)*(6*sig(iloop,2)+1)-(4*sig(iloop,1)+1)*(4*sig(iloop,2)+1));
	% end
	% diff=cell(numrand,3);
	
	
	
	% for iloop=1:numrand
		% p_temp(iloop)=p_temp(iloop)/((2*sig(iloop,1)+1)*(2*sig(iloop,2)+1));
	% end
	% order_p=ones(numrand,1);
	% for iloop=1:numrand
		% [max_num, max_idx]=max(p_temp(:));
		% order_p(iloop)=max_idx;
		% p_temp(max_idx)=0;
	% end
	
	
	% order_ppp=ones(3*numrand,2);
	% for iloop=1:3*numrand
		% [max_num, max_idx]=max(ppp_temp(:));		
		% [maxx,maxy]=ind2sub(size(ppp_temp),max_idx);
		% order_ppp(iloop,1)=maxx;
		% order_ppp(iloop,2)=maxy;
		% ppp_temp(maxx,maxy)=0;
	% end
	
	
	% for iloop1=1:numrand
		% diff{iloop1,1}=zeros((2*sig(iloop1,1)+1)*(2*sig(iloop1,2)+1),2);
		% iloop=1;
		% for iloop3=-1*sig(iloop1,1):sig(iloop1,1)
			% for iloop4=-1*sig(iloop1,2):sig(iloop1,2)
				% diff{iloop1,1}(iloop,1)=iloop3;
				% diff{iloop1,1}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end
	% end
	% %2\delta
	% for iloop1=1:numrand		
		% diff{iloop1,2}=zeros(((4*sig(iloop1,1)+1)*(4*sig(iloop1,2)+1))-((2*sig(iloop1,1)+1)*(2*sig(iloop1,2)+1)),2);
		% iloop=1;
		% for iloop3=-2*sig(iloop1,1):-1*sig(iloop1,1)-1
			% for iloop4=-2*sig(iloop1,2):2*sig(iloop1,2)
				% diff{iloop1,2}(iloop,1)=iloop3;
				% diff{iloop1,2}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end		
		% for iloop3=-1*sig(iloop1,1):sig(iloop1,1)
			% for iloop4=-2*sig(iloop1,2):-1*sig(iloop1,2)-1
				% diff{iloop1,2}(iloop,1)=iloop3;
				% diff{iloop1,2}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end
		% for iloop3=-1*sig(iloop1,1):sig(iloop1,1)
			% for iloop4=sig(iloop1,2)+1:2*sig(iloop1,2)
				% diff{iloop1,2}(iloop,1)=iloop3;
				% diff{iloop1,2}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end	
		% for iloop3=sig(iloop1,1)+1:2*sig(iloop1,1)
			% for iloop4=-2*sig(iloop1,2):2*sig(iloop1,2)
				% diff{iloop1,2}(iloop,1)=iloop3;
				% diff{iloop1,2}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end		
	% end
	% %3\delta
	% for iloop1=1:numrand		
		% %diff=[3 0;3 1; 3 2; 3 3; 2 3; 1 3;0 3; -1 3; -2 3; -3 3; -3 2; -3 1; -3 0; -3 -1; -3 -2;-3 -3; -2 -3; -1 -3; 0 -3; 1 -3; 2 -3; 3 -3; 3 -2; 3 -1];
		% diff{iloop1,3}=zeros(((6*sig(iloop1,1)+1)*(6*sig(iloop1,2)+1))-((4*sig(iloop1,1)+1)*(4*sig(iloop1,2)+1)),2);
		% iloop=1;
		% for iloop3=-3*sig(iloop1,1):-2*sig(iloop1,1)-1
			% for iloop4=-3*sig(iloop1,2):3*sig(iloop1,2)
				% diff{iloop1,3}(iloop,1)=iloop3;
				% diff{iloop1,3}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end		
		% for iloop3=-2*sig(iloop1,1):2*sig(iloop1,1)
			% for iloop4=-3*sig(iloop1,2):-2*sig(iloop1,2)-1
				% diff{iloop1,3}(iloop,1)=iloop3;
				% diff{iloop1,3}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end
		% for iloop3=-2*sig(iloop1,1):2*sig(iloop1,1)
			% for iloop4=2*sig(iloop1,2)+1:3*sig(iloop1,2)
				% diff{iloop1,3}(iloop,1)=iloop3;
				% diff{iloop1,3}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end		
		% for iloop3=2*sig(iloop1,1)+1:3*sig(iloop1,1)
			% for iloop4=-3*sig(iloop1,2):3*sig(iloop1,2)
				% diff{iloop1,3}(iloop,1)=iloop3;
				% diff{iloop1,3}(iloop,2)=iloop4;
				% iloop=iloop+1;
			% end
		% end
	% end
	
	% y_temp=y;
	% t=0;
	% tfinal=-1;
	% loc_follow_maximum_prev=[1,1];
	% loc_follow_maximum_next=[1,1];
	
	% for iloop1=1:3*numrand
		% for iloop3=1:size(diff{order_ppp(iloop1,1),order_ppp(iloop1,2)})
			% loc_follow_maximum_next(1)=mu(order_ppp(iloop1,1),1)+a+1+diff{order_ppp(iloop1,1),order_ppp(iloop1,2)}(iloop3,1);
			% loc_follow_maximum_next(2)=mu(order_ppp(iloop1,1),2)+a+1+diff{order_ppp(iloop1,1),order_ppp(iloop1,2)}(iloop3,2);
			% maxdis=sqrt((loc_follow_maximum_prev(1)-loc_follow_maximum_next(1))^2+(loc_follow_maximum_prev(2)-loc_follow_maximum_next(2))^2);		
			% t=t+maxdis;
			% if((loc_follow_maximum_next(1)>0) && (loc_follow_maximum_next(2)>0)&& (loc_follow_maximum_next(1)<2*a+2)&& (loc_follow_maximum_next(2)<2*a+2) &&(y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))~=0))
				% ptemp=y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2));
				% y_temp(loc_follow_maximum_next(1), loc_follow_maximum_next(2))=0;
				% y_temp=y_temp*(1-ptemp);
			% end	
			% loc_follow_maximum_prev(1)=loc_follow_maximum_next(1);
			% loc_follow_maximum_prev(2)=loc_follow_maximum_next(2);	
			% if((loc_follow_maximum_prev(2)==R(1,2)) && (loc_follow_maximum_prev(1)==R(1,1)))
				% tfinal=t;
			% end
		% end
	% end
	
	
end
end
i=3;