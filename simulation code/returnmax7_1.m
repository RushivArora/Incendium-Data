function [maxtemp, maxsum, poslist] = returnmax7_1(x,pos,n,lambda)
    temp=-1*ones(size(x,1),size(x,2));
    %lambda=0;
    ditinguishvalue=-2;
	if(pos(1)-n>=1)
		i1=pos(1)-n;
	else
		i1=1;
	end
	if(pos(1)+n<=size(x,1))
		i2=pos(1)+n;
	else
		i2=size(x,1);
	end
	if(pos(2)-n>=1)
		j1=pos(2)-n;
	else
		j1=1;
	end
	if(pos(2)+n<=size(x,2))
		j2=pos(2)+n;
	else
		j2=size(x,2);
	end
	temp(i1:i2,j1:j2)=x(i1:i2,j1:j2);
	temp(pos(1),pos(2))=0;
	[max_num, max_idx]=max(temp(:));
	[maxx,maxy]=ind2sub(size(temp),max_idx);
	maxtemp=zeros(2,1);
	maxtemp(1)=maxx(1);
	maxtemp(2)=maxy(1);
	sumdir=zeros(8,1);
	%sumdir(1)=sum(temp(pos(1):i2,pos(2)));
    prob1=0;
	for indexi=pos(1)+1:i2
		sumdir(6)=sumdir(6)+(indexi-pos(1)+1)*temp(indexi,pos(2));
        prob1=prob1+temp(indexi,pos(2));
	end
    sumdir(6)=(sumdir(6)+lambda*(1-prob1));
    if(pos(1)==i2 || prob1==0)
        sumdir(6)=ditinguishvalue;
    end
	%sumdir(2)=sum(temp(i1:pos(1),pos(2)));
    prob1=0;
	for indexi=i1:pos(1)-1
		sumdir(2)=sumdir(2)+(indexi-i1+1)*temp(indexi,pos(2));
        prob1=prob1+temp(indexi,pos(2));
	end
    sumdir(2)=(sumdir(2)+lambda*(1-prob1));
    if(pos(1)==i1  || prob1==0)
        sumdir(2)=ditinguishvalue;
    end
	%sumdir(3)=sum(temp(pos(1),pos(2):j2));
    prob1=0;
	for indexi=pos(2)+1:j2
		sumdir(3)=sumdir(3)+(indexi-pos(2)+1)*temp(pos(1),indexi);
        prob1=prob1+temp(pos(1),indexi);
	end
    sumdir(3)=(sumdir(3)+lambda*(1-prob1));
    if(pos(2)==j2  || prob1==0)
        sumdir(3)=ditinguishvalue;
    end
	%sumdir(4)=sum(temp(pos(1),j1:pos(2)));
    prob1=0;
	for indexi=j1:pos(2)-1
		sumdir(4)=sumdir(4)+(indexi-j1+1)*temp(pos(1),indexi);
        prob1=prob1+temp(pos(1),indexi);
	end
    sumdir(4)=(sumdir(4)+lambda*(1-prob1));
    if(pos(2)==j1  || prob1==0)
        sumdir(4)=ditinguishvalue;
    end
	mini1j1=min(pos(1)-i1,pos(2)-j1);
    prob1=0;
	for indexi=1:mini1j1
		sumdir(5)=sumdir(5)+indexi*temp(pos(1)-indexi,pos(2)-indexi);
        prob1=prob1+temp(pos(1)-indexi,pos(2)-indexi);
    end
    sumdir(5)=(sumdir(5)+lambda*(1-prob1));
    if(mini1j1==0  || prob1==0)
        sumdir(5)=ditinguishvalue;
    end
	mini2j2=min(i2-pos(1),j2-pos(2));
    prob1=0;
	for indexi=1:mini2j2
		sumdir(1)=sumdir(1)+indexi*temp(pos(1)+indexi,pos(2)+indexi);
        prob1=prob1+temp(pos(1)+indexi,pos(2)+indexi);
    end
    sumdir(1)=(sumdir(1)+lambda*(1-prob1));
    if(mini2j2==0 || prob1==0)
        sumdir(1)=ditinguishvalue;
    end
	mini1j2=min(pos(1)-i1,j2-pos(2));
    prob1=0;
	for indexi=1:mini1j2
		sumdir(7)=sumdir(7)+indexi*temp(pos(1)-indexi,pos(2)+indexi);
        prob1=prob1+temp(pos(1)-indexi,pos(2)+indexi);
    end
    sumdir(7)=(sumdir(7)+lambda*(1-prob1));
    if(mini1j2==0 || prob1==0)
        sumdir(7)=ditinguishvalue;
    end
	mini2j1=min(i2-pos(1),pos(2)-j1);
    prob1=0;
	for indexi=1:mini2j1
		sumdir(8)=sumdir(8)+indexi*temp(pos(1)+indexi,pos(2)-indexi);
        prob1=prob1+temp(pos(1)+indexi,pos(2)-indexi);
    end
    sumdir(8)=(sumdir(8)+lambda*(1-prob1));
    if(mini2j1==0 || prob1==0)
        sumdir(8)=ditinguishvalue;
    end
	[max_num, max_idx]=min(sumdir(:));
	[maxx1,maxy1]=ind2sub(size(sumdir),max_idx);
	if(maxx1(1)==6)
        maxtemp(1)=pos(1)+1;
        maxtemp(2)=pos(2);
        poslist=zeros(i2-pos(1),2);
        for indexi=pos(1)+1:i2
            poslist(indexi-pos(1),1)=indexi;
            poslist(indexi-pos(1),2)=pos(2);
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)+1;                
           end
        end
	elseif(maxx1(1)==2)
		maxtemp(1)=pos(1)-1;
		maxtemp(2)=pos(2);
        poslist=zeros(pos(1)-i1,2);
        for indexi=i1:pos(1)-1
            poslist(indexi-i1+1,1)=indexi;
            poslist(indexi-i1+1,2)=pos(2);
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)-1;                
           end
        end
	elseif(maxx1(1)==3)
		maxtemp(1)=pos(1);
		maxtemp(2)=pos(2)+1;
        poslist=zeros(j2-pos(2),2);
        for indexi=pos(2)+1:j2
            poslist(indexi-pos(2),1)=pos(1);
            poslist(indexi-pos(2),2)=indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(2)=maxtemp(2)+1;                
           end
        end
	elseif(maxx1(1)==4)
		maxtemp(1)=pos(1);
		maxtemp(2)=pos(2)-1;
        poslist=zeros(pos(2)-j1,2);
        for indexi=j1:pos(2)-1
            poslist(indexi-j1+1,1)=pos(1);
            poslist(indexi-j1+1,2)=indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(2)=maxtemp(2)-1;                
           end
        end
	elseif(maxx1(1)==5)
		maxtemp(1)=pos(1)-1;
		maxtemp(2)=pos(2)-1;
        mini1j1=min(pos(1)-i1,pos(2)-j1);
        poslist=zeros(mini1j1,2);
        for indexi=1:mini1j1
            poslist(indexi,1)=pos(1)-indexi;
            poslist(indexi,2)=pos(2)-indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)-1;
                maxtemp(2)=maxtemp(2)-1;
           end
        end
	elseif(maxx1(1)==1)
		maxtemp(1)=pos(1)+1;
		maxtemp(2)=pos(2)+1;
        mini2j2=min(i2-pos(1),j2-pos(2));
        poslist=zeros(mini2j2,2);
        for indexi=1:mini2j2
            poslist(indexi,1)=pos(1)+indexi;
            poslist(indexi,2)=pos(2)+indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)+1;
                maxtemp(2)=maxtemp(2)+1;
           end
        end
	elseif(maxx1(1)==7)
		maxtemp(1)=pos(1)-1;
		maxtemp(2)=pos(2)+1;
        mini1j2=min(pos(1)-i1,j2-pos(2));
        poslist=zeros(mini1j2,2);
        for indexi=1:mini1j2
            poslist(indexi,1)=pos(1)-indexi;
            poslist(indexi,2)=pos(2)+indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)-1;
                maxtemp(2)=maxtemp(2)+1;
           end
        end
	elseif(maxx1(1)==8)
		maxtemp(1)=pos(1)+1;
		maxtemp(2)=pos(2)-1;
        mini2j1=min(i2-pos(1),pos(2)-j1);
        poslist=zeros(mini2j1,2);
        for indexi=1:mini2j1
            poslist(indexi,1)=pos(1)+indexi;
            poslist(indexi,2)=pos(2)-indexi;
        end
        if(max_num~=ditinguishvalue)
           while(temp(maxtemp(1),maxtemp(2))<=0)
                maxtemp(1)=maxtemp(1)+1;
                maxtemp(2)=maxtemp(2)-1;
           end
        end
	end
	maxsum=sumdir(maxx1(1));
end