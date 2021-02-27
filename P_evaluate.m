function metric = P_evaluate(operation,functionvalue,parameter1)
% ��P_evaluate ���ض��ڽ⼯functionvalue�ڲ�ͬ����ָ���µ�����ֵ
%              ����������(CM),ת���������(IGD),������(C),�����(HV),���ָ��(SM),�������(GD)��ƽ��Hausdorff����(delta)
%              CM��GD�������۽��������
%              SM�������۽�ķֲ���
%              IGD,HV��delta�����ۺ����۽�������Ժͷֲ���
%              C�������۴˽⼯�������һ���⼯�ĺû�
%              ����CM,IGD,SM,GD,deltaֵԽСԽ��,C��HVֵԽ��Խ��
%
% ����:
%       ��operation:     �ַ���.����ָ������������,��ȡֵ����:
%                        1.'CM',��������������
%                        2.'IGD',��������ת���������
%                        3.'C',�������۸�����
%                        4.'HV',�������۳����
%                        5.'SM',�������ۼ��ָ��
%                        6.'GD',�������۴������
%                        7.'delta',��������ƽ��Hausdorff����
%       ��functionvalue: ��ά����.�����۵Ľ⼯,����ÿһ�ж�Ӧһ������,ÿһ�ж�ӦһάĿ�꺯��ֵ
%       ��parameter1:    һά�������ά����.��������:
%                        1.��operationȡֵΪ'CM'ʱ,Ϊ��ʵǰ�����ϵľ���ȡ���㼯��(��ά����),����ÿһ�ж�Ӧһ������,ÿһ�ж�ӦһάĿ�꺯��ֵ
%                        2.��operationȡֵΪ'IGD'ʱ,ͬ1
%                        3.��operationȡֵΪ'C'ʱ,Ϊ��һ���⼯(��ά����),����ÿһ�ж�Ӧһ������,ÿһ�ж�ӦһάĿ�꺯��ֵ
%                        4.��operationȡֵΪ'HV'ʱ,Ϊ�ο�������(һά����),����ÿһλ��Ӧһά����ֵ
%                        5.��operationȡֵΪ'SM'ʱ,ʡ��
%                        6.��operationȡֵΪ'GD'ʱ,ͬ1
%                        7.��operationȡֵΪ'delta'ʱ,ͬ1
%
% ���:
%       ��metric:        ��ֵ.��������:
%                        1.��operationȡֵΪ'CM'ʱ,Ϊ�⼯functionvalue�������ʵǰ����parameter1��������ֵ
%                        2.��operationȡֵΪ'IGD'ʱ,Ϊ�⼯functionvalue�������ʵǰ����parameter1��ת���������ֵ
%                        3.��operationȡֵΪ'C'ʱ,Ϊ�⼯functionvalue֧��⼯parameter1�и���İٷֱ�
%                        4.��operationȡֵΪ'HV'ʱ,Ϊ�⼯functionvalue����ڲο���parameter1�ĳ����ֵ
%                        5.��operationȡֵΪ'SM'ʱ,Ϊ�⼯functionvalue�Ľ���ָ��ֵ
%                        6.��operationȡֵΪ'GD'ʱ,Ϊ�⼯functionvalue�������ʵǰ����parameter1�Ĵ������
%                        7.��operationȡֵΪ'delta'ʱ,Ϊ�⼯functionvalue�������ʵǰ����parameter1��ƽ��Hausdorff����

    switch operation
        %������
        case 'CM'
            popnum=size(functionvalue,1);
            truenum=size(parameter1,1);
            dis=zeros(popnum,truenum);
            fmax=max(parameter1);
            fmin=min(parameter1);
            for i=1:popnum
                for j=1:truenum
                    dis(i,j)=sqrt(sum(((functionvalue(i,:)-parameter1(j,:))./(fmax-fmin)).^2));
                end
            end
            dis=min(dis,[],2);
            metric=mean(dis);
        %ת���������
        case 'IGD'
            popnum=size(functionvalue,1);
            truenum=size(parameter1,1);
            dis=zeros(truenum,popnum);
            for i=1:truenum
                for j=1:popnum
                    dis(i,j)=sqrt(sum((functionvalue(j,:)-parameter1(i,:)).^2));
                end
            end
            dis=min(dis,[],2);
            metric=mean(dis);
        %������
        case 'C'
            temp=zeros(1,size(parameter1,1));
            for i=1:size(parameter1,1)
                for j=1:size(functionvalue,1)
                    k=any(functionvalue(j,:)-parameter1(i,:)<0)-any(functionvalue(j,:)-parameter1(i,:)>0);
                    if k==1
                        temp(i)=1;
                        break
                    end
                end
            end
            metric=sum(temp)/size(parameter1,1);
        %�����
        case 'HV'
            functionvalue(sum(functionvalue-repmat(parameter1,size(functionvalue,1),1)<=0,2)<size(functionvalue,2),:)=[];
            if isempty(functionvalue)
                metric=0;
                return
            end
            if size(functionvalue,2)<5
                n = size(functionvalue,2);
                pl = sortrows(functionvalue,1);
                S(1,1) = {1};
                S(1,2) = {pl};
                for k =1:n-1
                    S_ = {};
                    j = size(S,1);
                    for l = 1:j
                        S_l = cell2mat(S(l,2));
                        S_temp = Slice(S_l,k,parameter1);
                        p = size(S_temp,1);
                        for q=1:p
                            cell_(1,1) = {cell2mat(S_temp(q,1))*cell2mat(S(l,1))};
                            cell_(1,2) = S_temp(q,2);
                            S_ = Add(cell_,S_);
                        end
                    end
                    S = S_;
                end
                vol = 0;
                k = size(S,1);
                for l = 1 : k
                    p = Head(cell2mat(S(l,2)));
                    num = abs(p(n)-parameter1(n));
                    vol = vol + cell2mat(S(l,1))*num;
                end
            else
                k=1000000;
                maxvalue=parameter1;
                minvalue=min(functionvalue,[],1);
                samples=repmat(minvalue,k,1)+rand(k,size(functionvalue,2)).*repmat((maxvalue-minvalue),k,1);
                m=zeros(1,k);
                for i=1:size(functionvalue,1)
                    m(sum(repmat(functionvalue(i,:),k,1)-samples<=0,2)==size(functionvalue,2))=1;
                end
                vol=prod(maxvalue-minvalue)*sum(m)/k;
            end
            metric=vol;
        %���ָ��
        case 'SM'
            popnum=size(functionvalue,1);
            dis=zeros(popnum);
            for i=1:popnum
                dis(i,i)=inf;
                for j=i+1:popnum
                    dis(i,j)=sum(abs(functionvalue(i,:)-functionvalue(j,:)));
                    dis(j,i)=dis(i,j);
                end
            end
            dis=min(dis,[],2);
            metric=sqrt(1/(popnum-1)*sum((mean(dis)-dis).^2));
        %�������
        case 'GD'
            metric=P_evaluate('IGD',parameter1,functionvalue);
        %ƽ��Hausdorff����
        case 'delta'
            metric=max(P_evaluate('IGD',functionvalue,parameter1),P_evaluate('GD',functionvalue,parameter1));
    end
end

%���������ڼ��㳬���(��ȷ)�ĸ�������
function [ S ] = Slice( pl , k ,refPoint)
    % ����ÿ����Ƭ�ϵĵ㼯
    p = Head(pl);
    pl = Tail(pl);
    ql = [];
    S = {};
    while ~isempty(pl)
        ql = Insert(p , k+1 ,ql);
        p_ = Head(pl);
        cell_(1,1) = {abs(p(k) - p_(k))};
        cell_(1,2) = {ql};
        S = Add(cell_,S);
        p = p_;
        pl = Tail(pl);
    end
    ql = Insert(p , k+1 , ql);
    cell_(1,1) = {abs(p(k) - refPoint(k))};
    cell_(1,2) = {ql};
    S = Add(cell_,S);
end

function [ ql ] = Insert( p,k,pl )
    % �㼯���������
    flag1 = 0;
    flag2 = 0;
    ql = [];
    hp = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = Append( hp , ql );
        pl = Tail(pl);
        hp = Head(pl);
    end
    ql = Append( p , ql );
    m = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k:m
            if p(i) < q(i)
                flag1 = 1;
            else
                if p(i) > q(i)
                    flag2 = 1;
                end
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = Append(Head(pl) , ql);
        end
        pl = Tail(pl);
    end  
end

function [ p ] = Head( pl )
    % ��ȡ�㼯�е�һ����
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end
end

function [ ql ] = Tail( pl )
    % ��ȡ�㼯�г���һ���������еĵ㼯��
    m = size(pl,1);
    if m == 0 || m == 1
        ql = [];
    else
        ql = pl(2:m,:);
    end
end

function [ S_ ] = Add( cell_ , S )
    % ���µ�Sֵ���벢����
    n = size(S,1);
    m = 0;
    for k=1:n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1)) + cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;     
end

function [ pl ] = Append( p , ql )
    % ���������׷���µ�
    m = size(ql,1);
    ql(m+1,:) = p;
    pl = ql;
end