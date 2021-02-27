function [K,G,para] = P_paraset(FUNC,MOEA,M)
%P_paraset ����ָ�����⡢ָ���㷨�Ĳ�������
%
% ����:
%       ��FUNC:  �ַ���.��������,������'DTLZ1'~'DTLZ7'�Լ�'WFG1'~'WFG9'����һ��
%       ��MOEA:  �ַ���.�㷨����,Ŀǰ֧��'GrEA','KnEA','HypE','MOEAD','PICEA'�Ĳ�������
%       ��M:     ��ֵ.�����ά��
%
% ���:
%       ��K:     ��ֵ��һά����.�������ز���(����)
%       ��G:     ��ֵ.�㷨�ĵ�������
%       ��para:  ��ֵ.�㷨����ز���

    D=str2double(FUNC(end));
    FUNC=FUNC(1:end-1);
    switch FUNC
        case 'DTLZ'
            K=[5 10 10 10 10 10 20];
            K=K(D);
            G=[1000 300 1000 300 300 1000 300];
            G=G(D);
            switch MOEA
                case 'GrEA'
                    para=[55 16 10 10 10  0 10  0 11
                          45 15 10  9  8  0  7  0  8
                          45 17 11 11 11  0 10  0 11
                          55 15 10  9  8  0  7  0  8
                          55 40 35 29 14  0 11  0 11
                          55 33 36 24 50  0 50  0 50
                          16 15  9  8  6  0  5  0  4];
                    para=para(D,M-1);
                case 'KnEA'
                    para=[.6 .6 .6 .2 .2  0 .1  0 .1
                          .6 .5 .5 .5 .5  0 .5  0 .5
                          .6 .6 .4 .3 .2  0 .1  0 .1
                          .6 .5 .5 .5 .5  0 .5  0 .5
                          .6 .6 .5 .5 .5  0 .3  0 .3
                          .6 .6 .5 .4 .4  0 .3  0 .3
                          .6 .6 .6 .6 .6  0 .5  0 .6];
                    para=para(D,M-1);
                case 'HypE'
                    para=10000;
                case 'MOEAD'
                    para=10;
                case 'PICEA'
                    para=M*100;
                case 'NSGAIII'
                    para=[99 12  8  6  4  0  3  0  3
                           0  0  0  0  2  0  2  0  2];
                    para=para(:,M-1);
            end
        case 'WFG'
            K=[4 4 6 8 10 0 7 0 9];
            K=K(M-1);
            K=[K,10];
            G=[1000 1000 300 300 300 300 300 300 300];
            G=G(D);
            switch MOEA
                case 'GrEA'
                    para=[45 10  8  9  9  0  7  0 10
                          45 12 11 11 11  0 11  0 11
                          55 22 18 18 18  0 16  0 22
                          45 15 10  9  9  0  7  0  8
                          45 15 10  9  9  0  7  0  8
                          45 15 10  9  9  0  7  0  8
                          45 15 10  9  9  0  7  0  8
                          45 15 10  9  9  0  7  0  8
                          45 15 10  9  9  0  7  0  8];
                    para=para(D,M-1);
                case 'KnEA'
                    para=[.6 .5 .5 .5 .5  0 .5  0 .5];
                    para=para(M-1);
                case 'HypE'
                    para=10000;
                case 'MOEAD'
                    para=10;
                case 'PICEA'
                    para=M*100;
                case 'NSGAIII'
                    para=[99 12  8  6  4  0  3  0  3
                           0  0  0  0  2  0  2  0  2];
                    para=para(:,M-1);
            end
    end
end

