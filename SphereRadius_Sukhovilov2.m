function [r, sigma_R, sigma_max] = SphereRadius_Sukhovilov2(S, sigma, sigma_m, r)

    % 2-й метод оценки радиуса c ошибками модели,
    % и ошибками измерения расстояний (расчет через собственные значения)
    % Построение матрицы S2
    C = 0.5 * S.^2;
    
    n = size(C,1);    
    options = optimset('Display', 'on');
    [r] = fzero(@(r) r_fun(r,n,C), r, options);

    % Сохраним результат для оценки статистики методом Монте_Карло
    % r2=[r2 r];

    % Вычисление СКО оценки радиуса
    [U,L]=eig((r .^2)*ones(n,n)-C);   % U-собственные вектора и L-собственные значения
    [~,ind]=sort(-abs(diag(L)));  % ф-ия 'sort' сортирует в порядке возрастания
                                % j-вектор первоначальных номеров отсортированных элементов
    % возмем 3 вектора, соответствующие трем наибольшим сообственным значениям
    U=U(:,[ind(1) ind(2) ind(3)]);
    summa=0;
    for i=1:3
      summa=summa+(sum(U(:,i))).^2;
    end 
    % summa

    kov=zeros(3,3);
    for i=1:3
     for j=1:3
      for k=1:n
       for l=1:n
        kov(i,j)=kov(i,j)+C(l,k)*(U(l,i)*U(k,i))*(U(l,j)*U(k,j));
       end 
      end
     end
    end

    % sumkov=2*sigma^2*(sum(sum(kov))); % сумма элементов ковариационной матрицы

    % Определение СКО оценки радиуса  - sigmar
    F=n-summa;
    sigma_R=sqrt((sigma_m^2)./F+(2*sigma^2*(sum(sum(kov))))/(4*r^2*F^2));

    % Верхняя граница СКО оценки радиуса
    sigma_max=sqrt(sigma_m^2/F+18*(sigma/F)^2);

end

function f = r_fun(r,n,C)
    rr=(r.^2)*ones(n,n);
    L=eig(rr-C);   %L-собственные значения
%     di=diag(L);
    [~,j]=sort(-abs(L));  % ф-ия 'sort' сортирует в порядке возрастания
                           % j-вектор первоначальных номеров отсортированных элементов
    f=L(j(1))+L(j(2))+L(j(3))-n*r^2;
end