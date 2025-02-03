function sumKahan = kahanSumMatrix(matrix)
    sumKahan = 0;
    compensation = 0;

    elements = matrix(:);

    for i = 1:length(elements)
        y = elements(i) - compensation;
        t = sumKahan + y;
        compensation = (t - sumKahan) - y;
        sumKahan = t;
    end
end
