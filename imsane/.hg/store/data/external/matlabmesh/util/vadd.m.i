         j   i       ��������$vs来:D�s�?HP��4            ufunction [ Mout ] = vadd( M, v )
%VADD adds v to each row of M

Mout = M + repmat(v, size(M,1), 1);

end
