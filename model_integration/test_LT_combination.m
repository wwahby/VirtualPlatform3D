
Ladd = 0;
Lnom = 1313;
Tnom = 87;

slack = 0.1;

L = Lnom + Ladd;
[L T N gf_L gf_T] = find_LT_combination(L,Tnom,slack)