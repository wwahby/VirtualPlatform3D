function [w_tsv h_tsv] = size_tsvs(Ach, tsv_area_ratio, n_tsvs, AR_tsv )
%Lch_gp - Length of one side of the chip (assumed to be square) in gate pitches

if(n_tsvs == 0)
    w_tsv = 0;
    h_tsv = 0;
else
    Area_for_tsvs = Ach*tsv_area_ratio;

    A_tsv = Area_for_tsvs/n_tsvs;
    w_tsv = sqrt(A_tsv);
    h_tsv = w_tsv*AR_tsv;
end